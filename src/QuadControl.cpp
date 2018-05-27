#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

namespace {
  template<class T>
  constexpr const T& clamp(const T& v, const T& lo, const T& hi)
  {
    return v < lo ? lo : hi < v ? hi : v;
  }
  
  float norm(float x, float y) {
    return sqrt(x * x + y * y);
  }
}

void QuadControl::Init()
{
  BaseController::Init();
  
  // variables needed for integral control
  integratedAltitudeError = 0;
  
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
  
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config + ".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
  
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);
  
  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);
  
  kpPQR = config->Get(_config + ".kpPQR", V3F());
  
  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);
  
  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);
  
  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float desired_total_thrust, V3F desired_moment)
{
  // Convert a desired 3-axis moment and collective thrust command to
  //   individual motor thrust commands
  // INPUTS:
  //   desCollectiveThrust: desired collective thrust [N]
  //   desMoment: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]
  
  // HINTS:
  // - you can access parts of desMoment via e.g. desMoment.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa
  
  float l_x = L / sqrt(2);
  float l_y = L / sqrt(2);
  
  float f_total = -desired_total_thrust;
  float f_x = desired_moment.x / l_x;
  float f_y = desired_moment.y / l_y;
  float f_z = desired_moment.z / kappa;
  
  cmd.desiredThrustsN[0] = -(f_total - f_x - f_y + f_z) / 4.0f;
  cmd.desiredThrustsN[1] = -(f_total + f_x - f_y - f_z) / 4.0f;
  cmd.desiredThrustsN[2] = -(f_total - f_x + f_y - f_z) / 4.0f;
  cmd.desiredThrustsN[3] = -(f_total + f_x + f_y + f_z) / 4.0f;
  
  for (int i = 0; i < 4; ++i) {
    cmd.desiredThrustsN[i] = clamp(cmd.desiredThrustsN[i], minMotorThrust, maxMotorThrust);
  }
  
  return cmd;
}

V3F QuadControl::BodyRateControl(V3F desired_pqr, V3F current_pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS:
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes
  
  // HINTS:
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)
  
  //return V3F();
  V3F error = desired_pqr - current_pqr;
  V3F output;
  output.x = error.x * kpPQR.x * Ixx;
  output.y = error.y * kpPQR.y * Iyy;
  output.z = error.z * kpPQR.z * Izz;
  
  float scale = 1.0f;
  {
    const float max_moment_x = 2.0f * maxMotorThrust * ArmLengthX();
    const float abs_moment_x = fabs(output.x);
    if (abs_moment_x > max_moment_x) scale = std::min(max_moment_x / abs_moment_x, scale);
  }
   
  {
    const float max_moment_y = 2.0f * maxMotorThrust * ArmLengthY();
    const float abs_moment_y = fabs(output.y);
    if (abs_moment_y > max_moment_y) scale = std::min(max_moment_y / abs_moment_y, scale);
  }
   
  {
    const float max_moment_z = 2.0f * maxMotorThrust * kappa;
    const float abs_moment_z = fabs(output.z);
    if (abs_moment_z > max_moment_z) scale = std::min(max_moment_z / abs_moment_z, scale);
  }
  
  return output * scale;
}

float QuadControl::ArmLengthX() const {
  return L / sqrt(2.0);
}

float QuadControl::ArmLengthY() const {
  return L / sqrt(2.0);
}

// returns a desired roll and pitch rate
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS:
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)
  
  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first
  
  //return V3F();
  
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  
  float collectiveAccel = -collThrustCmd / mass;
  float xyAccel = accelCmd.magXY();
  float absCollAccel = fabs(collectiveAccel);
  
  if (xyAccel > absCollAccel) {
    accelCmd *= absCollAccel / xyAccel;
  }
  
  const float target_roll = clamp(accelCmd.x / collectiveAccel, -maxTiltAngle, maxTiltAngle);
  const float target_pitch = clamp(accelCmd.y / collectiveAccel, -maxTiltAngle, maxTiltAngle);

  const float current_roll = R(0, 2);
  const float current_pitch = R(1, 2);
  
  const float roll_error = target_roll - current_roll;
  const float pitch_error = target_pitch - current_pitch;
  
  const float roll_dot = roll_error * kpBank;
  const float pitch_dot = pitch_error * kpBank;

  const float p = (roll_dot * R(1,0) - pitch_dot * R(0, 0)) / R(2, 2);
  const float q = (roll_dot * R(1,1) - pitch_dot * R(0, 1)) / R(2, 2);
  
  return V3F(p, q, /*r=*/0.0f);
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical
  //   acceleration feed-forward command
  // INPUTS:
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]
  
  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER
  
  //return 0;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  
  velZCmd = clamp(velZCmd, -maxAscentRate, maxDescentRate);
  const float positionError = posZCmd - posZ;
  const float velocityError = velZCmd - velZ;
  integratedAltitudeError += positionError * dt;
  
  const float output = positionError * kpPosZ +
      velocityError * kpVelZ +
      integratedAltitudeError * KiPosZ +
      accelZCmd;
  
  const float min_collective_thrust = 4.0f * minMotorThrust;
  const float max_collective_thrust = 4.0f * maxMotorThrust;
  const float output_thrust = (CONST_GRAVITY - output) / R(2, 2) * mass;
  return clamp(output_thrust, min_collective_thrust, max_collective_thrust);
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS:
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations.
  //     the Z component should be 0
  // HINTS:
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY
  
  const float error_pos_x = posCmd.x - pos.x;
  const float error_pos_y = posCmd.y - pos.y;
  
  {
    const float vel_xy = norm(velCmd.x, velCmd.y);
    if (vel_xy > maxSpeedXY) {
      const float scale = maxSpeedXY / vel_xy;
      velCmd.x *= scale;
      velCmd.y *= scale;
    }
  }
  
  const float error_vel_x = velCmd.x - vel.x;
  const float error_vel_y = velCmd.y - vel.y;
  
  float accel_x = kpPosXY * error_pos_x + kpVelXY * error_vel_x + accelCmdFF.x;
  float accel_y = kpPosXY * error_pos_y + kpVelXY * error_vel_y + accelCmdFF.y;
  
  {
    const float accel_xy = norm(accel_x, accel_y);
    if (accel_xy > maxAccelXY) {
      const float scale = maxAccelXY / accel_xy;
      accel_x *= scale;
      accel_y *= scale;
    }
  }
  
  return V3F(accel_x, accel_y, 0.0f);
}

// returns desired yaw rate
float QuadControl::YawControl(float desired_yaw, float current_yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS:
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS:
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b].
  //  - use the yaw control gain parameter kpYaw
  
  float error = desired_yaw - current_yaw;
  if (error > F_PI) {
    error -= 2.0f * F_PI;
  } else if (error < -F_PI) {
    error += 2.0f * F_PI;
  }
  const float result = error * kpYaw;
  return result;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);
  
  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);
  
  // reserve some thrust margin for angle control
  const float thrustMargin = (maxMotorThrust - minMotorThrust) / 10.0f;
  const float min_collective_thrust = (minMotorThrust + thrustMargin) * 4.0f;
  const float max_collective_thrust = (maxMotorThrust - thrustMargin) * 4.0f;
  collThrustCmd = clamp(collThrustCmd, min_collective_thrust, max_collective_thrust);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());
  
  V3F desMoment = BodyRateControl(desOmega, estOmega);
  
  return GenerateMotorCommands(collThrustCmd, desMoment);
}
