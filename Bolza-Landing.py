import time
import os
import logging

import krpc
import numpy as np

import screen
from UPG import UPG, Status
import utils
import Apollo

logging.basicConfig(filename='upg.log',
                    filemode='w',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

open("failed_fun.csv", 'w+').close()

os.system('cls' if os.name == 'nt' else 'clear')
# Establish krpc connection
conn = krpc.connect(name="UPG")
print("krpc connected!")
print("krpc version: " + str(conn.krpc.get_status().version))
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)
specific_impulse = conn.add_stream(getattr, vessel, 'specific_impulse')
max_thrust = conn.add_stream(getattr, vessel, 'available_thrust')
current_thrust = conn.add_stream(getattr, vessel, 'thrust')
mass = conn.add_stream(getattr, vessel, 'mass')
situation = conn.add_stream(getattr, vessel, 'situation')
ut = conn.add_stream(getattr, space_center, 'ut')
height = conn.add_stream(getattr, flight, 'surface_altitude')
mean_altitude = conn.add_stream(getattr, flight, 'mean_altitude')
vertical_speed = conn.add_stream(getattr, flight, 'vertical_speed')
speed = conn.add_stream(getattr, flight, 'speed')
velocity = conn.add_stream(vessel.velocity, body_frame)
position = conn.add_stream(vessel.position, body_frame)
ap_err = conn.add_stream(getattr, ap, "error")


# body constants
r_e = body.equatorial_radius
g_0 = body.surface_gravity

########## target and parameters ##########
lon = -50.1
lat = 0
print("landing target:")
print("longitude: " + str(lon) + " latitude: " + str(lat))

# Final touch down speed
final_velocity = 3.
# The target height above target that guidance ends.
# A vertical fixed speed descent will be used for the final touch down.
# UPG will be set to 300m + this value, then use terminal to correct it
final_height = 50.
upg_add_h = 300.
# Minimum throttle. This must be exactly equal to your engine's min throttle.
# Alternatively you can set to 0 for a coast phase
min_throttle = 0.1
# Coefficient for bolza landing. The larger the closer, but the solver will more likely to fail.
# You don't need to be too accurate with guidance, terminal guidance will make you there.
k = 5e2
# First burn time, in secondes. This should be fixed and have little effect on performance.
t_1 = 10.
# The time you want to burn at most.
# Does not guranteed to below this value, only used for generate initial guess.
# You will want to tune this value for better performance
burn_time = 600
# Set index to 0 for bolza landing, 1 for pinpoint landing.
mode = ['bolza', 'pinpoint'][1]
print("Using mode: {0}".format(mode))
# fall back to bolza mode when solver fails
fall_back = True
# start terminal guidance when the solver still fails.
early_terminal = True
# The distance difference where the guidance start.
# For instance, the guidance start to work when the target downrange is less than 
# predicted downrange to target minus this value.
PDI_margin = 1e3
# if your throttle saturate, and thrust error is always negative, 
# then your time is too short. Increase this value to lower thrust.
time_coef = 1.3

# Terminal guidance related
# When to swtich to guidance, by default this is set to 0.02 * body's radius
terminal_downrange = 2 * 1e3
# Time till end when the teriminal guidance kicks in. Defaults to 0.1 * total time
terminal_time = 0.05
# Set this to False, if you don't need to be very accurate to allow less fuel usage, as
# the terminal guidance will use the final solution of upg as the target, not the specified one.
precise_mode = True

# If you find your craft turning too slow, alter the following values.
# The time for craft to stop rotating. Since the torque is constant, the
# large these values is, the fast you will rotate.
st = 0.9
ap.stopping_time = (st, st, st)
# ap.deceleration_time = (3., 3., 3.)
# This controls the precision. Smaller values will lead to unstablity.
aa = 0.8
ap.attenuation_angle = (aa, aa, aa)
########### parameter ends ##############


# target vector, with final height added.
target_height = body.surface_height(lat, lon)
r_target = utils.vector(body.position_at_altitude(
    lat, lon, target_height + final_height + upg_add_h, body_frame))

# get final velocity and up vector in body frame
up = utils.normalize(r_target)
v_target = -(final_velocity + 5) * up * 0.

# draw the target
target_frame = utils.target_reference_frame(space_center, body, lat, lon)
line = conn.drawing.add_line((0., 0., 0.), (20., 0., 0.), target_frame)
line.thickness = 2

### for debug: draw upg target
line = conn.drawing.add_line(tuple(utils.swap_yz(r_target).flatten()), 
                             tuple(utils.swap_yz(r_target + 20 * up).flatten()), body_frame)
line.thickness = 5

# get initial values
v_0 = utils.vector(velocity())
r_0 = utils.vector(position())
# initial guess
pv = -utils.normalize(v_0)
pr = np.zeros((3, 1))
t_f = np.array([[burn_time]])

print("Initiating UPG")
upg = UPG(r_0=r_e, g_0=g_0, r_target=r_target, v_target=v_target, r_t=r_0, v_t=v_0,
          mass=mass(), max_thrust=max_thrust(), min_throttle=min_throttle,
          specific_impulse=specific_impulse(), k=k, t_0=ut(), t_1=t_1, p_r_0=pr, p_v_0=pv,
          t_f=t_f, t_2_bound=[5, burn_time / 4.])  # limit t2 not to be too large, and not to be too small.

print("Evaluating t2...")
upg.solve_t_2()
print("Optimal t2: {:.2f} s".format(upg.get_t_2))
print("Initiating UPG...")
# Wait a bit before clear terminal
time.sleep(0.5)
os.system('cls' if os.name == 'nt' else 'clear')
screen.background()

ap.reference_frame = body_frame
ap.engage()
guidance_interval = 0.1

while True:
    v_0 = utils.vector(velocity())
    r_0 = utils.vector(position())
    upg.update_state(r_t=r_0, v_t=v_0, mass=mass(),
                     max_thrust=max_thrust(), t_0=ut())

    if upg.status == Status.PDI:
        upg.update_start_time(ut())
        upg.solve(mode='soft')
        dr_guidance = utils.downrange(r_0, upg.r_final)
        dr_target = utils.downrange(r_0, r_target)
        # The direction of soft landing may be different from the bolza/pinpoint one.
        # Solve the problem again before start to give the correct direction.
        # if dr_guidance > dr_target - 30. * PDI_margin:
        #     upg.solve(mode=mode)
        # start descent when overshoot the target with margin of 1km.
        if dr_guidance > dr_target + PDI_margin:
            upg.status = Status.Powered_descent
            m_initial = mass()
            # set terminal conditions
            terminal_time *= upg.t_f
            logging.debug("Using terminal conditions: downrange=%f, time=%f", terminal_downrange, terminal_time)
            continue

    elif upg.status == Status.Powered_descent:
        # move to before solve to avoid stucks in solver
        vessel.control.throttle = upg.throttle
        
        upg.solve(mode=mode)
        with open("failed_fun.csv", 'a') as f:
            np.savetxt(f, upg.fun.flatten(), newline=' ', delimiter=",")
            f.write('\n')
        if not upg.convergence:
            pass
            #logging.info("Failed to converge, the results are:\n" + str(upg.fun))
            

        # correct target to the correct height as the distance to target could be large.
        # if k == 0.:
        #     r_target = utils.move_position2height(
        #         final_height, upg.r_final, body, body_frame)
        #     up = utils.normalize(r_target)
        #     v_target = -final_velocity * up
        #     upg.set_target(r_target, v_target)
        #
        # Update time again as the solver may take a while
        upg.update_state(r_t=r_0, v_t=v_0, mass=mass(),
                     max_thrust=max_thrust(), t_0=ut())

        dr_target = utils.downrange(r_0, r_target)
        dr_guidance = utils.downrange(r_0, upg.r_final)
        if dr_target < terminal_downrange or upg.t_f < terminal_time:
            upg.status = Status.Terminal
            break
        if upg.last_convergence > 10:
            if mode =='pinpoint' and fall_back:
                mode = 'bolza'
            if mode == 'bolza' and early_terminal:
                logging.debug('Can\'t converge, switch to terminal.')
        # logging.debug("landing error: dr %d, height %d", utils.downrange(upg.r_final, r_target), 
        #               np.linalg.norm(upg.r_final)-np.linalg.norm(r_target))

    ap.target_direction = upg.thrust_direction
    screen.update_upg(upg.status, dr_target, dr_guidance, upg.t_f, upg.t_1_go, upg.t_2_go, upg.v_go,
                      upg.convergence, upg.norm, upg.solver_status, upg.last_convergence, ap_err(),
                      np.linalg.norm(upg.r_final - r_target))
    time.sleep(guidance_interval)

# Terminal guidance: use apollo guidance algorithm.

if k != 0 or mode == 'pinpoint':
    r_target = utils.vector(body.position_at_altitude(
        lat, lon, target_height + final_height, body_frame))
# Change to target to upg's prediction
landing_err = utils.downrange(upg.r_final, r_target)
if k == 0 or np.arctan2((mean_altitude() - target_height), landing_err) -\
        np.arcsin(-vertical_speed() / speed()) < np.radians(5)  or (not precise_mode) or landing_err > 2e3:
    # make another landing frame with the predicted final location
    r_target = utils.move_position2height(
        final_height, upg.r_final, body, body_frame)
    up = utils.normalize(r_target)
    new_lat = body.latitude_at_position(tuple(utils.swap_yz(r_target.flatten())), body_frame)
    new_lon = body.longitude_at_position(tuple(utils.swap_yz(r_target.flatten())), body_frame)
    target_height = body.surface_height(new_lat, new_lon)
    target_frame = utils.target_reference_frame(space_center, body, new_lat, new_lon)
    line = conn.drawing.add_line((0., 0., 0.), (20., 0., 0.), target_frame)
    line.thickness = 3

a_f = 2.  # final acceleration N g(g with respect to body)
k_r = 8.0

# calculate t_go
t_f = time_coef * upg.t_f
t_f = t_f + ut()    # convert to ut
apdg = Apollo.APDG(r_target, final_velocity, a_f, t_f, ut(),
                   mass(), r_0, v_0, g_0, k_r)
while True:
    v_0 = utils.vector(velocity())
    r_0 = utils.vector(position())
    apdg.update(r_0, v_0, mass(), ut())
    apdg.compute()

    vessel.control.throttle = utils.thrust2throttle(
        apdg.thrust, max_thrust(), min_throttle)
    ap.target_direction = apdg.thrust_direction

    r_go = np.linalg.norm(r_0 - r_target)

    dr_target = utils.downrange(r_0, r_target, body, body_frame)

    screen.update_terminal(dr_target, apdg.t_go,
                           current_thrust() - apdg.thrust, ap_err())
    # End early enough to avoid singularity
    if apdg.t_go < 5 or r_go < 10 or height() < final_height:
        break

    time.sleep(guidance_interval)

# Final touch down with gravity turn trajectory
while True:
    v_0 = utils.vector(velocity())
    r_0 = utils.vector(position())
    acc_dir, acc = Apollo.gravity_turn(vertical_speed(), speed(
    ), mean_altitude(), final_velocity, target_height, g_0, r_0, v_0)
    ap.target_direction = acc_dir
    thrust = utils.thrust2throttle(
        mass() * acc, max_thrust(), min_throttle)
    vessel.control.throttle = thrust
    screen.print_error(current_thrust() - thrust, ap_err(), None)
    if situation() == space_center.VesselSituation.landed \
            or situation() == space_center.VesselSituation.splashed\
            or flight.vertical_speed > -0.5:
        break
ap.disengage()
vessel.control.throttle = 0.
upg.status = Status.Finished

# Update the final landing error
r_0 = utils.vector(position())
target_height = body.surface_height(lat, lon)
r_target = utils.vector(body.position_at_altitude(
    lat, lon, target_height, body_frame))
l_err = np.linalg.norm(r_0 - r_target)
screen.print_error(0, 0, l_err)
screen.print_status(upg.status)
screen.end()

print("Delta-v used: {:.2f}".format(specific_impulse()
      * 9.80665 * np.log(m_initial / mass())))
