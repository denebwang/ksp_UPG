import time
import krpc
import json
import numpy as np
from scipy import optimize
import UPG
from Utilities import \
    target_reference_frame, thrust2throttle, move_position2height
import Apollo

conn = krpc.connect()
print("krpc connected")
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)

# load params
with open("target.json") as f:
    params = json.load(f)
lon = params['lon']
lat = params['lat']
v = params['vf']
min_throttle = params['min_throttle']
min_thrust = 0.
k = params['Bolza_k']
t1 = params['t1']
rated_burn_time = params['rated_burn_time']

# target vector
r_f = np.asarray(body.position_at_altitude(
    lat, lon, body.surface_height(lat, lon) + 300, body_frame)).reshape(-1, 1)

# get final velocity and up vector in body frame
target_frame = target_reference_frame(space_center, body, lat, lon)

v_f = np.asarray(space_center.transform_velocity(
    (0., 0., 0.), (-v, 0., 0.), target_frame, body_frame)).reshape(-1, 1)
up = np.asarray(space_center.transform_velocity(
    (0., 0., 0.), (1, 0, 0), target_frame, body_frame)).reshape(-1, 1)
up /= np.linalg.norm(up)

# draw the target
line = conn.drawing.add_line((0., 0., 0.), (20., 0., 0.), target_frame)
line.thickness = 3

# create data streams
vacuum_specific_impulse = conn.add_stream(getattr, vessel,
                                          'vacuum_specific_impulse')
max_thrust = conn.add_stream(getattr, vessel, 'max_thrust')
current_thrust = conn.add_stream(getattr, vessel, 'thrust')
mass = conn.add_stream(getattr, vessel, 'mass')
ut = conn.add_stream(getattr, space_center, 'ut')
height = conn.add_stream(getattr, flight, 'surface_altitude')
speed = conn.add_stream(getattr, flight, 'speed')
velocity = conn.add_stream(vessel.velocity, body_frame)
position = conn.add_stream(vessel.position, body_frame)

# get initial values
thrust = max_thrust()
specific_impulse = vacuum_specific_impulse()
v_0 = np.asarray(velocity()).reshape(-1, 1)
r_0 = np.asarray(position()).reshape(-1, 1)
t_0 = ut()

# initial guess
pv = np.asarray(-v_0 / np.linalg.norm(v_0)).reshape(-1, 1)
pr = np.zeros((3, 1))
tf = np.array([rated_burn_time]).reshape(1, 1)

print("initiating UPG")
upg = UPG.UPG(body, k=k, r_f=r_f, v_f=v_f, r_t=r_0, v_t=v_0, mass=mass(),
              max_thrust=thrust, min_throttle=0.,
              specific_impulse=specific_impulse, t_0=t_0, t_1=t1,
              pv0=pv, pr0=pr, t_f=tf)

print("evaluating t2...")
upg.solve_t_2()
print("optimal t2: {:.2f}".format(upg.t_2))

# PDI determination
print("waiting for powered descent initial")
ap.reference_frame = body_frame
ap.engage()
guidance_interval = 1
while True:
    v_0 = np.asarray(velocity()).reshape(-1, 1)
    r_0 = np.asarray(position()).reshape(-1, 1)
    upg.update(r_t=r_0, v_t=v_0, mass=mass(),
               max_thrust=max_thrust(), t_0=ut())
    if upg.status == UPG.Status.PDI:
        upg.update_time()
        upg.solve(mode='soft')
        r_tf = upg.r_tf
        r_guidance = np.linalg.norm(r_0 - r_tf)
        r_togo = np.linalg.norm(r_0 - r_f)
        print("guidance distance: {0:.2f}m, target distance: {1:.2f}m"
              .format(r_guidance, r_togo))
        if r_guidance > r_togo + 1000.:
            print("PDI, UPG engaged")
            vessel.control.throttle = 1
            guidance_interval = 0.2
            upg.status = UPG.Status.Powered_descent
            continue
        print("{:.2f}m togo".format(r_togo + 1000. - r_guidance))
    elif upg.status == UPG.Status.Powered_descent:
        upg.solve(mode='bolza')
        tgo = upg.t_f
        t_1_go = upg.t_1_go
        t_2_go = upg.t_2_go
        distance = np.linalg.norm(r_0 - r_f)
        v_go = upg.v_go
        print("\ntgo: {:.2f}".format(tgo), end=" ")
        if t_1_go > 0:
            print("t1go: {:.2f}".format(t_1_go), end=" ")
        if t_2_go > 0:
            print("t2go: {:.2f}".format(t_2_go), end=" ")
        print("\ndistance: {:.2f}".format(distance))
        print("vgo: {:.2f}".format(v_go))
        # correct target to the correct height
        if k == 0.:
            r_f = move_position2height(300., upg.r_tf, body, body_frame)
            upg.set_target(r_f, v_f)

        vessel.control.throttle = upg.throttle

        if tgo < 20.:
            print("UPG disengaged")
            break
    else:
        raise Exception("UPG in wrong status!")
    thrust_direction = upg.thrust_direction
    ap.target_direction = tuple(thrust_direction.flatten())
    time.sleep(guidance_interval)
print("Switching to Apollo terminal guidance...")
# APDG terminal guidance
# change r_f to ground, cancel downrange to ensure a safe land
r_f_ori = np.asarray(body.position_at_altitude(
    lat, lon, body.surface_height(lat, lon), body_frame)).reshape(-1, 1)

r_f = move_position2height(0, upg.r_tf, body, body_frame)

a_f = 0 * up
v_f = -1 * up
# calculate tgo
sin_gamma = flight.vertical_speed / speed()
V_0 = speed()
h_0 = height()
a_t = optimize.root_scalar(Apollo.fun_at, args=(V_0, sin_gamma, h_0, upg.g0),
                           method='brentq', bracket=(0, thrust / mass())).root
t_f = Apollo.t_togo(V_0, sin_gamma, a_t, upg.g0)
print('tgo: {:.2f}'.format(t_f))
t_f = t_f + ut()
apdg = Apollo.APDG(r_f, v_f, a_f, t_f, t_0, mass(), r_0, v_0)
print("Terminal guidance engaged")

while True:
    v_0 = np.asarray(velocity()).reshape(-1, 1)
    r_0 = np.asarray(position()).reshape(-1, 1)
    apdg.update(r_0, v_0, mass(), ut())
    apdg.compute()
    print("\nt_go: {:.2f}".format(apdg.t_go))
    r_go = np.linalg.norm(r_f - r_0)
    print("distance to land: {:.2f}".format(r_go))
    vessel.control.throttle = thrust2throttle(
        apdg.thrust, max_thrust(), min_throttle)
    ap.target_direction = tuple(apdg.thrust_direction.flatten())
    if apdg.t_go < 0.5 or r_go < 2 or height() < 2 or speed() < 0.:
        print("APDG disengaged")
        vessel.control.throttle = 0.
        break

    time.sleep(0.1)

while True:
    if vessel.situation == space_center.VesselSituation.landed \
            or vessel.situation == space_center.VesselSituation.splashed:
        print("landed!")
        break
ap.disengage()
vessel.control.sas = True
print("Landing precision: {:.2f}m".format(np.linalg.norm(r_0 - r_f_ori)))
