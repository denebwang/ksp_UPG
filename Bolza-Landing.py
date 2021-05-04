import time
import krpc
import json
import numpy as np
from scipy import optimize
import UPG
from Utilities import target_reference_frame
from Apollo import *

conn = krpc.connect()
print("krpc connected")
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)

upg = UPG.UPG(conn)

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
r_f = np.asarray(body.position_at_altitude(lat, lon,
                                           body.surface_height(lat, lon) + 300, body_frame)).reshape(-1, 1)
# get final velocity and up vector in body frame
target_frame = target_reference_frame(space_center, body, lat, lon)
target_v = (-v, 0., 0.)
v_f = np.asarray(space_center.transform_velocity((0., 0., 0.), target_v, target_frame, body_frame)).reshape(-1, 1)
up = np.asarray(space_center.transform_velocity((0., 0., 0.), (1, 0, 0), target_frame, body_frame)).reshape(-1, 1)
up /= np.linalg.norm(up)

# normalize to dimensionless
r_f = r_f / upg.position_multiplier
v_f = v_f / upg.velocity_multiplier
t1 = t1 / upg.time_multiplier

# create data streams
vacuum_specific_impulse = conn.add_stream(getattr, vessel, 'vacuum_specific_impulse')
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
m0 = mass()
v_0 = np.asarray(velocity()).reshape(-1, 1) / upg.velocity_multiplier
r_0 = np.asarray(position()).reshape(-1, 1) / upg.position_multiplier
x = np.vstack((r_0, v_0))
t_0 = ut() / upg.time_multiplier
t_1 = t_0 + t1

# initial guess
pv = np.asarray(-v_0 / np.linalg.norm(v_0)).reshape(-1, 1)
pr = np.zeros((3, 1))
tf = np.array([(ut() + rated_burn_time) / upg.time_multiplier]).reshape(1, 1)
z0 = np.hstack((pv.T, pr.T, tf)).T

print("evaluating t2...")
t_2 = optimize.minimize_scalar(
    upg.target_function_t2_sl, bounds=(t_1, tf[0, 0]),
    args=(z0, r_f, v_f, x, t_0, t_1, thrust, min_throttle, specific_impulse, m0), method='Bounded').x
# convert to time interval
t2 = t_2 - t_0
print("optimal t2: {:.2f}".format(t2 * upg.time_multiplier))

# PDI determination
print("initiating UPG")
print("waiting for powered descent initial")
ap.reference_frame = body_frame
ap.engage()
while True:
    x, m0, thrust, t_0 = upg.update(velocity, position, mass, max_thrust, ut)
    t_1 = t_0 + t1
    t_2 = t_0 + t2
    sol = optimize.root(upg.target_function_sl, z0,
                        args=(r_f, v_f, x, t_0, t_1, t_2, thrust, min_thrust, specific_impulse, m0),
                        method='lm', jac=False)
    z = sol.x.reshape(-1, 1)
    z0 = z
    pv = z[0:3, :]
    unit_t = pv / np.linalg.norm(pv)
    ap.target_direction = tuple(unit_t.flatten())
    # check the PDI condition
    x_f = upg.x_f
    r_f_t = x_f[0:3, :]
    r_0 = x[0:3, :]
    r_guidance = np.linalg.norm(r_0 - r_f_t) * upg.position_multiplier
    r_togo = np.linalg.norm(r_0 - r_f) * upg.position_multiplier
    print("guidance distance: {0:.2f}m, target distance: {1:.2f}m".format(r_guidance, r_togo))
    if r_guidance > r_togo + 1000.:
        print("PDI")
        break
    time.sleep(1)

# UPG Guidance
t_1 = t_0 + t1
t_2 = t_0 + t2
vessel.control.throttle = 1
print("UPG engaged")
while True:
    x, m0, thrust, t_0 = upg.update(velocity, position, mass, max_thrust, ut)

    r_0 = x[0:3, :]
    sol = optimize.root(
        upg.target_function_bl, z0,
        args=(k, r_f, v_f, x, t_0, t_1, t_2, thrust, min_thrust, specific_impulse, m0), jac=False)
    z = sol.x.reshape(-1, 1)
    z0 = z
    pv = z[0:3, :]
    tf = z[6, 0]
    tgo = (tf - t_0) * upg.time_multiplier
    tgo_1 = (t_1 - t_0) * upg.time_multiplier
    tgo_2 = (t_2 - t_0) * upg.time_multiplier
    distance = np.linalg.norm((r_0 - r_f) * upg.position_multiplier)
    v_go = upg.v_go(specific_impulse, m0)
    print("tgo: {:.2f}".format(tgo))
    if tgo_1 > 0:
        print("tgo1: {:.2f}".format(tgo_1))
    if tgo_2 > 0:
        print("tgo2: {:.2f}".format(tgo_2))
    print("distance: {:.2f}".format(distance))
    print("vgo: {:.2f}".format(v_go))
    unit_t = pv / np.linalg.norm(pv)
    ap.target_direction = tuple(unit_t.flatten())
    if k == 0.:
        # correct target to the correct height
        x_f = upg.x_f
        r_f_t = tuple(upg.position_multiplier * x_f[0:3, :].flatten())
        lat = body.latitude_at_position(r_f_t, body_frame)
        lon = body.longitude_at_position(r_f_t, body_frame)
        r_f = np.asarray(body.position_at_altitude(
            lat, lon, body.surface_height(lat, lon) + 300, body_frame)).reshape(-1, 1) / upg.position_multiplier
    if t_1 < t_0 < t_2:
        vessel.control.throttle = min_thrust
    else:
        vessel.control.throttle = 1
    if tgo < 20.:
        print("UPG disengaged")
        break
    time.sleep(0.2)

# ADGP terminal guidance
# change r_f to ground, cancel downrange to ensure a safe land
r_f_ori = np.asarray(body.position_at_altitude(
    lat, lon, body.surface_height(lat, lon), body_frame)).reshape(-1, 1)
x_f = upg.x_f
r_f_t = tuple(upg.position_multiplier * x_f[0:3, :].flatten())
lat = body.latitude_at_position(r_f_t, body_frame)
lon = body.longitude_at_position(r_f_t, body_frame)
r_f = np.asarray(body.position_at_altitude(
    lat, lon, body.surface_height(lat, lon), body_frame)).reshape(-1, 1)

a_f = upg.g0 * up
v_f = np.zeros((3, 1))
# calculate tgo
sin_gamma = flight.vertical_speed / speed()
V_0 = speed()
h_0 = height()
a_t = optimize.root_scalar(fun_at, args=(V_0, sin_gamma, h_0, upg.g0),
                           method='brentq', bracket=(0, thrust / m0)).root
t_f = t_togo(V_0, sin_gamma, a_t, upg.g0)
print('tgo: {:.2f}'.format(t_f))
t_f = t_f + ut()
print("Terminal guidance engaged")
while True:
    time.sleep(0.1)

    v_0 = np.asarray(velocity()).reshape(-1, 1)
    r_0 = np.asarray(position()).reshape(-1, 1)
    m0 = mass()
    thrust = max_thrust()
    t_0 = ut()
    t_go = t_f - t_0
    print("t_go: {:.2f}".format(t_go))
    r_go = np.linalg.norm(r_f - r_0)
    print("distance to land: {:.2f}".format(r_go))
    a_t = (-6 / t_go) * (v_f - v_0) + (12 / t_go ** 2) * (r_f - r_0 - v_0 * t_go) + a_f
    a = np.linalg.norm(a_t)
    T = a * m0
    if T > thrust:
        print('Warning: insufficient thrust')
    throttle = (T - thrust * min_throttle) / (thrust * (1. - min_throttle))
    vessel.control.throttle = throttle
    unit_t = a_t / a
    ap.target_direction = tuple(unit_t.flatten())
    if tgo < 0.5 or r_go < 2 or height() < 2 or speed() < 0.:
        print("ADGP disengaged")
        vessel.control.throttle = 0.
        break
while True:
    if vessel.situation == space_center.VesselSituation.landed \
            or vessel.situation == space_center.VesselSituation.splashed:
        print("landed!")
        break
ap.disengage()
vessel.control.rcs = False
print("Landing precision: {:.2f}m".format(np.linalg.norm(r_0 - r_f_ori)))
