import time
import krpc
import json
import numpy as np
from scipy import optimize
import UPG

conn = krpc.connect()
print("krpc connected")
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)
# multipliers to normalize variables
upg = UPG.UPG(conn)
# soft landing prob with 3-arc thrust

# load params
with open("target.json") as f:
    params = json.load(f)
lon = params['lon']
lat = params['lat']
v = params['vf']
min_throttle = params['min_throttle']
# target vector
r_f = np.asarray(body.position_at_altitude(lat, lon,
                                           body.surface_height(lat, lon) + 250, body_frame)).reshape(-1, 1)
# get final velocity in body frame
temp1 = space_center.ReferenceFrame.create_relative(body_frame,
                                                    rotation=(0., np.sin(-lon / 2. * np.pi / 180), 0.,
                                                              np.cos(-lon / 2. * np.pi / 180)))
temp2 = space_center.ReferenceFrame.create_relative(temp1,
                                                    rotation=(0., 0., np.sin(lat / 2. * np.pi / 180),
                                                              np.cos(lat / 2. * np.pi / 180)))
height = body.surface_height(lat, lon) + body.equatorial_radius
target_frame = space_center.ReferenceFrame.create_relative(temp2, position=(height, 0., 0.))
del temp1, temp2
target_v = np.array([-v, 0., 0.])
v_f = np.asarray(space_center.transform_velocity((0., 0., 0.), target_v, target_frame, body_frame)).reshape(-1, 1)
t1 = 15.

# normalize to dimensionless
r_f = r_f / upg.position_multiplier
v_f = v_f / upg.velocity_multiplier

# create data streams
vacuum_specific_impulse = conn.add_stream(getattr, vessel, 'vacuum_specific_impulse')
max_thrust = conn.add_stream(getattr, vessel, 'max_thrust')
current_thrust = conn.add_stream(getattr, vessel, 'thrust')
mass = conn.add_stream(getattr, vessel, 'mass')
ut = conn.add_stream(getattr, space_center, 'ut')
height = conn.add_stream(getattr, flight, 'surface_altitude')

# get initial values
velocity = conn.add_stream(vessel.velocity, body_frame)
position = conn.add_stream(vessel.position, body_frame)
thrust = max_thrust()
m0 = mass()
specific_impulse = vacuum_specific_impulse()
v_0 = np.asarray(velocity()).reshape(-1, 1) / upg.velocity_multiplier
r_0 = np.asarray(position()).reshape(-1, 1) / upg.position_multiplier
t_0 = ut() / upg.time_multiplier
t_1 = t_0 + t1 / upg.time_multiplier

# initial guess

pv = np.asarray(-v_0 / np.linalg.norm(v_0)).reshape(-1, 1)
pr = np.zeros((3, 1))
tf = np.array([(ut() + 350.) / upg.time_multiplier]).reshape(1, 1)
x = np.vstack((r_0, v_0))

z0 = np.hstack((pv.T, pr.T, tf)).T


def fun(z, x_0, t_0, t_1, t_2, thrust, min_t, isp, m0):
    pv = z[0:3].reshape(-1, 1)
    pr = z[3:6].reshape(-1, 1)
    tf = z[6]
    # final position constraint
    x_f = upg.x_t(x_0, tf, t_0, t_1, t_2, pv, pr, thrust, min_t, isp, m0)
    r_tf = x_f[0:3, :]
    s1 = r_tf.T @ r_tf - r_f.T @ r_f
    # final velocity constraint
    v_tf = x_f[3:6, :]
    s2 = v_tf - v_f  # s2 includes 3 equations
    # transversality condition
    pvf, prf = upg.lambda_t(tf, t_0, pv, pr)

    s3 = prf.T @ v_tf - pvf.T @ r_tf + \
         np.linalg.norm(pvf) * thrust / (upg.mass_t(tf, t_0, specific_impulse, m0, thrust) * upg.g0) - 1
    s4 = r_tf[2, 0] * prf[0, 0] - r_tf[0, 0] * prf[2, 0]
    s5 = r_tf[2, 0] * prf[1, 0] - r_tf[1, 0] * prf[2, 0]
    return [s1[0, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0], s4, s5]


def fun_2(t2, z0, x_0, t_0, t_1, thrust, min_t, isp, m0):
    sol = optimize.root(fun, z0, args=(x_0, t_0, t_1, t2, thrust, min_t, isp, m0), method='lm', jac=False)
    z = sol.x.reshape(-1, 1)
    tf = z[6, 0]
    m_t1 = upg.mass_t(t_1, t_0, isp, m0, thrust)
    m_t2 = upg.mass_t(t2, t_1, isp, m_t1, thrust * min_t)
    m_tf = upg.mass_t(tf, t2, isp, m_t2, thrust)
    return m0 - m_tf


print("evaluating t2")
result = optimize.minimize_scalar(fun_2, bounds=(t_1, tf[0, 0]),
                                  args=(z0, x, t_0, t_1, thrust, min_throttle, specific_impulse, m0), method='Bounded')
t_2 = result.x
# convert to time interval
t2 = t_2 - t_0
print("optimal t2:" + str(t2 * upg.time_multiplier))
ap.reference_frame = body_frame
ap.engage()

# run once
print("initiating UPG")
sol = optimize.root(fun, z0, args=(x, t_0, t_1, t_2, thrust, min_throttle, specific_impulse, m0),
                    method='lm', jac=False)
z = sol.x.reshape(-1, 1)
z0 = z
pv = z[0:3, :]
pr = z[3:6, :]
tf = z[6, 0]
unit_t = pv / np.linalg.norm(pv)
ap.target_direction = tuple(unit_t.reshape(1, -1)[0])
print("pointing to correct direction")
ap.wait()



v_0 = np.asarray(velocity()).reshape(-1, 1) / upg.velocity_multiplier
r_0 = np.asarray(position()).reshape(-1, 1) / upg.position_multiplier
x = np.vstack((r_0, v_0))
m0 = mass()
thrust = max_thrust()
t_0 = ut() / upg.time_multiplier
t_1 = t_0 + t1 / upg.time_multiplier
t_2 = t_0 + t2

vessel.control.throttle = 1
while True:
    sol = optimize.root(fun, z0, args=(x, t_0, t_1, t_2, thrust, min_throttle, specific_impulse, m0),
                        method='lm', jac=False)
    z = sol.x.reshape(-1, 1)
    z0 = z
    pv = z[0:3, :]
    pr = z[3:6, :]
    tf = z[6, 0]
    tgo = (tf - t_0) * upg.time_multiplier
    tgo_1 = (t_1 - t_0) * upg.time_multiplier
    tgo_2 = (t_2 - t_0) * upg.time_multiplier
    distance = np.linalg.norm((r_0 - r_f) * upg.position_multiplier)

    print("tgo: {:.2f}".format(tgo))
    if tgo_1 > 0:
        print("tgo1: {:.2f}".format(tgo_1))
    if tgo_2 > 0:
        print("tgo2: {:.2f}".format(tgo_2))
    print("distance: {:.2f}".format(distance))
    unit_t = pv / np.linalg.norm(pv)
    ap.target_direction = tuple(unit_t.flatten())
    if t_1 < t_0 < t_2:
        vessel.control.throttle = 0.00001
    else:
        vessel.control.throttle = 1
    if tgo < 15. and flight.surface_altitude < 2000:
        print("UPG disengaged")
        break
    time.sleep(0.2)

    # correct target to the correct height
    x_f = upg.x_t(x, tf, t_0, t_1, t_2, pv, pr, thrust, min_throttle, specific_impulse, m0)
    r_f_t = tuple(upg.position_multiplier * x_f[0:3, :].flatten())
    lat = body.latitude_at_position(r_f_t, body_frame)
    lon = body.longitude_at_position(r_f_t, body_frame)
    r_f = np.asarray(body.position_at_altitude(
        lat, lon, body.surface_height(lat, lon) + 300, body_frame)).reshape(-1, 1) / upg.position_multiplier

    # get new values for next loop
    v_0 = np.asarray(velocity()).reshape(-1, 1) / upg.velocity_multiplier
    r_0 = np.asarray(position()).reshape(-1, 1) / upg.position_multiplier
    x = np.vstack((r_0, v_0))
    m0 = mass()
    thrust = max_thrust()
    t_0 = ut() / upg.time_multiplier

# vessel.control.throttle = 0.

# change r_f to ground
r_f = np.asarray(body.position_at_altitude(
    lat, lon, body.surface_height(lat, lon) + 10, body_frame)).reshape(-1, 1)
# ADGP
a_f = np.zeros((3, 1))
v_f = np.zeros((3, 1))
# calculate tgo
sin_gamma = flight.vertical_speed / flight.speed


def fun_at(a_t, v_0, s_gamma, h0, g0):
    return a_t ** 2 + ((v_0 ** 2 * s_gamma) / (2 * h0)) * a_t - \
           (v_0 ** 2 * g0 * (1 + s_gamma ** 2) / (4 * h0)) - g0 ** 2


V_0 = flight.speed
h_0 = height()
a_t = optimize.root_scalar(fun_at, args=(V_0, sin_gamma, h_0, upg.g0),
                           method='brentq', bracket=(0, thrust / m0)).root
t_f = 1.3 * (V_0 / 2) * (((1 + sin_gamma) / (a_t + upg.g0)) + ((1 - sin_gamma) / (a_t - upg.g0)))
print('tgo:' + str(t_f))
t_f = t_f + ut()
print("Terminal guidance engaged")
while True:
    time.sleep(0.05)

    v_0 = np.asarray(velocity()).reshape(-1, 1)
    r_0 = np.asarray(position()).reshape(-1, 1)
    m0 = mass()
    thrust = max_thrust()
    curr_t = current_thrust()
    t_0 = ut()
    t_go = t_f - t_0
    print("t_go: {:.2f}".format(t_go))
    r_go = np.linalg.norm(r_f - r_0)
    print("target distance: {:.2f}".format(r_go))
    a_t = (-6 / t_go) * (v_f - v_0) + (12 / t_go ** 2) * (r_f - r_0 - v_0 * t_go) + a_f
    a = np.linalg.norm(a_t)
    T = a * m0
    if T > thrust:
        print('Warning: insufficient thrust')
    throttle = (T - thrust * min_throttle) / (thrust * (1. - min_throttle))
    vessel.control.throttle = throttle
    unit_t = a_t / a
    ap.target_direction = tuple(unit_t.flatten())
    if tgo < 0.5 or r_go < 5 or height() < 5:
        print("ADGP disengaged")
        break
while True:
    if vessel.situation == space_center.VesselSituation.landed:
        break
ap.disengage()
vessel.control.throttle = 0.
vessel.control.rcs = False
print("Landing pricision: {:.2f}m".format(distance))

