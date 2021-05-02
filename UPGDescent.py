import time
import krpc
import json
import numpy as np
from scipy import optimize

conn = krpc.connect()
print("krpc connected")
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)
r0 = body.equatorial_radius
g0 = body.surface_gravity
ge = 9.80665  # standard gravity which is used to calculate exhaust velocity
# multipliers to normalize variables
time_multiplier = np.sqrt(r0 / g0)
position_multiplier = r0
velocity_multiplier = np.sqrt(r0 * g0)

# soft landing prob with 3-arc thrust

# load params
with open("target.json") as f:
    params = json.load(f)
lon = params['lon']
lat = params['lat']
v = params['vf']
min_throttle = params['min_throttle']
print(min_throttle)
# target vector
r_f = np.asarray(body.position_at_altitude(lat, lon,
                                           body.surface_height(lat, lon) + 150, body_frame)).reshape(-1, 1)

temp1 = space_center.ReferenceFrame.create_relative(body_frame,
                                                    rotation=(0., np.sin(-lon / 2. * np.pi / 180), 0.,
                                                              np.cos(-lon / 2. * np.pi / 180)))
temp2 = space_center.ReferenceFrame.create_relative(temp1,
                                                    rotation=(0., 0., np.sin(lat / 2. * np.pi / 180),
                                                              np.cos(lat / 2. * np.pi / 180)))
height = body.surface_height(lat, lon) + body.equatorial_radius
target_frame = space_center.ReferenceFrame.create_relative(temp2, position=(height, 0., 0.))
del temp1, temp2
'''line = conn.drawing.add_line((10,0,0),(0,0,0), target_frame)
line.color = (1,0,0)
line.thickness = 2'''
target_v = np.array([-v, 0., 0.])
v_f = np.asarray(space_center.transform_velocity((0., 0., 0.), target_v, target_frame, body_frame)).reshape(-1, 1)
t_1 = 10.

# normalize to dimensionless
r_f = r_f / position_multiplier
v_f = v_f / velocity_multiplier


def mass_t(t, t0, isp, m0, thrust):
    """
    calculate the expected mass at dimensionless time t
    :param t: dimensionless time
    :param isp: engine's specific impulse
    :param m0: initial mass
    :param thrust: engine thrust
    :return: mass at t
    """
    return m0 - (t - t0) * time_multiplier * thrust / (isp * ge)


def transfer_mat(t, t0):
    return np.block([[np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)],
                     [-np.sin(t - t0) * np.eye(3), np.cos(t - t0) * np.eye(3)]])


def gamma_mat(t, t0):
    return np.block([[np.sin(t - t0) * np.eye(3), -np.cos(t - t0) * np.eye(3)],
                     [np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)]])


def lambda_t(t, t0, pv0, pr0):
    """
    calculates the costate vectors at t(dimensionless)
    :param t: desired t
    :param t0: initial t
    :param pr0: initial Pr, 3x1
    :param pv0: initial Pv, 3x1
    :return: Pv,Pr at t
    """
    lambda_0 = np.vstack((pv0, -pr0))
    temp = transfer_mat(t, t0) @ lambda_0
    return temp[0:3, :], temp[3:6, :]


def x_t(x0, t, t0, t1, t2, pv0, pr0, thrust, min_t, isp, m0):

    def I(t, t0, pv0, pr0, thrust, isp, m0):
        i_c = thrust_integral('c', t, t0, pv0, pr0, thrust, isp, m0)
        i_s = thrust_integral('s', t, t0, pv0, pr0, thrust, isp, m0)
        return np.vstack((i_c, i_s))

    if t1 > t2:
        print("warning: t1>t2, t1="+str(t1)+" t2="+str(t2))
        t2 = t1
    if t0 < t1:
        # x at t1
        i_t1 = I(t1, t0, pv0, pr0, thrust, isp, m0)
        x_temp = transfer_mat(t1, t0) @ x0 + gamma_mat(t1, t0) @ i_t1
        pv_t1, pr_t1 = lambda_t(t1, t0, pv0, pr0)
        m_t1 = mass_t(t1, t0, isp, m0, thrust)
        # x at t2
        i_t2 = I(t2, t1, pv_t1, pr_t1, thrust * min_t, isp, m_t1)
        x_temp = transfer_mat(t2, t1) @ x_temp + gamma_mat(t2, t1) @ i_t2
        pv_t2, pr_t2 = lambda_t(t2, t1, pv_t1, pr_t1)
        m_t2 = mass_t(t2, t1, isp, m_t1, thrust * min_t)
        # x at t
        i_t = I(t, t2, pv_t2, pr_t2, thrust, isp, m_t2)
        return transfer_mat(t, t2) @ x_temp + gamma_mat(t, t2) @ i_t
    elif t0 > t1:
        if t0 < t2:
            # x at t2
            i_t2 = I(t2, t0, pv0, pr0, thrust * min_t, isp, m0)
            pv_t2, pr_t2 = lambda_t(t2, t0, pv0, pr0)
            m_t2 = mass_t(t2, t0, isp, m0, thrust * min_t)
            x_temp = transfer_mat(t2, t0) @ x0 + gamma_mat(t2, t0) @ i_t2
            # x at t
            i_t = I(t, t2, pv_t2, pr_t2, thrust, isp, m_t2)
            return transfer_mat(t, t2) @ x_temp + gamma_mat(t, t2) @ i_t
        elif t0 > t2:
            i_t = I(t, t0, pv0, pr0, thrust, isp, m0)
            return transfer_mat(t, t0) @ x0 + gamma_mat(t, t0) @ i_t


def thrust_integral(type_name, t, t0, pv0, pr0, thrust, isp, m0):
    """
    calculate thrust integral with milne's rule
    :param type_name: determins this is the cos integral or sin integral, should be 'c' or 's'
    :param t: end time
    :param t0: initial time
    :param pv0: initial pv
    :param pr0: initial pr
    :param thrust: engine thrust
    :param isp: engine isp
    :param m0: initial mass
    :return: the estimated thrust integral
    """

    def i(t):
        pv, _ = lambda_t(t, t0, pv0, pr0)
        pv = pv / np.linalg.norm(pv)
        temp = pv * thrust / (mass_t(t, t0, isp, m0, thrust) * g0)
        if type_name == 'c':
            return temp * np.cos(t - t0)
        elif type_name == 's':
            return temp * np.sin(t - t0)
        else:
            raise Exception()

    delta = (t - t0) / 4.
    return ((t - t0) / 90.) * (7 * i(t0) + 32 * i(t0 + delta) + 12 * i(t0 + 2 * delta) + 32 * i(t0 + 3 * delta) +
                               7 * i(t0 + 4 * delta))


# get initial values
specific_impulse = vessel.vacuum_specific_impulse
max_thrust = conn.add_stream(getattr, vessel, 'max_thrust')
mass = conn.add_stream(getattr, vessel, 'mass')
ut = conn.add_stream(getattr, space_center, 'ut')
velocity = conn.add_stream(vessel.velocity, body_frame)
position = conn.add_stream(vessel.position, body_frame)
thrust = max_thrust()
m0 = mass()
v_0 = np.asarray(velocity()).reshape(-1, 1) / velocity_multiplier
r_0 = np.asarray(position()).reshape(-1, 1) / position_multiplier
t_0 = ut() / time_multiplier
t_1 = t_0 + t_1 / time_multiplier

# initial guess

pv = np.asarray(-v_0 / np.linalg.norm(v_0)).reshape(-1, 1)
pr = np.zeros((3, 1))
tf = np.array([(ut() + 350.) / time_multiplier]).reshape(1, 1)
x = np.vstack((r_0, v_0))

z0 = np.hstack((pv.T, pr.T, tf)).T


def fun(z, x_0, t_0, t_1, t_2, thrust, min_t, isp, m0 ):
    pv = z[0:3].reshape(-1, 1)
    pr = z[3:6].reshape(-1, 1)
    tf = z[6]
    # final position constraint
    x_f = x_t(x_0, tf, t_0, t_1, t_2, pv, pr, thrust, min_t, isp, m0)
    r_tf = x_f[0:3, :]
    s1 = r_tf.T @ r_tf - r_f.T @ r_f
    # final velocity constraint
    v_tf = x_f[3:6, :]
    s2 = v_tf - v_f  # s2 includes 3 equations
    # transversality condition
    pvf, prf = lambda_t(tf, t_0, pv, pr)

    s3 = prf.T @ v_tf - pvf.T @ r_tf + \
         np.linalg.norm(pvf) * thrust / (mass_t(tf, t_0, specific_impulse, m0, thrust) * g0) - 1
    s4 = r_tf[2, 0] * prf[0, 0] - r_tf[0, 0] * prf[2, 0]
    s5 = r_tf[2, 0] * prf[1, 0] - r_tf[1, 0] * prf[2, 0]
    return [s1[0, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0], s4, s5]


def fun_2(t2, z0, x_0, t_0, t_1, thrust, min_t, isp, m0):
    sol = optimize.root(fun, z0, args=(x_0, t_0, t_1, t2, thrust, min_t, isp, m0), method='lm', jac=False)
    z = sol.x.reshape(-1, 1)
    tf = z[6, 0]
    m_t1 = mass_t(t_1, t_0, isp, m0, thrust)
    m_t2 = mass_t(t2, t_1, isp, m_t1, thrust * min_t)
    m_tf = mass_t(tf, t2, isp, m_t2, thrust)
    return m0-m_tf

print("evaluating t2")
result = optimize.minimize_scalar(fun_2, bounds=(t_1, tf[0, 0]),
                               args=(z0, x, t_0, t_1, thrust, min_throttle, specific_impulse, m0), method='Bounded')
t_2 = result.x
# check if t2 smaller than t1
if t_2 < t_1:
    t_2 = t_1
# convert to time interval
t2 = t_2 - t_0
print("optimal t2:"+str(t_2*time_multiplier))
ap.reference_frame = body_frame
ap.engage()

#run once
print("initiating UPG")
sol = optimize.root(fun, z0, args=(x, t_0, t_1, t_2, thrust, min_throttle, specific_impulse, m0),
                    method='lm', jac=False)
z = sol.x.reshape(-1, 1)
z0 = z
pv = z[0:3, :]
unit_t = pv / np.linalg.norm(pv)
ap.target_direction = tuple(unit_t.reshape(1, -1)[0])
print("pointing to correct position")
ap.wait()
v_0 = np.asarray(velocity()).reshape(-1, 1) / velocity_multiplier
r_0 = np.asarray(position()).reshape(-1, 1) / position_multiplier
x = np.vstack((r_0, v_0))
m0 = mass()
thrust = max_thrust()
t_0 = ut() / time_multiplier
t_1 = t_0 + 10 / time_multiplier
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
    tgo = (tf - t_0) * time_multiplier
    tgo_1 = (t_1 - t_0) * time_multiplier
    tgo_2 = (t_2 - t_0) * time_multiplier
    # distance = np.linalg.norm((r_0 - r_f) * r0)

    print("tgo:" + str(tgo))
    if tgo_1 > 0:
        print("tgo1:" + str(tgo_1))
    if tgo_2 > 0:
        print("tgo2:" + str(tgo_2))
    # print("distance:" + str(distance))
    unit_t = pv / np.linalg.norm(pv)
    ap.target_direction = tuple(unit_t.reshape(1, -1)[0])
    ap.wait()
    if t_1 < t_0 < t_2:
        vessel.control.throttle = 0.00001
    else:
        vessel.control.throttle = 1
    if tgo < 2. and flight.surface_altitude < 2000:  # TODO: disengage criterion not robust
        print("UPG disengaged")
        break
    time.sleep(0.3)
    # get new values for next loop
    v_0 = np.asarray(velocity()).reshape(-1, 1) / velocity_multiplier
    r_0 = np.asarray(position()).reshape(-1, 1) / position_multiplier
    x = np.vstack((r_0, v_0))
    m0 = mass()
    thrust = max_thrust()
    t_0 = ut() / time_multiplier

vessel.control.throttle = 0.
ap.disengage()
vessel.control.sas = True
vessel.control.sas_mode = space_center.SASMode.retrograde

