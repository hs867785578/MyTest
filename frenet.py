import numpy as np
from math import *
import csv
import matplotlib.pyplot as plt

def MyFmod(_X,  _Y):
    return _X - (int)(_X / _Y) * _Y

def NormalizeAngle(angle):	
    if angle> 2.0*pi or angle < -2.0*pi:
        a = MyFmod(angle + pi, 2.0 * pi)
        if a < 0.0 :
            a += 2.0 * pi
        return a - 2.0*pi
    else:
        return angle

def cartesian_to_frenet1D(rs, rx, ry, rtheta, x, y):
    s_condition = np.zeros(1)
    d_condition = np.zeros(1)
    
    dx = x - rx
    dy = y - ry
    
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    d_condition[0] = copysign(sqrt(dx * dx + dy * dy), cross_rd_nd)    
    
    s_condition[0] = rs
    
    return s_condition, d_condition

def cartesian_to_frenet2D(rs, rx, ry, rtheta, rkappa, x, y, v, theta):
    s_condition = np.zeros(2)
    d_condition = np.zeros(2)
    
    dx = x - rx
    dy = y - ry
    
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    d_condition[0] = copysign(sqrt(dx * dx + dy * dy), cross_rd_nd)
    
    delta_theta = theta - rtheta
    tan_delta_theta = tan(delta_theta)
    cos_delta_theta = cos(delta_theta)
    
    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    d_condition[1] = one_minus_kappa_r_d * tan_delta_theta
    
    
    s_condition[0] = rs
    s_condition[1] = v * cos_delta_theta / one_minus_kappa_r_d

    return s_condition, d_condition

def cartesian_to_frenet3D(rs, rx, ry, rtheta, rkappa, rdkappa, x, y, v, a, theta, kappa):
    s_condition = np.zeros(3)
    d_condition = np.zeros(3)
    
    dx = x - rx
    dy = y - ry
    
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    cross_rd_nd = cos_theta_r * dy - sin_theta_r * dx
    d_condition[0] = copysign(sqrt(dx * dx + dy * dy), cross_rd_nd)
    
    delta_theta = theta - rtheta
    tan_delta_theta = tan(delta_theta)
    cos_delta_theta = cos(delta_theta)
    
    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    d_condition[1] = one_minus_kappa_r_d * tan_delta_theta
    
    kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]
    
    d_condition[2] = (-kappa_r_d_prime * tan_delta_theta + 
      one_minus_kappa_r_d / cos_delta_theta / cos_delta_theta *
          (kappa * one_minus_kappa_r_d / cos_delta_theta - rkappa))
    
    s_condition[0] = rs
    s_condition[1] = v * cos_delta_theta / one_minus_kappa_r_d
    
    delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * kappa - rkappa
    s_condition[2] = ((a * cos_delta_theta -
                       s_condition[1] * s_condition[1] *
                       (d_condition[1] * delta_theta_prime - kappa_r_d_prime)) /
                          one_minus_kappa_r_d)
    return s_condition, d_condition


def frenet_to_cartesian1D(rs, rx, ry, rtheta, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]    

    return x, y

def frenet_to_cartesian2D(rs, rx, ry, rtheta, rkappa, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]

    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
    delta_theta = atan2(d_condition[1], one_minus_kappa_r_d)
    cos_delta_theta = cos(delta_theta)
    
    theta = NormalizeAngle(delta_theta + rtheta)    
    
    d_dot = d_condition[1] * s_condition[1]
    
    v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d * s_condition[1] * s_condition[1] + d_dot * d_dot)   

    return x, y, v, theta 

def frenet_to_cartesian3D(rs, rx, ry, rtheta, rkappa, rdkappa, s_condition, d_condition):
    if fabs(rs - s_condition[0])>= 1.0e-6:
        print("The reference point s and s_condition[0] don't match")
        
    cos_theta_r = cos(rtheta)
    sin_theta_r = sin(rtheta)
    
    x = rx - sin_theta_r * d_condition[0]
    y = ry + cos_theta_r * d_condition[0]

    one_minus_kappa_r_d = 1 - rkappa * d_condition[0]
    tan_delta_theta = d_condition[1] / one_minus_kappa_r_d
    delta_theta = atan2(d_condition[1], one_minus_kappa_r_d)
    cos_delta_theta = cos(delta_theta)
    
    theta = NormalizeAngle(delta_theta + rtheta)
    kappa_r_d_prime = rdkappa * d_condition[0] + rkappa * d_condition[1]
        
    kappa = ((((d_condition[2] + kappa_r_d_prime * tan_delta_theta) *
                 cos_delta_theta * cos_delta_theta) /
                    (one_minus_kappa_r_d) +
                rkappa) *
               cos_delta_theta / (one_minus_kappa_r_d))
    
    
    d_dot = d_condition[1] * s_condition[1]
    
    v = sqrt(one_minus_kappa_r_d * one_minus_kappa_r_d * s_condition[1] * s_condition[1] + d_dot * d_dot)
    
    delta_theta_prime = one_minus_kappa_r_d / cos_delta_theta * (kappa) - rkappa     
    a = (s_condition[2] * one_minus_kappa_r_d / cos_delta_theta +
           s_condition[1] * s_condition[1] / cos_delta_theta *
               (d_condition[1] * delta_theta_prime - kappa_r_d_prime))
    return x, y, v, a, theta, kappa 


waypoints_np   = None
waypoints_file = "./racetrack_waypoints.txt"
with open(waypoints_file) as waypoints_file_handle:
    waypoints = list(csv.reader(waypoints_file_handle, 
        delimiter=',',
        quoting=csv.QUOTE_NONNUMERIC))
    waypoints_np = np.array(waypoints)

x = waypoints_np[:, 0]
y = waypoints_np[:, 1]

plt.plot(x, y, 'r')
plt.show()

import scipy.signal.spline as spline

ds = 0.1# [m] distance of each interpolated points
sp = spline.qspline2d(x, y)
s = np.arange(0, sp.s[-1], ds)


left_bound = []
right_bound = []
new_x = []

for i_s in s:
    rx, ry = sp.calc_position(i_s)
    rtheta = sp.calc_yaw(i_s)

    new_x.append(np.array([rx, ry]))

    l_s_condition = np.array([i_s])
    l_d_condition = np.array([1.5])
    l_x, l_y = frenet_to_cartesian1D(i_s, rx, ry, rtheta, l_s_condition, l_d_condition)
    left_bound.append(np.array([l_x, l_y]))
    
    r_s_condition = np.array([i_s])
    r_d_condition = np.array([-1.5])
    r_x, r_y = frenet_to_cartesian1D(i_s, rx, ry, rtheta, r_s_condition, r_d_condition)   
    right_bound.append(np.array([r_x, r_y]))

left_bound = np.array(left_bound)
right_bound = np.array(right_bound)
new_x = np.array(new_x)
  
plt.plot(new_x[:, 0], new_x[:, 1], 'y')
plt.plot(left_bound[:,0],left_bound[:,1], 'b')
plt.plot(right_bound[:,0],right_bound[:,1], 'b')


ego = np.array([-126.87, -526.30])
plt.plot(ego[0], ego[1], 'o')


def polyfit(coeffs, t): 
  return coeffs[0] + coeffs[1] * t + coeffs[2] * t * t + coeffs[3] * t * t * t

def calc_cubic_poly_curve_coeffs(x_s, y_s, y_s_1d, x_e, y_e, y_e_1d):
    A = np.array([
        [1, x_s, pow(x_s, 2.0), pow(x_s, 3.0)],
        [0, 1,   2 * x_s,       3 * pow(x_s, 2.0)],
        [1, x_e, pow(x_e, 2.0), pow(x_e, 3.0)],
        [0, 1,   2 * x_e,       3 * pow(x_e, 2.0)]
    ])

    A = A.astype(np.float)

    b = np.array([[y_s], [y_s_1d], [y_e], [y_e_1d]])

    A_inv = np.linalg.inv(A)

    coffes = np.dot(A_inv, b)

    return coffes
  
planning_horizon = 200 * ds

targets = []
for i in range(11):
    x0 = ego_s
    y0 = ego_d
    dx0 = 0

    x1 = ego_s + planning_horizon
    y1 = -1.5 + i * 0.3
    dx1 = 0

    coffes = calc_cubic_poly_curve_coeffs(x0, y0, dx0, x1, y1, dx1)

    targets.append(coffes)