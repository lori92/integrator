import numpy as np
import math as mt
import scipy as sc
import matplotlib 
from matplotlib import pyplot as plt
import integrators
from integrators import bEuler as bEuler
from integrators import fEuler as fEuler
from integrators import CN as CN
from integrators import rk4 as rk4



######## integration of ODE dy/dt = f(y,t) ###########
def fun(y,t):
    return  2*t

######################################
x0 = np.array([1])
t0  = 0
f0 = 0
######################################
t_end = 8
N     = 401
dt    = (t_end - t0)/(N-1)
time = np.linspace(t0, t_end, N)
######################################
t_old = t0

y_old_bEuler = f0
y_old_CN     = f0
y_old_rk4    = f0

y_bEuler = np.zeros([N, np.max(np.size(f0))])
y_CN     = np.zeros([N, np.max(np.size(f0))])
y_rk4    = np.zeros([N, np.max(np.size(f0))])

y_bEuler[0][::] = f0
y_CN    [0][::] = f0
y_rk4   [0][::] = f0

#########################################
for i in range(1,N):
  y_new           = bEuler(y_old_bEuler, dt, time[i], fun )
  y_bEuler[i][::] = y_new
  y_old_bEuler    = y_new

  #y_new       = CN(y_old_CN, dt, time[i], time[i-1], fun )
  #y_CN[i][::] = y_new
  #y_old_CN    = y_new

  #y_new        = rk4(y_old_rk4, dt, time[i], fun )
  #y_rk4[i][::] = y_new
  #y_old_rk4    = y_new

######################################
fig, ax = plt.subplots()
ax.plot(time, y_bEuler, 'r', label="backwar Euler")
#ax.plot(time, y_CN,     "b", label="CN")
#ax.plot(time, y_rk4,    "k", label="rk4")
ax.plot(time,   time**2.,    "--", label="exact")

ax.set_xlabel("t")
ax.set_title("y(t)")
plt.show()







