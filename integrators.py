### time integrators #####
#####################################
def bEuler( y_old, dt, t, fun):
    ### y^n+1 = y^n + dt f(y^n, t^n)
    return y_old + dt*fun(y_old, t)
######################################
def fEuler( y_old, dt, t, fun):
    ### y^n+1 = y^n + dt f(y^n+1, t^n+1)
    n_iter = 0
    tol    = 1e-2
    value_new = y_old + dt * fun(y_old, t)
    value_old = y_old
    
    while ( abs(value_new - value_old)/abs(value_new) >= tol):
     value_old = value_new
     value_new  = y_old + dt * fun(value_old, t)    
     n_iter = n_iter + 1
     
    return value_new
######################################
def CN( y_old, dt, t, t_old, fun):
    ### y^n+1 = y^n + dt/2 * [f(y^n+1, t^n+1) + f(y^n, t^n)]
    n_iter = 0
    tol    = 1e-4
    value_new = y_old + dt * fun(y_old, t)
    value_old = y_old
    
    while ( abs(value_new - value_old)/abs(value_new) >= tol):
     value_old = value_new
     value_new  = y_old + (dt/2.) * ( fun(value_old, t)  + fun(y_old, t_old))  
     n_iter = n_iter + 1
     
    return value_new
def rk4 ( y_old, dt, t, fun):
    # RK-4 method
     k1 = dt * (fun( y_old , t))
     k2 = dt * (fun( y_old+k1/2, t+dt/2.))
     k3 = dt * (fun( y_old+k2/2, t+dt/2.))
     k4 = dt * (fun( y_old+k3  , t+dt   ))
     k = (k1+2*k2+2*k3+k4)/6.
     return  y_old + k    