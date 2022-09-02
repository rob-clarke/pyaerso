import math
import scipy.optimize
from mxs import vehicle, calc_state

def get_trim_error(state,input):
    vehicle.statevector = state
    derivative = vehicle.get_derivative(state,input)
    
    z_dot = derivative[2]
    
    u_dot = derivative[3]
    w_dot = derivative[5]
    
    orientation_dot = derivative[6:10]
    
    q_dot = derivative[11]
    
    # print(f"{z_dot=} {u_dot=} {w_dot=} {orientation_dot=} {q_dot=}")
    
    sum_error = z_dot**2 + u_dot**2 + w_dot**2 + q_dot**2 + \
        orientation_dot[0]**2 + orientation_dot[1]**2 + orientation_dot[2]**2 + orientation_dot[3]**2
    
    return sum_error, [z_dot,u_dot,w_dot,orientation_dot,q_dot]

def trim_trial(params,airspeed):
    [alpha,elevator,throttle] = params
    state = calc_state(alpha,airspeed)
    sum_error, _ = get_trim_error(state, [0, elevator, throttle])
    return sum_error

if __name__ == "__main__":
    airspeed = 22
    initial_state = calc_state(math.radians(10),airspeed)
    initial_input = [0,0,0]
    
    print(f"{initial_state=}")
    print(get_trim_error(initial_state, initial_input))

    x0 = [0.02, math.radians(-5), 0.5]
    bounds = [
        (0, math.radians(25)),
        (math.radians(-30), math.radians(30)),
        (0, 1)
    ]
    res = scipy.optimize.minimize(trim_trial,x0,args=(airspeed),bounds=bounds,options={'ftol': 1e-18, 'eps': 1e-9, 'gtol': 1e-9})
    print(res)
    [alpha,elevator,throttle] = res.x
    print(f"{math.degrees(alpha)=} {math.degrees(elevator)=} {throttle=}")
    
    trim_state = calc_state(res.x[0],airspeed)
    print( get_trim_error(trim_state, [0,*res.x[1:]]) )
