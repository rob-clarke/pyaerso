import math
import scipy.optimize
from mxs import calc_state

def get_trim_error(vehicle, state, input):
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



def get_trim_condition(vehicle, airspeed, debug=False):
    initial_state = calc_state(math.radians(10),airspeed)
    initial_input = [0,0,0,0,0,0]

    def trim_trial(params,airspeed):
        [alpha,elevator,throttle] = params
        state = calc_state(alpha,airspeed)
        sum_error, _ = get_trim_error(vehicle, state, [0, elevator, throttle, 0, 0, 0])
        return sum_error

    if debug:
        print(f"---- {airspeed:4.1f} m/s ----")
        print(f"{initial_state=}")
        print(f"initial_error={get_trim_error(vehicle, initial_state, initial_input)}")

    x0 = [0.02, math.radians(-5), 0.5]
    bounds = [
        (-math.radians(10), math.radians(25)),
        (math.radians(-30), math.radians(30)),
        (0, 1)
    ]
    res = scipy.optimize.minimize(
        trim_trial,
        x0,
        args=(airspeed),
        bounds=bounds,
        options={'ftol': 1e-18, 'eps': 1e-9, 'gtol': 1e-9}
    )
    
    [alpha,elevator,throttle] = res.x
    alpha = math.degrees(alpha)
    elevator = math.degrees(elevator)
    
    if debug:
        print(res.message)
        print(f"{alpha=} {elevator=} {throttle=}")
    
    trim_state = calc_state(res.x[0],airspeed)
    
    if debug:
        print( get_trim_error(vehicle, trim_state, [0,*res.x[1:],0,0,0]) )
    
    return (alpha,elevator,throttle)


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        airspeeds = [ float(v) for v in sys.argv[1].split(',') ]
    else:
        airspeeds = [15, 17.5, 20, 22.5, 25]
    
    from mxs import vehicle
    
    trim_states = {}
    
    for airspeed in airspeeds:
        trim_states[airspeed] = get_trim_condition(vehicle, airspeed, debug=True)
    
    print("----- Results -----")
    print("{ airspeed: (alpha(deg), elevator(deg), throttle) ... }")
    import pprint
    
    pprint.pprint(trim_states)
