import copy
import math
import numpy as np

from mxs import vehicle, calc_state
from calc_trim import get_trim_condition

# \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]

def linearize(vehicle, operating_point):
    # Want a model of the form:
    #  x_dot = a*x + b*u
    # Need to calculate a and b...
    a = np.zeros((5,5))
    op_derivative = np.array(vehicle.get_derivative(*operating_point))
    #                                      u, w,  q, theta, z
    for aj_index,state_index in enumerate([3, 5, 11, None, 2]):
        (state,input) = copy.deepcopy(operating_point)
        if state_index is None:
            # Theta is represented by a quaternion so needs special handling
            # Calculate the equivalent theta
            theta = math.asin(2*(state[9]*state[7]))
            if theta != 0:
                delta_xi = theta * 0.01
            else:
                delta_xi = -0.001
            state[9] = math.cos((theta+delta_xi)/2)
            state[7] = math.sin((theta+delta_xi)/2)
        else:
            if state[state_index] != 0:
                delta_xi = state[state_index] * 0.01
            else:
                delta_xi = -0.001
            state[state_index] += delta_xi
        
        derivative = vehicle.get_derivative(state,input)
        double_derivative = (np.array(derivative) - op_derivative) / delta_xi

        for ai_index, deriv_index in enumerate([3,5,11,None,2]):
            if deriv_index is None:
                # Compute the rate of change of theta from the rate of change of quaternion
                # https://ahrs.readthedocs.io/en/latest/filters/angular.html
                # qx = qz = 0
                # wx = wz = 0
                # dqw = 0.5 * (-wy*qy)
                # dqy = dqj = 0.5 * (wy*qw - wz*qx + wx*qz)
                #           = 0.5 * wy*qw
                
                # dqw
                wy = 2.0 * double_derivative[9] / -state[7]
                # dqj
                wy2 = 2.0 * double_derivative[7] / state[9]
                delta_wy = wy - wy2
                assert delta_wy < 0.00001
                a[ai_index,aj_index] = wy
            else:
                a[ai_index,aj_index] = double_derivative[deriv_index]
    
    return a

def print_eqn(a,b):
    print("\u250c   \u2510   \u250c                                              \u2510 \u250c   \u2510")
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('u\u0307',*a[0,:],'u'))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('w\u0307',*a[1,:],'w'))
    print("\u2502 {} \u2502 = \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('q\u0307',*a[2,:],'q'))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('\u03b8\u0307',*a[3,:],'\u03b8'))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('z\u0307',*a[4,:],'z'))
    print("\u2514   \u2518   \u2514                                              \u2518 \u2514   \u2518")

if __name__ == "__main__":
    np.set_printoptions(precision=3)
    
    airspeeds = np.linspace(12.5,25,100)
    eigenvalues = np.zeros((len(airspeeds),4), dtype='complex128')
    
    for i,airspeed in enumerate(airspeeds):
        (alpha,elevator,throttle) = get_trim_condition(airspeed)
        state = calc_state(math.radians(alpha),airspeed)
        a = linearize(vehicle, (state,[0,math.radians(elevator),throttle]))
        print(f"----- {airspeed:5.1f} m/s -----")
        print_eqn(a,a)
        vals, vecs = np.linalg.eig(a)
        eigenvalues[i,:] = vals[1:]
        for val,vec in zip(vals,vecs):
            print(f"{val}: {vec}")

    import matplotlib.pyplot as plt
    
    for eigval_index in range(4):
        plt.scatter(
            np.real(eigenvalues[:,eigval_index]),
            np.imag(eigenvalues[:,eigval_index]),
            c=airspeeds,
        )
    
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.minorticks_on()
    plt.grid(True,'both')
    plt.colorbar(label="Airspeed (m/s)")
    plt.show()
