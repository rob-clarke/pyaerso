import copy
import math
import numpy as np

from mxs import calc_state
from calc_trim import get_trim_condition

# \[position,velocity(body),attitude_quaternion(i,j,k,w),axis_rates(body)\]

def linearize(vehicle, operating_point, control_indices=[1,2]):
    # Want a model of the form:
    #  x_dot = a*x + b*u
    # Need to calculate a and b...
    a = np.zeros((5,5))
    b = np.zeros((5,len(control_indices)))
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
    
    for bj_index,input_index in enumerate(control_indices):
        (state,input) = copy.deepcopy(operating_point)
        
        if input[input_index] != 0:
            delta_ui = input[input_index] * 0.01
        else:
            delta_ui = -0.001
        input[input_index] += delta_ui
        
        derivative = vehicle.get_derivative(state,input)
        double_derivative = (np.array(derivative) - op_derivative) / delta_ui

        for bi_index, deriv_index in enumerate([3,5,11,None,2]):
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
                b[bi_index,bj_index] = wy
            else:
                b[bi_index,bj_index] = double_derivative[deriv_index]
    
    return a,b


def print_eqn(a,b):
    print("\u250c   \u2510   \u250c                                              \u2510 \u250c   \u2510   \u250c                   \u2510")
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} \u2502 \u250c   \u2510".format('u\u0307',*a[0,:],'u',*b[0,:]))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('w\u0307',*a[1,:],'w',*b[1,:],'\u03b7'))
    print("\u2502 {} \u2502 = \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502 + \u2502 {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502".format('q\u0307',*a[2,:],'q',*b[2,:],'\u03c4'))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} \u2502 \u2514   \u2518".format('\u03b8\u0307',*a[3,:],'\u03b8',*b[3,:]))
    print("\u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} \u2502 \u2502 {} \u2502   \u2502 {:8.2f} {:8.2f} \u2502".format('z\u0307',*a[4,:],'z',*b[4,:]))
    print("\u2514   \u2518   \u2514                                              \u2518 \u2514   \u2518   \u2514                   \u2518")

def pretty_eqn(a,b,control_symb=['\u03b7','\u03c4']):
    from prettyeqn import Equation, SymbolVector, Matrix, Symbol
    result = SymbolVector([
        'u\u0307',
        'w\u0307',
        'q\u0307',
        '\u03b8\u0307',
        'z\u0307',
    ])
    a = Matrix(a)
    state = SymbolVector([
        'u',
        'w',
        'q',
        '\u03b8',
        'z',
    ])
    b = Matrix(b)
    input = SymbolVector(control_symb)
    equation = Equation([result,Symbol('='),a,state,Symbol('+'),b,input])
    
    print(equation)
    

if __name__ == "__main__":
    vehicle_select = "mxs"
    # vehicle_select = "bixler"
    if vehicle_select == "mxs":
        from mxs import vehicle
        control_indices = [1,2]
        control_symb = [
            '\u03b7',
            '\u03c4',
            ]
        with_sweeps = False
    else:
        from bixler import vehicle
        control_indices = [1,2,4]
        control_symb = [
            '\u03b7',
            '\u03c4',
            '\u039b',
            ]
        with_sweeps = False

    np.set_printoptions(precision=3)
    
    if with_sweeps:
        airspeeds = [15]
        sweeps = np.linspace(-10,20,200)
        eigenvalues = np.zeros((len(sweeps),5), dtype='complex128')
    else:
        airspeeds = np.linspace(12,25,50)
        sweeps = [0]
        eigenvalues = np.zeros((len(airspeeds),5), dtype='complex128')
    
    for i,airspeed in enumerate(airspeeds):
        for j,sweep in enumerate(sweeps):
            (alpha,elevator,throttle) = get_trim_condition(vehicle,airspeed)
            state = calc_state(math.radians(alpha),airspeed)
            a,b = linearize(vehicle, (state,[0,math.radians(elevator),throttle,0,math.radians(sweep),0]), control_indices)
            print(f"----- {airspeed:5.1f} m/s -----")
            # print_eqn(a,b)
            pretty_eqn(a,b,control_symb)
            vals, vecs = np.linalg.eig(a)
            index = j if with_sweeps else i
            eigenvalues[index,:] = vals
            for val,vec in zip(vals,vecs):
                print(f"{val}: {vec}")

    import matplotlib.pyplot as plt
    
    for eigval_index in range(5):
        plt.scatter(
            np.real(eigenvalues[:,eigval_index]),
            np.imag(eigenvalues[:,eigval_index]),
            c=(sweeps if with_sweeps else airspeeds),
        )
    
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.minorticks_on()
    plt.grid(True,'both')
    if with_sweeps:
        plt.colorbar(label="Sweep (deg)")
    else:
        plt.colorbar(label="Airspeed (m/s)")
    plt.axis('equal')
    # plt.ylim([-20,20])
    plt.show()
