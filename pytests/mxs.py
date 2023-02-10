#!/usr/bin/env python3

import math

from numpy import mean, var

from pyaerso import AffectedBody, AeroBody, Body, Force, Torque

def calc_state(alpha,airspeed,combined=True):
    u = airspeed * math.cos(alpha)
    w = airspeed * math.sin(alpha)
    orientation = [0, math.sin(alpha/2), 0, math.cos(alpha/2)]
    if combined:
        return [0,0,0, u,0,w, *orientation, 0,0,0]
    else:
        return [0,0,0], [u,0,w], orientation, [0,0,0]

# Trim points airspeed(m/s): (alpha(deg),elev(deg),throttle)
trim_points = {
    15.0: ( 2.0574563762079516,  -6.8727513934777225, 0.24329486036084907),
    17.5: ( 1.0706920728883487,  -7.920398290900842,  0.32366463553969765),
    20.0: ( 0.4300108890771458,  -8.847841922284317,  0.43160143104358356),
    22.5: (-0.00919178945532153, -9.600216932754634,  0.5593081244658515 ),
    25.0: (-0.3232762564889445,  -9.890387987246884,  0.700473854936423  )
}


selected_trim_point = 17.5

mass = 1.221
inertia = [
    [0.019, 0.0,   0.0],
    [0.0,   0.09, 0.0],
    [0.0,   0.0,   0.121]]
position,velocity,attitude,rates = calc_state(
    math.radians(trim_points[selected_trim_point][0]),
    selected_trim_point,
    False
)

class WindModel:
    def get_wind(self,position):
        return [0,0,0]
    def step(self,dt):
        pass


class DensityModel:
    def get_density(self,position):
        return 1.225
        # import isacalc
        # [_,_,density,_,_] = isacalc.calculate_at_h(-position[2])
        # return density

def H(x,p,k=10):
    return 1.0 / ( 1.0 + math.exp(-2*k*(x-p)) )

def S(x,l,h,k=10):
    return H(x,l,k)*(1-H(x,h,k))

def deg2rad(deg):
    return deg/180.0 * math.pi

Sw = 0.263
cw = 0.24
S_t = 0.0825
x_t = -0.585
c_t = 0.165

def c_lta(alpha):
    a_stall = math.radians(10)
    a_decay = math.radians(40)
    return 3.5810*S(alpha,-a_stall,a_stall)*alpha \
        + 0.65 * S(alpha,a_stall,a_decay) \
        + (math.pi/2-0.25-0.8*alpha) * S(alpha,a_decay,math.pi) * (1-H(alpha,math.radians(145),3)) \
        - 0.65 * S(alpha,-a_decay,-a_stall) \
        + (-math.pi/2+0.25-0.8*alpha) * S(alpha,-math.pi,-a_decay) * H(alpha,math.radians(-145),3)

def c_dta(alpha):
    a_lim = math.radians(9)
    ao_lim = math.radians(30)
    am_lim = math.radians(150)
    return 2 * S(alpha,-a_lim,a_lim,12) * (0.5*alpha)**2 \
        + S(alpha,-ao_lim,ao_lim,10) * (1-S(alpha,-a_lim,a_lim,12)) * 1.5 * (abs(alpha) - 0.05) \
        + 1.05 * (1-S(alpha,-ao_lim,ao_lim,10)) * S(alpha,-am_lim,am_lim,5) \
        + 0.018

def poly_manifold(x,y,z,order,*args):    
    if order == 0:
        return args[0]
    if order == 1:
        return args[0] * x + args[1] * y + args[2] * z + args[3]
    if order == 2:
        return args[0]*x*x + args[1]*x*y + args[2]*x*z + args[3]*y*y + args[4]*y*z + args[5]*z*z \
            + args[6]*x + args[7]*y + args[8]*z \
            + args[9]
    if order == 3:
        return args[0]*x*x*x + args[1]*x*x*y + args[2]*x*x*z + args[3]*x*y*y + args[4]*x*y*z + args[5]*x*z*z + args[6]*y*y*x + args[7]*y*y*y + args[8]*y*y*z + args[9]*y*z*z + args[10]*z*z*z \
            + args[11]*x*x + args[12]*x*y + args[13]*y*y + args[14]*x*z + args[15]*y*z + args[16]*z*z \
            + args[17]*x + args[18]*y + args[19]*z \
            + args[20]

def clamp(x,u,l):
    return max(min(x,u),l)

class Combined:
    def __init__(self):
        self.time = 0

    def get_elevator_moment_coeff_delta(self,airstate,rates,input):
        coeffs = [
            -1.09419902e+00,
            -2.11376232e-01,
            2.68838823e-02,
            -2.19078045e+02,
            -1.85268630e-01,
            4.81678617e-03,
            2.19458992e+02,
            -4.51607587e-01,
            2.38686253e-03,
            -1.87506464e-04,
            1.22846953e-03,
            -4.61863302e-01,
            3.63517580e+00,
            5.27810084e-01,
            -1.51973289e-01,
            4.08504874e-03,
            -5.98825396e-02,
            2.59581816e+00,
            -1.73627947e-01,
            9.31570574e-01,
            -4.57208842e+00
        ]
        elev = input[1]
        throttle = input[2]
        aspd = clamp(airstate[2],10,22)
        return poly_manifold(elev,throttle,aspd,3,*coeffs)
    
    def get_thrust(self,airspeed,throttle):
        return 1.1576e1 * throttle**2 \
            + -4.8042e-1 * throttle * airspeed \
            + -1.1822e-2 * airspeed**2 \
            + 1.3490e1 * throttle \
            + 4.6518e-1 * airspeed \
            + -4.1090
    
    def c_l_curve(self,alpha):
        [cl_0, cl_alpha, pstall, nstall] = [ 0.14760915,  5.02005742,  0.21634284, -0.31512553]
        k_stall = 25
        k_decay = 10
        a_pdecay = math.radians(30)
        a_ndecay = math.radians(-35)
        linear_regime = S(alpha,nstall,pstall,k_stall) * (cl_0 + cl_alpha * alpha)
        
        post_stall_pos = H(alpha,pstall,k_stall) * (1-H(alpha,a_pdecay,k_decay)) * (alpha + 0.45)
        high_alpha_regime = H(alpha,a_pdecay,k_decay) * 0.8 * ((math.pi/2)-alpha) * (1-H(alpha,math.radians(145),3))
        
        post_stall_neg = H(alpha,a_ndecay,k_decay) * (1-H(alpha,nstall,k_stall)) * (alpha - 0.35)
        low_alpha_regime = (1-H(alpha,a_ndecay,k_decay)) * 0.8 * ((-math.pi/2)-alpha) * H(alpha,math.radians(-145),3)
        
        return linear_regime \
            + post_stall_pos + high_alpha_regime \
            + post_stall_neg + low_alpha_regime
    
    def c_d_curve(self,alpha):
        [cd_0,cd_alpha,alpha_cd0] = [-0.15458247, 1.1308512,  0.05247888]
        alpha_lim = math.radians(40)
        k = 2
        return S(alpha,-alpha_lim,alpha_lim,k) * (cd_alpha*(alpha-alpha_cd0)**2 + cd_0) \
            + (1.0-S(alpha,-alpha_lim,alpha_lim,k)) * 1.8
    
    def get_cm(self,alpha,c_lw):
        downwash_angle = 2 * c_lw / (math.pi * 4.54)
        alpha = alpha - downwash_angle + math.radians(-0.5)
        tail_lift = c_lta(alpha)
        tail_drag = c_lta(alpha)
        resolved_moment = x_t * S_t * (tail_lift * math.cos(alpha) + tail_drag * math.sin(alpha))
        effective_cm = resolved_moment / (Sw*cw)
        return effective_cm
    
    def get_effect(self,airstate,rates,input):
        alpha = airstate[0]
        c_l = self.c_l_curve(alpha)
        c_d = self.c_d_curve(alpha)
        c_m = self.get_cm(alpha,c_l)
        
        q = airstate[3]
        V = airstate[2]
        
        downwash_angle = (2 * c_l / (math.pi * 4.54))
        alpha_t = alpha - downwash_angle + math.radians(-0.75)
        alpha_q = math.atan2( (-rates[1]*x_t + V*math.sin(alpha_t)), (V*math.cos(alpha_t)) )
        m_t_aq = x_t * q * S_t * (
            math.cos(alpha_q) * c_lta(alpha_q) \
            + math.sin(alpha_q) * c_dta(alpha_q) \
            - math.cos(math.radians(alpha_t)) * c_lta(math.radians(alpha_t)) \
            - math.sin(math.radians(alpha_t)) * c_dta(math.radians(alpha_t))
        )
        
        # m_t_aq = S(rates[1],-2,2) * x_t * q * S_t * (c_lta(alpha_q) - c_lta(alpha)) \
        #        + (1-S(rates[1],-2,2)) * -0.4 * rates[1]
        
        dc_m_elev = self.get_elevator_moment_coeff_delta(airstate,rates,input)
        
        lift = q * Sw * c_l
        drag = q * Sw * c_d
        moment = q * Sw * cw * (c_m + dc_m_elev) + m_t_aq
        
        thrust = self.get_thrust(V,input[2])
        # print(f"t: {self.time:.3f} T: {thrust:.3f}, V: {V:.3f}, Alpha: {alpha:.3f}, q: {q:.3f},  c_m: {c_m:.3f}, alpha_q: {alpha_q:.3f}, c_m_eff: {c_m+dc_m_elev:.3f}, moment: {moment:.3f}, dc_m_elev: {dc_m_elev:.3f}, mtaq: {m_t_aq:.3f}")
        self.time += 0.01
        x = thrust - drag * math.cos(alpha) + lift * math.sin(alpha)
        z = -lift * math.cos(alpha) - drag * math.sin(alpha)
        
        # Force and moment fits are found relative to the loadcell point
        # Computed effect forces/moments should be supplied about CG
        #x_cg = 0.03 # m
        #z_cg = -0.016 # m Load cell reference point was 16mm below spar centreline
        
        x_cg = 0.03 # m
        z_cg = 0 # m
        moment = moment - z * x_cg + x * z_cg
        
        return (
            Force.body([x,0.0,z]),
            Torque.body([0,moment,0])
            )


def sensible(array,dp=4):
    elements = ", ".join([f"{v:.{dp}f}" for v in array])
    return f"[{elements}]"

body = Body(mass,inertia,position,velocity,attitude,rates)

# aerobody = AeroBody(body,WindModel(),("StandardDensity",[]))
# aerobody = AeroBody(body)
aerobody = AeroBody(body,None,DensityModel())

# vehicle = AffectedBody(aerobody,[Lift(),Drag(),Moment()])
vehicle = AffectedBody(aerobody,[Combined()])

def get_alpha_q(alpha,q,V=15):
    result = math.atan2( (-q*x_t + V*math.sin(alpha)), (V*math.cos(alpha)) )
    if result - alpha > 1:
        result = result - 2*math.pi
    return result

def get_cmq(alpha,q,V=15):
    dyn_press = 0.5 * 1.225 * V**2
    alpha_q = get_alpha_q(alpha,q,V)
    m_t_aq = x_t * dyn_press * S_t * (c_lta(alpha_q) - c_lta(alpha))
    return m_t_aq / (Sw * dyn_press * cw * q)

def get_cmqd(alpha, q, V=15):
    dyn_press = 0.5 * 1.225 * V**2
    alpha_q = get_alpha_q(alpha,q,V)
    m_t_aq = x_t * dyn_press * S_t * (
        math.cos(alpha_q) * c_lta(alpha_q) + math.sin(alpha_q) * c_dta(alpha_q)
        - math.cos(alpha) * c_lta(alpha) - math.sin(alpha) * c_dta(alpha)
    )
    return m_t_aq / (Sw * dyn_press * cw * q)

def get_cmqd_dw(alpha, q, V=15):
    dyn_press = 0.5 * 1.225 * V**2

    downwash_angle = 2 * Combined().c_l_curve(alpha) / (math.pi * 4.54) * math.cos(alpha)
    alpha_t = alpha - downwash_angle + math.radians(-0.5)

    alpha_q = get_alpha_q(alpha_t,q,V)
    m_t_aq = x_t * dyn_press * S_t * (
        math.cos(alpha_q) * c_lta(alpha_q) + math.sin(alpha_q) * c_dta(alpha_q)
        - math.cos(alpha_t) * c_lta(alpha_t) - math.sin(alpha_t) * c_dta(alpha_t)
    )
    return m_t_aq / (Sw * dyn_press * cw * q)
    

# print(vehicle.airstate)
if __name__ == "__main__":
    # import sys
    # import numpy as np
    # import matplotlib.pyplot as plt
    
    # alphas = np.linspace(-180,180,200)
    
    # # cds = np.zeros_like(alphas)
    # # for i,alpha in enumerate(alphas):
    # #     cds[i] = c_dta(math.radians(alpha))
    # # plt.plot(alphas,cds)

    # qs = np.linspace(1,60,3)
    # cmqs = np.zeros((len(qs),len(alphas),3))
    # c_ltas = np.zeros((len(qs),len(alphas),2))
    # alpha_qs = np.zeros((len(qs),len(alphas)))
    
    # for j,alpha in enumerate(alphas):
    #     for i,q in enumerate(qs):
    #         cmqs[i,j,0] = get_cmq(math.radians(alpha),math.radians(q))
    #         cmqs[i,j,1] = get_cmqd(math.radians(alpha),math.radians(q))
    #         cmqs[i,j,2] = get_cmqd_dw(math.radians(alpha),math.radians(q))
    #         alpha_q = get_alpha_q(math.radians(alpha),math.radians(q))
    #         c_ltas[i,j,0] = c_lta(alpha_q) - c_lta(math.radians(alpha))
    #         c_ltas[i,j,1] = math.cos(alpha_q) * c_lta(alpha_q) + math.sin(alpha_q) * c_dta(alpha_q) \
    #                         - math.cos(math.radians(alpha)) * c_lta(math.radians(alpha)) - math.sin(math.radians(alpha)) * c_dta(math.radians(alpha))
    #         alpha_qs[i,j] = alpha_q

    # for i,q in enumerate(qs):
    #     # plt.plot(alphas,cmqs[i,:,0])
    #     # plt.plot(alphas,cmqs[i,:,1])
    #     plt.plot(alphas,cmqs[i,:,2])
    # plt.legend(list(map(lambda v: str(v),qs)))
    # plt.grid(True,'both')
    
    # plt.figure()
    # for i,q in enumerate(qs):
    #     plt.plot(alphas,np.degrees(alpha_qs[i,:]))
    # plt.legend(list(map(lambda v: str(v),qs)))
    
    # plt.figure()
    # for i,q in enumerate(qs):
    #     # plt.plot(alphas,c_ltas[i,:,0])
    #     plt.plot(alphas,c_ltas[i,:,1])
    # plt.legend(list(map(lambda v: str(v),qs)))
    
    # plt.show()
    
    # sys.exit(0)
    
    scale = 1
    deltaT = 0.01/scale

    print(sensible(vehicle.statevector))

    import sys
    import time

    outfile = None
    if len(sys.argv) > 1:
        outfile = open(sys.argv[1],"w")
        outfile.write("time,x,y,z,u,v,w,qx,qy,qz,qw,p,q,r,alpha,elevator\n")

    samples = 1
    sample_times = []

    TRIM_ELEVATOR = trim_points[selected_trim_point][1]
    TRIM_THROTTLE = trim_points[selected_trim_point][2]

    def get_doublet_elevator_input(count):
        if count < 500*scale:
            return math.radians(TRIM_ELEVATOR), TRIM_THROTTLE
        if count < 550*scale:
            return math.radians(TRIM_ELEVATOR+5.0), TRIM_THROTTLE
        if count < 600*scale:
            return math.radians(TRIM_ELEVATOR-5.0), TRIM_THROTTLE
        
        return math.radians(TRIM_ELEVATOR), TRIM_THROTTLE

    def get_loop_input(count):
        throttle = TRIM_THROTTLE
        duration = 600
        if count < 500*scale:
            return math.radians(TRIM_ELEVATOR), throttle
        if count < (500+duration)*scale:
            return math.radians(TRIM_ELEVATOR+15.0), throttle

        if count < 1200*scale:
            return math.radians(TRIM_ELEVATOR), throttle
        if count < (1200+duration-100)*scale:
            return math.radians(TRIM_ELEVATOR+15.0), throttle

        if count < 1900*scale:
            return math.radians(TRIM_ELEVATOR), throttle
        if count < (1900+duration)*scale:
            return math.radians(TRIM_ELEVATOR+15.0), throttle

        return math.radians(TRIM_ELEVATOR), throttle

    def get_elevator_input(count):
        if count < 500*scale:
            return math.radians(TRIM_ELEVATOR), TRIM_THROTTLE
        if count < 600*scale:
            return math.radians(TRIM_ELEVATOR+2.0), TRIM_THROTTLE
        # if count < 700*scale:
        #     return math.radians(TRIM_ELEVATOR), 0.65
        # if count < 1100*scale:
        #     return math.radians(TRIM_ELEVATOR+5.0), 0.65

        return math.radians(TRIM_ELEVATOR), TRIM_THROTTLE
    
    def constant(level):
        class Constant:
            def __call__(self, t):
                return level
        return Constant
    
    def pulse(amplitude, t_start, period):
        class Pulse:
            def __call__(self, t):
                if t < t_start:
                    return 0
                if t < (t_start + period):
                    return amplitude
                return 0
    
    def doublet(amplitude, t_start, period):
        start_pulse = pulse(amplitude, t_start, period/2)
        end_pulse = pulse(amplitude, t_start + period/2, period/2)
        class Doublet:
            def __call__(self, t):
                return start_pulse(t) + end_pulse(t)
        return Doublet

    for i in range(samples):
        count = 0
        simtime = 0
        
        body = Body(mass,inertia,position,velocity,attitude,rates)
        # aerobody = AeroBody(body,None,DensityModel())
        aerobody = AeroBody(body,None,("StandardDensity",[]))
        vehicle = AffectedBody(aerobody,[Combined()])
        
        start = time.process_time()
        while count < 3000*scale:
            elevator, throttle = get_doublet_elevator_input(count)
            # elevator, throttle = get_loop_input(count)
            # elevator, throttle = get_elevator_input(count)
            vehicle.step(deltaT,[0,elevator,throttle])
            count += 1
            simtime += deltaT
            vehicle.statevector
            if outfile:
                outfile.write(f"{simtime},"+sensible(vehicle.statevector,10)[1:-1] + f",{vehicle.airstate[0]},{elevator}\n")
        end = time.process_time()
        sample_times.append(end-start)
        print(sensible(vehicle.statevector))

    print(f"Mean: {mean(sample_times)}\nVar: {var(sample_times)}")
    if outfile is not None:
        outfile.close()
