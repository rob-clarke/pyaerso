#!/usr/bin/env python3

import math

from numpy import mean, var

from pyaerso import AffectedBody, AeroBody, Body, Force, Torque

mass = 1.285
inertia = [
    [0.042,     0, 0    ],
    [    0, 0.029, 0    ],
    [    0,     0, 0.066],
]
Sw = 0.26
cw = 0.2

def calc_state(alpha,airspeed,combined=True):
    u = airspeed * math.cos(alpha)
    w = airspeed * math.sin(alpha)
    orientation = [0, math.sin(alpha/2), 0, math.cos(alpha/2)]
    if combined:
        return [0,0,0, u,0,w, *orientation, 0,0,0]
    else:
        return [0,0,0], [u,0,w], orientation, [0,0,0]

def interpolate(x_target, x_data, y_data):
    if x_target < x_data[0]:
        m = (y_data[0] - y_data[1])/(x_data[0] - x_data[1])
        return m * x_target + (y_data[0] - m * x_data[0])
    if x_target > x_data[-1]:
        m = (y_data[-2] - y_data[-1])/(x_data[-2] - x_data[-1])
        return m * x_target + (y_data[-1] - m * x_data[-1])
    #highidx = np.argmax( x_data < np.array(x_target) )
    highidx = 1
    for i in range(1,len(x_data)):
        if x_data[i] > x_target:
            highidx = i
            break
    lowidx = highidx - 1
    m = (y_data[lowidx] - y_data[highidx])/(x_data[lowidx] - x_data[highidx])
    return m * x_target + (y_data[lowidx] - m * x_data[lowidx])

class BixlerForces:
    def get_effect(self, airstate, rates, input):
        (
            C_L0, C_Lalpha, C_Lq,
            C_D0, C_Dalpha,
            C_Yv,
            C_l0, C_ldTipStbd, C_ldTipPort, C_lp,
            C_m0, C_malpha, C_melev, C_msweep, C_mwashout, C_mq,
            C_n0, C_nrudder, C_nr
        ) = self.get_coefficients(airstate, rates, input)
        
        alpha_rad = airstate[0]
        alpha = math.degrees(airstate[0])
        Q = airstate[3]
        elev = math.degrees(input[1])
        sweep = math.degrees(input[4])
        washout = math.degrees(input[5])
        
        C_D = C_D0 + C_Dalpha * alpha
        D = Q * Sw * C_D
        X = -D

        # Lift
        C_L = C_L0 + (C_Lalpha * alpha) + (C_Lq * rates[1])
        L = Q * Sw * C_L
        Z = -L

        l = 0.0
        
        # Pitching moment
        C_m = C_m0 + (C_malpha * alpha) + (C_melev * elev) + (C_msweep * sweep) + (C_mwashout * washout) + (C_mq * rates[1])
        m = Q * Sw * cw * C_m
        #print("0: {}, a: {}, e: {}, s: {}, w: {}, q: {}".format(C_m0,(C_malpha * self.alpha),(C_melev * self.elev),(C_msweep * self.sweep),(C_mwashout * self.washout),(C_mq * q)))
        #print("{},{},{},{},{},{}".format(C_m0,(C_malpha * self.alpha),(C_melev * self.elev),(C_msweep * self.sweep),(C_mwashout * self.washout),(C_mq * q)))
        #print("a: {}, e: {}, s: {}, w: {}, q: {}".format(self.alpha,self.elev,self.sweep,self.washout,q))
        #print("alpha: {}, u: {}, w: {}".format(self.alpha,self.velocity_b[0,0],self.velocity_b[2,0]))
        #print("Moment: {}".format(m))
        
        # Yawing moment
        n = 0.0
        
        x = 10*input[2] - D * math.cos(alpha_rad) + L * math.sin(alpha_rad)
        z = -L * math.cos(alpha_rad) - D * math.sin(alpha_rad)
        
        # Force and moment fits are found relative to the loadcell point
        # Computed effect forces/moments should be supplied about CG
        #x_cg = 0.03 # m
        #z_cg = -0.016 # m Load cell reference point was 16mm below spar centreline
        
        # x_cg = 0.03 # m
        # z_cg = -0.016 # m
        # moment = moment - z * x_cg + x * z_cg
        
        return (
            Force.body([x,0.0,z]),
            Torque.body([0,m,0])
            )
    
    def get_coefficients(self, airstate, rates, input):
        C_L0, C_Lalpha, C_Lq = self._get_coefficients_C_L(airstate, rates, input)
        C_D0, C_Dalpha = self._get_coefficients_C_D(airstate, rates, input)
        C_Yv = self._get_coefficients_C_Y(airstate, rates, input)
        C_l0, C_ldTipStbd, C_ldTipPort, C_lp = self._get_coefficients_C_l(airstate, rates, input)
        C_m0, C_malpha, C_melev, C_msweep, C_mwashout, C_mq = self._get_coefficients_C_m(airstate, rates, input)
        C_n0, C_nrudder, C_nr = self._get_coefficients_C_n(airstate, rates, input)
        
        return (
            C_L0, C_Lalpha, C_Lq,
            C_D0, C_Dalpha,
            C_Yv,
            C_l0, C_ldTipStbd, C_ldTipPort, C_lp,
            C_m0, C_malpha, C_melev, C_msweep, C_mwashout, C_mq,
            C_n0, C_nrudder, C_nr
            )

    def _get_coefficients_C_L(self, airstate, rates, input):
        # Definition of characteristic Aerodynamic Coefficient Arrays obtained from wind tunnel testing 
        # The tested airspeeds being [6 8 10 12 14] m/s
        airspeeds = [6, 8, 10, 12, 14]
        
        alpha = airstate[0]
        V = airstate[2]
        sweep = math.degrees(input[4])
        
        # Lift coefficient data
        CL0_m5to5 = [ 0.3417,  0.3771, 0.39, 0.41, 0.43 ]
        CL0_5to10 = [ 0.634, 0.6106, 0.604, 0.586, 0.593 ]
        CL0_over10 = [ 1.103, 1.077, 1.0296, 1.026, 1.0313 ]

        CLalpha_m5to5 = [ 0.085, 0.08117, 0.0791, 0.0778, 0.0753 ]
        CLalpha_5to10 = [ 0.0299, 0.0354, 0.037, 0.039, 0.039 ]
        CLalpha_over10 = [ -0.0170, -0.0113, -0.0056, -0.005, -0.0048 ] # Post-stall

        # Dynamic Pitching Lift coefficient
        speed_sample_dynamic = [6, 8, 10]
        sweep_sample_dynamic = [0, 10, 20, 30]

        CLq_6ms  = [0.01012, 0.01217, 0.01005, 0.010104]
        CLq_8ms  = [0.008201, 0.0081213, 0.008229, 0.008633]
        CLq_10ms = [0.00555, 0.0073, 0.007708, 0.006868]

        C_L0_samples = []
        C_Lalpha_samples = []
        
        if alpha < 5:
            C_L0_samples = CL0_m5to5
            C_Lalpha_samples = CLalpha_m5to5
        elif (alpha >= 5) and (alpha < 10):
            C_L0_samples = CL0_5to10
            C_Lalpha_samples = CLalpha_5to10
        elif (alpha >= 10):
            C_L0_samples = CL0_over10
            C_Lalpha_samples = CLalpha_over10

        C_L0     = interpolate(V, airspeeds, C_L0_samples)
        C_Lalpha = interpolate(V, airspeeds, C_Lalpha_samples)

        C_Lq = 0.0

        if rates[1] < 0:
            C_Lq = 0
        else:
            # scipy.interpolate.RectBivariateSpline?
            CLq6 = interpolate(sweep, sweep_sample_dynamic, CLq_6ms)
            CLq8 = interpolate(sweep, sweep_sample_dynamic, CLq_8ms)
            CLq10 = interpolate(sweep, sweep_sample_dynamic, CLq_10ms)
            CLq_speeds = [CLq6, CLq8, CLq10];		
            C_Lq = interpolate( V, speed_sample_dynamic, CLq_speeds)

        return (C_L0, C_Lalpha, C_Lq)

    def _get_coefficients_C_D(self, airstate, rates, input):
        # The tested airspeeds being [6 8 10 12 14] m/s
        airspeeds = [6, 8, 10, 12, 14]

        CD0_m5to5 =  [  0.0682 ,  0.075 ,  0.0765, 0.087 ,  0.084  ]
        CD0_5to10 =  [ -0.02125, -0.0054,  0.0259, 0.0394,  0.0456 ]
        CD0_over10 = [  0.0027 , -0.0167, -0.0715, -0.06 , -0.039  ]
        
        
        CDalpha_m5to0 =  [ 0.0006 , 0      , 0.0003 , 0.0003 , 0.0005  ] 
        CDalpha_0to5 =   [ 0.0067 , 0.00463, 0.00584, 0.00414, 0.00436 ]
        CDalpha_5to10 =  [ 0.02457, 0.02071, 0.01596, 0.01366, 0.01204 ]
        CDalpha_over10 = [ 0.02217, 0.02185, 0.0257 , 0.0236 , 0.0205  ] #Post-Stall

        C_D0_samples = []
        C_Dalpha_samples = []
        
        alpha = airstate[0]
        V = airstate[2]
        
        if alpha < 0:
            C_D0_samples = CD0_m5to5
            C_Dalpha_samples = CDalpha_m5to0
        elif (alpha >= 0) and (alpha < 5):
            C_D0_samples = CD0_m5to5
            C_Dalpha_samples = CDalpha_0to5
        elif (alpha >= 5) and (alpha < 10):
            C_D0_samples = CD0_5to10
            C_Dalpha_samples = CDalpha_5to10
        elif alpha >= 10:
            C_D0_samples =  CD0_over10
            C_Dalpha_samples = CDalpha_over10

        C_D0     = interpolate(V, airspeeds, C_D0_samples)
        C_Dalpha = interpolate(V, airspeeds, C_Dalpha_samples)

        return (C_D0, C_Dalpha)

    def _get_coefficients_C_Y(self, airstate, rates, input):
        C_Yv = -0.05 # Coeff defined relative to sideways velocity!
        return C_Yv

    def _get_coefficients_C_l(self, airstate, rates, input):
        # Wing Tip Port and Starboard effectiveness Rolling Moment Coefficient
        sweep_sample_wingtip = [-10, 0, 10, 20]
        alpha_sample_wingtip = [-5, 0, 5, 10, 15, 20, 25]

        CMXportwingtip_m10sweep = [0.0230, 0.02153, 0.022415, 0.016924, 0.01618, 0.011757, 0.010183]
        CMXportwingtip_0sweep = [0.026895, 0.02629, 0.02659, 0.02242, 0.018992, 0.01786, 0.01563]
        CMXportwingtip_10sweep = [0.02987, 0.03098, 0.0323, 0.03146, 0.0292, 0.02529, 0.01958]
        CMXportwingtip_20sweep = [0.0289, 0.0304, 0.03098, 0.0225, 0.03213, 0.02575, 0.0293]

        CMXstarboardwingtip_m10sweep = [-0.01822, -0.019, -0.020148, -0.02041, -0.02316, -0.02308, -0.018754]
        CMXstarboardwingtip_0sweep = [-0.02083, -0.02190, -0.02313, -0.02769, -0.0327, -0.026486, -0.02888]
        CMXstarboardwingtip_10sweep = [-0.02062, -0.0253, -0.03089, -0.03504, -0.03468, -0.0349, -0.03094]
        CMXstarboardwingtip_20sweep = [-0.02167, -0.02885, -0.02985, -0.02913, -0.03361, -0.027941, -0.027445]
        
        C_l0 = 0.0
        C_ldTipPort = 0.0
        C_ldTipStbd = 0.0
        C_lp = -0.4950 # Eye average from Stanford student projects. Units? Probably per radian...
    
        alpha = airstate[0]
        V = airstate[2]
        elev = math.degrees(input[1])
        sweep = math.degrees(input[4])
        washout = math.degrees(input[5])
    
        if (alpha + 100) >= 15:
            C_ldTipPort = 0.0
        else:
            CMXportwingtipM10 = interpolate(alpha,alpha_sample_wingtip,CMXportwingtip_m10sweep)
            CMXportwingtip0 = interpolate(alpha,alpha_sample_wingtip,CMXportwingtip_0sweep)
            CMXportwingtip10 = interpolate(alpha,alpha_sample_wingtip,CMXportwingtip_10sweep)
            CMXportwingtip20 = interpolate(alpha,alpha_sample_wingtip,CMXportwingtip_20sweep)
            
            CMXportwingtip_sweeps = [CMXportwingtipM10, CMXportwingtip0, CMXportwingtip10, CMXportwingtip20]

            C_ldTipPort = interpolate(sweep,sweep_sample_wingtip,CMXportwingtip_sweeps)
    
        if (alpha + 100) >= 15:
            C_ldTipStbd = 0.0
        else:
            CMXstarboardwingtipM10 = interpolate(alpha,alpha_sample_wingtip,CMXstarboardwingtip_m10sweep)
            CMXstarboardwingtip0 = interpolate(alpha,alpha_sample_wingtip,CMXstarboardwingtip_0sweep)
            CMXstarboardwingtip10 = interpolate(alpha,alpha_sample_wingtip,CMXstarboardwingtip_10sweep)
            CMXstarboardwingtip20 = interpolate(alpha,alpha_sample_wingtip,CMXstarboardwingtip_20sweep)
            CMXstarboardwingtip_sweeps = [CMXstarboardwingtipM10, CMXstarboardwingtip0, CMXstarboardwingtip10, CMXstarboardwingtip20]

            C_ldTipStbd = interpolate(sweep, sweep_sample_wingtip, CMXstarboardwingtip_sweeps)

        return (C_l0, C_ldTipStbd, C_ldTipPort, C_lp)

    def _get_coefficients_C_m(self, airstate, rates, input):
        # The tested airspeeds being [6 8 10 12 14] m/s
        airspeeds = [6, 8, 10, 12, 14]
        # Dynamic sample points
        speed_sample_dynamic = [6, 8, 10]
        sweep_sample_dynamic = [0, 10, 20, 30]
        
        # Pitching Moment Coefficients
        CMY0_m5to14 = [  0.01045,  0.0215,  0.02 ,  0.01 , -0.006  ]
        CMY0_over14 = [ -0.29075, -0.304 , -0.323, -0.307, -0.3065 ]

        CMYalpha_m5to0 =  [ -0.01185, -0.0115 , -0.01152, -0.01054, -0.01224 ]
        CMYalpha_0to14 =  [ -0.02151, -0.02325, -0.0245 , -0.02264, -0.02154 ]
        CMYalpha_over14 = [ 0       ,        0,        0,        0,        0 ]


        # Dynamic Pitching Moment Coefficient
        CMYq_6ms = [-0.0024608, -0.0033244, -0.003728, -0.0046806]
        CMYq_8ms = [-0.002354, -0.00269025, -0.0030326, -0.0038794]
        CMYq_10ms = [-0.002311, -0.0023768, -0.002771, -0.0029218]

        # Elevator effectiveness Pitching Moment Coefficient
        alpha_sample_elev = [-5, -2.5, 0, 2.5, 5]

        CMYelev_8p5to4p4 = [-0.00480, -0.00523, -0.006512, -0.004707, -0.00573]
        CMYelev_4p4to1p8 = [-0.01088, -0.01069, -0.011315, -0.01111, -0.01004]
        CMYelev_1p8tom1p2 = [-0.00533, -0.0085, -0.008993, -0.0099, -0.0103]
        CMYelev_m1p2tom5 = [-0.002184, -0.004632, -0.006184, -0.005978, -0.006237]

        CMYelev_overall = [-0.005355, -0.0068404, -0.0073926, -0.00745304, -0.0077259]

        # Wing Sweep effectiveness Pitching Moment Coefficient
        speed_sample_sweep = [6, 8, 10, 12, 14]
        alpha_sample_sweep = [-5, 0, 5, 10, 13, 15, 20, 25]

        # CMYsweep_8ms_5aoa = 0.011973;
        CMYsweep_6ms = [-0.00302, 0.000745, 0.003357, 0.00472, 0.005048, 0.00518, 0.00528, 0.00573]
        CMYsweep_8ms = [-0.00221, 0.00363, 0.00837, 0.01157, 0.01191, 0.0107, 0.0124, 0.0131]
        CMYsweep_10ms = [-0.00083, 0.00889, 0.01672, 0.02257, 0.02328, 0.02455, 0.0236, 0.0278]
        CMYsweep_12ms = [-0.00125, 0.0140, 0.02631, 0.03646, 0.03951]
        CMYsweep_14ms = [-0.00339, 0.01116, 0.02494]

        # Wing Tip Washout effectiveness Pitching Moment Coefficient
        sweep_sample_washout = [10, 20, 30]
        alpha_sample_washout = [-5, 0, 5, 10, 15, 20, 25]

        CMYwashout_0sweep_m30tom10deflec = [-0.00331, -0.003047, -0.00362, -0.00428, -0.00511, -0.001305, -0.000811]
        CMYwashout_0sweep_m10to5deflec = [0.01109, 0.009013, 0.01248, 0.01463]

        CMYwashout_10sweep = [0.00358, 0.00492, 0.0072, 0.009088, 0.00869, 0.009086, 0.01204]
        CMYwashout_20sweep = [0.01247, 0.0134, 0.01282, 0.01429, 0.01214, 0.01406, 0.01528]
        CMYwashout_30sweep = [0.01679, 0.01582, 0.01279, 0.01254, 0.01307, 0.01519, 0.01659]

        C_m0_samples = []
        C_malpha_samples = []

        alpha = airstate[0]
        V = airstate[2]
        elev = math.degrees(input[1])
        sweep = math.degrees(input[4])
        washout = math.degrees(input[5])

        if alpha < 0:
            C_m0_samples = CMY0_m5to14
            C_malpha_samples = CMYalpha_m5to0
        elif (alpha >= 0) and (alpha < 14):
            C_m0_samples = CMY0_m5to14
            C_malpha_samples = CMYalpha_0to14
        elif alpha >= 14:
            C_m0_samples = CMY0_over14
            C_malpha_samples = CMYalpha_over14
        
        C_m0 = interpolate(V, airspeeds, C_m0_samples)
        #print('ASI: {}, result: {}, speeds: {}, samples: {}'.format(self.airspeed,C_m0,airspeeds,C_m0_samples))
        C_malpha = interpolate(V, airspeeds, C_malpha_samples)
        
        C_melev_samples = []
        
        # Elevator effectiveness coefficient estimation
        if elev > 4.4:
            C_melev_samples = CMYelev_8p5to4p4
        elif (elev <= 4.4) and (elev > 1.8):
            C_melev_samples = CMYelev_4p4to1p8
        elif (elev <= 1.8) and (elev > -1.2):
            C_melev_samples = CMYelev_1p8tom1p2
        elif elev <= -1.2:
            C_melev_samples = CMYelev_m1p2tom5

        # Override the above to match JAVA_MODEL
        C_melev_samples = CMYelev_overall

        C_melev = interpolate(alpha,alpha_sample_elev,C_melev_samples)

        # Wing Sweep effictiveness coefficient estimation
        CMYsweep6 = interpolate(alpha, alpha_sample_sweep,CMYsweep_6ms)
        CMYsweep8 = interpolate(alpha, alpha_sample_sweep,CMYsweep_8ms)
        CMYsweep10 = interpolate(alpha, alpha_sample_sweep,CMYsweep_10ms)
        CMYsweep12 = interpolate(alpha, alpha_sample_sweep[0:5],CMYsweep_12ms)
        CMYsweep14 = interpolate(alpha, alpha_sample_sweep[0:3],CMYsweep_14ms)

        CMYsweep_speeds = [CMYsweep6, CMYsweep8, CMYsweep10, CMYsweep12, CMYsweep14]
        C_msweep = interpolate(V,speed_sample_sweep,CMYsweep_speeds)

        # Wing Tip Symetric Washout effictiveness coefficient estimation
        C_mwashout = 0.0
        
        if sweep < 5:
            if washout <= -10:
                C_mwashout = interpolate(alpha, alpha_sample_washout, CMYwashout_0sweep_m30tom10deflec)
            elif washout > -10:
                C_mwashout = interpolate(alpha, alpha_sample_washout[3:7], CMYwashout_0sweep_m10to5deflec)
        elif sweep >= 5:
            CMYwashout10 = interpolate(alpha, alpha_sample_washout, CMYwashout_10sweep)
            CMYwashout20 = interpolate(alpha, alpha_sample_washout, CMYwashout_20sweep)
            CMYwashout30 = interpolate(alpha, alpha_sample_washout, CMYwashout_30sweep)
            CMYwashout_sweeps = [CMYwashout10, CMYwashout20, CMYwashout30]

            C_mwashout = interpolate(sweep, sweep_sample_washout, CMYwashout_sweeps)

        # Dynamic Pitching Moment Coefficient estimation
        CMYq6 = interpolate(sweep,sweep_sample_dynamic,CMYq_6ms)
        CMYq8 = interpolate(sweep,sweep_sample_dynamic,CMYq_8ms)
        CMYq10 = interpolate(sweep,sweep_sample_dynamic,CMYq_10ms)
        CMYq_speeds = [CMYq6, CMYq8, CMYq10]
        C_mq = interpolate(V,speed_sample_dynamic,CMYq_speeds)

        return (C_m0, C_malpha, C_melev, C_msweep, C_mwashout, C_mq)

    def _get_coefficients_C_n(self, airstate, rates, input):
        # Rudder effectiveness coefficient estimation
        #C_n0 = 0.0020716
        C_n0 = 0.0 # From JAVA_MODEL

        # Rudder effectiveness Yawing Moment Coefficient
        alpha_sample_rudd = [-5, -2.5, 0, 2.5, 5]

        CMZrudd_20to10 = [0.0003325, 0.0002946, 0.0002567, 0.000308, 0.00036241]
        CMZrudd_10to0 = [0.0003961, 0.00038266, 0.00036924, 0.000399, 0.00042889]
        CMZrudd_0tom10 = [0.0003996, 0.0003719, 0.00034422, 0.0003965, 0.0004488]
        CMZrudd_m10tom20 = [0.00025022, 0.0002496, 0.00024895, 0.00031153, 0.0003741]

        alpha = airstate[0]
        rudder = input[3]

        C_nrudder_samples = []
        if rudder > 10:
            C_nrudder_samples = CMZrudd_20to10
        elif (rudder <= 10) and (rudder > 0):
            C_nrudder_samples = CMZrudd_10to0
        elif (rudder <= 0) and (rudder > -10):
            C_nrudder_samples = CMZrudd_0tom10
        elif rudder <= -10:
            C_nrudder_samples = CMZrudd_m10tom20

        C_nrudder = interpolate(alpha, alpha_sample_rudd, C_nrudder_samples)

        #C_nr = -0.05 # Eye average from Stanford student projects. Units? Probably per radian...
        C_nr = -0.002 # From JAVA_MODEL. Per degree

        return (C_n0, C_nrudder, C_nr)
        
position = [0,0,0]
velocity = [15,0,0]
attitude = [0,0,0,1]
rates = [0,0,0]

body = Body(mass,inertia,position,velocity,attitude,rates)
aerobody = AeroBody(body)
vehicle = AffectedBody(aerobody,[BixlerForces()])

def sensible(array,dp=4):
    elements = ", ".join([f"{v:.{dp}f}" for v in array])
    return f"[{elements}]"

if __name__ == "__main__":
    import sys
    deltaT = 0.01
    
    elevator = math.radians(-5)
    throttle = 0.0
    sweep = math.radians(20)
    washout = 0.0
    
    count = 0
    simtime = 0
    
    outfile = None
    if len(sys.argv) > 1:
        outfile = open(sys.argv[1],"w")
        outfile.write("time,x,y,z,u,v,w,qx,qy,qz,qw,p,q,r,alpha,elevator\n")
    
    while count < 2000:
        vehicle.step(deltaT,[0,elevator,throttle,0,sweep,washout])
        count += 1
        simtime += deltaT
        
        if outfile:
                outfile.write(f"{simtime},"+sensible(vehicle.statevector,10)[1:-1] + f",{vehicle.airstate[0]},{elevator}\n")
