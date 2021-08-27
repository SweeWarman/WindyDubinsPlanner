import numpy as np
from numpy import mod,sin,cos,pi,array,dot,sqrt
from DubinsParameters import RotateZ,dist,Angle

def AngleDiffCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi + angle2 - angle1,2*pi)

def AngleDiffCCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi - angle2 + angle1,2*pi)

def GenerateDubinsAirPath(airspeed,ps,xs,pe,xe,R,case,params):

    if (case <= 4):
        L        = params[0]
        cs       = params[1]
        lambda_s = params[2]
        ce       = params[3]
        lambda_e = params[4]
        z1       = params[5]
        q1       = params[6]
        z2       = params[7]
        z3       = params[8]
        q3       = params[9]
        vnu      = params[10]
        vnu2     = params[11]
    else:
        cs = params[0]
        ce = params[1]
        lambda_s = params[2]
        lambda_e = params[3]
        cm1 = params[4]
        cm2 = params[5]
        psi1 = params[6]
        psi2 = params[7]
        psi3 = params[8]

    X = []
    Y = []

    e1 = array([[1.0],[0.0],[0.0]])

    turn_rate = airspeed/R

    l = dist(cs,ce)
    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        # First turn segment
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Straight segment
        delta = l/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        diff2 = AngleDiffCW(vnu-pi/2,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = L/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        diff2 = AngleDiffCCW(vnu2 + pi,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu2 + pi - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate

    elif case == 3:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2, vnu + vnu2)

        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = L/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        diff2 = AngleDiffCW(vnu + vnu2 - pi,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + vnu2 - pi + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate

    elif case == 4:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,vnu + pi/2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Straight segment
        delta = l/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu  + pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        diff2 = AngleDiffCCW(vnu + pi/2,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate

    elif case == 5:
        # First turn segment
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,psi1)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(2*pi + xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        angle_ms = Angle(cs - cm1)
        diff2 = psi3
        delta  = diff2/100
        for i in range(101):
            pt = cm1 + R*dot(RotateZ(mod(2*pi + angle_ms - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        angle_em = Angle(cm1 - ce)
        diff3 = AngleDiffCW(angle_em,xe*pi/180 - pi/2)
        delta = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(2*pi + angle_em + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])


        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 6:
        # First turn segment
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,psi2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(2*pi + xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        angle_ms = Angle(cs - cm2)
        diff2 = mod(2*pi - psi3,2*pi)
        delta  = diff2/100
        for i in range(101):
            pt = cm2 + R*dot(RotateZ(mod(2*pi + angle_ms - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        angle_em = Angle(cm2 - ce)
        diff3 = AngleDiffCW(angle_em,xe*pi/180 - pi/2)
        delta = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(2*pi + angle_em + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 7:
        # First turn segment
        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,psi1)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        angle_ms = Angle(cs - cm1)
        diff2 = mod(2*pi - psi3,2*pi)
        delta  = diff2/100
        for i in range(101):
            pt = cm1 + R*dot(RotateZ(mod(angle_ms + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        angle_em = Angle(cm1 - ce)
        diff3 = AngleDiffCCW(angle_em,xe*pi/180 + pi/2)
        delta  = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(angle_em - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 8:
        # First turn segment
        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,psi2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Second turn segment
        angle_ms = Angle(cs - cm2)
        diff2 = psi3
        delta  = diff2/100
        for i in range(101):
            pt = cm2 + R*dot(RotateZ(mod(angle_ms + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        angle_em = Angle(cm2 - ce)
        diff3 = AngleDiffCCW(angle_em,xe*pi/180 + pi/2)
        delta  = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(angle_em - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    return (X,Y,timeTaken)

def GenerateDubinsGroundPath(turn_rate,airspeed,windspeed,winddirection,ps,xs,pe,xe,R,case,params):

    if (case <= 4):
        L        = params[0]
        cs       = params[1]
        lambda_s = params[2]
        ce       = params[3]
        lambda_e = params[4]
        z1       = params[5]
        q1       = params[6]
        z2       = params[7]
        z3       = params[8]
        q3       = params[9]
        vnu      = params[10]
        vnu2     = params[11]
    else:
        cs = params[0]
        ce = params[1]
        lambda_s = params[2]
        lambda_e = params[3]
        cm1 = params[4]
        cm2 = params[5]
        psi1 = params[6]
        psi2 = params[7]
        psi3 = params[8]

    X = []
    Y = []

    e1 = array([[1.0],[0.0],[0.0]])
    winddirection = winddirection*pi/180
    #turn_rate = turn_rate*pi/180
    l = dist(cs,ce)
    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        deltaT1 = 0
        # First turn segment
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0] + windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0] + windspeed*sin(winddirection)*deltaT1)
            deltaT1 += delta/turn_rate


        deltaT2 = 0
        # Straight segment
        delta = l/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0] + windspeed*cos(winddirection)*(deltaT1 + deltaT2))
            Y.append(pt[1,0] + windspeed*sin(winddirection)*(deltaT1 + deltaT2))
            #print windspeed*sin(winddirection)*(deltaT1 + deltaT2)
            deltaT2 += delta/airspeed



        deltaT3 = 0
        # Second turn segment
        diff2 = AngleDiffCW(vnu-pi/2,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0] + windspeed*cos(winddirection)*(deltaT1 + deltaT2 + deltaT3))
            Y.append(pt[1,0] + windspeed*sin(winddirection)*(deltaT1 + deltaT2 + deltaT3))
            deltaT3 +=  delta/turn_rate

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        deltaT1 = 0
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0] + windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0] + windspeed*sin(winddirection)*deltaT1)
            deltaT1 += delta/turn_rate

        deltaT2 = 0
        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = L/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0] + windspeed*cos(winddirection)*(deltaT1 + deltaT2))
            Y.append(pt[1,0] + windspeed*sin(winddirection)*(deltaT1 + deltaT2))
            deltaT2 +=  delta/airspeed

        deltaT3 = 0
        # Second turn segment
        diff2 = AngleDiffCCW(vnu2 + pi,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu2 + pi - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*(deltaT1 + deltaT2 + deltaT3))
            Y.append(pt[1,0]+windspeed*sin(winddirection)*(deltaT1 + deltaT2 + deltaT3))
            deltaT3 += delta/turn_rate

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate

    elif case == 3:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2, vnu + vnu2)

        deltaT1 = 0
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)

        deltaT2 = 0
        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = L/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*delta*i
            deltaT2 = deltaT1 + l/airspeed*i*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)

        deltaT3 = 0
        # Second turn segment
        diff2 = AngleDiffCW(vnu + vnu2 - pi,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + vnu2 - pi + delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate

    elif case == 4:

        deltaT1 = 0
        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,vnu + pi/2)
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)

        deltaT2 = 0
        # Straight segment
        delta = l/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu  + pi/2,2*pi)),e1) + q1*delta*i
            deltaT2 = deltaT1 + l/airspeed*i*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)

        deltaT3 = 0
        # Second turn segment
        diff2 = AngleDiffCCW(vnu + pi/2,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + pi/2 - delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate

    elif case == 5:
        # First turn segment
        deltaT1 = 0
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,psi1)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(2*pi + xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)

        deltaT2 = 0
        # Second turn segment
        angle_ms = Angle(cs - cm1)
        diff2 = psi3
        delta  = diff2/100
        for i in range(101):
            pt = cm1 + R*dot(RotateZ(mod(2*pi + angle_ms - delta*i,2*pi)),e1)
            deltaT2 = deltaT1 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)



        deltaT3 = 0
        angle_em = Angle(cm1 - ce)
        diff3 = AngleDiffCW(angle_em,xe*pi/180 - pi/2)
        delta = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(2*pi + angle_em + delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff3/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 6:
        deltaT1 = 0
        # First turn segment
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,psi2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(2*pi + xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)


        deltaT2 = 0
        # Second turn segment
        angle_ms = Angle(cs - cm2)
        diff2 = mod(2*pi - psi3,2*pi)
        delta  = diff2/100
        for i in range(101):
            pt = cm2 + R*dot(RotateZ(mod(2*pi + angle_ms - delta*i,2*pi)),e1)
            deltaT2 = deltaT1 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)

        deltaT3 = 0
        angle_em = Angle(cm2 - ce)
        diff3 = AngleDiffCW(angle_em,xe*pi/180 - pi/2)
        delta = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(2*pi + angle_em + delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff3/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 7:

        deltaT1 = 0
        # First turn segment
        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,psi1)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)


        deltaT2 = 0
        # Second turn segment
        angle_ms = Angle(cs - cm1)
        diff2 = mod(2*pi - psi3,2*pi)
        delta  = diff2/100
        for i in range(101):
            pt = cm1 + R*dot(RotateZ(mod(angle_ms + delta*i,2*pi)),e1)
            deltaT2 = deltaT1 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)

        deltaT3 = 0
        angle_em = Angle(cm1 - ce)
        diff3 = AngleDiffCCW(angle_em,xe*pi/180 + pi/2)
        delta  = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(angle_em - delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff3/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    elif case == 8:

        deltaT1 = 0
        # First turn segment
        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,psi2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            deltaT1 = i*diff1/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)

        deltaT2 = 0
        # Second turn segment
        angle_ms = Angle(cs - cm2)
        diff2 = psi3
        delta  = diff2/100
        for i in range(101):
            pt = cm2 + R*dot(RotateZ(mod(angle_ms + delta*i,2*pi)),e1)
            deltaT2 = deltaT1 + i*diff2/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)

        deltaT3 = 0
        angle_em = Angle(cm2 - ce)
        diff3 = AngleDiffCCW(angle_em,xe*pi/180 + pi/2)
        delta  = diff3/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(angle_em - delta*i,2*pi)),e1)
            deltaT3 = deltaT2 + i*diff3/turn_rate*0.01
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)

        timeTaken = diff1/turn_rate + diff2/turn_rate + diff3/turn_rate

    return (X,Y,timeTaken)
