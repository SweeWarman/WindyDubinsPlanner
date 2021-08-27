from Utils import *

def FindDubinsParameters2D(pStart,pEnd,hStart,hEnd,R,case):
    """
    :param pStart: start position - 3x1 array (x,y,z)
    :param pEnd: end position - 3x1 array
    :param hStart: initial heading (with respect to true north)
    :param hEnd: final heading (with respect to true north)
    :param R: turn radius
    :param case: case 1-4
    :return: (case,success,(L,cs,lambda_s,ce,lambda_e,z1,q1,z3,z3,q3),(ce,cs,lambda_s,lambda_e,cm1,cm2,psi1,psi2,psi3))

    L = length of straight segment
    cs = center of starting turn
    lambda_s = direction of starting turn (+ve clockwise)
    ce = center of ending turn
    lambda_e = direction of ending turn (+ve clockwise)
    z1,q1 = Hyperplan parameters
    z2 = z1
    z3,q3 = Hyperplan parameters
    cm1 = center of middle circle for CCC1 trajectories
    cm2 = center of middle circle for CCC2 trajectories
    psi1 = switching heading from the start to the middle circle in CCC1
    psi2 = switching heading from the start to the middle circle in CCC2
    """

    crs = pStart + R*dot(RotateZ(pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)]]))
    cls = pStart + R*dot(RotateZ(-pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)]]))
    cre = pEnd + R*dot(RotateZ(pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)]]))
    cle = pEnd + R*dot(RotateZ(-pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)]]))

    e1 = array([[1.0],[0.0]])

    z3 = pEnd
    q3 = dot(RotateZ(hEnd * pi/180),e1)
    paramsCSC = ()
    paramsCCC = ()
    success = False

    val = distH(pStart,pEnd)
    if case == 1:
        val = distH(crs,cre)
        if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cre - crs)
            L =  norm(crs - cre) + R*mod(2*pi + mod(vnu - pi/2,2*pi) - mod(hStart*pi/180 - pi/2,2*pi),2*pi) \
                  + R*mod(2*pi + mod(hEnd*pi/180 - pi/2,2*pi) - mod(vnu - pi/2,360),2*pi)
            cs = crs
            ce = cre
            lambda_s = 1
            lambda_e = 1
            q1 = (ce - cs)/norm(ce - cs)
            z1 = cs + R*dot(RotateZ(-pi/2),q1)
            z2 = ce + R*dot(RotateZ(-pi/2),q1)
            vnu1 = vnu
            vnu2 = 0
            success = True
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
    elif case == 2:
         val = distH(crs,cle)
         if (val > 2*R):
            l = norm(cle - crs)
            vnu = Angle(cle - crs)
            vnu2 = vnu - pi/2 + arcsin(2*R/l)
            L = sqrt(l**2 - 4*R**2) + R*mod(2*pi + mod(vnu2,2*pi) - mod(hStart*pi/180- pi/2,2*pi),2*pi) \
                 + R*mod(2*pi + mod(vnu2 + pi,2*pi) - mod(hEnd + pi/2,2*pi),2*pi)
            cs = crs
            ce = cle
            lambda_s = 1
            lambda_e = -1
            q1 = dot(RotateZ(vnu2 + pi/2),e1)
            z1 = cs + R*dot(RotateZ(vnu2),e1)
            z2 = ce + R*dot(RotateZ(vnu2 + pi),e1)
            vnu1 = vnu
            vnu2 = vnu2
            success = True
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
    elif case == 3:
         val = distH(cls,cre)
         if (val > 2*R):
            l = norm(cre - cls)
            vnu = Angle(cre - cls)
            vnu2 = arccos(2*R/l)
            L = sqrt(l**2 - 4*R**2) +R*mod(2*pi + mod(hStart*pi/180 + pi/2,2*pi) - mod(vnu + vnu2,2*pi),2*pi)  \
                 +R*mod(2*pi + mod(hEnd*pi/180 - pi/2,2*pi) - mod(vnu + vnu2 - pi,2*pi),2*pi)
            cs = cls
            ce = cre
            lambda_s = -1
            lambda_e = 1
            q1 = dot(RotateZ(vnu + vnu2 - pi/2),e1)
            z1 = cs + R*dot(RotateZ(vnu + vnu2),e1)
            z2 = ce + R*dot(RotateZ(vnu + vnu2 - pi),e1)
            vnu1 = vnu
            vnu2 = vnu2
            success = True
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
    elif case == 4:
         val = distH(cls,cle)
         if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cle - cls)
            L =  norm(cls - cle) + R*mod(2*pi + mod(vnu + pi/2,2*pi) + mod(hStart*pi/180 + pi/2,2*pi),2*pi) \
                  + R*mod(2*pi - mod(hEnd*pi/180 + pi/2,2*pi) + mod(vnu + pi/2,360),2*pi)
            cs = cls
            ce = cle
            lambda_s = -1
            lambda_e = -1
            q1 = (ce - cs)/norm(ce - cs)
            z1 = cs + R*dot(RotateZ(pi/2),q1)
            z2 = ce + R*dot(RotateZ(pi/2),q1)
            vnu1 = vnu
            vnu2 = 0
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            success = True
    elif case == 5 or case == 6:
        val = distH(crs,cre)
        if(val <= 4*R):
            VecSE = cre - crs
            dse = distH(cre,crs)
            if (dse <= 4*R):
                angleSE = Angle(VecSE)
                angleTheta = mod(2*pi + arccos(dse/(4*R)),2*pi)
                psi1 = mod(2*pi + angleSE - angleTheta,2*pi)
                psi2 = mod(2*pi + angleSE + angleTheta,2*pi)
                psi3 = mod(2*pi + pi - (2*angleTheta),2*pi)
                cm1 = crs + 2*R*array([[cos(psi1)],[sin(psi1)],[0.0]])
                cm2 = crs + 2*R*array([[cos(psi2)],[sin(psi2)],[0.0]])
                lambda_s = 1
                lambda_e = 1
                paramsCCC = (crs,cre,lambda_s,lambda_e,cm1,cm2,psi1,psi2,psi3)
                success = True
    elif case == 7 or case == 8:
        val = distH(cls,cle)
        if(val <= 4*R):
            VecSE = cle - cls
            dse = distH(cle,cls)
            if (dse <= 4*R):
                angleSE = Angle(VecSE)
                angleTheta = mod(2*pi + arccos(dse/(4*R)),2*pi)
                psi1 = mod(2*pi + angleSE - angleTheta,2*pi)
                psi2 = mod(2*pi + angleSE + angleTheta,2*pi)
                psi3 = mod(2*pi + pi - (2*angleTheta),2*pi)
                cm1 = cls + 2*R*array([[cos(psi1)],[sin(psi1)],[0.0]])
                cm2 = cls + 2*R*array([[cos(psi2)],[sin(psi2)],[0.0]])
                lambda_s = -1
                lambda_e = -1
                paramsCCC = (cls,cle,lambda_s,lambda_e,cm1,cm2,psi1,psi2,psi3)
                success = True

    return (case,success,paramsCSC,paramsCCC)

def GenerateDubinsAirPath2D(ps,xs,pe,xe,R,case,params):

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

    X = [ps[0,0]]
    Y = [ps[1,0]]

    e1 = array([[1.0],[0.0]])

    l = distH(cs,ce)
    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        # First turn segment
        delta = diff1/100
        for i in range(100):
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
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        length =  (diff1 + diff2,l)

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        delta = diff1/100
        for i in range(100):
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

        length =  (diff1 + diff2,l)

    elif case == 3:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2, vnu + vnu2)

        delta = diff1/100
        for i in range(100):
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

        length =  (diff1 + diff2,l)

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

        length =  (diff1 + diff2,l)

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


        length = (diff1 + diff2 + diff3,0)

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

        length = (diff1 + diff2 + diff3,0)

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

        length = (diff1 + diff2 + diff3,0)

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

        length = (diff1 + diff2 + diff3,0)

    return (X,Y,length)

def GenerateDubinsGroundPath2D(airspeed,windspeed,winddirection,ps,xs,pe,xe,R,case,params):

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

    e1 = array([[1.0],[0.0]])
    winddirection = winddirection*pi/180
    turn_rate = airspeed/R

    l = distH(cs,ce)
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

        length = (diff1 + diff2, l)

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

        length = (diff1 + diff2, l)

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

        length = (diff1 + diff2, l)

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

        length = (diff1 + diff2, l)

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

        length = (diff1 + diff2 + diff3, 0)

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

        length = (diff1 + diff2 + diff3, 0)

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

        length = (diff1 + diff2 + diff3, 0)

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

        length = (diff1 + diff2 + diff3, 0)

    return (X,Y,length)

def PlotPath2D(SP,EP,hStart,hEnd,turnR,case,status,paramsCSC,paramsCCC,showGroundPath=False,airspeed=2,windspeed=0,winddirection=0):

    # define the plot limits
    xmin = min([SP[0, 0], EP[0, 0]]) - (2 * turnR + 4)
    xmax = max([SP[0, 0], EP[0, 0]]) + (2 * turnR + 4)
    ymin = min([SP[1, 0], EP[1, 0]]) - (2 * turnR + 4)
    ymax = max([SP[1, 0], EP[1, 0]]) + (2 * turnR + 4)

    # define parameters to indicate initial and final orientation vectors
    arrow_start_delta = 1.0 * array([[cos(hStart * pi / 180)], [sin(hStart * pi / 180)]])
    arrow_end_delta = 1.0 * array([[cos(hEnd * pi / 180)], [sin(hEnd * pi / 180)]])

    if (case <= 4):
        plt.figure(1)
        plt.subplot(2, 2, case)
        params = paramsCSC
    else:
        plt.figure(2)
        plt.subplot(2, 2, case % 5 + 1)
        params = paramsCCC

    if case == 1:
        type = 'RSR'
    elif case == 2:
        type = 'RSL'
    elif case == 3:
        type = 'LSR'
    elif case == 4:
        type = 'LSL'
    elif case == 5:
        type = 'RLR1'
    elif case == 6:
        type = 'RLR2'
    elif case == 7:
        type = 'LRL2'
    elif case == 8:
        type = 'LRL1'

    if status:
        (X, Y, length1) = GenerateDubinsAirPath2D(SP, hStart, EP, hEnd, turnR, case, params)
        if showGroundPath:
            (Xg,Yg,length2) = GenerateDubinsGroundPath2D(airspeed,windspeed,winddirection,SP,hStart,EP,hEnd,turnR,case,params)

    else:
        (X, Y) = (SP[0, 0], EP[1, 0])
        if showGroundPath:
            (Xg, Yg) = (X,Y)

    plt.plot(SP[1, 0], SP[0, 0], 'ro')
    plt.plot(EP[1, 0], EP[0, 0], 'ro')
    plt.title(type)
    plt.arrow(SP[1, 0], SP[0, 0], arrow_start_delta[1, 0], arrow_start_delta[0, 0], head_width=0.5)
    plt.arrow(EP[1, 0], EP[0, 0], arrow_end_delta[1, 0], arrow_end_delta[0, 0], head_width=0.5)
    plt.ylim([xmin, xmax])
    plt.xlim([ymin, ymax])
    plt.plot(Y, X, 'b')
    if showGroundPath:
        plt.plot(Yg, Xg, 'g')



