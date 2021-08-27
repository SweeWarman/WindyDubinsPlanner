import numpy as np
from numpy import sin,cos,arctan2,\
                  arcsin,arccos,mod,pi,\
                  array,sqrt,dot
from numpy.linalg import norm

def dist(A,B):
    return np.sqrt(np.dot((A-B).T,(A-B)))[0,0]

def Angle(X):
    return arctan2(X[1,0],X[0,0])

def RotateZ(theta):
    return np.array([[cos(theta),-sin(theta),0],
                     [sin(theta), cos(theta),0],
                     [0,                   0,1]])

def FindDubinsParameters(pStart,pEnd,hStart,hEnd,R,case,deltaH_s,deltaH_t):
    """
    :param pStart: start position - 3x1 array (x,y,z)
    :param pEnd: end position - 3x1 array
    :param hStart: initial heading (with respect to true north)
    :param hEnd: final heading (with respect to true north)
    :param R: turn radius
    :param case: case 1-4
    :param deltaH_s: altitude loss per unit distance (m)
    :param deltaH_t: altitude loss per turn
    :return: (case,success,(L,cs,lambda_s,ce,lambda_e,z1,q1,z3,z3,q3),(ce,cs,lambda_s,lambda_e,cm1,cm2,psi1,psi2,psi3),(deltaX,finalTurns))

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

    # Find gliding angle for given airspeed (assuming you have a mapping between airspeed and glide angle).


    # Find bank angle/turn rate for given turn rate/bank angle

    crs = pStart + R*dot(RotateZ(pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)], [0]]))
    cls = pStart + R*dot(RotateZ(-pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)], [0]]))
    cre = pEnd + R*dot(RotateZ(pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)], [0]]))
    cle = pEnd + R*dot(RotateZ(-pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)], [0]]))

    e1 = array([[1.0],[0.0],[0.0]])

    z3 = pEnd
    q3 = dot(RotateZ(hEnd * pi/180),e1)
    paramsCSC = ()
    paramsCCC = ()
    paramsALT = ()
    success = False
    altitudeDiff = pStart[2,0] - pEnd[2,0]

    val = dist(pStart,pEnd)
    if case == 1:
        val = dist(crs,cre)
        if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cre - crs)
            l = norm(crs - cre)
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

            # Altitude loss
            diff1 = AngleDiffCW(hStart*pi/180 - pi/2,vnu - pi/2)
            diff2 = AngleDiffCW(vnu-pi/2,hEnd*pi/180 - pi/2)

            deltah_1 = deltaH_t/2*pi * diff1
            deltah_2 = deltaH_t/2*pi * diff2
            deltah_3 = deltaH_s  * l
            deltah_t = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltaH_t
            numTurns = int(np.floor(remainingAltDiff/deltaH_t))

            excessAlt = (remainingAltDiff/deltaH_t - numTurns)*deltaH_t
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            z2 = z2 + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns)



    elif case == 2:
         val = dist(crs,cle)
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

            diff1 = AngleDiffCW(hStart*pi/180 - pi/2,vnu2)
            diff2 = AngleDiffCCW(vnu2 + pi,hEnd*pi/180 + pi/2)

            deltah_1 = deltaH_t/2*pi * diff1
            deltah_2 = deltaH_t/2*pi * diff2
            deltah_3 = deltaH_s  * l
            deltah_t = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltaH_t
            numTurns = int(np.floor(remainingAltDiff/deltaH_t))

            excessAlt = (remainingAltDiff/deltaH_t - numTurns)*deltaH_t
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            z2 = z2 + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns)

    elif case == 3:
         val = dist(cls,cre)
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

            diff1 = AngleDiffCCW(hStart*pi/180 + pi/2, vnu + vnu2)
            diff2 = AngleDiffCW(vnu + vnu2 - pi,hEnd*pi/180 - pi/2)

            deltah_1 = deltaH_t/2*pi * diff1
            deltah_2 = deltaH_t/2*pi * diff2
            deltah_3 = deltaH_s  * l
            deltah_t = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltaH_t
            numTurns = int(np.floor(remainingAltDiff/deltaH_t))

            excessAlt = (remainingAltDiff/deltaH_t - numTurns)*deltaH_t
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            z2 = z2 + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns)

    elif case == 4:
         val = dist(cls,cle)
         if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cle - cls)
            l = norm(cls - cle)
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

            diff1 = AngleDiffCCW(hStart*pi/180 + pi/2,vnu + pi/2)
            diff2 = AngleDiffCCW(vnu + pi/2,hEnd*pi/180 + pi/2)

            deltah_1 = deltaH_t/2*pi * diff1
            deltah_2 = deltaH_t/2*pi * diff2
            deltah_3 = deltaH_s  * l
            deltah_t = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltaH_t
            numTurns = int(np.floor(remainingAltDiff/deltaH_t))

            excessAlt = (remainingAltDiff/deltaH_t - numTurns)*deltaH_t
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            z2 = z2 + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)],[0.0]])*extendedFinal
            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns)


    elif case == 5 or case == 6:
        val = dist(crs,cre)
        if(val <= 4*R):
            VecSE = cre - crs
            dse = dist(cre,crs)
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

        if case == 5:
            diff1 = AngleDiffCW(hStart*pi/180 - pi/2,psi1)
            diff2 = psi3
            angle_em = Angle(cm1 - cre)
            diff3 = AngleDiffCW(angle_em,hEnd*pi/180 - pi/2)

        else:
            diff1 = AngleDiffCW(hStart*pi/180 - pi/2,psi2)
            diff2 = mod(2*pi - psi3,2*pi)
            angle_em = Angle(cm2 - cre)
            diff3 = AngleDiffCW(angle_em,hEnd*pi/180 - pi/2)

        deltah_1 = deltaH_t/2*pi * diff1
        deltah_2 = deltaH_t/2*pi * diff2
        deltah_3 = deltaH_t/2*pi * diff3
        deltah_t = deltah_1 + deltah_2 + deltah_3

    elif case == 7 or case == 8:
        val = dist(cls,cle)
        if(val <= 4*R):
            VecSE = cle - cls
            dse = dist(cle,cls)
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

        if case == 7:
            diff1 = AngleDiffCCW(hStart*pi/180 + pi/2,psi1)
            diff2 = mod(2*pi - psi3,2*pi)
            angle_em = Angle(cm1 - cle)
            diff3 = AngleDiffCCW(angle_em,hEnd*pi/180 + pi/2)
        else:
            diff1 = AngleDiffCCW(hStart*pi/180 + pi/2,psi2)
            diff2 = psi3
            angle_em = Angle(cm2 - cle)
            diff3 = AngleDiffCCW(angle_em,hEnd*pi/180 + pi/2)

        deltah_1 = deltaH_t/2*pi * diff1
        deltah_2 = deltaH_t/2*pi * diff2
        deltah_3 = deltaH_t/2*pi * diff3
        deltah_t = deltah_1 + deltah_2 + deltah_3

    return (case,success,paramsCSC,paramsCCC,paramsALT)

def AngleDiffCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi + angle2 - angle1,2*pi)

def AngleDiffCCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi - angle2 + angle1,2*pi)

def GenerateDubinsAirPath(airspeed,ps,xs,pe,xe,R,case,params,paramsALT):

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
    Z = [ps[2,0]]

    e1 = array([[1.0],[0.0],[0.0]])

    deltaH_s      = paramsALT[0]
    deltaH_t      = paramsALT[1]
    extendedFinal = paramsALT[2]
    numTurns      = paramsALT[3]

    turn_rate = airspeed/R

    l = dist(cs,ce)
    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        # First turn segment
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            pt[2,0] = pt[2,0] - delta*R*deltaH_t
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(pt[2,0])

        # Straight segment
        delta = (l + extendedFinal)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*delta*i
            pt[2,0] = Z[-1] - delta*deltaH_s
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(pt[2,0])

        # Second turn segment
        diff2 = AngleDiffCW(vnu-pi/2,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            pt[2,0] = Z[-1] - delta*R*deltaH_t
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(pt[2,0])

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])

        # Straight segment
        L =  2*sqrt( l**2/4 - R**2) + extendedFinal
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
        L =  2*sqrt( l**2/4 - R**2) + extendedFinal
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
        delta = (l + extendedFinal)/100
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

    return (X,Y,Z,timeTaken)

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