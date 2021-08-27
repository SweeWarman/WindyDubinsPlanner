from Utils import *

def ComputeSturnSegment(R,altDropNeeded,deltaH_s,deltaH_t):

    diff = 1e3
    deltaL = 1
    while np.abs(diff) > 1e-1:
        if diff < 0:
            deltaL = deltaL*1.2
        elif diff > 0:
            deltaL = deltaL*0.8

        deltaX = 4*R*cos(0.5*(pi - deltaL/(2*R)))

        diff = deltaL*deltaH_t - deltaX*deltaH_s - altDropNeeded

    return deltaX

def FindDubinsParameters3D(pStart,pEnd,hStart,hEnd,R,case,altDiff,deltaH_s,deltaH_t):
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

    crs = pStart + R*dot(RotateZ(pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)]]))
    cls = pStart + R*dot(RotateZ(-pi/2),array([ [cos(hStart*pi/180)], [sin(hStart*pi/180)]]))
    cre = pEnd + R*dot(RotateZ(pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)]]))
    cle = pEnd + R*dot(RotateZ(-pi/2),array([ [cos(hEnd*pi/180)], [sin(hEnd*pi/180)]]))

    e1 = array([[1.0],[0.0]])

    z3 = pEnd
    q3 = dot(RotateZ(hEnd * pi/180),e1)
    paramsCSC = ()
    paramsCCC = ()
    paramsALT = ()
    success = False
    altitudeDiff = altDiff

    psi1 = 0
    psi2 = 0
    psi3 = 0

    cm1 = 0
    cm2 = 0

    val = distH(pStart,pEnd)
    if case == 1:
        val = distH(crs,cre)
        if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cre - crs)
            l = val
            L =  val + R*mod(2*pi + mod(vnu - pi/2,2*pi) - mod(hStart*pi/180 - pi/2,2*pi),2*pi) \
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


            deltah_1 = deltaH_t * diff1*R
            deltah_2 = deltaH_t * diff2*R
            deltah_3 = deltaH_s  * l
            deltah_T = deltah_1 + deltah_2 + deltah_3

            while(deltah_T > altitudeDiff):
                deltaH_s = deltaH_s*0.9
                deltah_1 = deltaH_t * diff1*R
                deltah_2 = deltaH_t * diff2*R
                deltah_3 = deltaH_s  * l
                deltah_T = deltah_1 + deltah_2 + deltah_3


            remainingAltDiff = altitudeDiff - deltah_T
            numTurns = int(np.floor(remainingAltDiff/(deltaH_t*2*pi*R)))

            excessAlt = (remainingAltDiff/(deltaH_t*2*pi*R) - numTurns)*(deltaH_t*2*pi*R)
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)]])*extendedFinal
            newL = distH(ce,cs)
            extension = newL - l
            q1 = (ce - cs)/norm(ce - cs)
            z1 = cs + R*dot(RotateZ(-pi/2),q1)
            z2 = ce + R*dot(RotateZ(-pi/2),q1)
            vnu1 = Angle(ce - cs)


            additionalAltDropNeeded = excessAlt/2 - extension*deltaH_s
            deltaX = ComputeSturnSegment(R,additionalAltDropNeeded,deltaH_s,deltaH_t)


            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns,deltaX)

    elif case == 2:
         val = distH(crs,cle)
         if (val > 2*R):
            l = val
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

            deltah_1 = deltaH_t* diff1*R
            deltah_2 = deltaH_t* diff2*R
            deltah_3 = deltaH_s  * l
            deltah_T = deltah_1 + deltah_2 + deltah_3


            while(deltah_T > altitudeDiff):
                deltaH_s = deltaH_s*0.9
                deltah_1 = deltaH_t * diff1*R
                deltah_2 = deltaH_t * diff2*R
                deltah_3 = deltaH_s  * l
                deltah_T = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltah_T
            nfac = deltaH_t*2*pi*R
            numTurns = int(np.floor(remainingAltDiff/nfac))

            excessAlt = (remainingAltDiff/nfac - numTurns)*nfac
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)]])*extendedFinal
            newl = distH(cs,ce)
            extension = sqrt(newl**2 - 4*R**2) - sqrt(l**2 - 4*R**2)

            vnu1 = Angle(ce - cs)
            vnu2 = vnu1 - pi/2 + arcsin(2*R/newl)
            q1 = dot(RotateZ(vnu2 + pi/2),e1)
            z1 = cs + R*dot(RotateZ(vnu2),e1)
            z2 = ce + R*dot(RotateZ(vnu2 + pi),e1)

            additionalAltDropNeeded = excessAlt/2 - extension*deltaH_s
            deltaX = ComputeSturnSegment(R,additionalAltDropNeeded,deltaH_s,deltaH_t)


            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns,deltaX)

    elif case == 3:
         val = distH(cls,cre)
         if (val > 2*R):
            l = val
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

            deltah_1 = deltaH_t * diff1 *R
            deltah_2 = deltaH_t * diff2 *R
            deltah_3 = deltaH_s  * l
            deltah_T = deltah_1 + deltah_2 + deltah_3

            while(deltah_T > altitudeDiff):
                deltaH_s = deltaH_s*0.9
                deltah_1 = deltaH_t * diff1*R
                deltah_2 = deltaH_t * diff2*R
                deltah_3 = deltaH_s  * l
                deltah_T = deltah_1 + deltah_2 + deltah_3


            remainingAltDiff = altitudeDiff - deltah_T
            nfac = deltaH_t*2*pi*R
            numTurns = int(np.floor(remainingAltDiff/nfac))

            excessAlt = (remainingAltDiff/nfac - numTurns)*nfac
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)]])*extendedFinal
            newl = distH(cs,ce)
            extension = sqrt(newl**2 - 4*R**2) - sqrt(l**2 - 4*R**2)
            vnu2 = arccos(2*R/newl)
            vnu1 = Angle(ce - cs)
            q1 = dot(RotateZ(vnu1 + vnu2 - pi/2),e1)
            z1 = cs + R*dot(RotateZ(vnu1 + vnu2),e1)
            z2 = ce + R*dot(RotateZ(vnu1 + vnu2 - pi),e1)

            additionalAltDropNeeded = excessAlt/2 - extension*deltaH_s
            deltaX = ComputeSturnSegment(R,additionalAltDropNeeded,deltaH_s,deltaH_t)

            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns,deltaX)

    elif case == 4:
         val = distH(cls,cle)
         if (val >= 3*R) or (val < 3*R):
            vnu = Angle(cle - cls)
            l = val
            L =  val + R*mod(2*pi + mod(vnu + pi/2,2*pi) + mod(hStart*pi/180 + pi/2,2*pi),2*pi) \
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

            deltah_1 = deltaH_t * diff1 * R
            deltah_2 = deltaH_t * diff2 * R
            deltah_3 = deltaH_s  * l
            deltah_T = deltah_1 + deltah_2 + deltah_3

            while(deltah_T > altitudeDiff):
                deltaH_s = deltaH_s*0.9
                deltah_1 = deltaH_t * diff1*R
                deltah_2 = deltaH_t * diff2*R
                deltah_3 = deltaH_s  * l
                deltah_T = deltah_1 + deltah_2 + deltah_3

            remainingAltDiff = altitudeDiff - deltah_T
            nfac = deltaH_t*2*pi*R
            numTurns = int(np.floor(remainingAltDiff/nfac))

            excessAlt = (remainingAltDiff/nfac - numTurns)*nfac
            extendedFinal = excessAlt/(2*deltaH_s)

            ce = ce + np.array([[cos(hEnd*np.pi/180 + pi)],[sin(hEnd*np.pi/180 + pi)]])*extendedFinal
            newL = distH(ce,cs)
            extension = newL - l
            q1 = (ce - cs)/norm(ce - cs)
            z1 = cs + R*dot(RotateZ(pi/2),q1)
            z2 = ce + R*dot(RotateZ(pi/2),q1)
            vnu1 = Angle(ce - cs)

            additionalAltDropNeeded = excessAlt/2 - extension*deltaH_s
            deltaX = ComputeSturnSegment(R,additionalAltDropNeeded,deltaH_s,deltaH_t)

            paramsCSC = (L,cs,lambda_s,ce,lambda_e,z1,q1,z2,z3,q3,vnu1,vnu2)
            paramsALT = (deltaH_s,deltaH_t,extendedFinal,numTurns,deltaX)

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
                cm1 = crs + 2*R*array([[cos(psi1)],[sin(psi1)]])
                cm2 = crs + 2*R*array([[cos(psi2)],[sin(psi2)]])
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
                cm1 = cls + 2*R*array([[cos(psi1)],[sin(psi1)]])
                cm2 = cls + 2*R*array([[cos(psi2)],[sin(psi2)]])
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

def GenerateDubinsAirPath3D(airspeed,ps,xs,pe,xe,R,startAlt,case,params,paramsALT):

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
    Z = [startAlt]

    e1 = array([[1.0],[0.0]])

    deltaH_s      = paramsALT[0]
    deltaH_t      = paramsALT[1]
    extendedFinal = paramsALT[2]
    numTurns      = paramsALT[3]
    deltaX        = paramsALT[4]

    turn_rate = airspeed/R

    l = distH(cs,ce)
    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        # First turn segment
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*R*deltaH_t)

        # Straight segment
        delta = (l - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)

        # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        angle = Angle(np.array([[X[-1] - X[-3]],[Y[-1]-Y[-3]]]))

        s1 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s2 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)


        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s3 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        # S-segment 3
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)


        # Second turn segment
        diff2 = AngleDiffCW(vnu-pi/2,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*R*deltaH_t)

        # spiral segment
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2 + delta*i,2*pi)),e1)
                X.append(pt[0,0])
                Y.append(pt[1,0])
                Z.append(Z[-1] - delta*R*deltaH_t)

        # extended final
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*R*deltaH_t)

        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = (L - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*deltaH_s)

        # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        angle = Angle(np.array([[X[-1] - X[-3]],[Y[-1]-Y[-3]]]))

        s1 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s2 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)


        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s3 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        # S-segment 3
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)

        # Second turn segment
        diff2 = AngleDiffCCW(vnu2 + pi,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu2 + pi - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)

        # spiral segment
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2 - delta*i,2*pi)),e1)
                X.append(pt[0,0])
                Y.append(pt[1,0])
                Z.append(Z[-1] - delta*R*deltaH_t)

        # extended final
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)


        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 3:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2, vnu + vnu2)

        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*R*deltaH_t)

        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = (L - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*deltaH_s)

        # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        angle = Angle(np.array([[X[-1] - X[-3]],[Y[-1]-Y[-3]]]))

        s1 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s2 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)


        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s3 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        # S-segment 3
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)

        # Second turn segment
        diff2 = AngleDiffCW(vnu + vnu2 - pi,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + vnu2 - pi + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)

        # spiral segment
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2 + delta*i,2*pi)),e1)
                X.append(pt[0,0])
                Y.append(pt[1,0])
                Z.append(Z[-1] - delta*R*deltaH_t)

        # extended final
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 4:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,vnu + pi/2)
        delta = diff1/100
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*R*deltaH_t)

        # Straight segment
        delta = (l-deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu  + pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append(Z[-1] - delta*deltaH_s)

         # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        angle = Angle(np.array([[X[-1] - X[-3]],[Y[-1]-Y[-3]]]))

        s1 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s2 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)

        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)
        s3 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        # S-segment 3
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)

        # Second turn segment
        diff2 = AngleDiffCCW(vnu + pi/2,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*R*deltaH_t)

        # spiral segment
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2 - delta*i,2*pi)),e1)
                X.append(pt[0,0])
                Y.append(pt[1,0])
                Z.append(Z[-1] - delta*R*deltaH_t)

        # extended final
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0])
            Y.append(pt[1,0])
            Z.append( Z[-1] - delta*deltaH_s)

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

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

def GenerateDubinsGroundPath3D(airspeed,windspeed,winddirection,ps,xs,pe,xe,R,startAlt,case,params,paramsALT):

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
    Z = [startAlt]

    e1 = array([[1.0],[0.0]])

    deltaH_s      = paramsALT[0]
    deltaH_t      = paramsALT[1]
    extendedFinal = paramsALT[2]
    numTurns      = paramsALT[3]
    deltaX        = paramsALT[4]

    turn_rate = airspeed/R

    l = distH(cs,ce)

    winddirection = winddirection*pi/180

    if case == 1:
        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu - pi/2)

        deltaT1 = 0
        # First turn segment
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT1)
            Z.append(Z[-1] - delta*R*deltaH_t)
            deltaT1 += delta/turn_rate

        deltaT2 = deltaT1

        # Straight segment
        delta = (l - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0] + windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0] + windspeed*sin(winddirection)*deltaT2)
            Z.append( Z[-1] - delta*deltaH_s)
            deltaT2 += delta/airspeed

        # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        pt1 = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + diff1,2*pi)),e1)
        pt2 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l-deltaX)

        angle = Angle(pt2 - pt1)

        deltaT3 = deltaT2
        s1 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT3)
            Z.append( Z[-1] - delta*R*deltaH_t)
            db1 =  delta/turn_rate
            deltaT3 += delta/turn_rate

        s2 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)

        deltaT4 = deltaT3
        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT4)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT4)
            Z.append( Z[-1] - delta*R*deltaH_t)
            db2 =  delta/turn_rate
            deltaT4 += delta/turn_rate

        s3 = cs + R*dot(RotateZ(mod(vnu - pi/2,2*pi)),e1) + q1*(l) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        deltaT5 = deltaT4
        # S-segment 3
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT5)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT5)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT5 +=  delta/turn_rate


        deltaT6 = deltaT5
        # Second turn segment
        diff2 = AngleDiffCW(vnu-pi/2,xe*pi/180 - pi/2)
        delta = diff2/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(vnu - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT6)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT6)
            Z.append(Z[-1] - delta*R*deltaH_t)
            deltaT6 +=  delta/turn_rate

        deltaT7 = deltaT6
        # spiral segment
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2 + delta*i,2*pi)),e1)
                X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT7)
                Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT7)
                Z.append(Z[-1] - delta*R*deltaH_t)
                deltaT7 += delta/turn_rate

        deltaT8 = deltaT7
        # extended final
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT8)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT8)
            Z.append( Z[-1] - delta*deltaH_s)
            deltaT8 += delta/airspeed

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 2:

        diff1 = AngleDiffCW(xs*pi/180 - pi/2,vnu2)

        deltaT1 = 0
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT1)
            Z.append(Z[-1] - delta*R*deltaH_t)
            deltaT1 += delta/airspeed

        # Straight segment
        deltaT2 = deltaT1
        L =  2*sqrt( l**2/4 - R**2)
        delta = (L - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT2)
            Z.append(Z[-1] - delta*deltaH_s)
            deltaT2 += delta/turn_rate

        # S-segment 1
        deltaT3 = deltaT2
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)


        pt1 = cs + R*dot(RotateZ(mod(xs*pi/180 - pi/2 + diff1,2*pi)),e1)
        pt2 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L-deltaX)
        angle = Angle(pt2 - pt1)

        s1 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT3)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT3 += delta/turn_rate

        s2 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)

        # S-segment 2
        alpha = (pi - 2*theta)
        delta = alpha/100
        deltaT4 = deltaT3
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT4)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT4)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT4 += delta/turn_rate

        s3 = cs + R*dot(RotateZ(mod(vnu2,2*pi)),e1) + q1*(L) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)

        # S-segment 3
        deltaT5 = deltaT4
        delta = gamma/100
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0]+ windspeed*cos(winddirection)*deltaT5)
            Y.append(pt[1,0]+ windspeed*sin(winddirection)*deltaT5)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT5 += delta/turn_rate

        # Second turn segment
        deltaT6 = deltaT5
        diff2 = AngleDiffCCW(vnu2 + pi,xe*pi/180 + pi/2)
        delta = diff2/100
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu2 + pi - delta*i,2*pi)),e1)
            X.append(pt[0,0] +windspeed*cos(winddirection)*deltaT6)
            Y.append(pt[1,0] +windspeed*sin(winddirection)*deltaT6)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT6 += delta/turn_rate

        # spiral segment
        deltaT7 = deltaT6
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2 - delta*i,2*pi)),e1)
                X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT7)
                Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT7)
                Z.append(Z[-1] - delta*R*deltaH_t)
                deltaT7 += delta/turn_rate

        # extended final
        deltaT8 = deltaT7
        delta = extendedFinal/100
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT8)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT8)
            Z.append( Z[-1] - delta*deltaH_s)
            deltaT8 += delta/turn_rate

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 3:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2, vnu + vnu2)

        deltaT1 = 0
        delta = diff1/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0] +windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0] +windspeed*sin(winddirection)*deltaT1)
            Z.append(Z[-1] - delta*R*deltaH_t)
            deltaT1 += delta/airspeed

        deltaT2 = deltaT1
        # Straight segment
        L =  2*sqrt( l**2/4 - R**2)
        delta = (L - deltaX)/100
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)
            Z.append(Z[-1] - delta*deltaH_s)

        # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)

        pt1 = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - diff1,2*pi)),e1)
        pt2 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L-deltaX)
        angle = Angle(pt2 - pt1)

        s1 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        deltaT3 = deltaT2
        delta = gamma/100
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT3 += delta/turn_rate

        # S-segment 2
        s2 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)
        deltaT4 = deltaT3
        alpha = (pi - 2*theta)
        delta = alpha/100
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT4)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT4)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT4 += delta/turn_rate

        # S-segment 3
        s3 = cs + R*dot(RotateZ(mod(vnu + vnu2,2*pi)),e1) + q1*(L) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)
        delta = gamma/100
        deltaT5 = deltaT4
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT5)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT5)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT5 += delta/turn_rate

        # Second turn segment
        diff2 = AngleDiffCW(vnu + vnu2 - pi,xe*pi/180 - pi/2)
        delta = diff2/100
        deltaT6 = deltaT5
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + vnu2 - pi + delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT6)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT6)
            Z.append( Z[-1] - delta*deltaH_s)

        # spiral segment
        deltaT7 = deltaT6
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2 + delta*i,2*pi)),e1)
                X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT7)
                Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT7)
                Z.append(Z[-1] - delta*R*deltaH_t)
                deltaT7 += delta/turn_rate

        # extended final
        delta = extendedFinal/100
        deltaT8 = deltaT7
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 - pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT8)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT8)
            Z.append( Z[-1] - delta*deltaH_s)
            deltaT8 += delta/airspeed

        timeTaken = diff1/turn_rate + L/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

    elif case == 4:

        diff1 = AngleDiffCCW(xs*pi/180 + pi/2,vnu + pi/2)
        delta = diff1/100
        deltaT1 = 0
        for i in range(101):
            pt = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT1)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT1)
            Z.append(Z[-1] - delta*R*deltaH_t)
            deltaT1 += delta/turn_rate

        # Straight segment
        delta = (l-deltaX)/100
        deltaT2 = deltaT1
        for i in range(100):
            pt = cs + R*dot(RotateZ(mod(vnu  + pi/2,2*pi)),e1) + q1*delta*i
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT2)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT2)
            Z.append(Z[-1] - delta*deltaH_s)
            deltaT2 += delta/airspeed

         # S-segment 1
        theta = arccos(deltaX/(4*R))
        gamma = (pi - 2*theta)/2
        h = 2*R*sin(theta)


        pt1 = cs + R*dot(RotateZ(mod(xs*pi/180 + pi/2 - diff1,2*pi)),e1)
        pt2 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l-deltaX)
        angle = Angle(pt2 - pt1)

        s1 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l - deltaX) +  (R)*dot(RotateZ(mod(angle-pi/2,2*pi)),e1)
        delta = gamma/100
        deltaT3 = deltaT2
        for i in range(100):
            pt = s1 +  R*dot(RotateZ(mod(angle + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT3)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT3)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT3 += delta/turn_rate

        # S-segment 2
        s2 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l - deltaX + deltaX/2) + (h-R)*dot(RotateZ(mod(angle + pi/2,2*pi)),e1)
        alpha = (pi - 2*theta)
        delta = alpha/100
        deltaT4 = deltaT3
        for i in range(100):
            pt = s2 +  R*dot(RotateZ(mod(angle - pi/2 - alpha/2 + delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT4)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT4)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT4 += delta/turn_rate

        # S-segment 3
        s3 = cs + R*dot(RotateZ(mod(vnu + pi/2,2*pi)),e1) + q1*(l) +  R*dot(RotateZ(mod(angle - pi/2,2*pi)),e1)
        delta = gamma/100
        deltaT5 = deltaT4
        for i in range(100):
            pt = s3 +  R*dot(RotateZ(mod(angle + pi/2 + gamma - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT5)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT5)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT5 += delta/turn_rate

        # Second turn segment
        diff2 = AngleDiffCCW(vnu + pi/2,xe*pi/180 + pi/2)
        delta = diff2/100
        deltaT6 = deltaT5
        for i in range(101):
            pt = ce + R*dot(RotateZ(mod(vnu + pi/2 - delta*i,2*pi)),e1)
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT6)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT6)
            Z.append( Z[-1] - delta*R*deltaH_t)
            deltaT6 += delta/turn_rate

        # spiral segment
        deltaT7 = deltaT6
        if numTurns > 0:
            diff3 = numTurns*2*pi
            delta = diff3/(100*numTurns)
            for i in range(100*numTurns+ 1):
                pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2 - delta*i,2*pi)),e1)
                X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT7)
                Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT7)
                Z.append(Z[-1] - delta*R*deltaH_t)
                deltaT7 += delta/turn_rate

        # extended final
        delta = extendedFinal/100
        deltaT8 = deltaT7
        for i in range(100):
            pt = ce + R*dot(RotateZ(mod(xe*pi/180 + pi/2,2*pi)),e1) + q3*delta*i
            X.append(pt[0,0]+windspeed*cos(winddirection)*deltaT8)
            Y.append(pt[1,0]+windspeed*sin(winddirection)*deltaT8)
            Z.append( Z[-1] - delta*deltaH_s)
            deltaT8 += delta/airspeed

        timeTaken = diff1/turn_rate + l/airspeed + diff2/turn_rate + (2*gamma + alpha)/turn_rate + extendedFinal/airspeed

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



def PlotPath3D(SP,EP,hStart,hEnd,turnR,startAlt,case,status,paramsCSC,paramsCCC,paramsALT,showGroundPath=False,airspeed=2,windspeed=0,winddirection=0):
    fig3 = plt.figure(3)
    ax = fig3.gca(projection='3d')

    if(case <= 4):
        params = paramsCSC
    else:
        params = paramsCCC

    if status:
        (X,Y,Z,timeTaken) = GenerateDubinsAirPath3D(airspeed,SP,hStart,EP,hEnd,turnR,startAlt,case,params,paramsALT)
        (Xg,Yg,Zg,timeTaken) = GenerateDubinsGroundPath3D(airspeed,windspeed,winddirection,SP,hStart,EP,hEnd,turnR,startAlt,case,params,paramsALT)
    else:
        (X,Y,Z) = (SP[0,0],SP[1,0],startAlt)
        (Xg,Yg,Zg) = (X,Y,0)

    ax.plot(X,Y,Z)
    ax.plot(Xg,Yg,Zg,'g--')
