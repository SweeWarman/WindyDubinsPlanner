from DubinsPath2D import *


start_pos = array([[2.0],[5.0]])
end_pos   = array([[13.0],[6.0]])

R = 2
hStart = 50
hEnd = 90
turn_rate = 1
airspeed = 2
windspeed = 0.4
winddirection = 135
success = False
end_pos_virtual = []

alpha = 1
end_pos_virtual = end_pos
error = 1e3
input_case = 2

(case,status,paramsCSC,paramsCCC) = FindDubinsParameters2D(start_pos,end_pos,hStart,hEnd,R,input_case)

if status and not success:
    if(case <= 4):
        params = paramsCSC
    else:
        params = paramsCCC

    (X,Y,timeTaken) = GenerateDubinsAirPath2D(airspeed,start_pos,hStart,end_pos,hEnd,R,case,params)
    (Xg,Yg,timeTaken) = GenerateDubinsGroundPath2D(turn_rate,airspeed,windspeed,winddirection,start_pos,hStart,end_pos,hEnd,R,case,params)

    shift =  np.array([[-windspeed*cos(winddirection*pi/180)*timeTaken],[-windspeed*sin(winddirection*pi/180)*timeTaken]])

    shiftDist = dist(np.zeros((2,1)),shift)

    while shiftDist > 1e-1:

        # Shifted destination
        end_pos_virtual = end_pos_virtual + shift
        (case,status,paramsCSC,paramsCCC) = FindDubinsParameters2D(start_pos,end_pos_virtual,hStart,hEnd,R,case)

        if status is True:
            if(case <= 4):
                params = paramsCSC
            else:
                params = paramsCCC

            (Xg,Yg,timeTaken) = GenerateDubinsGroundPath2D(turn_rate,airspeed,windspeed,winddirection,start_pos,hStart,end_pos_virtual,hEnd,R,case,params)
            shift = np.array([[end_pos[0,0] - Xg[-1]],[end_pos[1,0] - Yg[-1]]])
            shiftDist = dist(np.zeros((2,1)),shift)
        else:
            break

    if shiftDist < 1e-1:
        success = True
        print end_pos_virtual
    else:
        success = False
