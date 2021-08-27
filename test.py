from DubinsPath3D import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

start_pos = array([[2.0],[5.0]])
end_pos   = array([[13.0], [6.0]])
start_alt = 60.0
#end_pos   = array([[18.76], [0.239], [0.0]])

R = 2
hStart = 50
hEnd = 90
turn_rate = 1
airspeed = 2
windspeed = 0.4
winddirection = 135

xmin = min([start_pos[0,0],end_pos[0,0]])-(2*R + 4)
xmax = max([start_pos[0,0],end_pos[0,0]])+(2*R + 4)
ymin = min([start_pos[1,0],end_pos[1,0]])-(2*R + 4)
ymax = max([start_pos[1,0],end_pos[1,0]])+(2*R + 4)

arrow_start_delta =  1.0*array([[cos(hStart*pi/180)],[sin(hStart*pi/180)],[0]] )
arrow_end_delta =1.0*array([[cos(hEnd*pi/180)],[sin(hEnd*pi/180)],[0]])
fig3 = plt.figure(3)
ax = fig3.gca(projection='3d')

for i in range(1,9):
    (case,status,paramsCSC,paramsCCC,paramsALT) = FindDubinsParameters3D(start_pos,end_pos,hStart,hEnd,R,i,start_alt,1.0,1.2)

    if(case <= 4):
        plt.figure(1)
        plt.subplot(2,2,i)
        params = paramsCSC
    else:
        plt.figure(2)
        plt.subplot(2,2,i%5 + 1)
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
        (X,Y,Z,timeTaken) = GenerateDubinsAirPath3D(airspeed,start_pos,hStart,end_pos,hEnd,R,start_alt,case,params,paramsALT)
        (Xg,Yg,Zg,timeTaken) = GenerateDubinsGroundPath3D(turn_rate,airspeed,windspeed,winddirection,start_pos,hStart,end_pos,hEnd,R,start_alt,case,params,paramsALT)
        #print type + ":"
        #print (windspeed*cos(winddirection*pi/180)*timeTaken,windspeed*sin(winddirection*pi/180)*timeTaken)
    else:
        (X,Y) = (start_pos[0,0],start_pos[1,0])
        (Xg,Yg) = (X,Y)

    plt.plot(start_pos[1,0],start_pos[0,0],'ro')
    plt.plot(end_pos[1,0],end_pos[0,0],'ro')
    plt.title(type)
    plt.arrow(start_pos[1,0],start_pos[0,0],arrow_start_delta[1,0],arrow_start_delta[0,0],head_width = 0.5)
    plt.arrow(end_pos[1,0],end_pos[0,0],arrow_end_delta[1,0],arrow_end_delta[0,0],head_width=0.5)
    plt.ylim([xmin,xmax])
    plt.xlim([ymin,ymax])
    plt.plot(Y,X,'b')
    plt.plot(Yg,Xg,'g')

    if(case == 1):
        print X[-1],Y[-1],Z[-1]
        ax.plot(X,Y,Z)
        ax.plot(Xg,Yg,Zg,'g--')

plt.show()
