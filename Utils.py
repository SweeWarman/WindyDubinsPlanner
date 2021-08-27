import numpy as np
from numpy import sin,cos,arctan2,\
                  arcsin,arccos,mod,pi,\
                  array,sqrt,dot
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

def dist(A,B):
    return np.sqrt(np.dot((A-B).T,(A-B)))[0,0]

def distH(A,B):
    return np.sqrt( (A[0,0] - B[0,0])**2 + (A[1,0] - B[1,0])**2)

def Angle(X):
    return arctan2(X[1,0],X[0,0])

def RotateZ(theta):
    return np.array([[cos(theta),-sin(theta)],
                     [sin(theta), cos(theta)]])

def AngleDiffCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi + angle2 - angle1,2*pi)

def AngleDiffCCW(angle1,angle2):
    angle1 = mod(angle1,2*pi)
    angle2 = mod(angle2,2*pi)
    return mod(2*pi - angle2 + angle1,2*pi)


