# -*- coding: utf-8 -*-
# basicICP.py

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
"""
Created on Mon Jan 23 18:08:58 2017

@author: agnivsen

A basic Iterative Closes Point matching betwwen two set of point cloud, 
                        or one point cloud and one point + normal pair.
                        
Basic point to plane matching has been done using a Least squares approach and a Gauss-Newton approach
Point to point matching has been done using Gauss-Newton only
"""

import numpy as np
import re
import transformations as transform



def read_file_original(file_path):
    a = []
    with open(file_path) as f:
        content = f.readlines()
        for line in content:
            x = float(re.split('\s+', line)[0])
            y = float(re.split('\s+', line)[1])
            z = float(re.split('\s+', line)[2])
            
            b = np.array([x,y,z])
            a.append(b)
            
    data = np.array(a)
    return data
    
    
def read_file_deformed(file_path):
    a = []
    with open(file_path) as f:
        content = f.readlines()
        for line in content:
            x = float(re.split('\s+', line)[0])
            y = float(re.split('\s+', line)[1])
            z = float(re.split('\s+', line)[2])
            
            nx = float(re.split('\s+', line)[3])
            ny = float(re.split('\s+', line)[4])
            nz = float(re.split('\s+', line)[5])
            
            b = np.array([x,y,z,nx,ny,nz])
            a.append(b)
                 
            
    data = np.array(a)
    return data
    

def icp_point_to_plane(source_points, dest_points,loop):
    """
    Point to plane matching using least squares
    
    source_points:  nx3 matrix of n 3D points
    dest_points: nx6 matrix of n 3D points + 3 normal vectors, which have been obtained by some rigid deformation of 'source_points'
    """
    
    A = []
    b = []
    
    for i in range (0,dest_points.shape[0]-1):
        
        #print dest_points[i][3],dest_points[i][4],dest_points[i][5]
        dx = dest_points[i][0]
        dy = dest_points[i][1]
        dz = dest_points[i][2]
        nx = dest_points[i][3]
        ny = dest_points[i][4]
        nz = dest_points[i][5]
        
        sx = source_points[i][0]
        sy = source_points[i][1]
        sz = source_points[i][2]
        
        _a1 = (nz*sy) - (ny*sz)
        _a2 = (nx*sz) - (nz*sx)
        _a3 = (ny*sx) - (nx*sy)
        
        _a = np.array([_a1, _a2, _a3, nx, ny, nz])
        
        _b = (nx*dx) + (ny*dy) + (nz*dz) - (nx*sx) - (ny*sy) - (nz*sz)
        
        A.append(_a)
        b.append(_b)
        
        
    
    A1 = np.array(A)
    b1 = np.array(b)
    

    A_ = np.linalg.pinv(A1)    
    
    tr = np.dot(A_,b)
    
    print str(tr[0])+','+str(tr[1])+','+str(tr[2])+','+str(tr[3])+','+str(tr[4])+','+str(tr[5])
    
    R = transform.euler_matrix(tr[0],tr[1],tr[2])
    R[0,3] = tr[3]
    R[1,3] = tr[4]
    R[2,3] = tr[5]
    
    source_transformed = []
    
    for i in range (0,dest_points.shape[0]-1):
        ss = np.array([(source_points[i][0]),(source_points[i][1]),(source_points[i][2]),(1)])
        p = np.dot(R,ss)
        source_transformed.append(p)
        
    source_points = np.array(source_transformed)
    
    
    loop = loop + 1
    
    if(loop < 3):   #although this should converge in one step (which it does), you might want to reiterate over and over, just for the fun of it!
    
        icp_point_to_plane(source_points,dest_points,loop)
                
        
def icp_point_to_point_lm(source_points, dest_points,initial,loop):
    """
    Point to point matching using Gauss-Newton
    
    source_points:  nx3 matrix of n 3D points
    dest_points: nx3 matrix of n 3D points, which have been obtained by some rigid deformation of 'source_points'
    initial: 1x6 matrix, denoting alpha, beta, gamma (the Euler angles for rotation and tx, ty, tz (the translation along three axis). 
                this is the initial estimate of the transformation between 'source_points' and 'dest_points'
    loop: start with zero, to keep track of the number of times it loops, just a very crude way to control the recursion            
                
    """
    
    J = []
    e = []
    
    for i in range (0,dest_points.shape[0]-1):
        
        #print dest_points[i][3],dest_points[i][4],dest_points[i][5]
        dx = dest_points[i][0]
        dy = dest_points[i][1]
        dz = dest_points[i][2]
        
        sx = source_points[i][0]
        sy = source_points[i][1]
        sz = source_points[i][2]
        
        alpha = initial[0][0]
        beta = initial[1][0]
        gamma = initial[2][0]
        tx = initial[3][0]        
        ty = initial[4][0]
        tz = initial[5][0]
        #print alpha
        
        a1 = (-2*beta*sx*sy) - (2*gamma*sx*sz) + (2*alpha*((sy*sy) + (sz*sz))) + (2*((sz*dy) - (sy*dz))) + 2*((sy*tz) - (sz*ty))
        a2 = (-2*alpha*sx*sy) - (2*gamma*sy*sz) + (2*beta*((sx*sx) + (sz*sz))) + (2*((sx*dz) - (sz*dx))) + 2*((sz*tx) - (sx*tz))
        a3 = (-2*alpha*sx*sz) - (2*beta*sy*sz) + (2*gamma*((sx*sx) + (sy*sy))) + (2*((sy*dx) - (sx*dy))) + 2*((sx*ty) - (sy*tx))
        a4 = 2*(sx - (gamma*sy) + (beta*sz) +tx -dx)
        a5 = 2*(sy - (alpha*sz) + (gamma*sx) +ty -dy)
        a6 = 2*(sz - (beta*sx) + (alpha*sy) +tz -dz)
        
        _residual = (a4*a4/4)+(a5*a5/4)+(a6*a6/4)
        
        _J = np.array([a1, a2, a3, a4, a5, a6])
        _e = np.array([_residual])
        
        J.append(_J)
        e.append(_e)
        
    jacobian = np.array(J)
    residual = np.array(e)
    
    update = -np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(jacobian),jacobian)),np.transpose(jacobian)),residual)
    
    #print update, initial
    
    initial = initial + update
    
    print np.transpose(initial)
    
    loop = loop + 1
    
    if(loop < 50):  # here lies the control variable, control the number of iteration from here
    
        icp_point_to_point_lm(source_points,dest_points,initial, loop)
        
        
def icp_point_to_plane_lm(source_points, dest_points,initial,loop):
    """
    Point to plane matching using Gauss Newton
    
    source_points:  nx3 matrix of n 3D points
    dest_points: nx6 matrix of n 3D points + 3 normal vectors, which have been obtained by some rigid deformation of 'source_points'
    initial: 1x6 matrix, denoting alpha, beta, gamma (the Euler angles for rotation and tx, ty, tz (the translation along three axis). 
                this is the initial estimate of the transformation between 'source_points' and 'dest_points'
    loop: start with zero, to keep track of the number of times it loops, just a very crude way to control the recursion            
                
    """
    
    J = []
    e = []
    
    for i in range (0,dest_points.shape[0]-1):
        
        #print dest_points[i][3],dest_points[i][4],dest_points[i][5]
        dx = dest_points[i][0]
        dy = dest_points[i][1]
        dz = dest_points[i][2]
        nx = dest_points[i][3]
        ny = dest_points[i][4]
        nz = dest_points[i][5]
        
        sx = source_points[i][0]
        sy = source_points[i][1]
        sz = source_points[i][2]
        
        alpha = initial[0][0]
        beta = initial[1][0]
        gamma = initial[2][0]
        tx = initial[3][0]        
        ty = initial[4][0]
        tz = initial[5][0]
        
        a1 = (nz*sy) - (ny*sz)
        a2 = (nx*sz) - (nz*sx)
        a3 = (ny*sx) - (nx*sy)
        a4 = nx
        a5 = ny
        a6 = nz
        
        _residual = (alpha*a1) + (beta*a2) + (gamma*a3) + (nx*tx) + (ny*ty) + (nz*tz) - (((nx*dx) + (ny*dy) + (nz*dz)) - ((nx*sx) + (ny*sy) + (nz*sz)))
    
        _J = np.array([a1, a2, a3, a4, a5, a6])
        _e = np.array([_residual])
        
        J.append(_J)
        e.append(_e)
        
    jacobian = np.array(J)
    residual = np.array(e)
    
    update = -np.dot(np.dot(np.linalg.pinv(np.dot(np.transpose(jacobian),jacobian)),np.transpose(jacobian)),residual)
    
    #print update, initial
    
    initial = initial + update
    
    print np.transpose(initial)
    
    loop = loop + 1
    
    if(loop < 50):  # here lies the control variable, control the number of iteration from here
    
        icp_point_to_point_lm(source_points,dest_points,initial, loop)


fileOriginal = '/home/asengupt/Documents/icp/data/original.xyz'
deformed = '/home/asengupt/Documents/icp/data/deformed.xyz'

source_points = read_file_original(fileOriginal)
dest_points_et_normal = read_file_deformed(deformed)

initial = np.array([[0.01], [0.05], [0.01], [0.001], [0.001], [0.001]])

# **********************************************************************
# Uncomment one of the following lines to trigger the respective module:
# **********************************************************************


#icp_point_to_plane(source_points,dest_points_et_normal,0)

#icp_point_to_point_lm(source_points,dest_points_et_normal,initial,0)

icp_point_to_plane_lm(source_points,dest_points_et_normal,initial,0)




