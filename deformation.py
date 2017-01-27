# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re
import transformations as tr
import numpy as np


alpha = 0.0872665
beta = 0.122173
gamma = 0

tx = 0.00
ty = 0.00
tz = 0.00

R = tr.euler_matrix(alpha,beta,gamma)

R[0,3] = tx
R[1,3] = ty
R[2,3] = tz

print R

with open('/home/asengupt/Downloads/sample/bunny2.xyz') as f:
    content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    for line in content:
        x = float(re.split('\s+', line)[0])
        y = float(re.split('\s+', line)[1])
        z = float(re.split('\s+', line)[2])
        
        print x,y,z
        
        b = np.array([(x),(y),(z),(1)])
        
        P = np.dot(R,b)
        
        print P
        
        if ((float(re.split('\s+', line)[3]) != 0) and (float(re.split('\s+', line)[4]) != 0) and (float(re.split('\s+', line)[5]) != 0)):
        
            with open("/home/asengupt/Downloads/sample/source/original.xyz", "a") as myfile:
                """myfile.write(str(x)+' '+str(y)+' '+str(z)+' '+str(230)+' '+str(0)+' '+str(10)+'\n')
                myfile.write(str(P[0])+' '+str(P[1])+' '+str(P[2])+' '+str(10)+' '+str(0)+' '+str(225)+'\n')"""
                
                #myfile.write(str(x)+' '+str(y)+' '+str(z)+' '+re.split('\s+', line)[3]+' '+re.split('\s+', line)[4]+''+re.split('\s+', line)[5]+'\n')
                #myfile.write(str(P[0])+' '+str(P[1])+' '+str(P[2])+' '+re.split('\s+', line)[3]+' '+re.split('\s+', line)[4]+''+re.split('\s+', line)[5]+'\n')
        
        
                #myfile.write(str(x)+' '+str(y)+' '+str(z)+' '+re.split('\s+', line)[3]+' '+re.split('\s+', line)[4]+' '+re.split('\s+', line)[5]+'\n')
                myfile.write(str(P[0])+' '+str(P[1])+' '+str(P[2])+'\n')