# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import re

with open('/home/asengupt/Documents/points.xyz') as f:
    content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    for line in content:
        x = float(re.split('\s+', line)[0])
        y = float(re.split('\s+', line)[1])
        z = float(re.split('\s+', line)[2])
        
        print x,y,z