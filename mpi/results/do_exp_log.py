#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:57:35 2019

@author: fuyuan
"""
import re
import os, sys
import numpy

#some_list = ['abc-123', 'def-456', 'ghi-789', 'abc-456']
#matching = [s for s in some_list if "abc" in s]
#print(matching)

f = open(sys.argv[1], 'r')
jobid=sys.argv[1].split(".")[1]
t_total=[]
t_get_spread_f=[]
t_phase_2 = []
temp = []
t_eqlb = []
t_stream=[]
t_bounce=[]
t_rho_u=[]
t_spread_vel=[]
t_tail=[]

for line in f.readlines():
    temp = line.split(",")
    for s in temp:
        if "Seconds" in s:
            t_total.append(s.split(":")[2])
        if "T_get_Sprd_F" in s:
            t_get_spread_f.append(s.split("=")[1])
        if "T_stream" in s:
            t_phase_2.append(s.split("=")[1])
        if "T_4_1" in s:
            t_eqlb.append(s.split("=")[1])
        if "T_4_2" in s:
            t_stream.append(s.split("=")[1])
        if "T4_3" in s:
            t_bounce.append(s.split("=")[1])
        if "T4_4" in s:
            t_rho_u.append(s.split("=")[1])
        if "T_Sprd_Vel" in s:
            t_spread_vel.append(s.split("=")[1])
        if "T_tail" in s:
            t_tail.append(s.split("=")[1])


#print(t_stream_lines)
#print(t_get_spread_f_lines)
result=(jobid,
      max(t_total),
      max(t_get_spread_f),
      max(t_phase_2),
      max(t_eqlb),
      max(t_stream),
      max(t_bounce),
      max(t_rho_u),
      max(t_spread_vel),
      max(t_tail)
      )
print(' '.join(result))
