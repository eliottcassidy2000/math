#!/usr/bin/env python3
"""Checksum and readout inherited from the six-vertex physical window."""
Q,H,A,B,F,E=432,45,129,258,15,30
K=(135,137,139)

def p(x):
    s=0
    while x:
        if not x&1:
            m=x//2;s+=m+(m&1)
        x//=2
    return s

def value(x,y,j):
    return Q*x*y+A*x+B*y+F*(x&1)+E*(y&1)+H*(p(x)+p(y))+K[j]

