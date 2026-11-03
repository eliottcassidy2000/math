#!/usr/bin/env python3
"""Exact EGF series check with an explicit missing-boundary hostile."""
from fractions import Fraction as Q
from math import factorial
import json
from audit_carrier import reduced_count
L=12
A=[Q(0) for _ in range(L+1)];B=A.copy();D=A.copy()
for a in range(L+1):
 for b in range((L-a)//2+1):A[a+2*b]+=Q(1,factorial(a)*2**b*factorial(b))
for h in range(L//3+1):
 for j in range((L-3*h)//2+1):B[3*h+2*j]+=Q(1,2**h*factorial(h)*factorial(h+2*j))
for a in range(L+1):
 for j in range(1,(L-a)//2+1):D[a+2*j]+=Q(1,factorial(a)*factorial(2*j))
raw=[sum(A[j]*B[k-j] for j in range(k+1)) for k in range(L+1)]
correct=[(raw[k]-D[k])*factorial(k) for k in range(L+1)]
expected=[reduced_count(k) for k in range(L+1)]
if correct!=expected:raise RuntimeError('EGF coefficient disagreement')
if raw[2]*2!=3 or correct[2]!=2:raise RuntimeError('k0 positive interior hostile did not trigger')
print(json.dumps({'orders':[0,L],'counts':expected,'formula':'exp(x+x^2/2) sum_{h,j>=0} x^(3h+2j)/(2^h h! (h+2j)!) - exp(x)(cosh(x)-1)','checks':14,'missing_boundary_hostile':{'without_subtraction_c2':3,'actual_c2':2}},indent=2))
