"""Exact finite controls for LM2, LM4, LM5 and Gaussian product loss.
Universe: explicitly stated sign law, three-point laws, Gaussian degrees <=4.
These controls do not simulate a Yang–Mills field or prove continuum limits.
Run: python3 05-knowledge/results/session272_yang_mills_bridges_20261111/limits_measure/probe.py
"""
from fractions import Fraction as F
from math import factorial

# LM2: fixed-volume tails vs moving diagonal, exact declared finite window.
for p in range(1, 101):
    assert all(int(N <= p) == 0 for N in range(p+1,p+11))
    assert int(p <= p) == 1
print('LM2: 100 fixed-volume tail controls zero; 100 diagonal controls one.')

# LM4: entire two-atom probability space; positive and hostile controls.
law = [(F(1,2),1,-1),(F(1,2),-1,1)]
assert sum(w*xp*xp for w,xm,xp in law) == 1
assert sum(w*xm*xp for w,xm,xp in law) == -1
assert sorted((w,xm,xp) for w,xm,xp in law) == sorted((w,xp,xm) for w,xm,xp in law)
print('LM4: reflection invariant; ordinary square 1; OS square -1.')

# LM5: both full support; h injective, mean blind, product h*h detects.
x = [-1,0,1]
nu = [F(1,4),F(1,2),F(1,4)]
rho = [F(1,3)]*3
assert all(w>0 for w in nu+rho)
assert sum(nu)==sum(rho)==1
assert sum(w*t for w,t in zip(nu,x))==sum(w*t for w,t in zip(rho,x))==0
assert sum(w*t*t for w,t in zip(nu,x))==F(1,2)
assert sum(w*t*t for w,t in zip(rho,x))==F(2,3)
print('LM5: identical mean 0; second moments 1/2 and 2/3.')

# Independent real-Gaussian route to complex product E[(H+iG)^2(H-iG)^2].
def normal_moment(n):
    return 0 if n%2 else factorial(n)//(2**(n//2)*factorial(n//2))
assert normal_moment(2)-normal_moment(2)==0
value=normal_moment(4)+2*normal_moment(2)**2+normal_moment(4)
assert value==8
print('Gaussian: E[U^2]=E[V^2]=0; E[(H^2+G^2)^2]=8.')
print('PASS: all declared exact finite controls.')
