"""Exact controls for explicit induced unit-distance four-path amplifier."""
from fractions import Fraction as Q
from itertools import product
import json

def check(b):
    if not b: raise RuntimeError('check failed')

def arch(j):
    t=Q(1,2)-Q(1,20*(j+1))
    c=(1-t*t)/(1+t*t); s=2*t/(1+t*t)
    rad=1-(1-c)**2
    return c,s,rad

# Rational inequalities establish the all-pair geometric proof's hypotheses.
arches=[arch(j) for j in range(1,101)]
for c,s,a in arches:
    check(Q(3,5)<c<Q(2,3)); check(c*c+s*s==1)
    check(a==1-(1-c)**2); check(a>Q(21,25))
check(len({c for c,s,a in arches})==100)
# Exact cross-incidence check independent of circle-intersection argument.
# For low x_i and high y_j, squared-distance minus 1 = A+B*sqrt(rad_j).
# Nonzero norm proves no equality. Reflection gives the z_i case as well.
for i,(ci,si,ai) in enumerate(arches):
    for j,(cj,sj,aj) in enumerate(arches):
        A=(ci-1)**2+(si-sj)**2+aj-1
        B=-2*(si-sj)
        if i==j: check(A==0 and B==0)
        else: check(A*A-B*B*aj != 0)
        check((ci-(2-cj))**2+(si-sj)**2 < 1)
        check((ci-cj)**2+(si-sj)**2 < 1)
# Kernel multiplication, independent of closed forms.
channels={}
for k in range(2,9):
    M=[[int(a!=b) for b in range(k)] for a in range(k)]
    A=[[int(a==b) for b in range(k)] for a in range(k)]
    for _ in range(4):
        A=[[sum(A[i][l]*M[l][j] for l in range(k)) for j in range(k)] for i in range(k)]
    e=((k-1)**4+k-1)//k; d=((k-1)**4-1)//k
    check(A[0][0]==e and A[0][1]==d and e-d==1)
    # Literal independent colouring of 3 internal vertices.
    for b,w in [(0,e),(1,d)]:
        count=sum(all(x!=y for x,y in zip((0,)+p,p+(b,))) for p in product(range(k),repeat=3))
        check(count==w)
    channels[k]=[e,d]
# Full parallel graph on two paths, independently enumerate internal colourings.
for k in (2,3,4):
    for b in (0,1):
        good=0
        for p in product(range(k),repeat=6):
            good+=all(all(x!=y for x,y in zip((0,)+p[l:l+3],p[l:l+3]+(b,))) for l in (0,3))
        check(good==channels[k][b]**2)
check(channels[6]==[105,104])
print(json.dumps({'status':'PASS','algebraic_arches':100,'exact_cross_incidence_pairs':10000,'palettes':channels,'six_colour_channels_r3':[105**3,104**3], 'scope':'rational hypotheses checked; exact radical-norm cross-incidence checked; all-r incidence proved in PROOF.md; matrix and literal colour counts agree'},sort_keys=True,indent=2))
