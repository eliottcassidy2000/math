#!/usr/bin/env python3
from fractions import Fraction
P=[(17,13*11),(19*7,17*5),(17,19),(23,17),(29*7*5,23*3),(23,29),(31,23),(37*7**3*3,31*2),(31,37),(41,31),(43*2,41*7),(41,43),(13,41)]
if not all(Fraction(a,b).numerator==a and Fraction(a,b).denominator==b for a,b in P): raise AssertionError('fractions must be reduced')
def run(m):
    v=13*11**m*2**17*3**5*5; steps=0
    while True:
        for a,b in P:
            if v%b==0:
                v=(v//b)*a;steps+=1;break
        else:return v,steps
        if steps>1_000_000:raise RuntimeError('step cap')
def valuation(v,p):
    a=0
    while v%p==0:v//=p;a+=1
    return a,v
h=[1,5,17]
for _ in range(10):h.append(3*h[-1]+h[-2]+h[-3])
for m in range(5):
    v,steps=run(m);out=[]
    for prime in [2,3,5]:
        e,v=valuation(v,prime);out.append(e)
    if out!=[h[m+2],h[m+1],h[m]] or v!=13:raise RuntimeError((m,out,v))
    print('m=',m,'steps=',steps,'A,B,C=',out)
print('fractions:',', '.join(f'{a}/{b}' for a,b in P))
