#!/usr/bin/env python3
"""Exact bounds and small examples for the edge-rooted satellite grid."""
from fractions import Fraction as F
from math import factorial, prod
from functools import cache

w=tuple(F(x,1000) for x in (1000,1000,1000,1000,1000,1000,1000,1071,1149,1233,1324,1421,1525,1637,1758,1887,2025,2174,2334,2505,2689,2886,3098,3326,3570,3832,4113,4415,4740,5088,5461,5862,6293,6755,7251,7784,8356,8970,9629,10336,11095,11910,12786,13725,14732,15816,16978,18225,19565,21003,22546,24203,25984,27891,29939,32142))
wa=F(33844,1000)
e0=F(27182818,10**7)
e1=F(271828183,10**8)
partial=sum(F(1,factorial(i)) for i in range(14))
assert e0<partial and partial+F(1,13*factorial(13))<e1

# The exact layer bounds also used in the native affine-threshold note.
f=1
for d in range(1,3252):
    f*=d
    assert f*e1.numerator**d<((d+5)*e1.denominator)**d
    if d<=59:
        assert f*(500*e1.numerator)**d<((510*d+1087)*e1.denominator)**d

from math import comb

def insertion(k,L,B,c,bits):
    N=1568
    r=N-k
    assert L>=k and B>=2*k and c>=1
    b=r-110
    E=27*r+1
    mean=F(k)+F(27*k*(k-1),L-1)
    share=1-mean/(c*k+1)
    assert share>0
    num=27**(27*N)*factorial(r)*b**E*comb(L,k)*(factorial(B)//factorial(B-2*k))*share.numerator
    den=28**(26*N)*(b+1)**E*share.denominator
    assert num**c > (den**c * factorial(c)**k) << (c*bits)

insertion(211,1568,487,5,21386)
insertion(128,434,381,9,20445)
insertion(131,448,371,9,20443)
print('Floors: 21386 cyclic halves; 20445 monotone and selected native; 20443 near-regular rows')

def ceil(x):
    return -(-x.numerator//x.denominator)

def upper(mate):
    m=56
    scale=10**6
    f=ceil((m*sum(wa/t for t in w)+5*wa)/e0*scale)
    for i,t in enumerate(w):
        backwards=[t/u for u in w[:i]]
        rowbase=sum(backwards)+(m-1)*sum(t/u for u in w[i+1:])+t/wa
        extra=sum(max(v-1,0) for v in backwards+[t/wa])
        ordinary=rowbase+F(1413,500)+F(1087,500)*t+F(extra,50)
        first=last=0
        if mate:
            first=-1
            if i:
                v=t/w[i-1]
                first+=v+(v-1)/50
            last=1-(t/w[i+1] if i+1<m else 0)
        mid={s:ceil((ordinary+s)/e0*scale) for s in (27,28)}
        def blocks(a,b):
            return (ceil((ordinary+a+first)/e0*scale)*
                    ceil((ordinary+b+last)/e0*scale)*
                    mid[27]**(28-(a==27)-(b==27))*mid[28]**(28-(a==28)-(b==28)))
        f*=max(blocks(a,b) for a in (27,28) for b in (27,28))
    assert 10*f < 9*(scale**(m*m+1)<<25628)

upper(False)
upper(True)
print('Near-regular rows: both upper endpoints certified')

def native_insertion():
    from pathlib import Path
    from runpy import run_path
    old=run_path(str(Path(__file__).with_name('spread-native-deficits.py')))
    m,q,E=old['m'],old['q'],old['E']
    alpha,beta=old['alpha'],old['beta']
    assert m==1547 and q==78 and E==52929
    k,L,B,c=88,381,300,14
    r=m-k
    b=r-238
    share=1-(F(k)+F(56*k*(k-1),L-1))/(c*k+1)
    assert share>0
    num=(22**(22*m)*77**(77*2*q)*alpha.numerator**(2*q)*beta.numerator*
         factorial(r)*b**E*comb(L,k)*(factorial(B)//factorial(B-2*k))*share.numerator)
    den=(23**(21*m)*78**(76*2*q)*factorial(q)**2*alpha.denominator**(2*q)*beta.denominator*
         (b+1)**E*share.denominator)
    a,b=num**c,den**c*factorial(c)**k
    bits=(a.bit_length()-b.bit_length())//c
    while a<=b<<(c*bits):
        bits-=1
    assert bits>=19503
    print('Uniform original-native floor:',bits)

if __name__=='__main__':
    native_insertion()
