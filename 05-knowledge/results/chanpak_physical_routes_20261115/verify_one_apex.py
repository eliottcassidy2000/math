#!/usr/bin/env python3
"""Reproduce explicit calibrated marked routes and the one-site bound."""
from random import Random

def need(condition):
    if not condition:
        raise RuntimeError('route theorem control failed')

def p(n):
    s=0
    while n:
        n,b=divmod(n,2)
        if not b:s+=n+(n&1)
    return s

def w(a,b,mask):
    L=bin(a)[2:] if a else ''
    R=''.join('0' if d=='1' else '1' for d in (bin(b)[2:] if b else '')[::-1])
    return int(L+f'{1158^mask:011b}'+R,2)

def verify_route(a,b,slices,M):
    vals=[1+p(w(a,b,s)) for s in slices]
    need(all(j-i==M for i,j in zip(vals,vals[1:])))
    need(all(w(a,b,s)&1==0 for s in slices))
    return vals

for t in range(400):
    a,b=9*t,10*t
    f=verify_route(a,b,[0,32,220,252],2+32*t)
    if t==1:need(f==[208063,208097,208131,208165])
    a,b=5*t+1,2*t
    verify_route(a,b,[32,0,252,220],30+96*t)

rng=Random(20261115)
n=0
for m in range(7,29):
    k=m//2
    G=(1<<(k-2))-1 if m&1 else ((1<<(k-2))-2)//2
    for _ in range(300):
        z=(1<<(m-1))|(rng.randrange(1<<(m-2))<<1)
        site=[i for i in range(m) if abs(p(z^(1<<i))-p(z))<=2*G]
        need(len(site)<=1)
        if site:
            y=z^(1<<site[0])
            next_sites=[i for i in range(m) if abs(p(y^(1<<i))-p(y))<=2*G]
            need(next_sites==site)
        n+=1
print('800 exact physical routes; '+str(n)+' marked-word small-component controls')
