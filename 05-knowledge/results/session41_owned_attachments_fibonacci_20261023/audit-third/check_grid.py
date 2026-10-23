#!/usr/bin/env python3
import json
from math import isqrt

def need(p,m):
    if not p:raise ValueError(m)

def floor_alpha(m):
    return (3*m-isqrt(5*m*m)-1)//2 if m else 0

def S(n):
    out=0;steps=0
    while n:
        k=floor_alpha(n)
        out+=k*(n+1)-3*k*(k+1)//2
        n=k;steps+=1
    return out,steps

w='C'
while len(w)<20000:w=''.join({'C':'CA','A':'B','B':'CA'}[c] for c in w)
counts={c:0 for c in 'CAB'};totals={c:0 for c in 'CAB'}
for h,c in enumerate(w[:10000],1):
    counts[c]+=1
    need(counts['C']==floor_alpha(h+2),'prefix C')
    need(counts['A']==floor_alpha(h+1),'prefix A')
    for z in counts:totals[z]+=counts[z]
    need(totals['C']==S(h+2)[0],'grid C')
    need(totals['A']==S(h+1)[0],'grid A')
    need(totals['B']==h*(h+1)//2-S(h+2)[0]-S(h+1)[0],'grid B')
    if h<=200:
        need(S(h)[0]==sum(floor_alpha(j) for j in range(1,h+1)),'direct floorsum')
# Irrational root bracketing independently proves each computed floor.
for m in list(range(1,10001))+[10**k+73 for k in (20,50,100,500)]:
    q=floor_alpha(m);a=3*m-2*q;b=a-2
    need(a>0 and a*a>5*m*m,'lower floor bracket')
    need(b<=0 or b*b<5*m*m,'upper floor bracket')
N=10**100
answer,steps=S(N)
print(json.dumps({'status':'PASS','prefix_grid_range':10000,'direct_sum_range':200,'large_horizon':str(N),'large_floor_sum':str(answer),'large_recursion_steps':steps},indent=2))
