#!/usr/bin/env python3
"""
lrc14_crux_anatomy — mac-mini-2026-06-17-S4

Anatomy of the LAST inequality: for covering 13-sets S=A u {w}, w==0 mod14,
M(S)=max_tau min_v||v tau|| is attained at a binding-pair sum-crossing
tau*=num/D, D=f+w (f=flank core runner), M=k/D, and we need M>=1/14 <=> D<=14k.
But k=(f*num) mod D is achievable for MANY num; the REAL constraint is that the
SAME num keeps ALL OTHER runners >=k/D (others-clear). So the honest crux is:
"some pair-crossing clears everyone at level >=1/14."

This script dissects: the binding pair (does it always contain w?), D, k, the
inequality D<=14k + margin, the others-clear slack, and the flank/gcd structure.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from functools import reduce
import random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))

def anatomy(S):
    S=sorted(set(S)); Mv,tau=M(S)
    D=tau.denominator; num=tau.numerator
    binders=[v for v in S if nrm(v*tau)==Mv]
    k=Mv.numerator if Mv.denominator==D else (Mv*D)  # M=k/D
    # others-clear slack: smallest distance among NON-binders minus M
    nonb=[float(nrm(v*tau)) for v in S if v not in binders]
    slack_others=(min(nonb)-float(Mv)) if nonb else 0.0
    return dict(M=Mv,tau=tau,D=D,num=num,binders=binders,k=int(Mv*D),
                ineq_ok=(D<=14*int(Mv*D)),margin=14*int(Mv*D)-D,
                has_w=any(v%14==0 for v in binders),slack_others=slack_others)

print("="*78)
print("Binding-pair anatomy of covering 13-sets: M=k/D, D=f+w, need D<=14k")
print("="*78)
# structured covering families {1..13}\{j} u {84m} (j in covering-drops) + searched
fams=[]
for j in [1,2,3,4,5,6,7,12]:
    for m in [1,2,3,5]:
        S=sorted([v for v in range(1,14) if v!=j]+[84*m])
        if len(set(S))==13 and covering(S): fams.append((f"drop{j} u{84*m}",S))
worst=(99,None)
for name,S in fams:
    a=anatomy(S)
    bp=a['binders']; pair = (bp[0],bp[-1]) if len(bp)>=2 else bp
    print(f"  {name:14s}: M={a['M']}={float(a['M']):.5f} D={a['D']} k={a['k']} binders={bp} "
          f"D<=14k:{a['ineq_ok']} margin={a['margin']} w-in-pair:{a['has_w']} others-slack={a['slack_others']:.4f}")
    rat=a['D']/(14*a['k'])
    if rat<99 and rat>worst[0] if False else (a['margin']<worst[0]):
        worst=(a['margin'],name)

print("\n"+"="*78)
print("Does the binding pair ALWAYS contain the parked runner w? (random covering sets)")
print("="*78)
random.seed(4); ntot=0; w_in=0; pair_cnt={}
samples=[]
# build random covering 13-sets: take {1..13}\{j} u {big mult of lcm-ish}
for _ in range(400):
    j=random.choice([1,2,3,4,5,6,7,12])
    m=random.randint(1,8)
    w=84*m*random.choice([1,1,1,2])
    S=sorted(set([v for v in range(1,14) if v!=j]+[w]))
    if len(S)!=13 or not covering(S): continue
    a=anatomy(S); ntot+=1
    if a['has_w']: w_in+=1
    pair_cnt[len(a['binders'])]=pair_cnt.get(len(a['binders']),0)+1
    if not a['ineq_ok']:
        print("  !!! D>14k (would break LRC):",S,a)
    samples.append(a)
print(f"  covering sets sampled: {ntot}")
print(f"  binding pair contains w (==0 mod14): {w_in}/{ntot} = {w_in/max(ntot,1):.1%}")
print(f"  #binders distribution: {dict(sorted(pair_cnt.items()))}")
print(f"  D<=14k held: {sum(1 for a in samples if a['ineq_ok'])}/{ntot}")
mn=min(samples,key=lambda a:a['margin'])
print(f"  smallest margin 14k-D = {mn['margin']} at D={mn['D']} k={mn['k']} (M={mn['M']})")

print("\n"+"="*78)
print("The flank f = D - w (when w in pair): which core runner flanks, and gcd(f,w)")
print("="*78)
for name,S in fams[:8]:
    a=anatomy(S)
    if a['has_w']:
        w=[v for v in a['binders'] if v%14==0][0]; f=a['D']-w
        print(f"  {name:14s}: w={w} flank f={f} gcd(f,w)={gcd(f,w)} k={a['k']} "
              f"(f*num mod D = {(f*a['num'])%a['D']}, =k? {(f*a['num'])%a['D']==a['k'] or a['D']-(f*a['num'])%a['D']==a['k']})")
print("\nKEY: M(S)=k/D with D=f+w; need D<=14k. k=(f*num mod D); num also makes ALL")
print("others >= k/D (others-clear). The honest crux = a clearing crossing exists.")
