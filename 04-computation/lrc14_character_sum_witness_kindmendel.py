#!/usr/bin/env python3
"""Beurling-Selberg/character-sum route to HYP-2864: the witness COUNT
N(S,D)=#{a coprime mod D: ||s a/D||>=1/14 for all s} has main term (6/7)^13 * phi(D).
If N>0 for covering S (error < main, via covering controlling resonances), LRC(14) hard case closes.
kind-mendel-2026-06-22-S5."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)
def is_covering(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def phi(D): return sum(1 for a in range(1,D+1) if gcd(a,D)==1)
def N_witness(S,D):
    "count multipliers a coprime to D with ||s a/D||>=1/14 for all s"
    return sum(1 for a in range(1,D) if gcd(a,D)==1
               and all(nrm(F(s*a,D))>=F(1,14) for s in S))

main=lambda D:(6/7)**13*phi(D)
print("Main term (6/7)^13 = %.4f"%((6/7)**13))
print("\n=== covering sets: witness count N vs main term (6/7)^13*phi(D) ===")
random.seed(1)
cov=[]
while len(cov)<5:
    S=sorted(random.sample(range(1,500),13))
    if gall(S)==1 and is_covering(S): cov.append(S)
cov.append([1,2,3,4,5,6,7,8,9,10,11,13,84])  # loosest
for S in cov:
    row=[]
    for D in [41,71,83,89,97]:
        N=N_witness(S,D); row.append(f"D={D}:N={N}(main{main(D):.1f})")
    tag='loosest' if S[-1]==84 else ''
    print(f"  {str(S[:4])+'..' :16s}{tag:8s} "+"  ".join(row))
print("\n=== contrast: the TIGHT AP {1..13} (non-covering, M=1/14) needs D|14 ===")
S=list(range(1,14))
for D in [14,28,41,42,83,89]:
    print(f"  AP: D={D}: N={N_witness(S,D)} (main {main(D):.1f})  [witness only when 14|D]")
print("\n=== why covering helps: N>0 robustly for covering; track growth with D ===")
S=cov[0]
for D in range(15,100,8):
    N=N_witness(S,D)
    print(f"  D={D:3d}: N={N:3d}  main={main(D):5.1f}  ratio N/main={N/main(D) if main(D)>0 else 0:.2f}")
