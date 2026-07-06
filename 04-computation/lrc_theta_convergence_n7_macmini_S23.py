#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S23 -- working the sole open piece: VALIDATE opus's theta identity
safe(S,beta)=SUM_{a in L(S)} prod hhat(a_i) and PIN ITS CONVERGENCE RATE at the feasible
n=7 scale, to ground the Beurling-Selberg N-scaling (my S22 HYP-4512).

For each n=7 family (6 runners), compute:
 (1) safe(S, 2/13) EXACTLY (arc-measure complement);
 (2) the theta-sum truncated to ALL relations with |a_i| <= N, for N=1..4;
 (3) the convergence trunc_N -> safe, and how large N must be to certify the sign.

The Beurling-Selberg majorant bound safe >= INT prod(1-g+) is a truncated theta with
the majorant's coefficients; its band-limit N must be large enough that the truncation
error < the window slack.  This measures that error directly.
"""
from fractions import Fraction as F
from math import sin, pi, gcd
from itertools import product, combinations

def hhat(m, beta):
    if m == 0: return 1 - 2*beta
    return -sin(2*pi*m*beta)/(pi*m)

def arcs_danger(v, beta):
    out=[]; j=0
    while F(j)-F(beta*10**9,10**9)/1 < v:   # crude; use float below instead
        break
    return out

def safe_exact(S, beta):
    """|{t: ||v_i t||>=beta all i}| via exact rational arcs."""
    b=F(beta).limit_denominator(10**6)
    iv=[]
    for v in S:
        j=0
        while F(j)-b < v:
            lo=max((F(j)-b)/v, F(0)); hi=min((F(j)+b)/v, F(1))
            if lo<hi: iv.append((lo,hi))
            j+=1
    iv.sort(); danger=F(0); ch=None; cl=None
    for lo,hi in iv:
        if ch is None: cl,ch=lo,hi
        elif lo<=ch: ch=max(ch,hi)
        else: danger+=ch-cl; cl,ch=lo,hi
    if ch is not None: danger+=ch-cl
    return float(F(1)-danger)

def theta_trunc(S, beta, N):
    """SUM over relations a with sum a_i v_i=0, |a_i|<=N, of prod hhat(a_i). Full (not sparse)."""
    n=len(S); total=0.0
    for a in product(range(-N,N+1), repeat=n):
        if sum(a[i]*S[i] for i in range(n))!=0: continue
        term=1.0
        for i in range(n): term*=hhat(a[i],beta)
        total+=term
    return total

beta=2/13
fams={
 "AP {1..6} (safe=0, tiles)": [1,2,3,4,5,6],
 "n=7 gap {1,5,6,11,16,17} (safe=0, tiles)": [1,5,6,11,16,17],
 "generic {1,2,3,4,5,20} (loose)": [1,2,3,4,5,20],
 "loose {1,2,4,8,9,15}": [1,2,4,8,9,15],
}
print(f"n=7, beta=2/13={beta:.4f}: theta-sum truncation |a_i|<=N -> exact safe\n")
print(f"  {'family':42s} {'safe(exact)':>11} {'N=1':>8} {'N=2':>8} {'N=3':>8}")
for name,S in fams.items():
    se=safe_exact(S,beta)
    ths=[theta_trunc(S,beta,N) for N in (1,2,3)]
    print(f"  {name:42s} {se:>11.5f} " + " ".join(f"{t:>8.4f}" for t in ths))
print("\nREADING: does trunc_N -> safe, and how big must N be to get the SIGN right (safe=0 vs >0)?")
print("The truncation error at N ~ the tail the Beurling-Selberg majorant must control;")
print("if N=1,2,3 already track the sign at n=7 (wide window 1/91), the route is sound and")
print("the n=13 obstruction is purely the WIDTH (N~2k^2=288), not a failure of the identity.")
