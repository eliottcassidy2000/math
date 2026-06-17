#!/usr/bin/env python3
"""
lrc14_gap_M_exact — mac-mini-2026-06-16-S3

THE ACTUAL LRC QUANTITY, exactly. LRC(14) for a 13-set S says
    M(S) := max_{tau} min_{v in S} ||v*tau||  >=  1/14.
A COUNTEREXAMPLE is any primitive 13-set with M(S) < 1/14 (the closed danger
bands cover the circle; no point is >=1/14 from all runners). This is sharper
than the open lonely measure L(S): L=0 is necessary but NOT sufficient for a
counterexample (need the CLOSED lonely set empty, i.e. M<1/14).

EXACT computation: g(tau)=min_i ||v_i tau|| is a lower envelope of tent (triangle)
waves; its maximum is attained at a vertex of the envelope, i.e. at
  tau = k/(v_a+v_b) or k/(v_a-v_b)  (crossing of two tents)  or  (2k+1)/(2v_i) (a tent peak).
We enumerate ALL such candidate tau in (0,1/2] (g is even, 1-periodic), evaluate g
exactly (fractions), and take the max = M(S). This is rigorous (the max IS at such
a vertex). Denominators are bounded by 2*max(S), so it is exact and fast.
"""
from fractions import Fraction as F
from math import gcd

def frac_norm(x):  # ||x|| = distance to nearest integer, exact
    r = x - int(x)            # in (-1,1)
    if r < 0: r += 1          # [0,1)
    return r if r <= F(1,2) else 1 - r

def g_exact(S, tau):
    return min(frac_norm(v*tau) for v in S)

def candidate_taus(S):
    S=sorted(set(S)); cands=set()
    mx=max(S)
    # tent peaks (2k+1)/(2v)
    for v in S:
        k=0
        while True:
            t=F(2*k+1,2*v)
            if t>F(1,2): break
            cands.add(t); k+=1
    # tent crossings k/(v_a +/- v_b)
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d<=0: continue
                k=1
                while True:
                    t=F(k,d)
                    if t>F(1,2): break
                    cands.add(t); k+=1
    cands.add(F(1,2))
    return cands

def M_exact(S):
    best=F(0); argbest=None
    for t in candidate_taus(S):
        val=g_exact(S,t)
        if val>best: best=val; argbest=t
    return best, argbest

def report(S, name=""):
    S=sorted(set(S))
    M,tau=M_exact(S)
    rel = "<" if M<F(1,14) else ("=" if M==F(1,14) else ">")
    verdict = "  *** M < 1/14 : COUNTEREXAMPLE TO LRC(14) ***" if M<F(1,14) else ""
    print(f"  {name or S}: M = {M} = {float(M):.6f}  {rel} 1/14   (lonely tau*={tau})" + verdict)
    return M

print("1/14 =", float(F(1,14)))
print("="*72); print("Sanity: the tight AP and known configs"); print("="*72)
report(list(range(1,14)),                         "{1..13} tight AP")
report([1,2,3,4,5,6,7,8,9,10,11,13,24],           "{1..11,13,24} sporadic tight (L=0)")
report([1,2,3,4,5,6,7,8,9,10,11,13,36],           "{1..11,13,36} inf-L=1/1260")
report([1,2,3,4,5,7,8,9,10,11,12,13,98],          "{1..13}\\{6} u 98  (my restricted inf)")
report(list(range(1,14))+[],                      "")  # noop placeholder

print("\n"+"="*72)
print("The L=0 (tight) configs: is M=1/14 (equality, LRC holds) or <1/14 (counterexample)?")
print("="*72)
# Enumerate single-perturbation tight configs (L=0): from S2, 12->24 was the only one.
# Also test small dilations / sporadics. Compute M for each L=0 config we can find.
def L0_single_perts(maxw=200):
    base=list(range(1,14)); out=[]
    # reuse the measure quickly via float lonely measure
    def Lf(S):
        arcs=[]
        for v in set(S):
            inv=1.0/(14*v)
            for k in range(v+1):
                lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
                if lo<hi: arcs.append((lo,hi))
        arcs.sort(); tot=0.0; cl,ch=arcs[0]
        for lo,hi in arcs[1:]:
            if lo<=ch: ch=max(ch,hi)
            else: tot+=ch-cl; cl,ch=lo,hi
        return 1.0-(tot+ch-cl)
    for e in base:
        for w in range(14,maxw+1):
            if w in base and w!=e: continue
            S=[x for x in base if x!=e]+[w]
            if len(set(S))!=13: continue
            if Lf(S)<1e-12: out.append(tuple(sorted(S)))
    return out
tights=L0_single_perts()
print(f"L=0 single-perturbation configs found (w<=200): {len(tights)}")
anyCE=False
for S in tights:
    M=report(list(S))
    if M<F(1,14): anyCE=True
print("\nANY counterexample (M<1/14) among tight configs:", anyCE)
print("(expected NONE: tight configs satisfy LRC with EQUALITY M=1/14; a counterexample would be historic)")
