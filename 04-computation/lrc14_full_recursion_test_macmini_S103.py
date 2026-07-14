#!/usr/bin/env python3
"""mac-mini-S103b: RIGOROUSLY test the full safe-peel recursion (skeptical, given S97 overclaim history).
Peel safe elements (M preserved) until STUCK; examine the TERMINAL binding core: is it non-covering
(sieve THM-366 => M>=1/14), or small/settled, or a STUCK covering core (the genuine gap)?"""
from fractions import Fraction as F
import numpy as np, random
def M_t0(S,dens=800_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def is_safe_peel(S,v):
    core=[u for u in S if u!=v]
    if len(core)<2: return True,0.5
    mu0,t0=M_t0(core,dens=400_000)
    r=(v*t0)%1.0; vn=min(r,1-r)
    return (vn>=mu0-2e-4), mu0   # v safe at core tight point => M(S)=M(core)
def covering_sub(S): return all(any(u%q==0 for u in S) for q in range(2,15))
def recurse(S):
    """peel safe elements; return (terminal_set, M_preserved). Stops when no safe peel."""
    S=list(S)
    while len(S)>2:
        peeled=False
        for v in sorted(S,reverse=True):   # prefer peeling large (outliers) first
            safe,_=is_safe_peel(S,v)
            if safe: S=[u for u in S if u!=v]; peeled=True; break
        if not peeled: break
    return S
def covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
rng=random.Random(55)
tot=0; term_noncov=0; term_small=0; stuck_cov=[]; term_Ms=[]
for _ in range(20000):
    S=set()
    for q in range(2,15):
        if any(u%q==0 for u in S): continue
        S.add(q*rng.randint(1,20))
    while len(S)<13: S.add(rng.randint(1,350))
    S=sorted(set(list(S)[:13]))
    if len(S)!=13 or 1 not in S or not covering(S): continue
    if max(S)/min([u for u in S if u>1]) < 8: continue   # spread/unsafe stratum
    tot+=1
    T=recurse(S); Mt,_=M_t0(T); term_Ms.append((len(T),Mt))
    if len(T)<=5: term_small+=1
    elif not covering_sub(T): term_noncov+=1
    else: stuck_cov.append((len(T),round(Mt,4),sorted(T)))
    if tot>=500: break
print(f"spread covering families: {tot}")
print(f"  recursion terminated at SMALL set (<=5 elts, LRC-trivial): {term_small}")
print(f"  terminated at NON-COVERING core (sieve THM-366 => M>=1/14): {term_noncov}")
print(f"  STUCK at a COVERING binding core (the genuine gap): {len(stuck_cov)}")
if stuck_cov:
    print("  stuck cores (size, M, set):")
    for sz,M,T in stuck_cov[:8]: print(f"     size={sz} M={M}: {T}")
    szs=[s[0] for s in stuck_cov]; print(f"  stuck-core sizes: min={min(szs)} median={sorted(szs)[len(szs)//2]} max={max(szs)}")
print()
if not stuck_cov:
    print("=> the safe-peel recursion ALWAYS terminates at a small or NON-COVERING core (sieve) => M>=1/14.")
    print("   THM-751 (safe branch) + THM-366 sieve would close the loose stratum ELEMENTARILY. VERIFY-then-prove.")
else:
    print("=> STUCK covering binding cores exist = the irreducible density-floor residual. THM-751 reduces")
    print("   the loose stratum to these small all-binding covering cores (opus's floor / finite check).")
