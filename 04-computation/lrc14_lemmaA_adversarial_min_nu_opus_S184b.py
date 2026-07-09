"""
lrc14_lemmaA_adversarial_min_nu_opus_S184b.py  (opus-2026-07-09-S184)

Scope Lemma A (nu(E) >= nuConsec(13) = 477/1078 = 0.44249). ADVERSARIAL search: hill-climb to MINIMIZE
nu(E)=meas{x: maxgap{frac(e_i x)}>1/7} over primitive 13-clusters. If the global min is the consecutive
cluster (nu=0.4425) and nothing beats it, Lemma A holds; the min-nu-vs-spread profile scopes the "hard
core" (where nu is near 0.4425) vs the easy tail (nu>>0.4425). Answers: is consecutive THE global min,
and how small is the core the finite check must cover?
"""
import numpy as np, random
NGc=24001                      # coarse grid for search
Xc=(np.arange(1,NGc))/NGc
NUCONSEC=477/1078

def nu(E, X):
    E=np.array(sorted(set(E)),dtype=float)
    P=np.sort((np.outer(E,X))%1.0,axis=0)
    mg=np.maximum(np.diff(P,axis=0).max(axis=0),(P[0]+1.0)-P[-1])
    return float((mg>1.0/7).mean())

# reference
print(f"nuConsec(13) = 477/1078 = {NUCONSEC:.5f}")
print(f"nu(consecutive {{0..12}}) coarse = {nu(list(range(13)),Xc):.5f}")
print("="*80)
print("ADVERSARIAL min-nu hill-climb (primitive 13-clusters, minimize nu):")
rng=random.Random(1)
globalmin=1.0; argmin=None
for trial in range(120):
    # random primitive 13-cluster, spread in [12, 60]
    sp=rng.randint(12,60)
    E=sorted(rng.sample(range(0,sp+1),13))
    cur=nu(E,Xc)
    for _ in range(400):
        i=rng.randrange(13); old=E[i]; cand=rng.randint(0,max(E)+3)
        if cand in E: continue
        E2=sorted(set(E)-{old}|{cand})
        if len(E2)!=13: continue
        v=nu(E2,Xc)
        if v<cur-1e-9: E,cur=E2,v
    if cur<globalmin: globalmin=cur; argmin=E[:]
print(f"  global min nu found = {globalmin:.5f}")
print(f"  minimizer = {argmin}")
diffs=[argmin[i+1]-argmin[i] for i in range(12)] if argmin else []
g=0
if argmin:
    import math
    g=math.gcd(*[a for a in argmin if a>0]) if any(a>0 for a in argmin) else 1
print(f"  minimizer diffs = {diffs}  gcd={g}  (consecutive => all 1s)")
print(f"  is min >= nuConsec? {globalmin >= NUCONSEC - 1e-4}  (margin {globalmin-NUCONSEC:+.5f})")
print("-"*80)
print("min nu by SPREAD (adversarial per spread, 30 climbs each):")
for sp in [12,13,14,16,18,20,24,30,40]:
    best=1.0
    for t in range(30):
        E=sorted(rng.sample(range(0,sp+1),13))
        if max(E)!=sp: E[-1]=sp; E=sorted(set(E))
        if len(E)!=13: continue
        cur=nu(E,Xc)
        for _ in range(200):
            i=rng.randrange(13); cand=rng.randint(0,sp)
            if cand in E: continue
            E2=sorted(set(E)-{E[i]}|{cand})
            if len(E2)==13 and max(E2)==sp and nu(E2,Xc)<cur: E,cur=E2,nu(E2,Xc)
        best=min(best,cur)
    print(f"  spread {sp:>3}: min nu ~ {best:.4f}  ({'>= nuConsec' if best>=NUCONSEC-1e-3 else 'BELOW!'})")
print("="*80)
print("READING: if global min = consecutive (all-1s diffs, nu=0.4425) and min-nu>=nuConsec at every")
print("spread, Lemma A holds; the core (nu near 0.4425) is the small-spread near-consecutive band.")
