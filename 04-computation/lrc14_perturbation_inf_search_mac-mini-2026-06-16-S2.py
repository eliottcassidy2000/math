#!/usr/bin/env python3
"""
lrc14_perturbation_inf_search — mac-mini-2026-06-16-S2

INDEPENDENT stress-test of kind-pasteur HYP-2561 / THM-522: is inf_S L(S) = 1/1260,
attained at the minimal perturbation {1..13} with 12->36? The tight AP {1..13} has
L=0 exactly; extremizers are its minimal perturbations. kind-pasteur searched single
+ 2-element perturbations with stranger <=72. We extend: single perturbations to
w<=600, AND 2-element perturbations (float-screened to w<=140, exact-confirmed for
any candidate below 1/1260), AND a coarse 3-element probe. Report anything <= 1/1260.

If nothing beats 1/1260, this independently corroborates HYP-2561 over a much larger
search. Uses fast FLOAT interval-union to screen, exact Fraction to confirm.
(Interval union, NOT grid counting: grid destroys the thin lonely arcs.)
"""
from fractions import Fraction as F
from itertools import combinations

# ---------- fast float lonely measure (interval union) ----------
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort()
    tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot

def danger(v):
    out=[]
    for k in range(v+1):
        lo=max(F(14*k-1,14*v),F(0)); hi=min(F(14*k+1,14*v),F(1))
        if lo<hi: out.append((lo,hi))
    return out
def L_exact(S):
    iv=[]
    for v in set(S): iv+=danger(v)
    iv.sort(); t=F(0); cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else: t+=ch-cl; cl,ch=lo,hi
    return F(1)-(t+ch-cl)

BASE=list(range(1,14))
TARGET=1/1260.0
EPS=1e-9
hits=[]   # (Lfloat, S) with L < 1/1260 + eps

print("="*74)
print("Tight AP {1..13}: L =", float(L_exact(BASE)), " (exactly 0 -> extremizers are perturbations)")
print("Reference inf candidate 1/1260 =", TARGET)
print("="*74)

# ---- single-element perturbations: replace one element by w ----
print("\n[1] SINGLE perturbations  e -> w,  e in 1..13,  w in 14..600")
best_single=(1.0,None)
for e in BASE:
    for w in range(14,601):
        if w in BASE and w!=e: continue
        S=[x for x in BASE if x!=e]+[w]
        if len(set(S))!=13: continue
        Lf=L_float(S)
        if Lf<best_single[0]: best_single=(Lf,(e,w,tuple(sorted(S))))
        if Lf<TARGET+EPS: hits.append((Lf,tuple(sorted(S)),f"single {e}->{w}"))
print("   best single:", best_single[0], "config", best_single[1])

# ---- two-element perturbations: replace e1,e2 by w1,w2 (float screen) ----
print("\n[2] TWO-element perturbations  {e1,e2}->{w1,w2},  w in 14..140 (float screen)")
best_double=(1.0,None)
WR=range(14,141)
for e1,e2 in combinations(BASE,2):
    rest=[x for x in BASE if x not in (e1,e2)]
    for w1,w2 in combinations(WR,2):
        if w1 in rest or w2 in rest: continue
        S=rest+[w1,w2]
        Lf=L_float(S)
        if Lf<best_double[0]: best_double=(Lf,(e1,e2,w1,w2,tuple(sorted(S))))
        if Lf<TARGET+EPS: hits.append((Lf,tuple(sorted(S)),f"double {e1},{e2}->{w1},{w2}"))
print("   best double:", best_double[0], "config", best_double[1])

# ---- coarse three-element probe (structured: multiples) ----
print("\n[3] THREE-element probe (replace 3 by multiples k*e, k in 2..6)")
best_triple=(1.0,None)
import itertools
for trip in combinations(BASE,3):
    rest=[x for x in BASE if x not in trip]
    for ks in itertools.product(range(2,7),repeat=3):
        repl=[trip[i]*ks[i] for i in range(3)]
        S=rest+repl
        if len(set(S))!=13: continue
        Lf=L_float(S)
        if Lf<best_triple[0]: best_triple=(Lf,(trip,ks,tuple(sorted(S))))
        if Lf<TARGET+EPS: hits.append((Lf,tuple(sorted(S)),f"triple {trip}*{ks}"))
print("   best triple(multiples):", best_triple[0], "config", best_triple[1])

# ---- exact-confirm all sub-target hits ----
print("\n" + "="*74)
print("EXACT confirmation of configs with L < 1/1260 + eps:")
print("="*74)
if not hits:
    print("  NONE found below 1/1260 across single(w<=600)+double(w<=140)+triple-multiples.")
    print("  => independently CORROBORATES kind-pasteur HYP-2561 (inf = 1/1260) over a")
    print("     much larger search space.")
else:
    seen=set(); strictly_below=[]
    for Lf,S,desc in sorted(hits):
        if S in seen: continue
        seen.add(S)
        Le=L_exact(list(S))
        rel = "<" if Le<F(1,1260) else ("=" if Le==F(1,1260) else ">")
        print(f"  L={Le} = {float(Le):.9f} {rel} 1/1260   {desc}  S={list(S)}")
        if Le<F(1,1260): strictly_below.append((Le,S,desc))
    if strictly_below:
        print("\n  *** FOUND configs STRICTLY BELOW 1/1260 — would REFUTE HYP-2561: ***")
        for Le,S,desc in sorted(strictly_below)[:5]:
            print(f"     L={Le}={float(Le):.9f}  {desc}  S={list(S)}")
    else:
        print("\n  All hits EQUAL 1/1260 (ties), none strictly below => HYP-2561 corroborated;")
        print("  the tight locus (minimizers) may be larger than one config.")
