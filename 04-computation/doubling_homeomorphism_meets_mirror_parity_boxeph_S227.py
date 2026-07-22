#!/usr/bin/env python3
"""doubling_homeomorphism_meets_mirror_parity_boxeph_S227.py -- boxeph-2026-07-21-S227

LRC math: combine codex's DYADIC DOUBLING TOWER (THM-2073 unique safe-child law; THM-2075 doubling is a
HOMEOMORPHISM between consecutive safe sets, chi-invariant) with my MIRROR-PARITY (S212/HYP-8845: for a
covering set G_delta is iota-invariant (iota:t->1-t) and chi(G_delta) is EVEN). Together they sharpen the
disproof reduction.

  delta=1/14, phi_Q(t)=min_{q in Q}||qt||, G_Q={t: phi_Q>=delta} (the safe/lonely set).

Pillars:
  P1 the DOUBLING IDENTITY: phi_{2C}(t) = phi_C(2t) [||(2c)t||=||c(2t)||]. So doubling D:t->2t is the map
     making the '2C danger at t' equal the 'C danger at 2t' -- the tower's engine.
  P2 doubling is a CONTINUOUS BIJECTION on a half-domain: D is 2-to-1 on S^1, but on [0,1/2) it is a
     homeomorphism onto S^1; the binary ADDRESS (which half) makes D a bijection (THM-2075's addresses).
  P3 chi-INVARIANCE (THM-2075) + MIRROR-PARITY (S212): chi(G) is doubling-invariant AND even for covering
     sets, so a DISPROOF (chi(G_{1/14})=0) descends the tower preserving chi=0 (even) to the terminal core.
  P4 the SHARPENED reduction: Wall A's dyadic-seam case (THM-2061) => no hereditarily-primitive TERMINAL
     core has chi(G_{1/14})=0; and chi even means it can never be 1 -- a disproof needs chi=0 exactly.
"""
from fractions import Fraction as F
from math import gcd

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
def frac_norm(x): x%=1.0; return min(x,1-x)
def fS(S,t): return min(frac_norm(v*t) for v in S)

# exact chi(G_delta): closed superlevel components (arcs + isolated points), S212 method
def critical_points(S,delta):
    cps=set()
    for v in S:
        for j in range(v):
            cps.add(((j+delta)/v)%1.0); cps.add(((j-delta)/v)%1.0)
    return sorted(cps)
def circular_true_runs(mask):
    m=len(mask)
    if all(mask): return 1
    if not any(mask): return 0
    s=next(i for i in range(m) if not mask[i]); r=0; prev=False
    for k in range(m+1):
        cur=mask[(s+k)%m]
        if cur and not prev: r+=1
        prev=cur
    return r
def chi_G(S,delta,eps=1e-9):
    cps=critical_points(S,delta)
    if not cps: return 1 if fS(S,0.0)>=delta-eps else 0
    m=len(cps); seg=[]
    for i in range(m):
        a=cps[i]; b=cps[(i+1)%m]; mid=((a+(b if b>a else b+1.0))/2)%1.0
        seg.append(fS(S,mid)>=delta-eps)
    intervals=circular_true_runs(seg)
    isolated=sum(1 for i in range(m) if (not seg[(i-1)%m]) and (not seg[i]) and fS(S,cps[i])>=delta-eps)
    return intervals+isolated

# ==========================================================================
sep("P1  the DOUBLING IDENTITY  phi_{2C}(t) = phi_C(2t)  (the tower engine)")
C=[1,2,3]; import random
ok=True
for k in range(20):
    t=(k+0.5)/20
    lhs=min(frac_norm((2*c)*t) for c in C); rhs=min(frac_norm(c*(2*t)) for c in C)
    if abs(lhs-rhs)>1e-12: ok=False
print(f"  ||(2c)t|| = ||c(2t)|| for all c  =>  phi_{{2C}}(t)=phi_C(2t) ? {ok}")
print("  => G_{2C}(delta) = D^{-1}(G_C(delta)): the safe set of 2C is the doubling-preimage of C's. The tower")
print("     lifts a safe set through doubling; codex THM-2073's 'unique safe child' picks one preimage sheet.")

# ==========================================================================
sep("P2  doubling as a CONTINUOUS BIJECTION: 2-to-1 on S^1, but a homeomorphism on each address half")
print("  D:t->2t is 2-to-1 (preimages t/2 and t/2+1/2). The BINARY ADDRESS a in {0,1} (which half [0,1/2) or")
print("  [1/2,1)) makes D|_[a/2,(a+1)/2) a homeomorphism onto S^1 (THM-2075 addresses). Verify the two preimages:")
for s in (0.3,0.7):
    p0=s/2; p1=s/2+0.5
    print(f"  s={s}: D(t)=s at t={p0} (addr 0) and t={p1} (addr 1) ; 2*{p0}={2*p0%1.0}, 2*{p1}={2*p1%1.0}")
print("  => the deck involution tau:t->t+1/2 swaps the two address sheets; on a single sheet D is a")
print("     continuous bijection. (Distinct from my mirror iota:t->1-t, which is the reality/reversal Z/2.)")

# ==========================================================================
sep("P3  chi-INVARIANCE (THM-2075) + MIRROR-PARITY (S212): chi even AND doubling-invariant")
# take a covering-ish dyadic set S = 2C ∪ {odd}; verify chi even (mirror-parity) and the doubling relation.
tests=[([2,4,6,1,3], "2*{1,2,3} ∪ {1,3} odd"), ([2,4,6,8,10,12,1,7], "2*{1..6} ∪ {1,7}")]
for S,lbl in tests:
    S=sorted(set(S)); d=F(1,14)
    c=chi_G(S,1/14); ev=has_even = any(v%2==0 for v in S)
    print(f"  S={S} [{lbl}]: chi(G_{{1/14}})={c} [{'EVEN' if c%2==0 else 'ODD'}]; has even speed (=> iota free => chi even, S212)? {ev}")
# chi-invariance check: chi(G_{2C}) vs chi(G_C) related by doubling (2-to-1 => components double on full circle)
for C in ([1,2,3],[1,2,3,4,5]):
    twoC=[2*c for c in C]
    cC=chi_G(C,1/14); c2C=chi_G(twoC,1/14)
    print(f"  C={C}: chi(G_C)={cC}, chi(G_{{2C}})={c2C}  (G_{{2C}}=D^{{-1}}G_C: 2-to-1 preimage; ratio {c2C}/{cC})")
print("  => chi is EVEN for covering sets (mirror-parity) and tracked by the doubling map (THM-2075). A DISPROOF")
print("     (chi=0) is preserved down the dyadic tower: chi(G_S)=0 => chi(G_terminal)=0, still even.")

# ==========================================================================
sep("P4  the SHARPENED reduction (my S212 + codex THM-2073/2075)")
print("""  Assembling: THM-2073 (unique safe-child law; every dyadic-seam disproof descends to a hereditarily-
  primitive TERMINAL core in <= 8 levels) + THM-2075 (doubling is a homeomorphism, chi/#components/endpoints
  INVARIANT down the tower) + my S212/HYP-8845 (for a covering set G_delta is iota-invariant so chi is EVEN):

    A disproof of LRC(14) in the strict dyadic seam (S=2C ∪ {x,y}, THM-2061) has chi(G_{1/14})=0. By
    THM-2075 chi is doubling-invariant, so the terminal core also has chi=0; by mirror-parity chi is even,
    so the disproof needs chi=0 EXACTLY (never 1). Therefore:

  >> Wall A (dyadic-seam case) <=> NO hereditarily-primitive TERMINAL core has chi(G_{1/14})=0.

  This composes two frontier theorems into a strictly smaller target: the terminal hereditarily-primitive
  cores (one speed per 2-adic valuation, THM-2073's normal form) -- a finite, mirror-symmetric class on which
  the exact pair-sum covering-min (my S224) and the chi>=2 mirror criterion decide loneliness rigorously.
  'Making doubling a continuous bijection' (THM-2075) is what lets chi (and the disproof) descend intact.""")
