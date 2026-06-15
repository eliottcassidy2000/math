#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Improving the Paley-Zygmund signed-short-vector route to the LRC(14) singular series.
kind-pasteur-2026-06-14-S2.  Builds on THM-501 (L(S), kind-pasteur S6) + THM-503
(structure + almost-Sidon loose class, mac-mini S1).

THE OBJECT (THM-501):
  L(S) = (6/7)^13 + Σ_{exact relations Σ t_v v = 0, t_v≠0} (6/7)^{13−|T|}(−1)^{|T|} ∏ s(t_v),
  s(t) = sin(πt/7)/(πt).   L>0 ⟹ loose.  Open prize: inf L>0 over the dilated cores
  d·{1..12}∪{r}, which carry abundant |T|=3 relations where THM-503's ABSOLUTE bound
  (almost-Sidon) fails.

THE PALEY-ZYGMUND IMPROVEMENT (HYP-2510):  the |T|=3 correction is a SIGNED sum
  Σ_3 = (6/7)^{10}·(−1)^3·Σ_{exact triples} ∏ s(t_v).
THM-503 bounds |Σ_3| ≤ A_3 := (6/7)^{10} Σ|∏s| (triangle inequality = no cancellation).
If the signs sign(∏s) = ∏ sign(sin(πt_v/7)) behave quasi-randomly, the TRUE |Σ_3|
should show SQUARE-ROOT cancellation ≈ E_3 := (6/7)^{10} √(Σ(∏s)²) ≪ A_3, extending
the proved-loose class.  THE CRUX (two-sided, MISTAKE-062): do the STRUCTURED cores'
|T|=3 relations have ALIGNED signs (no cancellation, PZ fails exactly on the hard
cores) or CANCELLING signs (PZ wins)?  Compute and let the data decide.

OUTPUT per config: baseline (6/7)^13; the |T|=3 signed sum Σ_3, absolute A_3,
L2-bound E_3; the cancellation ratio |Σ_3|/A_3 (1 = aligned/no cancel, →0 = strong);
and the full small-relation estimate of L vs the deficit-density ground truth.
"""

import sys, itertools, math
from math import gcd, sin, pi, sqrt
from functools import reduce
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

BASE = (6/7)**13

def s(t):
    if t % 7 == 0:  # 7-vanishing (THM-503 (1))
        return 0.0
    return sin(pi*t/7.0)/(pi*t)

def relations_T3(S, Tmax=6):
    """7-primitive exact relations t_a*a+t_b*b+t_c*c=0, support exactly 3,
    |t_*|<=Tmax, gcd of coeffs =1 (primitive), no t_v≡0 mod7. Returns list of
    (triple_of_speeds, (t_a,t_b,t_c))."""
    rels = []
    n = len(S)
    for ia, ib, ic in itertools.combinations(range(n), 3):
        a, b, c = S[ia], S[ib], S[ic]
        seen = set()
        for ta in range(-Tmax, Tmax+1):
            if ta == 0 or ta % 7 == 0: continue
            for tb in range(-Tmax, Tmax+1):
                if tb == 0 or tb % 7 == 0: continue
                num = -(ta*a + tb*b)
                if num % c != 0: continue
                tc = num // c
                if tc == 0 or tc % 7 == 0 or abs(tc) > Tmax: continue
                g = reduce(gcd, (abs(ta), abs(tb), abs(tc)))
                key = (ta//g, tb//g, tc//g)
                # only primitive, and dedup +/- (the term and its negation both appear;
                # keep both since the singular series sums over ALL nonzero (t_v), but
                # primitive reps with a canonical sign avoid double-count of g-multiples)
                if g != 1: continue
                if key in seen or tuple(-x for x in key) in seen: continue
                seen.add(key)
                rels.append(((a,b,c), key))
    return rels

def T3_sums(S, Tmax=6):
    """signed Σ∏s, absolute Σ|∏s|, L2 sqrt(Σ(∏s)^2) over |T|=3 exact relations.
    s is EVEN (s(−t)=s(t)), so a relation and its negation give the SAME ∏s and
    ADD (no ±cancellation). The cancellation question is whether sign(∏s) varies
    ACROSS DIFFERENT primitive relations (the t_v mod 14 sign structure: s(t)>0 on
    t mod14∈{1..6}, <0 on {8..13}). We sum one representative per ± class.
    cancellation ratio |signed|/abs: 1 = all same sign (NO cancel, almost-Sidon
    bound tight, PZ fails); →0 = strong √-cancellation (PZ improves)."""
    signed = 0.0; abss = 0.0; sq = 0.0; npair = 0
    for (sp, key) in relations_T3(S, Tmax):
        p = s(key[0])*s(key[1])*s(key[2])
        if p == 0: continue
        signed += p
        abss += abs(p)
        sq += p*p
        npair += 1
    return signed, abss, sqrt(sq), npair

def L_estimate(S, Tmax=6, maxT=4):
    """estimate L(S) by summing exact relations up to support maxT, coeffs<=Tmax."""
    L = BASE
    n = len(S)
    for k in range(2, maxT+1):
        for combo in itertools.combinations(range(n), k):
            sp = [S[i] for i in combo]
            # enumerate (t_v) with Σ t_v sp = 0, |t|<=Tmax, all nonzero, 7-primitive
            ranges = [range(-Tmax, Tmax+1)]*(k-1)
            for pre in itertools.product(*ranges):
                if any(t==0 or t%7==0 for t in pre): continue
                num = -sum(pre[j]*sp[j] for j in range(k-1))
                last = sp[-1]
                if num % last != 0: continue
                tl = num//last
                if tl==0 or tl%7==0 or abs(tl)>Tmax: continue
                tv = pre+(tl,)
                prod = 1.0
                for t in tv: prod *= s(t)
                if prod == 0: continue
                L += (6/7)**(13-k) * ((-1)**k) * prod
    return L

def main():
    print(f"baseline (6/7)^13 = {BASE:.5f}\n", flush=True)
    configs = {
        "generic-coprime {14,17,19,23,29,31,37,41,43,47,53,59,61}":
            [14,17,19,23,29,31,37,41,43,47,53,59,61],
        "evader 7*{1..12}+611": sorted([7*k for k in range(1,13)]+[611]),
        "evader 7*{1..12}+702": sorted([7*k for k in range(1,13)]+[702]),
        "d=3 core 3*{1..12}+182": sorted([3*k for k in range(1,13)]+[182]),
        "near-tight {1..12}+14": sorted(list(range(1,13))+[14]),
    }
    print("=== KEY TEST: |T|=3 signed sum vs absolute — is there sign cancellation? ===", flush=True)
    print("    NOTE: |t|<=6 forces s(t)>0 (pi*t/7 in (0,pi)); sign variation needs |t|>=8.", flush=True)
    print("    ratio=|signed|/abs: 1 => aligned (PZ FAILS); <<1 => cancellation (PZ IMPROVES).", flush=True)
    for Tmax in (6, 13):
        print(f"   --- Tmax={Tmax} ---", flush=True)
        for name, S in configs.items():
            signed, abss, l2, npair = T3_sums(S, Tmax=Tmax)
            ratio = abs(signed)/abss if abss > 0 else 0.0
            # split dominant (all small-coeff, |t|<=6) vs tail
            print(f"     {name[:38]:38s}: #rel={npair:4d}  |signed|={abs(signed):.4f}  "
                  f"abs={abss:.4f}  ratio={ratio:.3f}", flush=True)

    print("\n=== full small-relation L estimate (support<=4, coeff<=6) vs baseline ===", flush=True)
    for name, S in configs.items():
        L4 = L_estimate(S, Tmax=6, maxT=4)
        L3 = L_estimate(S, Tmax=6, maxT=3)
        L2 = L_estimate(S, Tmax=6, maxT=2)
        print(f"   {name[:42]:42s}: L(|T|<=2)={L2:.4f}  L(<=3)={L3:.4f}  L(<=4)={L4:.4f}", flush=True)

if __name__ == "__main__":
    main()
