#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the resonance expansion of the covering deficit, and "resonances die =>
witness at large q".  kth-pasteur-2026-06-14-S6.  Advances HYP-2489 (the deficit
as a circle-method / character sum) and THM-497 D (the over-correlated regime).

EXACT EXPANSION (additive characters mod q).  Deficit over the full group
  D(q,S) = #{ a in Z/q : v*a mod q not in B_q for all v in S },  B_q = +-{0..h}, h=floor(q/14).
Writing 1 - 1_{B_q}(va) = (1-beta) - sum_{t!=0} chat(t) e_q(t v a), beta=(2h+1)/q,
chat(t) = (1/q) sum_{x=-h}^h e_q(-tx) = (1/q) Dirichlet_h(t) (real), and summing the
product over a:
  D(q,S)/q = (1-beta)^13
           + sum_{T subset S, T nonempty} (1-beta)^{13-|T|} (-1)^{|T|}
                 sum_{(t_v)_{v in T}, t_v!=0, sum_{v in T} t_v v ≡ 0 (mod q)}  prod_{v in T} chat(t_v).
So the DEVIATION from the independence main term (1-beta)^13 is a sum over ADDITIVE
RESONANCES sum_{v in T} t_v v ≡ 0 (mod q) of the speeds, weighted by Dirichlet
kernels.  KEY: for fixed S, a relation sum t_v v = m != 0 (integer) only resonates
at q | m, finitely many q; so as q grows the resonances vanish and D(q,S)/q ->
(1-beta)^13 ~ (6/7)^13 ~ 0.1348 > 0  =>  a witness EXISTS at every large enough q
=> the config is loose.  The threshold is governed by the speed magnitudes / the
integer linear relations among the speeds (the Sidon/B_h structure, THM-446).

VERIFY (exact integer):
 (A) D/q -> (1-beta)^13 as q grows, for evaders and generic configs (resonances die).
 (B) the deviation is concentrated at q that DIVIDE small speed-relations (resonant q).
 (C) threshold scaling: ladder height (first witness q) vs max(speed) under dilation.
"""

import sys, itertools
from math import gcd, sin, pi
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def band(q, n=14):
    h = q // n
    return set(r % q for r in range(-h, h + 1)), h


def deficit_full(S, q, n=14):
    B, h = band(q, n)
    return sum(1 for a in range(q) if all((v * a) % q not in B for v in S))


def deficit_units(S, q, n=14):
    B, h = band(q, n)
    return sum(1 for a in range(q) if gcd(a, q) == 1 and all((v * a) % q not in B for v in S))


def main_term(q, n=14):
    h = q // n
    beta = (2 * h + 1) / q
    return (1 - beta) ** 13, beta


def five_evaders():
    base = [7 * k for k in range(1, 13)]
    return {r: sorted(base + [r]) for r in (611, 702, 793, 962, 1053)}


def first_witness(S, Hmax=4000, n=14):
    for q in range(2, Hmax + 1):
        if deficit_full(S, q, n) > 0:
            return q
    return None


def part_A():
    print("=== (A) D(q,S)/q -> (1-beta)^13 ~ 0.1348 as q grows (resonances die => witness) ===", flush=True)
    cfgs = {'evader r=611': five_evaders()[611],
            'generic A': sorted([3,11,17,22,29,34,41,46,53,58,14,67,73]),
            'dilated 7*{1..12}+13 (=evader-ish)': sorted([7*k for k in range(1,13)]+[13])}
    for name, S in cfgs.items():
        if not any(v % 14 == 0 for v in S):
            # ensure a multiple of 14 present for C'(14) relevance (else still a valid speed set)
            pass
        row = []
        for q in (50, 100, 200, 500, 1000, 2000):
            D = deficit_full(S, q)
            mt, beta = main_term(q)
            row.append(f"q={q}:D/q={D/q:.3f}(mt {mt:.3f})")
        print(f"   {name}: " + "  ".join(row), flush=True)
        print(f"      first witness shell (ladder height) = {first_witness(S)}", flush=True)


def part_B():
    print("\n=== (B) the deviation lives at RESONANT q (q | a small speed-relation) ===", flush=True)
    S = five_evaders()[611]
    mt_inf = (6/7)**13
    print(f"   evader r=611, speeds {S}", flush=True)
    print(f"   q where D/q is FAR below main term (large negative deviation = strong resonance):", flush=True)
    devs = []
    for q in range(28, 400):
        D = deficit_full(S, q); mt, _ = main_term(q)
        devs.append((D/q - mt, q, D))
    devs.sort()
    for dev, q, D in devs[:8]:
        # which small relations resonate: gcd-structure / does q divide pairwise diffs?
        small_div = [q for v in S for w in S if v!=w and (v-w)%q==0]
        print(f"      q={q}: D={D} (D/q={D/q:.3f}, dev={dev:+.3f}); q divides a speed-difference: {len(small_div)>0}", flush=True)


def part_C():
    print("\n=== (C) threshold scaling: ladder height vs max(speed) under dilation+shift ===", flush=True)
    print("   (relation sum t_v v = m resonates only at q|m, m ~ speed magnitude; so the", flush=True)
    print("    first-witness shell should NOT blow up with pure dilation — M is dilation-invariant)", flush=True)
    base = sorted([7*k for k in range(1,13)]+[611])
    for c in (1, 2, 3, 5):
        S = sorted(set(c*v for v in base))
        h = first_witness(S)
        print(f"      {c}*evader (max={max(S)}): ladder height = {h}  (dilation-invariant M => same loose/tight)", flush=True)
    # genuinely larger non-dilated configs
    print("   non-dilated growth (evader family r climbing):", flush=True)
    for r in (611, 1093, 2000, 5000):
        S = sorted([7*k for k in range(1,13)]+[r])
        if len(set(S))!=13: continue
        h = first_witness(S, Hmax=200)
        print(f"      r={r} (max={max(S)}): ladder height = {h}", flush=True)


if __name__ == "__main__":
    part_A(); part_B(); part_C()
