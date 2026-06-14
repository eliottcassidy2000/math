#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The LRC singular series L(S) = lim_{q->inf} D(q,S)/q, and the reframe C'(14) <=> L(S)>0.
kind-pasteur-2026-06-14-S6.  Advances HYP-2489 with the concrete circle-method object.

D(q,S) = #{a in Z/q : v a mod q not in B_q for all v in S}, B_q=+-{0..floor(q/14)}.
Additive-character expansion =>  D(q,S)/q = (6/7-ish main) + sum over additive
resonances of the speeds.  As q->inf the resonances from NON-zero relations die
(q|m), leaving the EXACT integer relations sum t_v v = 0; the Dirichlet coefficient
chat(t)=(1/q)sum_{x=-h}^h e_q(-tx) -> s(t)=sin(pi t/7)/(pi t).  So the SINGULAR SERIES

   L(S) = lim_q D(q,S)/q   exists and is controlled by the speed relation lattice.

TESTS (exact deficit, large-q averaging to estimate L):
 (A) L(S) exists (D/q stabilizes): average D(q,S)/q over a high-q window for many S;
     report the limit + its stability (std over the window).
 (B) the L-SPECTRUM by structure: generic/Sidon (L~0.135), structured/AP (L lower),
     near-tight (L~0).  Is L(S) > 0 for ALL primitive multiple-of-14 S (=C'(14))?
     -> min L over a large sample; is it bounded away from 0?
 (C) the threshold q*(S) (ladder height) vs L(S): small L => high threshold.
 (D) Euler-product test (HYP-2503): does L(S) ~ prod over small primes of a local
     density?  Probe via the per-prime residue structure of S.
"""

import sys, itertools, random
from math import gcd, sin, pi
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def band_set(q, n=14):
    h = q // n
    return set(r % q for r in range(-h, h + 1))


def deficit(S, q, n=14):
    B = band_set(q, n)
    return sum(1 for a in range(q) if all((v * a) % q not in B for v in S))


def L_estimate(S, qlo=1500, qhi=1560, n=14):
    """average D(q,S)/q over a high-q window -> estimate of the limit L(S)."""
    vals = [deficit(S, q, n) / q for q in range(qlo, qhi)]
    m = sum(vals) / len(vals)
    var = sum((x - m) ** 2 for x in vals) / len(vals)
    return m, var ** 0.5


def first_witness(S, Hmax=200, n=14):
    for q in range(2, Hmax + 1):
        if deficit(S, q, n) > 0:
            return q
    return None


def primitive_mult14(S):
    from functools import reduce
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)
def dominant(S, n=14):
    s = sorted(S); return s[-1] > (n - 1) * s[-2]


def main():
    print("=== (A)+(B) the singular series L(S) = lim D(q,S)/q, by structure ===", flush=True)
    base = [7 * k for k in range(1, 13)]
    cfgs = {
        'generic/Sidon-ish': sorted([1,2,4,8,16,32,64,5,11,23,47,95,14]),
        'evader r=611 (AP core)': sorted(base + [611]),
        'evader r=702': sorted(base + [702]),
        'AP 1..13 *? -> 14*{1..13} (non-prim tight x14)': sorted([14*k for k in range(1,14)]),
        'near-regular-ish + mult14': sorted([1,2,3,4,5,6,7,14,21,28,9,11,13]),
    }
    print("   config: L(S) (limit of D/q) +- window-std   ladder-height   prim-mult14? dominant?", flush=True)
    for name, S in cfgs.items():
        if len(set(S)) != 13:
            print(f"   {name}: |S|={len(set(S))} != 13, skip", flush=True); continue
        L, sd = L_estimate(S)
        h = first_witness(S)
        print(f"   {name[:34]:34s}: L={L:.4f} +-{sd:.4f}   h={h}   prim14={primitive_mult14(S)} dom={dominant(S)}", flush=True)

    print("\n=== (B) min L(S) over a large primitive-non-dominant multiple-of-14 sample (is L>0 always?) ===", flush=True)
    rng = random.Random(20260614)
    Ls = []
    tries = 0
    while len(Ls) < 120 and tries < 30000:
        tries += 1
        S = sorted(rng.sample(range(1, 100), 13))
        if not primitive_mult14(S) or dominant(S):
            continue
        L, sd = L_estimate(S, qlo=1500, qhi=1530)
        Ls.append((L, S))
    Ls.sort()
    print(f"   {len(Ls)} configs: L ranges [{Ls[0][0]:.4f}, {Ls[-1][0]:.4f}], median {Ls[len(Ls)//2][0]:.4f}", flush=True)
    print(f"   MIN L config: L={Ls[0][0]:.4f}  S={Ls[0][1]}  (L>0 => loose; is min bounded away from 0?)", flush=True)
    print(f"   #configs with L < 0.02 (near-tight, the hard ones): {sum(1 for L,_ in Ls if L<0.02)}", flush=True)

    print("\n=== (C) threshold q*(S) vs L(S): small singular series => high ladder ===", flush=True)
    for L, S in Ls[:3] + Ls[-3:]:
        h = first_witness(S)
        print(f"   L={L:.4f}: ladder height = {h}  S={S}", flush=True)


if __name__ == "__main__":
    main()
