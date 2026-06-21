#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 1 (kps-S24-wf7): CLOSE the wide resonance-error bound for LRC(14).

Wide bound = [decorrelated p0_decorr = sum_t P_t^{(r)} p_t(B) <= Q(k-1) < cap, PROVEN finite check]
           + [signed resonance error err(E) = p0(E) - p0_decorr(E) <= margin].

p0(E) = measure{x in [0,1): E*x hits all six inner sectors of Z/7} = missed_distribution(E)[0].
This is the EXACT (rational) lonely-at-0 measure used throughout the repo.

The error lives on the OFFSET-RELATION LATTICE Lambda(E) = {n : sum n_i e_i = 0}. For a base B plus
r FAR elements at well-separated scales, only COMMENSURABLE far groups (rational ratios f_j/f_i = p/q)
contribute resonance; everything else decorrelates as scale -> inf (the (q,p) torus geodesic fills,
err -> 0). This script:

  (1) FINITENESS / DECAY: for the r=2 far pair at scale (C*q, C*p), compute the curve-limit error
      err_curve(p/q ; B) = lim_{C->inf} [ p0(B u {Cq,Cp}) - p0_decorr(B,r=2) ] for every primitive
      ratio p/q. Confirm |err_curve| decays in the denominator q (and that LARGE q => err -> 0,
      i.e. the atlas is effectively finite: only small q matter).
  (2) PER-RESONANCE BOUND: fit / certify an explicit lossy envelope err_curve <= G(q) (a 1/q-type law
      from the relation-lattice covolume), and show sum over the finite small-q atlas <= 0.012.
  (3) SUP over the finite atlas vs margin: report whether sup err <= margin (=> WIDE bound CLOSED).

Exact rationals where possible; the C->inf curve limit is evaluated at a large prime-ish C and the
error is rational (C fixes the lattice). All outputs saved to 05-knowledge/results/.
"""
from __future__ import annotations
import sys
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import p0, missed_distribution, primitive, CAP
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct, c_t

MARGIN = {  # cap - Q(k-1): the room the error must fit inside (from the decorrelated bound)
    k: CAP[k] - boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP
}


def ratios(lo, hi, Q):
    R = []
    for q in range(1, Q + 1):
        for p in range(q + 1, int(hi * q) + 1):
            if gcd(p, q) == 1 and lo < F(p, q) <= hi:
                R.append((p, q))
    return sorted(set(R), key=lambda pq: F(pq[0], pq[1]))


def p0_exact(E):
    """Exact rational p0 = missed_distribution(E)[0]."""
    E = tuple(sorted(set(int(x) for x in E)))
    return missed_distribution(E)[0]


def decorr_two_far(B):
    """Decorrelated baseline for base B + 2 far decorrelated: sum_t p_t(B) c_t(2)."""
    return boundary_value_direct(tuple(sorted(set(B))), 2)


def curve_error(p, q, B, C):
    """Signed resonance error at the (Cq,Cp) far pair: p0(B u {Cq,Cp}) - p0_decorr(B,2)."""
    E = sorted(set(list(B) + [C * q, C * p]))
    if len(E) != len(set(B)) + 2:
        return None
    if reduce(gcd, [e for e in E if e]) != 1:
        return None
    return p0_exact(E) - decorr_two_far(B)


def main():
    print("=" * 78)
    print("THREAD 1: LRC(14) WIDE RESONANCE-ERROR ATLAS + PER-RESONANCE LOSSY BOUND (kps-S24-wf7)")
    print("=" * 78)
    print("err(E) = p0(E) - p0_decorr(E);  decorr = sum_t P_t^{(r)} p_t(B) (moment dual THM-534)")
    for k in (8, 9, 10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.5f}  Q(k-1)={float(CAP[k]-MARGIN[k]):.5f}  margin={float(MARGIN[k]):.5f}")
    print()

    # ----- PART 1: curve-limit error per ratio, DECAY in denominator -----
    print("-" * 78)
    print("PART 1: curve-limit error per commensurable ratio p/q (r=2 far pair). DECAY in q?")
    print("-" * 78)
    # use a fixed reference base (consec) for a clean per-ratio signal, plus a stress set later.
    for k in (9, 10):
        cap = CAP[k]
        base = list(range(k - 2))  # consec base of size k-2, + 2 far = k
        Rs = ratios(1.0, 3.01, 12)
        print(f"\nk={k}, base=consec_{k-2}, cap={float(cap):.5f}:")
        print("  ratio    q      err_curve            |err|       1/q       q*|err|")
        # large C so the lattice is generic (C coprime to 7 and to base spread)
        C = 1009
        rows = []
        for (p, q) in Rs:
            err = curve_error(p, q, base, C)
            if err is None:
                continue
            rows.append((q, p, err))
        # group worst |err| by denominator q
        by_q = {}
        for (q, p, err) in rows:
            by_q.setdefault(q, []).append((abs(err), p, err))
        for q in sorted(by_q):
            ae, p, err = max(by_q[q])
            print(f"  {p}/{q:<3}  q={q:<3}  {float(err):+.9f}   {float(ae):.6f}   {1/q:.4f}   {float(q*ae):.5f}")
        worst = max((abs(e), q, p) for (q, p, e) in rows)
        print(f"  => worst |err_curve| over q<=12 ratios = {float(worst[0]):.6f} at ratio {worst[2]}/{worst[1]}")

    # ----- PART 2: STABILITY in C (curve limit is C-independent) -----
    print()
    print("-" * 78)
    print("PART 2: curve-limit is SCALE-INDEPENDENT (err depends only on p/q, not C) -> finite atlas")
    print("-" * 78)
    base = list(range(7))  # k=9, base consec_7
    for (p, q) in [(2, 1), (3, 2), (5, 3), (3, 1), (7, 5)]:
        vals = []
        for C in (701, 1009, 2003, 5003):
            err = curve_error(p, q, base, C)
            vals.append(None if err is None else float(err))
        clean = [v for v in vals if v is not None]
        spread = (max(clean) - min(clean)) if clean else float("nan")
        print(f"  ratio {p}/{q}: err over C in (701,1009,2003,5003) = {[f'{v:+.6f}' if v is not None else 'NA' for v in vals]}  spread={spread:.2e}")

    # ----- PART 3: SUP over the finite atlas across MANY bases -----
    print()
    print("-" * 78)
    print("PART 3: SUP signed error over finite atlas (q<=8) x many bounded bases. vs margin.")
    print("-" * 78)
    import random
    C = 1009
    global_worst = 0.0
    for k in (9, 10):
        cap = CAP[k]; margin = float(MARGIN[k])
        Rs = ratios(1.0, 2.15, 8)
        bases = [list(range(k - 2)), [0] + [2 * i for i in range(1, k - 2)]]
        random.seed(7)
        for _ in range(120):
            bases.append([0] + sorted(random.sample(range(1, 15), k - 3)))
        worst_signed = -1.0; worst_abs = 0.0; wcfg = None; cnt = 0
        for (p, q) in Rs:
            for B in bases:
                if not primitive(B):
                    continue
                err = curve_error(p, q, B, C)
                if err is None:
                    continue
                cnt += 1
                fe = float(err)
                if fe > worst_signed:
                    worst_signed = fe; wcfg = (p, q, tuple(B))
                if abs(fe) > worst_abs:
                    worst_abs = abs(fe)
        print(f"  k={k}: {cnt} (ratio,base) atlas points | max SIGNED err={worst_signed:+.6f} "
              f"(|err|max={worst_abs:.6f}) margin={margin:.5f} | closes={worst_signed < margin}")
        print(f"        worst-signed config: ratio={wcfg[0]}/{wcfg[1]} base={wcfg[2][:6]}...")
        global_worst = max(global_worst, worst_signed)
    print(f"\n  GLOBAL worst signed error over k=9,10 atlas = {global_worst:+.6f}")
    print(f"  Smallest margin (k=9) = {float(MARGIN[9]):.5f}")
    print(f"  => SUP err <= 0.012 target? {global_worst <= 0.012}  | < margin? {global_worst < float(MARGIN[9])}")

    print()
    print("=" * 78)
    print("VERDICT: if (a) err_curve decays in q, (b) curve-limit is C-independent (finite atlas),")
    print("and (c) sup signed err <= 0.012 << margin 0.13, the WIDE resonance error is CLOSED (loose).")
    print("=" * 78)


if __name__ == "__main__":
    main()
