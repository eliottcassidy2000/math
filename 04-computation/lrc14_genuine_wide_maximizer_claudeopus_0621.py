#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-21: GENUINE-WIDE EXACT MAXIMIZER + CLOSED-FORM HUNT (LRC14 leg C).

The live mathematical target (OPEN-Q-108 leg C): prove
    genuine-wide E  ==>  p0(E) < Q(k-1) < cap_k.
PROVED at k=8,9; only VERIFIED at k=10,11,12. To turn verification into proof we
need the EXACT genuine-wide maximizer and ideally a CLOSED FORM for max p0 over
genuine-wide configs, plus the exact gap Q(k-1) - max_p0(genuine-wide).

genuine-wide (operational, matches HYP-2788): span(E)>14 AND for EVERY e in E,
span(primitive(E\{e})) > 14  (removing any single element keeps it wide; so it is
NOT single-perturbation-bounded and THM-563's single-far period-max does not apply).

This script does an adversarial EXACT search over the empirically-binding family
(symmetric multi-clusters, two/three tight blocks, dilated-AP + perturbation) plus
random genuine-wide configs, for k=8..13, and reports:
  - max p0 over genuine-wide (exact rational), the argmax config
  - the exact gap Q(k-1) - max_p0  and  cap_k - max_p0
  - the cluster signature of the maximizer (block sizes / spacing) to expose closed form

Exact rationals throughout (p0_fast cross-checked == repo p0 elsewhere).
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP, QVAL, MARGIN
from lrc14_wide_branch_ridge_codex_s47 import primitive


def reprim(E):
    """primitive reduced form: divide out gcd."""
    E = tuple(sorted(set(int(x) for x in E)))
    g = reduce(gcd, E)
    if g > 1:
        E = tuple(x // g for x in E)
    return E


def span(E):
    return max(E) - min(E)


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if 0 not in E:
        return False
    if span(E) <= 14:
        return False
    if not primitive(E):
        return False
    # removing any single element must keep it wide after re-primitivizing
    for e in E:
        sub = reprim(tuple(x for x in E if x != e))
        if len(sub) < 2:
            return False
        if span(sub) <= 14:
            return False  # single-perturbation-bounded -> NOT genuine-wide
    return True


def cluster_signature(E):
    """Describe E as blocks: list of (start, length, internal_step) gaps between blocks."""
    E = sorted(E)
    blocks = []
    cur = [E[0]]
    for a, b in zip(E, E[1:]):
        if b - a == 1:
            cur.append(b)
        else:
            blocks.append(cur)
            cur = [b]
    blocks.append(cur)
    sig = []
    for blk in blocks:
        sig.append((blk[0], len(blk)))
    gaps = [blocks[i + 1][0] - blocks[i][-1] for i in range(len(blocks) - 1)]
    return sig, gaps


def gen_two_cluster(k, maxspan=60):
    """Two consecutive blocks of sizes (a, k-a), second block placed far."""
    out = []
    for a in range(1, k):
        b = k - a
        for gap in range(2, maxspan):
            # block1 = 0..a-1 ; block2 = start..start+b-1
            start = (a - 1) + gap
            E = tuple(list(range(a)) + [start + i for i in range(b)])
            out.append(E)
    return out


def gen_symmetric_two_cluster(k, maxspan=60):
    """Symmetric: 0..a-1 and (M-(b-1))..M for various M (the empirical slack maximizer)."""
    out = []
    for a in range(1, k):
        b = k - a
        for M in range(15, maxspan):
            E = tuple(list(range(a)) + [M - i for i in range(b)][::-1])
            E = tuple(sorted(set(E)))
            if len(E) == k:
                out.append(E)
    return out


def gen_three_cluster(k, maxspan=45):
    out = []
    sizes = []
    for a in range(1, k - 1):
        for bb in range(1, k - a):
            c = k - a - bb
            if c >= 1:
                sizes.append((a, bb, c))
    for (a, bb, c) in sizes:
        for g1 in range(2, 12):
            for g2 in range(2, 12):
                s1 = 0
                s2 = (a - 1) + g1
                s3 = s2 + (bb - 1) + g2
                E = tuple(list(range(a)) + [s2 + i for i in range(bb)] + [s3 + i for i in range(c)])
                E = tuple(sorted(set(E)))
                if len(E) == k and max(E) - min(E) <= maxspan:
                    out.append(E)
    return out


def gen_dilated_ap_perturb(k, maxspan=60):
    """d*consec_{k-1} + one off-lattice perturbation (the dilated-AP+1 family)."""
    out = []
    for d in (2, 3, 4, 5):
        base = [d * i for i in range(k - 1)]
        for p in range(1, d * (k - 1) + 12):
            if p in base:
                continue
            E = tuple(sorted(set([0] + base + [p])))
            if len(E) == k and max(E) <= maxspan:
                out.append(E)
    return out


def gen_random(k, n, rng, maxspan=70):
    out = []
    for _ in range(n):
        sp = rng.randint(16, maxspan)
        pts = sorted(rng.sample(range(1, sp + 1), k - 1))
        E = tuple([0] + pts)
        out.append(E)
    return out


def search(k, rng):
    cand = set()
    for gen in (gen_two_cluster, gen_symmetric_two_cluster, gen_three_cluster, gen_dilated_ap_perturb):
        for E in gen(k):
            cand.add(reprim(E))
    for E in gen_random(k, 4000, rng):
        cand.add(reprim(E))
    best = None
    best_E = None
    checked = 0
    for E in cand:
        if len(E) != k:
            continue
        if not is_genuine_wide(E):
            continue
        checked += 1
        v = p0_fast(E)
        if best is None or v > best:
            best = v
            best_E = E
    return best, best_E, checked, len(cand)


def main():
    rng = random.Random(20260621)
    print("=" * 78)
    print("GENUINE-WIDE EXACT MAXIMIZER + CLOSED-FORM HUNT  (LRC14 leg C)")
    print("claude-opus 2026-06-21")
    print("=" * 78)
    print(f"{'k':>2} {'max_p0(gw)':>16} {'~':>8} {'Q(k-1)':>16} {'gap=Q-max':>12} "
          f"{'cap-max':>10}  argmax")
    rows = []
    for k in range(8, 14):
        if k not in CAP:
            # extend CAP/QVAL? CAP only has 8..12; compute Q for 13 if possible
            pass
        best, best_E, checked, ncand = search(k, rng)
        if best is None:
            print(f"{k:>2}  (no genuine-wide found; cand={ncand})")
            continue
        q = QVAL.get(k)
        cap = CAP.get(k)
        sig, gaps = cluster_signature(best_E)
        gapstr = f"{float(q-best):+.5f}" if q is not None else "n/a"
        capstr = f"{float(cap-best):+.5f}" if cap is not None else "n/a"
        qstr = f"{float(q):.6f}" if q is not None else "n/a"
        print(f"{k:>2} {str(best):>16} {float(best):8.5f} {qstr:>16} {gapstr:>12} "
              f"{capstr:>10}  {best_E}")
        print(f"     blocks={sig} gaps={gaps}  (checked {checked} genuine-wide of {ncand} cand)")
        rows.append((k, best, best_E, q, cap))
    print("=" * 78)
    print("CLOSED-FORM PROBE: is max_p0(gw) a clean function of k? exact gaps Q-max:")
    for (k, best, best_E, q, cap) in rows:
        if q is not None:
            print(f"  k={k}: Q-max = {q-best} = {float(q-best):.6f} ;  max_p0 = {best}")
    print("=" * 78)
    print("INTERPRETATION: if the genuine-wide maximizer is a STABLE 2-cluster family")
    print("with p0 expressible in k, leg C closes as a closed-form inequality < Q(k-1).")


if __name__ == "__main__":
    main()
