#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C / FRONTIER (kps-wf9): ADVERSARIAL stress-test of the gK8 wide-closure conjecture
(HYP-2809, claude-opus-2026-06-22). gK8 = (10,0,0,1,0,0,10) on the miss-distribution q:
    L_yK8(E) = 10*q0 + q3 + 10*q6,   q_t = meas{exactly t of the 6 inner sectors missed}.
CONJECTURE: max over ALL wide E of L_yK8 <= 10*cap_k  (=> p0=q0 <= cap_k since q3,q6>=0).
opus VERIFIED it on a few hand-picked binding families with ~20-40% slack. THIS script
stress-tests it ADVERSARIALLY:
  (1) the SHARP genuine-wide maximizer (HYP-2805 dilated doublet {0,2..14,15,16} at k=10) and
      its k-analogues -- the configs that EDGE OUT the consec doublet and where the 0.16 margin
      FAILS. Does gK8 still hold there, and with how much slack?
  (2) a BROAD exhaustive/random genuine-wide search per k: maximize L_yK8 and report the argmax
      and slack 10*cap - max L_yK8. If any config violates (slack<0), gK8 is REFUTED.
  (3) the doublet family E_f = base U {f,f+1} swept over f -- does L_yK8 stay bounded as f->inf
      (the frozen plateau of L_yK8)?
Exact rationals. A clean PASS (positive slack everywhere) supports gK8 as the wide closure;
a violation kills it and refocuses on the doublet P/R route (THM-564).
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

from lrc14_threadA_regime_dichotomy_kpswf8 import CAP
from lrc14_wide_branch_ridge_codex_s47 import primitive
ALL_INNER = 0b1111110


def miss_dist(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return [F(0)] * 7
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l; den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e; x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps); q = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi; mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        missed = 6 - bin(mask & ALL_INNER).count("1")
        q[missed] += F(hi - lo, d)
    return q


def LyK8(E):
    q = miss_dist(E)
    return 10 * q[0] + q[3] + 10 * q[6], q[0]


def reprim(S):
    S = tuple(sorted(set(int(x) for x in S)))
    g = reduce(gcd, S)
    return tuple(x // g for x in S) if g > 1 else S


def span(S):
    return max(S) - min(S)


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if 0 not in E or span(E) <= 14 or not primitive(E):
        return False
    n_far = sum(1 for e in E if e > 14)
    if n_far < 2:
        return False
    for e in E:
        sub = reprim([x for x in E if x != e])
        if len(sub) < 2 or span(sub) <= 14:
            return False
    return True


def main():
    print("=" * 92)
    print("ADVERSARIAL gK8 stress-test  L_yK8 = 10q0+q3+10q6 <= 10*cap_k ?  (kps-wf9)")
    print("=" * 92)

    # (1) SHARP dilated maximizers (HYP-2805) where the 0.16 margin FAILS
    print("\n[1] SHARP genuine-wide maximizers (HYP-2805 dilated doublets, margin<0.16 cases):")
    sharp = {
        8:  (0, 10, 11, 12, 13, 14, 15, 16),
        9:  (0, 4, 6, 8, 10, 12, 14, 15, 16),
        10: (0, 2, 4, 6, 8, 10, 12, 14, 15, 16),
        12: (0, 2, 4, 5, 6, 7, 8, 10, 12, 14, 20, 21),
    }
    for k, E in sharp.items():
        if k not in CAP:
            continue
        ly, p0 = LyK8(E)
        cap = CAP[k]; bound = 10 * cap
        slack = bound - ly
        print(f"  k={k} E={E}")
        print(f"     p0={float(p0):.6f} (cap-p0={float(cap-p0):+.6f})  L_yK8={float(ly):.6f}  "
              f"10cap={float(bound):.6f}  slack={float(slack):+.6f}  {'OK' if slack>=0 else 'VIOLATION!'}")

    # (2) BROAD genuine-wide search maximizing L_yK8
    print("\n[2] BROAD genuine-wide search: maximize L_yK8 (adversarial). VIOLATION if slack<0.")
    rng = random.Random(20260621)
    for k in (8, 9, 10):
        if k not in CAP:
            continue
        cap = CAP[k]; bound = 10 * cap
        best_ly = F(-1); best_E = None; n_checked = 0
        worst_slack = None; worst_E = None
        cand = set()
        # structured: dilated bases d*consec + adjacent far doublet/triple near 15
        for d in (1, 2, 3, 4):
            base = [d * i for i in range(k - 1)]
            base = [b for b in base if b <= 14]
            for f0 in range(15, 26):
                for nf in (2, 3):
                    far = [f0 + i for i in range(nf)]
                    E = tuple(sorted(set([0] + base + far)))
                    if len(E) == k:
                        cand.add(reprim(E) if not primitive(E) else E)
        # random bounded base (any subset) + adjacent far pairs
        for _ in range(6000):
            sz = rng.randint(2, k - 2)
            base = [0] + rng.sample(range(1, 15), sz - 1)
            nf = k - sz
            if nf < 2:
                continue
            f0 = rng.randint(15, 24)
            far = sorted(set(rng.sample(range(15, 30), nf)))
            if len(far) != nf:
                continue
            E = tuple(sorted(set(base + far)))
            if len(E) == k:
                cand.add(E)
        for E in cand:
            if len(E) != k or not is_genuine_wide(E):
                continue
            n_checked += 1
            ly, p0 = LyK8(E)
            slack = bound - ly
            if ly > best_ly:
                best_ly, best_E = ly, E
            if worst_slack is None or slack < worst_slack:
                worst_slack, worst_E = slack, E
        if best_E is None:
            print(f"  k={k}: no genuine-wide candidates found in search")
            continue
        ly, p0 = LyK8(best_E)
        print(f"  k={k}: max L_yK8 = {float(best_ly):.6f} at E={best_E}")
        print(f"        10cap={float(bound):.6f}  min slack={float(worst_slack):+.6f}  "
              f"({'ALL OK' if worst_slack>=0 else 'VIOLATION!'}); checked {n_checked} genuine-wide")

    # (3) doublet family frozen plateau of L_yK8
    print("\n[3] doublet L_yK8 as f->inf (consec base): does L_yK8 stay < 10cap at all f?")
    for k in (9, 10):
        if k not in CAP:
            continue
        cap = CAP[k]; bound = 10 * cap
        base = list(range(k - 1))
        vals = []
        for f in list(range(15, 30)) + [50, 100, 200, 400]:
            E = tuple(sorted(set(base + [f, f + 1])))
            ly, p0 = LyK8(E)
            vals.append((f, float(ly), float(bound - ly)))
        worst = max(vals, key=lambda t: t[1])
        print(f"  k={k}: max L_yK8 over sampled f = {worst[1]:.6f} at f={worst[0]}  10cap={float(bound):.6f}  "
              f"slack={worst[2]:+.6f}")

    print("\n" + "=" * 92)
    print("If slack >= 0 EVERYWHERE (esp. on the sharp dilated maximizers where the 0.16 margin")
    print("fails): gK8 holds as a candidate wide-closure certificate with the noted slack. Any")
    print("negative slack REFUTES gK8 and the closure must use the doublet P/R route (THM-564).")


if __name__ == "__main__":
    main()
