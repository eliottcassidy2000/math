#!/usr/bin/env python3
"""
PART 2 of adversarial verification: stress-test the angle's structural claims.

kind-pasteur-2026-06-18-S5-wf

(A) Scale-invariance L1: mu(g*E) == mu(E) for integer g>0, BOTH thresholds. The
    angle leans HEAVILY on this (single run <=> consecutive). Verify exactly.
(B) The angle's own SELF-GAP: "the 2/7 window-minimum keeps decreasing with
    spread, no clean plateau through spread 20". Reproduce the k=9 exhaustive
    spread18/19/20 numbers (0.198/0.252/0.205) -> confirms NO proven spread bound.
(C) HUNT for mu(E) = 0 (any k<=13, any spread). L2 says mu>0 per shape; if we
    EVER find mu=0 the whole edifice collapses. Search consecutive, perforated,
    common-factor, and structured (AP-of-runs) families.
(D) Adversarial: can we drive mu(E) BELOW the canon search-min c0(k) at any k?
    i.e. is the bounded-shape minimum the angle quotes actually a minimum, or can
    a cleverer shape beat it? This tests whether the floor (even if unproven) is
    at least where they say it is.
"""

from fractions import Fraction as F
from itertools import combinations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

# re-import the rigorous engine
import importlib.util
spec = importlib.util.spec_from_file_location(
    "verifmod",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "lrc14_Bk_verify_subtorusrelation_kps-S5-wf.py"))
vm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vm)
mu_rig = vm.mu_rig


def primitivize(E):
    from math import gcd
    g = 0
    for e in E:
        g = gcd(g, e)
    if g == 0:
        return tuple(E)
    return tuple(e // g for e in E)


if __name__ == "__main__":
    print("="*90)
    print("(A) SCALE-INVARIANCE L1: mu(g*E) == mu(E), both thresholds")
    print("="*90)
    bad = 0
    for thr in (F(2,7), F(1,7)):
        for k in range(4, 9):
            base = list(range(k))
            for g in (2,3,5,7):
                Eg = [g*e for e in base]
                v0 = mu_rig(base, thr)
                vg = mu_rig(Eg, thr)
                ok = (v0 == vg)
                bad += (not ok)
                tag = "OK" if ok else "*** FAIL ***"
                print(f"  thr={thr} k={k} g={g}: mu(base)={v0} mu(gE)={vg} {tag}")
        print()
    # Also test a NON-run dilation: mu(g*E)==mu(E) for arbitrary E
    print("  --- arbitrary E dilation ---")
    for E in [(0,2,3,4,5,6,8), (0,1,4,9,16), (0,3,5,11)]:
        for g in (2,3,5):
            v0 = mu_rig(list(E)); vg = mu_rig([g*e for e in E])
            ok = v0==vg; bad += (not ok)
            print(f"  E={E} g={g}: {v0} vs {vg} {'OK' if ok else '*** FAIL ***'}")
    print(f"\n  Scale-invariance failures: {bad}")

    print("\n" + "="*90)
    print("(B) k=9 EXHAUSTIVE bounded-spread minimum at spreads 18,19,20 (2/7)")
    print("    angle claims downward envelope 0.198 / 0.252 / 0.205 -> no plateau")
    print("="*90)
    for W in (17, 18, 19, 20):
        best = None; bestE = None
        # E = {0} ∪ 7 chosen from 1..W-1 ∪ {W}, gcd must be 1; exhaustive over
        # 7-subsets of interior with both endpoints fixed (0 and W).
        interior = list(range(1, W))
        cnt = 0
        for mid in combinations(interior, 7):
            E = (0,) + mid + (W,)
            # require gcd 1 (primitive) else it's a dilation handled by scale-inv
            from math import gcd
            g = 0
            for e in E: g = gcd(g, e)
            if g != 1:
                continue
            cnt += 1
            v = mu_rig(list(E))
            if best is None or v < best:
                best = v; bestE = E
        print(f"  spread {W}: min mu = {best} ~= {float(best):.5f}  at {bestE}  "
              f"({cnt} primitive shapes)")

    print("\n" + "="*90)
    print("(C) HUNT for mu(E)=0  (L2 says impossible: mu>=5/(7 maxE)>0)")
    print("="*90)
    found_zero = False
    # consecutive
    for k in range(3, 14):
        v = mu_rig(list(range(k)))
        if v == 0:
            found_zero = True
            print(f"  *** mu=0 consecutive k={k} ***")
    # perforated near-AP families k up to 10, modest spread
    import random
    random.seed(7)
    trials = 0
    minfound = None; minE = None
    for k in range(5, 11):
        for W in range(k, k+14):
            interior = list(range(1, W))
            if len(interior) < k-2:
                continue
            # random sample of shapes (full exhaustion too big for big W)
            n_samp = 4000
            for _ in range(n_samp):
                mid = tuple(sorted(random.sample(interior, k-2)))
                E = (0,) + mid + (W,)
                v = mu_rig(list(E))
                trials += 1
                if v == 0:
                    found_zero = True
                    print(f"  *** mu=0 found: E={E} ***")
                if minfound is None or v < minfound:
                    minfound = v; minE = E
    print(f"  sampled {trials} structured shapes; any mu=0? {found_zero}")
    print(f"  global min mu found in sample: {minfound} ~= {float(minfound):.5f} at {minE}")

    # structured AP-of-runs with large shared scale (the 'dangerous' subtori)
    print("\n  --- AP-of-runs (shared large scale q) ---")
    for q in (97, 101):
        for runlen in (2,3):
            for nruns in (3,4):
                E = []
                for r in range(nruns):
                    for s in range(runlen):
                        E.append(r*q + s)
                E = sorted(set(E))
                if len(E) > 13:
                    continue
                v = mu_rig(E)
                print(f"  q={q} runlen={runlen} nruns={nruns} k={len(E)}: mu={v} ~={float(v):.5f}")
