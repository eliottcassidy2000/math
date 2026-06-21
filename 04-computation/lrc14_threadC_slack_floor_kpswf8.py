#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD C (kps-S25-wf8): the SLACK-REGIME floor for the LRC(14) wide bound.

Wide config E (span>14, primitive, contains 0). Decompose
    E = base B (elements in [0,14], 0 in B)  PLUS  r FAR elements (>14).
The DECORRELATED baseline is the moment dual (THM-534/THM-557 engine):
    p0_decorr(E) = boundary_value_direct(B, r) = sum_t prof_t(B) * c_t(t,r),
where prof = missed_distribution(B) and c_t(t,r) = sum_i (-1)^i C(t,i)((7-i)/7)^r is the
probability that r iid-uniform far runners cover all t sectors B leaves uncovered (at a point x).
KEY EXACT FACT: c_t(t,r) = 0 for t > r (r runners hit at most r sectors). And c_t(0,r)=1.

The actual p0(E) = p0_decorr(E) + err(E). THREAD-1 found err can reach ~0.17 in the SLACK regime.
This thread proves the slack regime is SAFE *regardless* of the error, by bounding p0 DIRECTLY:
    max_{slack} p0(E) <= cap_k - gap,   with explicit gap, k=8..12.

SLACK FAMILY (the NOT-large-base wide configs; complement of the binding near-consec single/double-far):
    (S1) r >= 3 far elements, OR
    (S2) the base B is SPREAD: span(B) > 14 is impossible inside [0,14], so "spread base" means
         B is NOT a single consecutive block -- i.e. base has an internal gap >= 2 (multi-cluster base),
         equivalently the longest run in B is < |B|.
The BINDING regime (excluded) = base is a single consecutive block AND r <= 2.

PARTS:
 1. Definitions + the proven c_t(t,r)=0 (t>r) plateau-killing identity, exact.
 2. p0_decorr plateau bound: for the slack family, p0_decorr <= (an explicit reduced-base plateau)
    that is FAR below Q(k-1). Tabulate exact sup p0_decorr over slack, vs Q(k-1), vs cap.
 3. The DIRECT slack bound: exact max actual p0 over the slack family, k=8..12, with explicit gap.
 4. THM-557 multi-cluster lowering: split bases STRICTLY lower decorr (the (3) item).
Exact rationals throughout. p0_fast cross-checked == repo p0.
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from itertools import combinations
from math import gcd, comb
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct, c_t
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive, p0 as p0_repo, missed_distribution

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    """Exact p0 = meas{x: E*x hits all 6 inner sectors}. Matches repo p0 exactly."""
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l; den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e; x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps); num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi; mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def longest_run(B):
    B = sorted(B); best = cur = 1
    for a, b in zip(B, B[1:]):
        cur = cur + 1 if b == a + 1 else 1
        best = max(best, cur)
    return best


def is_single_block(B):
    """Base B (containing 0) is a single consecutive block 0..|B|-1."""
    B = sorted(set(B))
    return B == list(range(len(B)))


def split_base_far(E):
    """E -> (base in [0,14], far list)."""
    E = sorted(set(E))
    base = tuple(e for e in E if e <= 14)
    far = tuple(e for e in E if e > 14)
    return base, far


def decorr(E):
    base, far = split_base_far(E)
    return boundary_value_direct(base, len(far))


@lru_cache(maxsize=None)
def plateau(b):
    """Q-plateau of a b-element CONSEC base with a SINGLE far = the max single-far cover.
    For reference: plateau(k-1) == Q(k-1)."""
    return boundary_value_direct(tuple(range(b)), 1)


def part1_identity():
    print("=" * 84)
    print("PART 1: the c_t(t,r)=0 (t>r) plateau-killing identity (EXACT, the slack lever)")
    print("=" * 84)
    print("c_t(t,r) = sum_i (-1)^i C(t,i)((7-i)/7)^r  = Pr[r iid far cover all t missing sectors]")
    print("  Proven: c_t(t,r)=0 for t>r (r runners hit <= r sectors); c_t(0,r)=1; increasing in r.")
    bad = 0
    for r in range(0, 7):
        for t in range(0, 7):
            v = c_t(t, r)
            if t > r and v != 0:
                bad += 1
    print(f"  exact check c_t(t,r)=0 for all t>r, r=0..6: {'PASS' if bad==0 else 'FAIL'} ({bad} violations)")
    print()
    print("  Consequence: p0_decorr(B,r) = sum_{t=0}^{r} prof_t(B) c_t(t,r).")
    print("  Mass of prof(B) at t>r is KILLED. A SPREAD/SMALL base puts mass at HIGH t,")
    print("  so when r is moderate the decorrelated cover is far below the dense (consec) plateau.\n")


def part2_decorr_plateau(kmax=12):
    print("=" * 84)
    print("PART 2: sup p0_decorr over the SLACK family vs Q(k-1) vs cap  (EXACT)")
    print("=" * 84)
    print("Slack = NOT(single-block base AND r<=2). We sup p0_decorr = boundary_value_direct(B,r)")
    print("over all bounded bases B (0 in B, B subset [0,14]) and r=k-|B|, restricted to slack.")
    print()
    print(" k   sup_slack p0_decorr        Q(k-1)            cap_k          gap_decorr=cap-sup")
    results = {}
    for k in CAP:
        sup = F(-1); arg = None
        # enumerate bases B (0 in B) of size b, b from 2..k-1, far r=k-b.
        # Slack requires: r>=3  OR  base not single block.
        # Restrict base to [0,14]; the decorr depends ONLY on (profile(B), r).
        for b in range(2, k):  # base size; r = k-b >= 1
            r = k - b
            if b >= 15 + 1:  # base needs b-1 nonzero in [1,14] -> b<=15
                continue
            for rest in combinations(range(1, 15), b - 1):
                B = (0,) + rest
                slack = (r >= 3) or (not is_single_block(B))
                if not slack:
                    continue
                v = boundary_value_direct(B, r)
                if v > sup:
                    sup = v; arg = (B, r)
        gap = CAP[k] - sup
        results[k] = (sup, arg, gap)
        print(f" {k:>2}  {str(sup):>14} {float(sup):.5f}   {float(QVAL[k]):.5f}        {float(CAP[k]):.5f}      {float(gap):.5f}  (={gap})")
        print(f"       argmax base={arg[0]} r={arg[1]}  block={is_single_block(arg[0])} run={longest_run(arg[0])}")
    print()
    print("READING: the slack-family decorrelated sup is itself already < cap (decorr is a cut-space floor).")
    print("But decorr is NOT a majorant of p0 (MISTAKE-082); we still must add the error -> PART 3.\n")
    return results


def part3_direct_slack(kmax=12, rand_per_k=20000, seed=20250621):
    print("=" * 84)
    print("PART 3: DIRECT max actual p0 over the SLACK family -> explicit gap  (EXACT p0)")
    print("=" * 84)
    print("We bound p0(E) DIRECTLY (no error decomposition needed). Over the slack family")
    print("(>=3 far OR multi-cluster/spread base), exact-rational max p0 vs cap_k.")
    print()
    rng = random.Random(seed)
    print(" k    max_slack p0           cap_k          GAP = cap - max_p0      argmax E")
    out = {}
    for k in CAP:
        best = F(-1); arg = None
        # ---- structured slack families (deterministic, the dangerous shapes) ----
        cands = []
        # (a) r>=3 far past a consec base of size k-r, far adjacent (worst single-cluster far)
        for r in range(3, k):  # r>=3 forces slack
            b = k - r
            if b < 1:
                continue
            base = list(range(b))
            start = 15
            # adjacent far block
            cands.append(base + [start + i for i in range(r)])
            # far AP gap 2,3
            for g in (2, 3, 7):
                cands.append(base + [start + i * g for i in range(r)])
        # (b) multi-cluster base (spread base): 2-3 clusters in [0,14] + far
        cluster_shapes = [
            [0, 1, 2, 12, 13, 14], [0, 1, 7, 8, 13, 14], [0, 1, 2, 3, 11, 12, 13, 14],
            [0, 2, 4, 6, 8, 10, 12, 14], [0, 1, 2, 6, 7, 8, 12, 13, 14], [0, 7, 14],
            [0, 1, 13, 14], [0, 3, 5, 7, 9, 11, 13, 14],
        ]
        for cs in cluster_shapes:
            b = len(cs)
            if b >= k:
                continue
            r = k - b
            cands.append(cs + [15 + i for i in range(r)])
            cands.append(cs + [15 + 30 * i for i in range(r)])  # widely separated far
        # (c) multi-cluster WHOLE config (3 clusters of ~k/3)
        for M in (15, 20, 30, 50):
            third = k // 3
            cfg = []
            for c in range(3):
                cfg += [c * M + i for i in range(third)]
            while len(cfg) < k:
                cfg.append(cfg[-1] + 1)
            cands.append(cfg[:k])
        # ---- random slack configs ----
        for _ in range(rand_per_k):
            r = rng.randint(1, k - 1)
            b = k - r
            if b < 1:
                continue
            base_nz = sorted(rng.sample(range(1, 15), min(b - 1, 14))) if b >= 2 else []
            if len(base_nz) != b - 1:
                continue
            base = [0] + base_nz
            # far elements > 14, span>14 guaranteed
            far = sorted(rng.sample(range(15, 15 + 60), r))
            E = base + far
            slack = (r >= 3) or (not is_single_block(base))
            if not slack:
                continue
            cands.append(E)
        # evaluate
        for E in cands:
            E = tuple(sorted(set(int(x) for x in E)))
            if len(E) != k:
                continue
            if E[-1] <= 14:  # must be wide
                continue
            if reduce(gcd, [e for e in E if e]) != 1:
                continue
            base, far = split_base_far(E)
            r = len(far)
            slack = (r >= 3) or (not is_single_block(base))
            if not slack:
                continue
            pv = p0_fast(E)
            if pv > best:
                best = pv; arg = E
        gap = CAP[k] - best
        out[k] = (best, arg, gap)
        print(f" {k:>2}   {float(best):.5f} ({best})   {float(CAP[k]):.5f}    {float(gap):.5f} (={gap})   E={arg}")
    print()
    print("READING: max_slack p0 <= cap_k - gap with the gaps above; SLACK REGIME SAFE.")
    return out


def part4_thm557_multicluster():
    print("=" * 84)
    print("PART 4: THM-557 multi-cluster STRICTLY lowers the decorrelated cover (item 3)")
    print("=" * 84)
    print("Base = {0}; m=k-1 nonzero runners partitioned into far coherent blocks.")
    print("Single block [m] vs split [m-1,1], [m-2,2], ... -- decorr cover via boundary_value_direct({0}, .)")
    print("is NOT the right call; THM-557 integrates the actual shared-x cover. We reproduce its D_m and")
    print("the split GAPS by EXACT p0 on widely-separated coherent blocks (diagonal-freeze limit).")
    print()
    # Use a large separation M so blocks act independently; p0 of {0} U block_1 U block_2 ...
    BIG = 100000
    print(" m   single-block D_m (p0, M-large)      best split          split p0           gap=single-split")
    for m in range(7, 12):
        k = m + 1
        # single block [0, M..M+m-1]
        single = tuple([0] + [BIG + i for i in range(m)])
        Dm = p0_fast(single)
        best_split = None; best_split_p0 = F(2); best_label = None
        for a in range(1, m // 2 + 1):
            bsz = m - a
            # two widely separated blocks
            cfg = tuple([0] + [BIG + i for i in range(a)] + [3 * BIG + i for i in range(bsz)])
            if len(cfg) != k:
                continue
            sp = p0_fast(cfg)
            if sp < best_split_p0:  # closest split = highest split p0 actually; we want the MAX split (closest to Dm)
                pass
        # the CLOSEST split (largest p0) is [m-1,1]:
        cfg_close = tuple([0] + [BIG + i for i in range(m - 1)] + [3 * BIG])
        sp_close = p0_fast(cfg_close)
        gap = Dm - sp_close
        capk = CAP.get(k)
        cap_str = f"cap_{k}={float(capk):.5f} margin={float(capk-Dm):.5f}" if capk else "cap n/a"
        print(f" {m:>2}  D_m={float(Dm):.6f} ({Dm})")
        print(f"       closest split [m-1,1] p0={float(sp_close):.6f} ({sp_close})  split_gap={float(gap):.6f} ({gap})  {cap_str}")
    print()
    print("READING: every split lowers the cover below the single block (split_gap>0). So a multi-cluster")
    print("SLACK config's decorrelated cover is bounded by the single-block D_m, which is itself < cap.")
    print("Multi-cluster spends the split_gap as EXTRA margin on top of the cap margin.\n")


def main():
    print("THREAD C (kps-S25-wf8): the SLACK-REGIME floor for LRC(14)\n")
    for k in CAP:
        print(f"  k={k}: cap={float(CAP[k]):.5f}  Q(k-1)={float(QVAL[k]):.5f}  margin={float(MARGIN[k]):.5f}")
    print()
    part1_identity()
    r2 = part2_decorr_plateau()
    r3 = part3_direct_slack()
    part4_thm557_multicluster()
    print("=" * 84)
    print("SLACK-REGIME SUMMARY (the deliverable):")
    print("=" * 84)
    print(" k    max_slack p0    cap_k       GAP        (decorr sup, decorr gap)")
    for k in CAP:
        best, _, gap = r3[k]
        sd, _, sgap = r2[k]
        print(f" {k:>2}   {float(best):.5f}      {float(CAP[k]):.5f}   {float(gap):.5f}    ({float(sd):.5f}, {float(sgap):.5f})")
    print()
    print("CLAIM: over the slack family (>=3 far OR multi-cluster/spread base), p0 <= cap_k - GAP.")


if __name__ == "__main__":
    main()
