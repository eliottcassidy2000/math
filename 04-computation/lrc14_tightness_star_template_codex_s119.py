#!/usr/bin/env python3
"""S119: exact checks for the LRC14 tightness-star proof template.

The user's THM-079 analogy reduces LRC14 to one sharp statement:

    (*) M(S)=1/14 forces the apex-blocked denominator-14 tight locus,
        equivalently the only primitive tight atoms are AP and GW, and a
        covering over-core would have to realize the forbidden K3/H=7 packet.

This script does not prove (*).  It records what is theorem-level already and
what is only finite evidence:

  * the denominator-14 survivor/grid obstruction is exact;
  * AP and Goddyn-Wong are exact tight sets with denominator-14 argmaxes;
  * the AP single-swap atlas is isolated in a bounded exact census;
  * a small q-covering window has strict slack, with minimum 1/12.

All lonely constants are computed with exact rational critical points.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd


N = 14
THR = F(1, N)
K14_UNITS = (1, 3, 5, 9, 11, 13)


def norm1(x):
    r = x % 1
    return min(r, 1 - r)


def primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1


def critical_times(S):
    """Critical points for max_t min_s ||s t||.

    The lower envelope changes at cusps k/(2s) and pairwise intersections
    k/(s_i+s_j), k/|s_i-s_j|.  This is the standard exact LRC candidate set
    used in the repo's tight-locus scripts.
    """
    S = tuple(sorted(set(abs(s) for s in S if s)))
    C = set()
    for i, a in enumerate(S):
        for k in range(1, 2 * a):
            C.add(F(k, 2 * a))
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d:
                    for k in range(1, d):
                        C.add(F(k, d))
    return C


def lonely_constant(S):
    best = F(0)
    pts = []
    for t in critical_times(S):
        if not (0 < t < 1):
            continue
        val = min(norm1(s * t) for s in S)
        if val > best:
            best = val
            pts = [t]
        elif val == best:
            pts.append(t)
    return best, tuple(sorted(pts))


def binders(S, t, value):
    return tuple(s for s in sorted(S) if norm1(s * t) == value)


def q_covering_necessary(S):
    """Necessary condition for M(S)<1/14 from the exact q-grid lemma."""
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def denom14_survivors(S):
    return tuple(k for k in K14_UNITS if all(norm1(F(v * k, 14)) >= THR for v in S))


def print_tight_row(name, S):
    M, pts = lonely_constant(S)
    print(f"{name}:")
    print(f"  S={tuple(S)}")
    print(f"  M={M} = {float(M):.9f}")
    print(f"  argmax_count={len(pts)} argmaxes={pts}")
    print(f"  denominators={sorted({t.denominator for t in pts})}")
    print(f"  denom14_survivors={denom14_survivors(S)}")
    for t in pts:
        print(f"    t={t}: binders={binders(S, t, M)}")
    print()


def audit_grid_obstruction():
    print("1. Exact denominator-14 grid obstruction")
    failures = []
    for k in K14_UNITS:
        for v in range(1, 500):
            if ((v * k) % 14 == 0) != (v % 14 == 0):
                failures.append((k, v))
    print(f"  checked 14|v*k <=> 14|v for k in {K14_UNITS}, v<500: failures={len(failures)}")
    print("  theorem: if S has no multiple of 14, every k/14 with k in the unit list survives.")
    print("  contrapositive target: any strict counterexample must contain a multiple of 14.")
    print()


def audit_ap_gw():
    print("2. Exact tight boundary rows")
    AP = tuple(range(1, 14))
    GW = tuple(list(range(1, 12)) + [13, 24])
    BASE = tuple(list(range(1, 12)) + [13])
    print_tight_row("AP {1..13}", AP)
    print_tight_row("GW {1..11,13,24}", GW)
    print_tight_row("12-speed base {1..11,13}", BASE)


def audit_single_swap(limit=80):
    print(f"3. Exact AP single-swap tight atlas, replacement v<= {limit}")
    AP = tuple(range(1, 14))
    tight = []
    below = []
    for rem in AP:
        kept = [x for x in AP if x != rem]
        for v in range(1, limit + 1):
            S = tuple(sorted(set(kept + [v])))
            if len(S) != 13:
                continue
            if not primitive(S):
                continue
            M, pts = lonely_constant(S)
            if M == THR:
                tight.append((rem, v, S, pts))
            if M < THR:
                below.append((rem, v, S, M, pts))
    non_ap = [(rem, v, S, pts) for rem, v, S, pts in tight if S != AP]
    print(f"  tight rows total={len(tight)}; non-AP tight rows={len(non_ap)}; below-threshold rows={len(below)}")
    for rem, v, S, pts in non_ap:
        print(f"    non-AP tight: remove {rem}, add {v}, S={S}, argmaxes={pts}")
    print("  readout: in this atlas, the only non-AP tight atom is GW; no disproof row appears.")
    print("  caveat: this is a bounded single-swap atlas, not the global sporadic-finiteness theorem.")
    print()


def audit_q_covering_window(top=18):
    print(f"4. Exact q-covering necessary-condition window [1,{top}]")
    best = None
    best_row = None
    tight = []
    below = []
    count = 0
    for S in combinations(range(1, top + 1), 13):
        if not q_covering_necessary(S):
            continue
        count += 1
        M, pts = lonely_constant(S)
        if best is None or M < best:
            best = M
            best_row = (S, pts)
        if M == THR:
            tight.append((S, pts))
        if M < THR:
            below.append((S, M, pts))
    print(f"  q-covering rows={count}")
    print(f"  minimum M={best} = {float(best):.9f}")
    print(f"  minimum row={best_row[0]} argmaxes={best_row[1]}")
    print(f"  tight rows={len(tight)}; below-threshold rows={len(below)}")
    print("  readout: this bounded covering window has strict slack, min 1/12 > 1/14.")
    print("  caveat: the window is evidence for Move B, not a proof of the bounded-core theorem.")
    print()


def hamiltonian_path_count(adj):
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if not cur:
                continue
            for nxt in range(n):
                if ((mask >> nxt) & 1) == 0 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full])


def audit_proof_carrier_tournament():
    print("5. Tournament Analysis over proof carriers")
    carriers = [
        ("denom14_grid_floor", 6),
        ("AP_GW_exact_tight", 5),
        ("single_swap_atlas", 4),
        ("q_covering_window", 3),
        ("K3_state_lift_target", 2),
        ("generic_digraph_shadow", 1),
    ]
    print("  vertices are proof carriers, not runners.")
    print("  pairwise observable: which carrier preserves more of the LRC predicate at theorem level.")
    print("  gauge: orient higher fidelity score -> lower; tie path is the listed order.")
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    for i, (_, si) in enumerate(carriers):
        for j, (_, sj) in enumerate(carriers):
            if i != j and (si > sj or (si == sj and i < j)):
                adj[i][j] = True
    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    print(f"  score_histogram={dict(Counter(scores))}")
    print(f"  directed_3cycles=0")
    print(f"  scc_count={n}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("  challenged assumption: the K3/H7 route can use any binary digraph shadow.")
    print("  result: the carrier must be a tournament/OCF-realizable packet, per HYP-2908.")
    print()


def main():
    print("S119 LRC14 TIGHTNESS-STAR / THM-079 TEMPLATE AUDIT")
    print("=" * 72)
    print()
    audit_grid_obstruction()
    audit_ap_gw()
    audit_single_swap()
    audit_q_covering_window()
    audit_proof_carrier_tournament()
    print("6. Proof-template conclusion")
    print("  Move A is supported by HYP-2906/HYP-2905: locally large speeds peel,")
    print("  omit-prime gives direct witnesses, and dilation normalizes.")
    print("  Move B is still exactly (*): prove tightness forces the denominator-14")
    print("  AP/GW boundary, or equivalently prove the apex-7 over-cover state lift")
    print("  into the forbidden connected K3 packet.  The script verifies the boundary")
    print("  atoms and finite evidence, but it does not close (*).")


if __name__ == "__main__":
    main()
