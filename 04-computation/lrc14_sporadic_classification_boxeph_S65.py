#!/usr/bin/env python3
"""
THE SPORADIC CLASSIFICATION FOR ARBITRARY CLUSTERS (boxeph-2026-07-17-S65)
Owner directive: classify the sporadic survivors at general carry vectors
(LEM-035's named open).

THE ORBIT REDUCTION (LEM-036(A)).  At a class-0 crossing x = k/(7m) of a
7-full owner e = 7m in ANY cluster E, write k = mj + r (j in Z_7, r in
{0..m-1}).  The left-section vector of the other six runners is
    v_f(j, r) = phi_f j + o_f(r)   (mod 7),
    phi_f = f mod 7,   o_f(r) = floor(f r / m) - [m | f r],
so the survivor set of column r,  S_r = {j : 0 not in v and {1..5} <= v},
depends ONLY on the pair (phi, o(r)) in Z_7^6 x Z_7^6.  Runners with 7 | f
are CONSTANT coordinates (phi_f = 0): their section never moves with j.

THE FULL-COLUMN THEOREM (LEM-036(B)).  |S_r| = 6 iff phi is a permutation
of {1..6} and o(r) = t phi for some t (then v_f = phi_f (j + t) and
S_r = Z_7 minus {-t}: the SHIFTED MULTIPLICATION PERMUTATION).  Proof:
each moving coordinate hits 0 at the unique j_f = -o_f phi_f^{-1}; |S_r| = 6
forces all six zero-hits to coincide (o = t phi) and then coverage forces
phi injective nonzero (a non-injective phi has value-set of size <= 5 and
(j+t)-scaled sets equal {1..5} for at most one j).  |S_r| = 7 is impossible
whenever a moving coordinate exists (its zero-hit kills one j).
COROLLARY: LEM-035's clean columns are t = 0; the rescued column r = M/6 is
t = 0 after boundary cancellation; EVERYTHING ELSE (1 <= |S_r| <= 5) is a
SPORADIC (partial) column.  The classification is the grading by |S_r|.

CONSTANT-COORDINATE LAWS (LEM-036(C)): a constant coordinate with o_f = 0
kills the whole column; otherwise it contributes a fixed section, and the
five moving coordinates need only cover {1..5} minus what constants cover
-- this is why two-large (phi = (1,2,3,4,5,0)) has only partial columns yet
several of them.

Referee: 11 clusters (8 family + balanced + near-AP + two-large), every
7-full owner, every column: S_r by direct mod-7 orbit check == the exact
S64 survivor formula (census-verified there); full-column theorem both
directions; |S_r| = 7 never; the complete sporadic tables with mechanism
tags; the count law sum_r |S_r| = #survivors.
"""

import sys
from math import gcd

sys.path.insert(0, '04-computation')
from lrc14_family70_crt_boxeph_S64 import survivor_formula, leftvec_formula
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints

FIVE = frozenset(range(1, 6))

CLUSTERS = ([([1, 2, 3, 4, 5, 6, 7 * M], f"family{7 * M}") for M in
             range(8, 16)] +
            [([12, 15, 20, 21, 28, 30, 35], "balanced"),
             ([8, 9, 10, 12, 14, 15, 18], "near-AP"),
             ([1, 2, 3, 4, 5, 56, 84], "two-large")])


def column_data(E, e, m, r):
    """(phi, o) for column r; and S_r via the orbit."""
    others = [f for f in E if f != e]
    phi = [f % 7 for f in others]
    o = [((f * r // m) - (1 if (f * r) % m == 0 else 0)) % 7 for f in others]
    S = []
    for j in range(7):
        v = [(phi[i] * j + o[i]) % 7 for i in range(6)]
        if 0 not in v and FIVE <= set(v):
            S.append(j)
    return phi, o, S


def is_shifted_perm(phi, o):
    """o == t*phi for some t, and phi a permutation of {1..6}?"""
    if sorted(phi) != [1, 2, 3, 4, 5, 6]:
        return None
    inv = {x: pow(x, 5, 7) for x in range(1, 7)}
    ts = {(o[i] * inv[phi[i]]) % 7 for i in range(6)}
    return ts.pop() if len(ts) == 1 else None


if __name__ == "__main__":
    print("THE SPORADIC CLASSIFICATION -- LEM-036 referee (boxeph S65)")
    print("=" * 78)
    grand = {"full": 0, "sporadic": 0, "dead": 0}
    for E, name in CLUSTERS:
        if not endpoints(E, 0):
            print(f"  [{name}] R_0 empty; skipped")
            continue
        for e in [f for f in E if f % 7 == 0]:
            m = e // 7
            hist = [0] * 8
            spor = []
            full_cols = []
            nsurv_orbit = 0
            for r in range(m):
                phi, o, S = column_data(E, e, m, r)
                # orbit vs exact S64 formula
                S_direct = [j for j in range(7)
                            if (m * j + r) < 7 * m and (m * j + r) > 0
                            and survivor_formula(E, e, m, m * j + r)]
                assert S == S_direct or (r == 0 and sorted(S) == sorted(
                    S_direct + ([0] if 0 in S else []))), (name, e, r, S,
                                                          S_direct)
                # k = 0 (r=0, j=0) is the origin point, excluded from k-range
                if r == 0 and 0 in S:
                    S = [j for j in S if j != 0]
                nsurv_orbit += len(S)
                hist[len(S)] += 1
                t = is_shifted_perm(phi, o)
                assert len(S) < 7, (name, e, r, "seven-full column!")
                if len(S) == 6:
                    assert t is not None and (-t) % 7 not in S, (name, e, r)
                    full_cols.append((r, t))
                    grand["full"] += 1
                elif t is not None and len(S) != 6:
                    # shifted-perm form must give exactly 6 (both directions)
                    raise AssertionError((name, e, r, "shifted-perm not full"))
                elif 0 < len(S) < 6:
                    consts = [(f, o[i]) for i, f in
                              enumerate([f for f in E if f != e])
                              if f % 7 == 0]
                    spor.append((r, S, consts))
                    grand["sporadic"] += 1
                else:
                    grand["dead"] += 1
            surv_census = [k for k in range(1, e)
                           if survivor_formula(E, e, m, k)]
            assert nsurv_orbit == len(surv_census), (name, e)
            print(f"  [{name}] e={e} (m={m}): columns by |S_r| "
                  f"{{{', '.join(f'{i}:{hist[i]}' for i in range(7) if hist[i])}}}; "
                  f"full columns {full_cols if full_cols else 'NONE'}; "
                  f"survivors {nsurv_orbit} (= census {len(surv_census)})")
            for r, S, consts in spor:
                phi, o, _ = column_data(E, e, m, r)
                print(f"      SPORADIC r={r}: S_r = {S}; o = {o}"
                      + (f"; constant coords {consts}" if consts else ""))
    print(f"  GRAND: {grand['full']} full columns, {grand['sporadic']} "
          f"sporadic columns, {grand['dead']} dead columns")
    print()
    print("  full-column theorem verified BOTH directions on every column;")
    print("  |S_r| = 7 never occurs; count law exact on every owner.")
    print("=" * 78)
    print("done")
