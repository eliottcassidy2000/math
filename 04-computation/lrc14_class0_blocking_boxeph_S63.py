#!/usr/bin/env python3
"""
THE ADJACENT-CLASS SIGN LAW AND THE CLASS-0 LATTICE (boxeph-2026-07-17-S63)
Owner directive: the N_0 structural check (LEM-033's named open).

Resolution has TWO LAYERS plus a universal sign law:

(A) THE ADJACENT-CLASS SIGN LAW (universal, any cluster, any s): an endpoint
    whose attributed owner has boundary class c == s (mod 7) is ALWAYS an exit
    (sign -1): just after it, that owner sits in section s, so x+ is out of
    R_s.  Dually c == s+1 is ALWAYS an entry (sign +1): just before, the owner
    sits in section s.  Hence per owner N_s <= 0 <= N_{s+1} (signed).

(B) THE CLASS-0 COINCIDENCE LATTICE + ATTRIBUTION LEMMA (Layer 1): a class-0
    crossing of a 7-full owner e = 7m is x = k/e; runner f is co-boundary
    there iff m | fk.  In particular every f = mt in E is co-boundary at
    EVERY k (index tk).  The endpoint machinery attributes min(owners), so:
    if E contains mt with t <= 6 (hence mt < e), owner e can NEVER be
    attributed a class-0 endpoint: N_0^(e) = 0 IDENTICALLY, forced by the
    lattice + convention.  Balanced: 21 (12,15 < 21), 28 (12,20), 35 (15,20,
    30); near-AP: 14 (8,10,12).  Two-large: 56 (m=8), 84 (m=12) have NO small
    co-multiples -> attribution permits endpoints (at k coprime to the small-
    owner overlap).

(C) THE BOUNDARY CENSUS (Layer 2, convention-independent): is x- in R_0 at
    all?  Blocking lemmas: (B1) if some co-lattice runner f = m t' (m | fk)
    has boundary index fk/m == 1 (mod 7), then f sits in section 0 just
    before x -> both sides dead.  (B2) if 7 | k and #{f in E : m does not
    divide f} < 5, occupancy of sections 1..5 is impossible just before
    (all m-multiples sit at the top of section 6) -> dead.  Everything else:
    exact left-neighborhood evaluation, mechanism classified (some runner in
    section 0 / occupancy gap / SURVIVOR = true boundary of R_0).

(D) POSITIVE CONTROLS: two-large's survivors must coincide with actual
    endpoints() entries (exits, sole-owner k); family-70 {1..6,70} as a fresh
    attribution-free geometry.

Little statements: sign law across 6 clusters x 7 sections; sigma_e
decomposition into forced-sign classes; adjacent-class share of all endpoints.

Pure Python + Fractions (exact).  Reuses S25 endpoints, S26 owner_data.
"""

import sys
from bisect import bisect_left
from fractions import Fraction as Fr
from math import gcd, lcm

sys.path.insert(0, '04-computation')
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints
from lrc14_general_resonance_law_boxeph_S26 import owner_data

SIX = frozenset(range(1, 7))

CLUSTERS = [([12, 15, 20, 21, 28, 30, 35], "balanced"),
            ([1, 2, 3, 4, 5, 36, 60], "two-owner"),
            ([1, 2, 3, 4, 5, 6, 60], "family60"),
            ([8, 9, 10, 12, 14, 15, 18], "near-AP"),
            ([1, 2, 3, 4, 5, 56, 84], "two-large"),
            ([1, 2, 3, 4, 5, 6, 70], "family70")]


def section(f, x):
    return int((f * x % 1) * 7)


def occupancy(E, x):
    return set(section(f, x) for f in E)


def in_R0(E, x):
    occ = occupancy(E, x)
    return occ == SIX


# ---------------------------------------------------------------- PART A

def part_A():
    print("=" * 78)
    print("PART A -- the adjacent-class sign law (all clusters, all s)")
    tot_s = tot_s1 = 0
    tot_pts = 0
    viol = 0
    for E, name in CLUSTERS:
        shares = []
        for s in range(7):
            pts = endpoints(E, s)
            if not pts:
                continue
            adj = 0
            for p, sg, o in pts:
                j = p * 7 * o
                assert j.denominator == 1
                c = int(j) % 7
                if c == s:
                    tot_s += 1
                    adj += 1
                    if sg != -1:
                        viol += 1
                elif c == (s + 1) % 7:
                    tot_s1 += 1
                    adj += 1
                    if sg != +1:
                        viol += 1
            tot_pts += len(pts)
            shares.append(adj / len(pts))
        if shares:
            print(f"  [{name}] adjacent-class share of endpoints by s: " +
                  " ".join(f"{sh:.2f}" for sh in shares))
    print(f"  TOTALS: {tot_s} class-s endpoints (all exits), {tot_s1} "
          f"class-(s+1) endpoints (all entries), of {tot_pts} endpoints; "
          f"violations = {viol}")
    assert viol == 0
    # N_s <= 0 <= N_{s+1} per owner (signed), s = 0 clusters
    for E, name in CLUSTERS:
        P, M, data = owner_data(E, 0)
        if M == 0:
            continue
        for e, d in sorted(data.items()):
            assert d["N"][0] <= 0 <= d["N"][1], (name, e, d["N"])
    print("  N_s <= 0 <= N_{s+1} per owner: verified at s = 0, all clusters")


# ---------------------------------------------------------------- PART B

def part_B():
    print("=" * 78)
    print("PART B -- the coincidence lattice + attribution lemma (Layer 1)")
    for E, name in CLUSTERS:
        P, M, data = owner_data(E, 0)
        if M == 0:
            print(f"  [{name}] R_0 empty; skipped")
            continue
        for e in sorted(data):
            if e % 7 != 0:
                continue
            m = e // 7
            small = [m * t for t in range(1, 7) if m * t in E]
            # referee the lattice: co-owners at x = k/e are {f : m | f k}
            for k in (1, 2, 3, e - 1):
                x = Fr(k, e)
                own = [f for f in E if (x * 7 * f).denominator == 1]
                pred = [f for f in E if (f * k) % m == 0]
                assert own == pred, (name, e, k)
            n0 = data[e]["N"][0]
            tag = ("attribution-forced N_0 = 0 (small co-multiples "
                   f"{small})" if small else "attribution-free")
            print(f"  [{name}] owner {e} = 7*{m}: N_0 = {n0:+d}; {tag}")
            if small:
                assert n0 == 0
    print("  lattice rule own = {f : m | fk} exact; forced owners all N_0 = 0")


# ---------------------------------------------------------------- PART C

def part_C():
    print("=" * 78)
    print("PART C -- boundary census at class-0 crossings of 7-full owners")
    print("  (Layer 2: is x- in R_0 at all?  mechanism per crossing)")
    results = {}
    for E, name in CLUSTERS:
        pts = endpoints(E, 0)
        if not pts:
            continue
        ep_pos = {p for p, sg, o in pts}
        bps = sorted(set(Fr(k, 7 * f) for f in E for k in range(7 * f + 1)))
        sevens = [e for e in E if e % 7 == 0]
        if not sevens:
            continue
        for e in sevens:
            m = e // 7
            nonmult = [f for f in E if f % m != 0]
            cnt = {"B1": 0, "zero-other": 0, "occupancy": 0, "survivor": 0}
            surv = []
            for k in range(e):
                x = Fr(k, e)
                if k == 0:
                    # left side wraps to 1-: all runners at top of section 6
                    cnt["occupancy"] += 1
                    continue
                i = bisect_left(bps, x)
                assert bps[i] == x
                mid = (bps[i - 1] + x) / 2
                occ = occupancy(E, mid)
                if occ == SIX:
                    cnt["survivor"] += 1
                    surv.append(k)
                    continue
                zeros = [f for f in E if section(f, mid) == 0]
                if zeros:
                    # B1-form: a co-lattice runner with index fk/m == 1 mod 7?
                    if any((f * k) % m == 0 and (f * k // m) % 7 == 1
                           for f in zeros):
                        cnt["B1"] += 1
                    else:
                        cnt["zero-other"] += 1
                else:
                    cnt["occupancy"] += 1
            tag = f"non-multiples of {m}: {len(nonmult)}"
            print(f"  [{name}] e={e} (m={m}, {tag}): {e} crossings -> "
                  f"B1 {cnt['B1']}, other-zero {cnt['zero-other']}, "
                  f"occupancy {cnt['occupancy']}, SURVIVORS {cnt['survivor']}"
                  + (f" at k = {surv}" if surv else ""))
            results[(name, e)] = (cnt, surv, ep_pos, m, E)
    return results


# ---------------------------------------------------------------- PART D

def part_D(results):
    print("=" * 78)
    print("PART D -- survivors vs actual endpoints (Layer 1 x Layer 2)")
    for (name, e), (cnt, surv, ep_pos, m, E) in sorted(results.items()):
        if not surv:
            continue
        P, M, data = owner_data(E, 0)
        n_owned = 0
        for k in surv:
            x = Fr(k, e)
            assert x in ep_pos, (name, e, k, "survivor not an endpoint!")
            own = [f for f in E if (x * 7 * f).denominator == 1]
            if min(own) == e:
                n_owned += 1
        n0 = data[e]["N"][0]
        print(f"  [{name}] e={e}: {len(surv)} survivors, all are endpoints "
              f"(exits); {n_owned} attributed to {e}; owner N_0 = {n0:+d} "
              f"(= -{n_owned}: {n0 == -n_owned})")
        assert n0 == -n_owned


# ---------------------------------------------------------------- PART E

def part_E():
    print("=" * 78)
    print("PART E -- sigma_e decomposition: forced classes vs mixed (s = 0)")
    for E, name in CLUSTERS:
        P, M, data = owner_data(E, 0)
        if M == 0:
            continue
        row = []
        for e in sorted(data):
            N = data[e]["N"]
            sig = sum(N)
            row.append(f"{e}:sig={sig:+d}(N0={N[0]:+d},N1={N[1]:+d})")
        print(f"  [{name}] " + "  ".join(row))


if __name__ == "__main__":
    print("THE ADJACENT-CLASS SIGN LAW AND THE CLASS-0 LATTICE (boxeph S63)")
    part_A()
    part_B()
    res = part_C()
    part_D(res)
    part_E()
    print("=" * 78)
    print("done")
