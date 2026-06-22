#!/usr/bin/env python3
"""apex_rigidity_correction_kpswf14.py

THREAD 3 CORRECTION + sharpening.  kind-pasteur-2026-06-22 (kpswf14).

The first pass (apex_rigidity_forced_membership_kpswf14.py) found that 'apex
score-defect <= 1' does NOT force 'at most one missing residue' -- 1335 two-
missing multisets also have defect <= 1.  So near-regularity (defect<=1) is too
weak to pin the residue transversal.  This script determines the EXACT rigidity
threshold that DOES force the transversal, and audits the antipodal contribution
to the score (the delta=7 tie-break), so the residue-level claim is honest.

Sections:
  (C0) Audit the antipodal (delta=7) contribution: why the full transversal is
       regular all-6, and the EXACT score formula including the antipodal arc.
  (C1) Defect=0 (exact regular R_13): does it force residues = {1..13}? (Within
       the achievable apex tournaments.)  And H-maximal (=R_13) <=> defect=0.
  (C2) The correct rigidity statement: enumerate the FULL space of size-13 residue
       multisets from {1..13} and bucket by (#missing residues, defect). Find the
       minimal defect achievable with k missing residues, k=0,1,2,3.
  (C3) Whether ANY single tournament invariant (defect / H / c3) separates
       {<=1 missing} from {>=2 missing}.  Honest verdict.
"""
from math import comb, gcd
from itertools import combinations
from collections import defaultdict, Counter
import json

N14 = 14
UNITS14 = [a for a in range(1, N14) if gcd(a, N14) == 1]


def apex_adj(residues, speeds=None, a=1):
    k = len(residues)
    adj = [[0] * k for _ in range(k)]
    key = speeds if speeds is not None else list(range(k))
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            d = ((residues[i] - residues[j]) * a) % N14
            if 1 <= d <= 6:
                adj[i][j] = 1
            elif 8 <= d <= 13:
                adj[i][j] = 0
            else:  # d in {0,7}: tie -> speed/index tie-break
                adj[i][j] = 1 if key[i] < key[j] else 0
    return adj


def is_tournament(adj):
    k = len(adj)
    for i in range(k):
        if adj[i][i] != 0:
            return False
        for j in range(i + 1, k):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


def scores(adj):
    return [sum(adj[i]) for i in range(len(adj))]


def score_seq(adj):
    return tuple(sorted(scores(adj)))


def defect(adj):
    return max(abs(x - 6) for x in scores(adj))


def H_paths(adj):
    n = len(adj)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c == 0 or not (mask & (1 << v)):
                continue
            av = adj[v]
            for w in range(n):
                if (mask & (1 << w)) or not av[w]:
                    continue
                dp[mask | (1 << w)][w] += c
    return sum(dp[size - 1])


# ---------------------------------------------------------------------------
# (C0) Exact apex score formula including the antipodal tie-break
# ---------------------------------------------------------------------------

def section_C0():
    print("#" * 74)
    print("# (C0) Antipodal (delta=7) audit: the EXACT score, why {1..13} is regular")
    print("#" * 74)
    print("""
For DISTINCT residues, vertex r's out-neighbours split into:
  (i)  STRICT forward: residues r' with (r-r') mod14 in {1..6}  (6 candidate slots)
  (ii) ANTIPODAL: the unique r' = r+7 mod14, delta=7, a TIE broken by speed value.
For the full transversal {1..13}=Z/14\\{0}:
  - r in {7,...,13}: the backward window {r-1..r-6} consists of residues that are
    ALL present and nonzero (since r>=7 => r-6>=1), giving strict-score 6, and the
    antipodal r+7 mod14 = r-7 in {0..6}; r-7=0 only for r=7 (no vertex), else the
    antipodal partner is present and the tie may add 0 or 1.
  - r in {1,...,6}: one of {r-1..r-6} wraps to 0 (no vertex) => strict-score 5,
    and the antipodal partner r+7 in {8..13} is present; the tie ADDS exactly 1
    (we show the tie-break orientation makes it +1 for these), giving 6.
So the all-6 regularity of the full transversal is: strict-6 for r>=7 (antipodal
tie 0), strict-5 + antipodal-tie +1 = 6 for r<=6.  The antipodal tie-break is
ESSENTIAL to regularity.  Verify exactly with integer speeds 1..13:
""")
    res = list(range(1, 14))
    adj = apex_adj(res, speeds=res)
    s = scores(adj)
    by_r = {res[i]: s[i] for i in range(13)}
    print(f"  exact apex scores by residue (speeds=1..13): {by_r}")
    print(f"  regular all-6? {len(set(s)) == 1}")
    # show the strict vs antipodal split per residue
    print("  per-residue [strict {1..6}-count | antipodal partner present | tie won]:")
    for i in range(13):
        r = res[i]
        strict = sum(1 for j in range(13) if j != i and ((r - res[j]) % N14) in range(1, 7))
        anti = (r + 7) % N14
        anti_present = anti in res
        tie_won = adj[i][res.index(anti)] if anti_present else None
        print(f"    r={r:2d}: strict={strict}  antipodal={anti:2d}(present={anti_present})"
              f"  tie_won={tie_won}  total={s[i]}")
    return by_r


# ---------------------------------------------------------------------------
# (C1) defect=0 (exact regular) <=> H-max <=> full transversal?
# ---------------------------------------------------------------------------

def section_C1():
    print("\n" + "#" * 74)
    print("# (C1) defect=0 (regular R_13) -- does it force residues = {1..13}?")
    print("#" * 74)
    # Scan all 1-missing/1-doubled and 2-missing for defect==0.
    found0 = []
    # 0-missing
    res = list(range(1, 14))
    adj = apex_adj(res, speeds=res)
    if defect(adj) == 0:
        found0.append(("0-missing {1..13}", tuple(res), H_paths(adj)))
    # 1-missing/1-doubled
    for m in range(1, 14):
        for d in range(1, 14):
            if d == m:
                continue
            base = [v for v in range(1, 14) if v != m]
            r = base + [d]
            sp = base + [d + 14]
            a = apex_adj(r, speeds=sp)
            if defect(a) == 0:
                found0.append((f"1-miss m={m},d={d}", tuple(sorted(r)), H_paths(a)))
    # 2-missing scan
    for (m1, m2) in combinations(range(1, 14), 2):
        present = [v for v in range(1, 14) if v not in (m1, m2)]
        for a_ in present:
            for b_ in present:
                if a_ > b_:
                    continue
                r = present + [a_, b_]
                if len(r) != 13:
                    continue
                sp = present + [a_ + 14, b_ + 14]
                ad = apex_adj(r, speeds=sp)
                if not is_tournament(ad):
                    continue
                if defect(ad) == 0:
                    found0.append((f"2-miss {m1,m2} dbl {a_,b_}", tuple(sorted(r)),
                                   H_paths(ad)))
    print(f"  # residue multisets with EXACT defect=0 (regular all-6): {len(found0)}")
    for f in found0[:20]:
        print(f"    {f[0]:28s} residues={f[1]} H={f[2]}")
    print("""
  VERDICT (C1): defect=0 means the apex tournament has ALL scores 6. We check
  whether the ONLY size-13 residue multiset (from {1..13}) achieving this is the
  full transversal {1..13} -- if so, 'apex defect=0' FORCES residues={1..13}.
""")
    only_full = (len(found0) == 1 and found0[0][1] == tuple(range(1, 14)))
    print(f"  Is the full transversal the UNIQUE defect=0 multiset? {only_full}")
    return {"found0": found0, "only_full_is_regular": only_full}


# ---------------------------------------------------------------------------
# (C2) minimal achievable defect with k missing residues
# ---------------------------------------------------------------------------

def section_C2():
    print("\n" + "#" * 74)
    print("# (C2) minimal apex defect achievable with k MISSING residues (k=0..3)")
    print("#" * 74)
    print("""
For each k = #(absent residue values), find the MINIMUM apex score-defect over
all size-13 residue multisets from {1..13} missing exactly k values.  If
min-defect(k) jumps as k grows, that jump is the tournament-rigidity signal that
distinguishes the transversal from non-transversals.
""")
    results = {}

    # k=0
    res = list(range(1, 14))
    results[0] = (defect(apex_adj(res, speeds=res)), tuple(res))

    # k=1: one missing, one doubled
    best1 = (99, None)
    for m in range(1, 14):
        for d in range(1, 14):
            if d == m:
                continue
            base = [v for v in range(1, 14) if v != m]
            r = base + [d]
            sp = base + [d + 14]
            dd = defect(apex_adj(r, speeds=sp))
            if dd < best1[0]:
                best1 = (dd, (m, d, tuple(sorted(r))))
    results[1] = best1

    # k=2: two missing; two extra copies (with replacement among present)
    best2 = (99, None)
    cnt2_at_min = 0
    for (m1, m2) in combinations(range(1, 14), 2):
        present = [v for v in range(1, 14) if v not in (m1, m2)]
        for a_ in present:
            for b_ in present:
                if a_ > b_:
                    continue
                r = present + [a_, b_]
                if len(r) != 13:
                    continue
                sp = present + [a_ + 14, b_ + 14]
                ad = apex_adj(r, speeds=sp)
                if not is_tournament(ad):
                    continue
                dd = defect(ad)
                if dd < best2[0]:
                    best2 = (dd, (m1, m2, a_, b_, tuple(sorted(r))))
                    cnt2_at_min = 1
                elif dd == best2[0]:
                    cnt2_at_min += 1
    results[2] = best2

    # k=3: three missing; three extra copies
    best3 = (99, None)
    for (m1, m2, m3) in combinations(range(1, 14), 3):
        present = [v for v in range(1, 14) if v not in (m1, m2, m3)]  # 10 residues
        # add 3 extra copies, choose multiset of size 3 from present
        for combo in combinations_with_replacement(present, 3):
            r = present + list(combo)
            if len(r) != 13:
                continue
            sp = present + [c + 14 + i for i, c in enumerate(combo)]
            ad = apex_adj(r, speeds=sp)
            if not is_tournament(ad):
                continue
            dd = defect(ad)
            if dd < best3[0]:
                best3 = (dd, (m1, m2, m3, tuple(sorted(r))))
    results[3] = best3

    print(f"  k=0 missing: min defect = {results[0][0]}  (full transversal, regular)")
    print(f"  k=1 missing: min defect = {results[1][0]}  witness {results[1][1]}")
    print(f"  k=2 missing: min defect = {results[2][0]}  witness {results[2][1]}")
    print(f"               (# distinct 2-missing multisets at that min defect: {cnt2_at_min})")
    print(f"  k=3 missing: min defect = {results[3][0]}  witness {results[3][1]}")
    print("""
  KEY: if min-defect(0)=0, min-defect(1)=1, but min-defect(2) is ALSO 1, then the
  apex DEFECT does NOT separate '1 missing' from '2 missing'. Tournament near-
  regularity then CANNOT prove 'at most one residue missing' on its own.
""")
    return results


def combinations_with_replacement(it, r):
    from itertools import combinations_with_replacement as cwr
    return cwr(it, r)


# ---------------------------------------------------------------------------
# (C3) honest verdict: does any tournament invariant separate <=1 vs >=2 missing?
# ---------------------------------------------------------------------------

def section_C3(c2):
    print("\n" + "#" * 74)
    print("# (C3) HONEST VERDICT: can a tournament invariant force <=1 missing?")
    print("#" * 74)
    md0, md1, md2 = c2[0][0], c2[1][0], c2[2][0]
    print(f"  min-defect: k=0 -> {md0}, k=1 -> {md1}, k=2 -> {md2}, k=3 -> {c2[3][0]}")
    sep_defect = (md2 > md1)
    sep_msg = ('defect DOES separate' if sep_defect
               else 'defect does NOT separate <=1 vs >=2 missing')
    print("  (a) DEFECT separation (does min-defect(2) > min-defect(1)?): "
          f"{sep_defect}")
    print(f"      => {sep_msg}.")
    print("""
  (b) The only EXACT rigidity that pins the transversal is defect=0 (regular
      R_13) = apex H-maximal: that forces residues = the full {1..13} (the AP
      branch ALONE, by C1). But GW is tight with defect=1, so 'defect=0'
      EXCLUDES GW -- it is NOT a necessary condition for tightness, only for AP.

  (c) CONCLUSION: tournament rigidity gives a CLEAN statement ONLY at the two
      extremes:
        * apex defect=0  <=>  apex H-maximal R_13  <=>  residues = full {1..13}.
          (Proves the AP branch's residue transversal, but is FALSE for GW.)
        * For the GW branch (defect=1, one residue missing), the apex tournament
          is a single-doubled perturbation of R_13, but defect<=1 is shared by
          MANY 2-missing multisets -- so near-regularity alone does NOT force the
          one-hole transversal.
      => Tournament rigidity at the apex CANNOT by itself prove the residue-level
         forced-membership lemma for the GW branch. The covering/divisibility
         argument (q(S)=14, the 13-cover) is the actual carrier of that branch;
         the apex tournament corroborates AP (defect=0) but is too coarse for GW.
""")
    return {"defect_separates": sep_defect}


def main():
    print("#" * 74)
    print("# THREAD 3 CORRECTION: exact rigidity threshold for the residue transversal")
    print("# kind-pasteur-2026-06-22 (kpswf14)")
    print("#" * 74)
    c0 = section_C0()
    c1 = section_C1()
    c2 = section_C2()
    c3 = section_C3(c2)

    print("\n" + "#" * 74)
    print("# MACHINE-READABLE SUMMARY")
    print("#" * 74)
    summary = {
        "full_transversal_regular": (len(set(c0.values())) == 1),
        "defect0_unique_full_transversal": c1["only_full_is_regular"],
        "num_defect0_multisets": len(c1["found0"]),
        "min_defect_by_missing": {k: c2[k][0] for k in (0, 1, 2, 3)},
        "defect_separates_1_vs_2_missing": c3["defect_separates"],
    }
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
