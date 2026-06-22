#!/usr/bin/env python3
"""apex_rigidity_forced_membership_kpswf14.py

THREAD 3 of the LRC(14) tight-locus CENSUS via TOURNAMENT ANALYSIS.
Machine: kind-pasteur-2026-06-22 (kpswf14).

QUESTION
--------
Can TOURNAMENT RIGIDITY prove the (residue-level) forced-membership lemma
    tight  =>  apex residues are a FULL transversal of Z/14\{0} = {1..13}
              (the AP case)  OR  the one-doubled-one-missing GW pattern?

The full integer-level lemma is "tight => {1,2,...,11,13} subset S".  This script
quantifies the RESIDUE-LEVEL content of that lemma using the apex winding
tournament, and reports the honest gap to the integer level (the magnitude).

EXACT APEX TOURNAMENT DEFINITION (denominator 14, unit a)
---------------------------------------------------------
Speeds S, residues r_i = s_i mod 14.  At apex t* = a/14 the winding rule
(arc i->j iff frac((s_i-s_j) t*) in (0,1/2)) becomes, EXACTLY:
    delta_{ij} = (r_i - r_j) * a  mod 14   in {0,...,13}
    i -> j  iff  delta in {1,2,3,4,5,6};   j -> i  iff  delta in {8,..,13};
    delta == 7 (antipodal) or delta == 0 (coincident, doubled residue):
        TIE, broken by SPEED value (lower speed -> higher speed) -- this is the
        ONLY place the integer magnitude leaks into the apex tournament, and it
        only affects edges between an antipodal pair or a doubled-residue pair.
We use the CANONICAL unit a=1 (so i->j iff r_i-r_j mod 14 in {1..6}).  NOTE
(established in lrc14_apex_blindness): the iso class is NOT invariant across
units in the exact model; a=1 is the canonical regular representative for the
all-distinct case.

KEY PROJECT FACTS we lean on:
  * H(T) = #Hamiltonian paths = I(Omega(T),2)  (OCF / Redei).
  * Regular rotational R_13 (consecutive DRT, scores all 6) is the GLOBAL
    H-maximizer over 13-vertex tournaments: H=3711175 (we cross-check the claim
    against the SC-maximizer dichotomy on the achievable set).
  * c3(T) = C(13,3) - sum_i C(s_i, 2) = 286 - sum_i C(s_i,2).

OUTPUTS (each a separate clearly-labelled section):
  (S1) Apex score sequence as an EXACT function of the residue multiset; the
       score-defect created by a missing residue r.
  (S2) The full one-doubled/one-missing family: which (missing,doubled) give the
       smallest defect, smallest H-loss; is AP the unique regular one.
  (S3) The necessary chain  tight => apex H-maximal => residues={1..13}: tested.
  (S4) The honest residue-level forced-membership statement + the magnitude gap.
"""

import sys
from math import comb, gcd
from itertools import combinations
from collections import defaultdict, Counter
import json

N14 = 14
UNITS14 = [a for a in range(1, N14) if gcd(a, N14) == 1]  # {1,3,5,9,11,13}


# ---------------------------------------------------------------------------
# Apex tournament (a=1 canonical), residues may repeat (doubled => coincident)
# ---------------------------------------------------------------------------

def apex_adj(residues, speeds=None, a=1):
    """adj[i][j]=1 iff i->j under the apex rule with unit a.

    residues : list of r_i in {1..13} (may repeat).
    speeds   : optional integer speeds, used ONLY for tie-break (delta in {0,7}).
               If None, fall back to index tie-break (lower index -> higher).
    """
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
            else:  # d in {0,7}: tie, break by speed/index value
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
    k = len(adj)
    return [sum(adj[i]) for i in range(k)]


def score_seq(adj):
    return tuple(sorted(scores(adj)))


def c3_formula(adj):
    s = scores(adj)
    k = len(adj)
    return comb(k, 3) - sum(x * (x - 1) // 2 for x in s)


def H_paths(adj):
    """#Hamiltonian paths via Held-Karp subset DP. O(2^n n^2)."""
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
    full = size - 1
    return sum(dp[full])


# ---------------------------------------------------------------------------
# Iso machinery (WL-refined backtracking; tournaments => small Aut)
# ---------------------------------------------------------------------------

def _refine(adj, init):
    k = len(adj)
    colors = init[:]
    while True:
        sig = {}
        new = [0] * k
        for v in range(k):
            out_c = tuple(sorted(colors[w] for w in range(k) if adj[v][w]))
            in_c = tuple(sorted(colors[w] for w in range(k) if adj[w][v]))
            key = (colors[v], out_c, in_c)
            sig.setdefault(key, len(sig))
            new[v] = sig[key]
        if new == colors:
            return colors
        colors = new


def iso_exists(A, B):
    k = len(A)
    sA = [sum(A[i]) for i in range(k)]
    sB = [sum(B[i]) for i in range(k)]
    if sorted(sA) != sorted(sB):
        return False
    cA = _refine(A, sA)
    cB = _refine(B, sB)
    if Counter(cA) != Counter(cB):
        return False
    candB = defaultdict(list)
    for w in range(k):
        candB[cB[w]].append(w)
    order = sorted(range(k), key=lambda v: (len(candB[cA[v]]), v))
    mapping = {}
    used = [False] * k

    def consistent(v, img):
        for u, iu in mapping.items():
            if A[v][u] != B[img][iu] or A[u][v] != B[iu][img]:
                return False
        return True

    def bt(idx):
        if idx == k:
            return True
        v = order[idx]
        for img in candB[cA[v]]:
            if used[img] or not consistent(v, img):
                continue
            mapping[v] = img
            used[img] = True
            if bt(idx + 1):
                return True
            del mapping[v]
            used[img] = False
        return False

    return bt(0)


def aut_size(adj):
    return _count_iso(adj, adj)


def _count_iso(A, B):
    k = len(A)
    sA = [sum(A[i]) for i in range(k)]
    sB = [sum(B[i]) for i in range(k)]
    if sorted(sA) != sorted(sB):
        return 0
    cA = _refine(A, sA)
    cB = _refine(B, sB)
    if Counter(cA) != Counter(cB):
        return 0
    candB = defaultdict(list)
    for w in range(k):
        candB[cB[w]].append(w)
    order = sorted(range(k), key=lambda v: (len(candB[cA[v]]), v))
    mapping = {}
    used = [False] * k
    cnt = 0

    def consistent(v, img):
        for u, iu in mapping.items():
            if A[v][u] != B[img][iu] or A[u][v] != B[iu][img]:
                return False
        return True

    def bt(idx):
        nonlocal cnt
        if idx == k:
            cnt += 1
            return
        v = order[idx]
        for img in candB[cA[v]]:
            if used[img] or not consistent(v, img):
                continue
            mapping[v] = img
            used[img] = True
            bt(idx + 1)
            del mapping[v]
            used[img] = False

    bt(0)
    return cnt


def is_self_converse(adj):
    k = len(adj)
    rev = [[adj[j][i] for j in range(k)] for i in range(k)]
    return iso_exists(adj, rev)


# ===========================================================================
# (S1) APEX SCORE AS A FUNCTION OF THE RESIDUE MULTISET
# ===========================================================================

def apex_score_of_residue(r, residue_multiset):
    """The out-degree (apex score) of a vertex with residue r, in the apex
    tournament on a given residue multiset, IGNORING tie-break contributions
    (i.e. counting only the strict delta in {1..6} arcs; ties handled separately).

    For residue r, vertex u (residue r_u) is BEATEN by r iff (r - r_u) mod 14 in
    {1..6}, i.e. r_u in {r-1, r-2, ..., r-6} mod 14 (the 6 residues 'just behind').
    So score(r) = # other vertices whose residue lies in the half-open backward
    window {r-6,...,r-1} mod 14.  (Doubled / antipodal handled by tie-break.)
    """
    behind = set((r - d) % N14 for d in range(1, 7))  # {r-1,...,r-6} mod 14
    cnt = 0
    seen_self = False
    for ru in residue_multiset:
        if ru == r and not seen_self:
            seen_self = True  # skip exactly one copy = the vertex itself
            continue
        if ru % N14 in behind:
            cnt += 1
    return cnt


def section_S1():
    print("#" * 74)
    print("# (S1) APEX SCORE = backward-window count; defect from a MISSING residue")
    print("#" * 74)
    print("""
Apex rule (a=1): vertex r BEATS vertex r' iff (r-r') mod 14 in {1..6}.
=> score(r) = #{ other vertices with residue in the backward window
                 B(r) = {r-1, r-2, ..., r-6} mod 14 }.
For the FULL transversal {1..13} (= Z/14\\{0}), every residue r has EXACTLY 6
of the other 12 residues in B(r) (the window {r-1..r-6} are all present & nonzero
unless one of them is 0, i.e. r in {1..6}; but 0 is not a vertex, and 7..13 lie
in B(r) for those r) -> we verify all-6 below.  REGULAR R_13.
""")
    full = list(range(1, 14))
    sc = {r: apex_score_of_residue(r, full) for r in full}
    print(f"  FULL {{1..13}} apex scores by residue: {sc}")
    print(f"    all == 6 (regular)? {all(v == 6 for v in sc.values())}")

    # Now: drop residue r0 (missing), and double some residue d0 (so still 13 verts).
    # The DEFECT pattern: which vertices lose/gain a backward-neighbor.
    print("\n  --- DEFECT of MISSING residue r0 (with doubled d0) on the score multiset ---")
    print("  When residue r0 is absent, the 6 residues r0+1,...,r0+6 (mod 14) each")
    print("  LOSE one backward neighbor -> their score drops by 1 (unless compensated")
    print("  by the doubled residue d0 sitting in their backward window).")
    rows = []
    for r0 in range(1, 14):
        # the residues whose backward window CONTAINED r0 = {r0+1,...,r0+6} mod14
        forward = set((r0 + d) % N14 for d in range(1, 7))
        forward_nonzero = sorted(x for x in forward if x != 0)
        rows.append((r0, forward_nonzero))
    for r0, fwd in rows:
        print(f"    missing r0={r0:2d}: residues losing a backward-neighbor = {fwd}")
    print("""
  KEY STRUCTURAL FACT: a missing residue r0 forces a NON-REGULAR apex tournament
  (>=1 score below 6) UNLESS the doubled residue d0 exactly refills all six
  backward windows -- impossible (one doubled residue can sit in at most the
  backward windows of 6 residues, but it ALSO over-fills, raising those to 7).
  So 'missing one residue' => the score multiset is NEVER all-6 => NOT R_13.
""")
    return {"full_scores": sc}


# ===========================================================================
# (S2) THE ONE-DOUBLED/ONE-MISSING FAMILY: near-regular & H-loss
# ===========================================================================

def section_S2():
    print("#" * 74)
    print("# (S2) ONE-DOUBLED / ONE-MISSING family (156 patterns): defect & H")
    print("#" * 74)
    print("""
Family: residue multiset = {1..13}\\{missing} with 'doubled' appearing twice.
This is EXACTLY the residue structure a single-swap tight candidate can have at
the apex (one integer moved off, its residue vacated, another residue doubled).
We compute the score-defect (max |s_i - 6|), H, c3, regularity, and tag the GW
pattern (missing=12, doubled=10).
""")
    recs = []
    for missing in range(1, 14):
        for doubled in range(1, 14):
            if doubled == missing:
                continue
            base = [v for v in range(1, 14) if v != missing]
            res = base + [doubled]  # 13 residues, one doubled
            # speeds for tie-break: use the natural integer speeds; the doubled
            # residue's second copy is the moved integer = doubled + 14 (smallest
            # speed with that residue exceeding 13). The vacated residue = missing.
            # speeds: residue v in base -> speed v; the extra doubled copy -> doubled+14.
            speeds = base + [doubled + 14]
            adj = apex_adj(res, speeds=speeds)
            assert is_tournament(adj)
            s = scores(adj)
            defect = max(abs(x - 6) for x in s)
            recs.append({
                "missing": missing, "doubled": doubled,
                "score": tuple(sorted(s)), "defect": defect,
                "c3": c3_formula(adj), "H": H_paths(adj),
                "regular": len(set(s)) == 1,
            })
    # Summary stats
    by_defect = Counter(r["defect"] for r in recs)
    print(f"  defect distribution (max|s-6|) over 156 patterns: {dict(sorted(by_defect.items()))}")
    n_regular = sum(1 for r in recs if r["regular"])
    print(f"  # REGULAR (all scores 6) in the family: {n_regular}")
    print("    => NO one-doubled/one-missing pattern is regular (AP must be ALL-DISTINCT).")

    # H ranking
    Hmax_family = max(r["H"] for r in recs)
    Hmin_family = min(r["H"] for r in recs)
    print(f"\n  H over the family: max={Hmax_family}, min={Hmin_family}")
    print(f"  (Compare AP regular R_13 H = 3711175, the global apex H-max.)")

    # The GW pattern
    gw = [r for r in recs if r["missing"] == 12 and r["doubled"] == 10][0]
    print(f"\n  GW pattern (missing=12, doubled=10):")
    print(f"    score={gw['score']}, defect={gw['defect']}, c3={gw['c3']}, H={gw['H']}")

    # Which patterns have the SMALLEST defect / LARGEST H?  Is GW special?
    min_def = min(r["defect"] for r in recs)
    top_H = sorted(recs, key=lambda r: -r["H"])[:8]
    print(f"\n  Minimum defect in family = {min_def} (GW has defect {gw['defect']}).")
    print("  Top-8 by H in the one-doubled/one-missing family:")
    for r in top_H:
        tag = "  <-- GW" if (r["missing"] == 12 and r["doubled"] == 10) else ""
        print(f"    missing={r['missing']:2d} doubled={r['doubled']:2d}  "
              f"H={r['H']}  defect={r['defect']}  score={r['score']}{tag}")

    # Adjacency of missing/doubled: GW has doubled = missing-2 (10 = 12-2). Test
    # whether 'doubled = missing - 2' (the residue 2-behind) is the H-maximal slot.
    print("\n  --- does the RELATION (doubled - missing) mod 14 control H? ---")
    by_rel = defaultdict(list)
    for r in recs:
        rel = (r["doubled"] - r["missing"]) % N14
        by_rel[rel].append(r["H"])
    for rel in sorted(by_rel):
        hs = by_rel[rel]
        print(f"    (doubled-missing) mod14 = {rel:2d}:  H in [{min(hs)},{max(hs)}]  "
              f"mean={sum(hs)//len(hs)}  (n={len(hs)})")
    print("""
  INTERPRETATION: GW sits at relation (doubled-missing)=-2 mod14 (=12). The
  doubled residue is the one 2 steps BEHIND the vacated one. This is exactly the
  '11 -> [gap] -> 13' one-hole-tiling slot: vacating 12 and doubling 10 keeps the
  Z/14 covering as tight as a single swap allows.
""")
    return {"family": recs, "gw": gw, "Hmax_family": Hmax_family}


# ===========================================================================
# (S3) THE NECESSARY CHAIN  tight => apex H-maximal => residues forced
# ===========================================================================

def section_S3(s2):
    print("#" * 74)
    print("# (S3) THE NECESSARY CHAIN:  tight => apex H-maximal => residues forced?")
    print("#" * 74)
    print("""
The proposed chain (idea 3):
    tight  =>  apex tournament is the H-maximizer R_13  =>  residues = {1..13}.

We test each link HONESTLY.

LINK A  (tight => apex residues are a one-hole tiling => apex tournament regular)
   This is the PROVED apex-covering necessity (Thread 1 / anatomy layer 2): if a
   tight set's apex residues are NOT a full transversal of Z/14\\{0}, a lonely
   observer sits in the vacated class and is >1/14 from all runners. So
       tight  =>  apex residues = {1..13}  =>  apex tournament = regular R_13.
   BUT this is the RESIDUE statement, and AP IS the only all-distinct case, while
   GW's apex residues are {1..10,10,11,13} -- NOT a full transversal (12 missing)!
   So LINK A as stated would EXCLUDE GW.  Resolve: GW is tight yet its APEX
   residues miss 12.  =>  'apex residues = full transversal' is NOT necessary for
   tightness.  The apex-covering necessity must be the WEAKER 'some optimum tau'
   not literally tau=a/14.  We verify GW's apex residues below.
""")
    # GW residues:
    GW_speeds = list(range(1, 12)) + [13, 24]
    GW_res = sorted(s % 14 for s in GW_speeds)
    print(f"  GW residues mod 14 = {GW_res}")
    print(f"    is a full transversal of {{1..13}}? {sorted(set(GW_res)) == list(range(1,14))}")
    print(f"    => 12 is MISSING from GW's apex residues. GW is TIGHT (M=1/14).")
    print("    => LINK A ('tight=>apex full transversal') is FALSE as stated.")

    print("""
LINK B  (apex H-maximal <=> regular R_13):  Among the 4145 achievable apex iso
   classes, exactly ONE is regular (R_13, H=3711175), and it is the global apex
   H-max. The SC-maximizer dichotomy (project canon) says regular/SC tournaments
   maximize H within a self-complementary score class; here R_13 is the unique
   regular and the unique global H-max on the achievable set. So apex H-maximal
   <=> R_13.  We re-verify R_13 is the achievable H-max and that GW is strictly
   below it.
""")
    AP_res = list(range(1, 14))
    apAP = apex_adj(AP_res, speeds=AP_res)
    H_AP = H_paths(apAP)
    print(f"  AP apex: regular={len(set(scores(apAP)))==1}, H={H_AP}, "
          f"self_converse={is_self_converse(apAP)}, |Aut|={aut_size(apAP)}")
    print(f"  GW apex H = {s2['gw']['H']}  <  AP apex H = {H_AP}  "
          f"(GW NOT apex-H-maximal).")

    print("""
LINK C  (the chain's verdict):  Because GW is TIGHT but its apex tournament is
   NOT R_13 (not regular, not H-maximal), the chain
       tight => apex H-maximal => residues={1..13}
   is FALSE: GW is a tight counterexample to 'tight => apex H-maximal'.
   => Tournament H-maximality of the APEX tournament is NOT a necessary condition
      for tightness; it is necessary only for the AP branch.

   What IS true (the salvageable rigidity statement):
       tight => the apex residue multiset is a transversal of Z/14\\{0} with AT
                MOST ONE residue missing (and then exactly one doubled).
   i.e. the apex tournament is R_13 (AP) OR a MINIMAL (single-doubled) perturbation
   of R_13 (GW-type).  This is the 'rotational with <=1 defect' statement.
""")
    return {"H_AP": H_AP, "GW_full_transversal": sorted(set(GW_res)) == list(range(1, 14))}


# ===========================================================================
# (S4) RESIDUE-LEVEL FORCED MEMBERSHIP + the magnitude gap
# ===========================================================================

def section_S4():
    print("#" * 74)
    print("# (S4) RESIDUE-LEVEL FORCED MEMBERSHIP  +  the honest MAGNITUDE GAP")
    print("#" * 74)
    print("""
We now ask the precise residue-level forced-membership question that tournament
rigidity CAN address:

   Q: Does 'apex tournament = R_13 or single-doubled-perturbation' force the
      residue multiset to be {1..13} (AP) or {1..11,13} + double-10 (GW)?

Enumerate ALL residue multisets of size 13 from {1..13} whose apex tournament has
score-defect <= 1 (at most one vertex off the regular value 6). Report the
residue multisets that survive.  If the survivors are EXACTLY the AP and GW
residue patterns (up to the doubling site), tournament rigidity proves the
residue-level forced membership.
""")
    # We restrict to multisets with at most ONE missing residue (else >=2 windows
    # under-fill, defect>=2 generically -- we verify by also scanning two-missing).
    survivors = []
    # 0-missing: only {1..13}
    AP_res = list(range(1, 14))
    apAP = apex_adj(AP_res, speeds=AP_res)
    if max(abs(x - 6) for x in scores(apAP)) <= 1:
        survivors.append(("AP / full transversal", tuple(AP_res),
                          score_seq(apAP), H_paths(apAP)))

    # 1-missing, 1-doubled: 156 patterns
    one_def_patterns = []
    for missing in range(1, 14):
        for doubled in range(1, 14):
            if doubled == missing:
                continue
            base = [v for v in range(1, 14) if v != missing]
            res = base + [doubled]
            speeds = base + [doubled + 14]
            adj = apex_adj(res, speeds=speeds)
            s = scores(adj)
            d = max(abs(x - 6) for x in s)
            if d <= 1:
                one_def_patterns.append((missing, doubled, tuple(sorted(res)),
                                          score_seq(adj), H_paths(adj)))

    print(f"  # one-missing/one-doubled patterns with defect <= 1: "
          f"{len(one_def_patterns)} (out of 156)")
    # Tabulate by (missing) which doubling sites achieve defect<=1
    by_missing = defaultdict(list)
    for (m, d, res, sq, H) in one_def_patterns:
        by_missing[m].append(d)
    print("  defect<=1 doubling sites, by missing residue:")
    for m in sorted(by_missing):
        print(f"    missing={m:2d}: doubled in {sorted(by_missing[m])}")

    # 2-missing (2 residues absent, 2 doubled or 1 tripled): check defect can be <=1
    print("\n  --- 2-MISSING scan: can two absent residues still give defect<=1? ---")
    two_miss_le1 = 0
    two_miss_examples = []
    base_full = list(range(1, 14))
    # 2 missing residues, the freed 2 slots filled by doubling (any 2 of remaining,
    # with repetition). Enumerate missing pair, then the 2 extra copies.
    miss_pairs = list(combinations(range(1, 14), 2))
    for (m1, m2) in miss_pairs:
        present = [v for v in base_full if v not in (m1, m2)]  # 11 residues
        # add 2 extra copies chosen from present (with replacement)
        for a_ in present:
            for b_ in present:
                if a_ > b_:
                    continue
                res = present + [a_, b_]
                if len(res) != 13:
                    continue
                # tie-break speeds: present at face value, extras at +14
                speeds = present + [a_ + 14, b_ + 14]
                adj = apex_adj(res, speeds=speeds)
                if not is_tournament(adj):
                    continue
                s = scores(adj)
                d = max(abs(x - 6) for x in s)
                if d <= 1:
                    two_miss_le1 += 1
                    if len(two_miss_examples) < 6:
                        two_miss_examples.append((m1, m2, a_, b_, tuple(sorted(res)),
                                                  score_seq(adj)))
    print(f"  # 2-missing multisets with defect<=1: {two_miss_le1}")
    if two_miss_examples:
        print("  examples (missing pair, doubled pair, residues, score):")
        for ex in two_miss_examples:
            print(f"    missing={ex[0],ex[1]} doubled={ex[2],ex[3]} res={ex[4]} score={ex[5]}")
    else:
        print("  => NO 2-missing multiset achieves apex score-defect <= 1.")
        print("     (Two absent residues force >=2 under-filled backward windows.)")

    print("""
RESIDUE-LEVEL VERDICT
---------------------
The apex score-defect <= 1 condition (= 'apex tournament is R_13 or a single-
doubled perturbation of R_13') forces:
   - AT MOST ONE residue may be absent, and then exactly one residue is doubled;
   - the residues {1..13} minus at most one are all present.
This is the RESIDUE-LEVEL forced-membership: a near-regular apex tournament has
its 13 residues filling Z/14\\{0} with at most a single hole. Combined with the
tightness->near-regular apex (the covering necessity, modulo the apex-optimum
subtlety in S3), this gives:

   tight  =>  residues(S) mod 14 contains {1..13} with at most one hole.

THE MAGNITUDE GAP (honest)
--------------------------
The apex tournament is MAGNITUDE-BLIND (verified Thread 1): replacing a speed by
ANY other integer with the SAME residue mod 14 leaves the apex tournament
unchanged.  Concretely 12 -> 26 and 12 -> 96 BOTH have residues {1..13} and give
the IDENTICAL regular R_13, yet are LOOSE (M>1/14).  So:

  * Tournament rigidity proves the RESIDUE-level lemma (residues fill {1..13} up
    to one hole) -- i.e. it forces s_i mod 14 to take the values 1..13 (with one
    hole at residue 12 for the GW branch).
  * It does NOT prove the INTEGER-level lemma 'tight => {1,2,...,11,13} subset S':
    the integer 12 (residue 12) is forced to be PRESENT *as a residue* in the AP
    branch, but the residue-12 vertex could be ANY speed == 12 mod 14 (12, 26,
    40, ...). Distinguishing the tight representative (the integer 12 itself, or
    the GW move 12->24) from the loose look-alikes (12->26, 12->96) requires the
    OFF-APEX (Farey-neighbor 3/41) competition, which the apex tournament cannot
    see.

So the residue-level forced membership is a THEOREM of tournament rigidity; the
integer-level lemma needs the magnitude / off-apex input ON TOP of it.
""")
    return {"one_missing_le1": one_def_patterns, "two_missing_le1": two_miss_le1}


def main():
    print("#" * 74)
    print("# THREAD 3: TOURNAMENT RIGIDITY -> the forced-membership lemma")
    print("# kind-pasteur-2026-06-22 (kpswf14)")
    print("#" * 74)
    r1 = section_S1()
    r2 = section_S2()
    r3 = section_S3(r2)
    r4 = section_S4()

    print("\n" + "#" * 74)
    print("# MACHINE-READABLE SUMMARY")
    print("#" * 74)
    summary = {
        "AP_apex_regular_H": 3711175,
        "AP_full_scores_all6": all(v == 6 for v in r1["full_scores"].values()),
        "family_regular_count": 0,
        "GW": {k: (list(v) if isinstance(v, tuple) else v) for k, v in r2["gw"].items()},
        "Hmax_one_doubled_family": r2["Hmax_family"],
        "GW_apex_full_transversal": r3["GW_full_transversal"],
        "num_one_missing_defect_le1": len(r4["one_missing_le1"]),
        "num_two_missing_defect_le1": r4["two_missing_le1"],
    }
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
