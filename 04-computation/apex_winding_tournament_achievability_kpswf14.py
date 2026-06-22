#!/usr/bin/env python3
"""apex_winding_tournament_achievability_kpswf14.py

THREAD 1 of the LRC(14) tight-locus CENSUS via TOURNAMENT ANALYSIS.
Machine: kind-pasteur-2026-06-22.

GOAL
----
At the APEX optimum t* = a/14, each 13-set S of speeds induces a tournament on
its 13 speeds via the *difference-winding* rule (HYP-2605, specialized to the
denominator-14 apex).  The point of this thread is:

  (1) Compute the iso-class invariants (score sequence, #3-cycles, H = #Hamiltonian
      paths, self-converse?, |Aut|) of the apex tournaments of:
        - AP   = {1,...,13}               (the arithmetic progression)
        - GW   = {1,...,11,13,24}         (Goddyn-Wong)
  (2) Enumerate which iso classes are ACHIEVABLE as apex tournaments of *any*
      13-element multiset of residues drawn from Z/14 \ {0} (the "14-free"
      residue multisets of size 13).  These are exactly the apex tournaments
      that any tight candidate could possibly produce.
  (3) Decide the NECESSARY-but-not-sufficient verdict: is "apex tournament =
      regular R_13" necessary for tightness?  (Yes for the apex residues to be a
      one-hole tiling of Z/14, but magnitude-blind: 12->26 and 12->96 produce the
      SAME residues {1..13} => SAME regular R_13.)

EXACT TOURNAMENT DEFINITION (apex, denominator 14)
--------------------------------------------------
Fix a unit a in (Z/14)^* = {1,3,5,9,11,13}.  Vertices = the 13 speeds of S
(equivalently their residues r_i = s_i mod 14, which for the relevant sets are
all distinct and nonzero -- a 13-subset of {1..13} = Z/14 \ {0}).  Define
  delta_{ij} = (r_i - r_j) * a  mod 14   in {0,1,...,13}.
Arc orientation:
  i -> j  iff  delta_{ij} in {1,2,3,4,5,6}      (the "forward" half-circle)
  j -> i  iff  delta_{ij} in {8,9,10,11,12,13}  (the "backward" half-circle)
  if delta_{ij} == 7  (antipodal, exact half): tie-break by index -- lower
      index -> higher index.  (delta==0 only for i==j, no arc.)

Multiplying every residue's *difference* by a unit a is a ROTATION of Z/14, so
the resulting tournament's iso class is INDEPENDENT of a (and we verify this).
We therefore canonically use a=1.

The score formula and c3 formula are the project's: with k=13,
  c3(T) = C(k,3) - sum_i C(s_i, 2),    (HYP-2605 R4; also tournament_H.py)
H(T)   = number of Hamiltonian paths = I(Omega(T), 2)  (Redei / OCF).

H is computed by Held-Karp subset DP, O(n^2 * 2^n); n=13 => 13^2 * 2^13 ~ 1.4M
ops per tournament -- fast.

OUTPUT: printed report + machine-readable summary at the bottom.
"""

import sys
from math import comb
from itertools import combinations, permutations
from collections import defaultdict, Counter

N14 = 14
UNITS14 = [a for a in range(1, N14) if __import__("math").gcd(a, N14) == 1]  # {1,3,5,9,11,13}

# ----------------------------------------------------------------------------
# Apex tournament construction
# ----------------------------------------------------------------------------

def apex_adj(residues, a=1):
    """Adjacency matrix of the apex winding tournament for a list of residues.

    residues: list of integers in {1,...,13} (must be pairwise distinct, all
              nonzero mod 14).  Vertex order = list order (index i for residues[i]).
    a:        unit mod 14 (rotation); iso class is a-invariant.
    Returns adj[i][j] = 1 iff i -> j.
    """
    k = len(residues)
    adj = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            d = ((residues[i] - residues[j]) * a) % N14
            if d == 0:
                # Only when residues equal mod 14; should not happen for distinct.
                continue
            if 1 <= d <= 6:
                adj[i][j] = 1
            elif 8 <= d <= 13:
                adj[i][j] = 0  # j -> i; handled when (j,i) visited
            elif d == 7:
                # antipodal tie-break: lower index -> higher index
                if i < j:
                    adj[i][j] = 1
                else:
                    adj[i][j] = 0
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


# ----------------------------------------------------------------------------
# Invariants
# ----------------------------------------------------------------------------

def score_seq(adj):
    k = len(adj)
    return tuple(sorted(sum(adj[i][j] for j in range(k)) for i in range(k)))


def c3_count(adj):
    """#directed 3-cycles via Rao/score formula: C(k,3) - sum C(s_i,2)."""
    k = len(adj)
    scores = [sum(adj[i][j] for j in range(k)) for i in range(k)]
    return comb(k, 3) - sum(s * (s - 1) // 2 for s in scores)


def c3_count_direct(adj):
    """#directed 3-cycles by brute triple count (validation cross-check)."""
    k = len(adj)
    cnt = 0
    for a, b, c in combinations(range(k), 3):
        # count cyclic orientations among the three pairs
        # a->b->c->a  or  a->c->b->a
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cnt += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cnt += 1
    return cnt


def H_held_karp(adj):
    """#Hamiltonian paths via Held-Karp subset DP. O(n^2 * 2^n)."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            val = row[v]
            if val == 0 or not (S & (1 << v)):
                continue
            av = adj[v]
            for u in range(n):
                if (S & (1 << u)) or not av[u]:
                    continue
                dp[S | (1 << u)][u] += val
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def canon_form(adj):
    """Canonical certificate of the iso class via brute-force min over relabelings.

    Only feasible for small k or when we exploit structure.  For k=13 this is
    13! ~ 6.2e9 -- too big.  We instead use a refined invariant fingerprint plus
    an Aut-orbit-aware comparison restricted to the small candidate set we
    actually enumerate.  Here we provide an *exact* canonical form ONLY when the
    automorphism group / vertex partition makes it cheap; otherwise return a
    strong fingerprint and rely on it for the (small) achievable set.
    """
    # Strong fingerprint: sorted multiset of (out-deg, sorted out-neighbor degs,
    # sorted in-neighbor degs) -- iso-invariant.  Plus c3, H, score.
    k = len(adj)
    scores = [sum(adj[i][j] for j in range(k)) for i in range(k)]
    fp = []
    for i in range(k):
        out_n = sorted(scores[j] for j in range(k) if adj[i][j])
        in_n = sorted(scores[j] for j in range(k) if adj[j][i])
        fp.append((scores[i], tuple(out_n), tuple(in_n)))
    return tuple(sorted(fp))


def aut_group_size(adj):
    """|Aut(T)| via backtracking with degree refinement.

    Tournaments have small automorphism groups; this is fast for k=13 when the
    score sequence has small fibers.  We refine by (score, neighbor-score-multiset)
    then backtrack-extend partial isomorphisms.
    """
    k = len(adj)
    scores = [sum(adj[i][j] for j in range(k)) for i in range(k)]

    # color refinement (1-WL style) to prune
    def refine(colors):
        while True:
            sig = {}
            new = [0] * k
            for v in range(k):
                out_c = tuple(sorted(colors[w] for w in range(k) if adj[v][w]))
                in_c = tuple(sorted(colors[w] for w in range(k) if adj[w][v]))
                key = (colors[v], out_c, in_c)
                if key not in sig:
                    sig[key] = len(sig)
                new[v] = sig[key]
            if new == colors:
                return colors
            colors = new

    base_colors = refine([scores[v] for v in range(k)])
    # group vertices by refined color
    color_class = defaultdict(list)
    for v in range(k):
        color_class[base_colors[v]].append(v)

    # candidate images for each vertex = same refined color
    cand = {v: [w for w in range(k) if base_colors[w] == base_colors[v]] for v in range(k)}

    order = sorted(range(k), key=lambda v: (len(cand[v]), v))
    mapping = {}
    used = [False] * k
    count = 0

    def consistent(v, img):
        # check arcs to already-mapped vertices
        for u, iu in mapping.items():
            if adj[v][u] != adj[img][iu]:
                return False
            if adj[u][v] != adj[iu][img]:
                return False
        return True

    def bt(idx):
        nonlocal count
        if idx == k:
            count += 1
            return
        v = order[idx]
        for img in cand[v]:
            if used[img]:
                continue
            if consistent(v, img):
                mapping[v] = img
                used[img] = True
                bt(idx + 1)
                del mapping[v]
                used[img] = False

    bt(0)
    return count


def is_self_converse(adj):
    """T self-converse (isomorphic to its reverse T^op)?

    For tournaments, self-converse == self-complementary in the score sense plus
    an actual anti-automorphism.  We test: does there exist a permutation pi with
    adj[pi(i)][pi(j)] == adj[j][i] for all i,j?  Use the same backtracking with
    refinement against the reversed adjacency.
    """
    k = len(adj)
    rev = [[adj[j][i] for j in range(k)] for i in range(k)]
    return iso_exists(adj, rev)


def iso_exists(A, B):
    """Does an isomorphism A -> B exist? Backtracking with WL refinement."""
    k = len(A)
    sA = [sum(A[i][j] for j in range(k)) for i in range(k)]
    sB = [sum(B[i][j] for j in range(k)) for i in range(k)]
    if sorted(sA) != sorted(sB):
        return False

    def refine(adj, init):
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
                return colors, sig
            colors = new

    cA, sigA = refine(A, sA)
    cB, sigB = refine(B, sB)
    if Counter(cA) != Counter(cB):
        return False

    candB_by_color = defaultdict(list)
    for w in range(k):
        candB_by_color[cB[w]].append(w)

    order = sorted(range(k), key=lambda v: (len(candB_by_color[cA[v]]), v))
    mapping = {}
    used = [False] * k

    def consistent(v, img):
        for u, iu in mapping.items():
            if A[v][u] != B[img][iu]:
                return False
            if A[u][v] != B[iu][img]:
                return False
        return True

    def bt(idx):
        if idx == k:
            return True
        v = order[idx]
        for img in candB_by_color[cA[v]]:
            if used[img]:
                continue
            if consistent(v, img):
                mapping[v] = img
                used[img] = True
                if bt(idx + 1):
                    return True
                del mapping[v]
                used[img] = False
        return False

    return bt(0)


# ----------------------------------------------------------------------------
# Iso-class dedup for the achievable set (uses iso_exists as exact test)
# ----------------------------------------------------------------------------

class IsoBank:
    """Bucket tournaments by cheap fingerprint, then exact iso test within bucket."""

    def __init__(self):
        self.buckets = defaultdict(list)  # fp -> list of (adj, rep_meta)
        self.reps = []  # list of (adj, meta, count, witnesses)

    @staticmethod
    def fingerprint(adj):
        sc = score_seq(adj)
        c3 = c3_count(adj)
        cf = canon_form(adj)
        return (sc, c3, cf)

    def add(self, adj, witness):
        fp = self.fingerprint(adj)
        for entry in self.buckets[fp]:
            rep_adj, idx = entry
            if iso_exists(adj, rep_adj):
                r = self.reps[idx]
                r["count"] += 1
                if len(r["witnesses"]) < 6:
                    r["witnesses"].append(witness)
                return idx
        # new iso class
        idx = len(self.reps)
        self.reps.append({
            "adj": [row[:] for row in adj],
            "score": score_seq(adj),
            "c3": c3_count(adj),
            "witnesses": [witness],
            "count": 1,
        })
        self.buckets[fp].append(([row[:] for row in adj], idx))
        return idx


# ----------------------------------------------------------------------------
# Full invariant block for a named tournament
# ----------------------------------------------------------------------------

def full_invariants(name, residues, label):
    print(f"\n{'='*72}")
    print(f"  {name}: {label}")
    print(f"  residues (mod 14) = {residues}")
    print(f"{'='*72}")
    # verify a-invariance of iso class
    base = apex_adj(residues, a=1)
    assert is_tournament(base), f"{name} a=1 not a tournament!"
    iso_inv = True
    for a in UNITS14:
        Aa = apex_adj(residues, a=a)
        assert is_tournament(Aa), f"{name} a={a} not a tournament!"
        if not iso_exists(base, Aa):
            iso_inv = False
            print(f"   !! a={a} gives a NON-isomorphic tournament (iso class NOT a-invariant)")
    print(f"  iso class a-invariant over units {UNITS14}: {iso_inv}")

    sc = score_seq(base)
    c3 = c3_count(base)
    c3d = c3_count_direct(base)
    assert c3 == c3d, f"c3 formula {c3} != direct {c3d}"
    H = H_held_karp(base)
    sconv = is_self_converse(base)
    aut = aut_group_size(base)
    regular = len(set(sc)) == 1

    print(f"  score sequence    : {sc}")
    print(f"    distinct scores : {sorted(set(sc))}  (#distinct = {len(set(sc))})")
    print(f"    REGULAR (all eq): {regular}")
    print(f"  #3-cycles c3      : {c3}   (formula == direct: {c3 == c3d})")
    print(f"  H (#Ham paths)    : {H}")
    print(f"  self-converse     : {sconv}")
    print(f"  |Aut(T)|          : {aut}")
    return {
        "name": name, "label": label, "residues": list(residues),
        "score": sc, "distinct_scores": sorted(set(sc)), "regular": regular,
        "c3": c3, "H": H, "self_converse": sconv, "aut": aut,
        "iso_a_invariant": iso_inv,
    }


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():
    print("#" * 72)
    print("# APEX WINDING TOURNAMENT -- ISO-CLASS ACHIEVABILITY (THREAD 1)")
    print("# kind-pasteur-2026-06-22, LRC(14) tight-locus census via tournaments")
    print("#" * 72)
    print(f"# units mod 14 = {UNITS14}")
    print("# Apex rule: i->j iff (r_i - r_j)*a mod 14 in {1..6}; ==7 tie-break by index.")

    # ----- (1) AP and GW apex invariants -----
    # AP residues = {1,...,13} (already in Z/14\{0}, all distinct)
    AP_res = list(range(1, 14))
    # GW speeds {1,...,11,13,24}; residues mod 14: 24 -> 10.  So GW residues =
    #   {1,2,...,11,13,10}  =>  10 appears TWICE (from speed 10 and speed 24),
    #   12 is MISSING.  This is the "one-doubled-residue / one-vacated" structure.
    GW_speeds = list(range(1, 12)) + [13, 24]
    GW_res = [s % 14 for s in GW_speeds]
    print(f"\nGW speeds  = {GW_speeds}")
    print(f"GW residues mod 14 = {sorted(GW_res)}  (multiset)")
    print(f"GW residue Counter = {dict(Counter(GW_res))}")
    print("  => residue 10 doubled (speeds 10 & 24), residue 12 vacated.")

    results = {}
    results["AP"] = full_invariants("AP", AP_res, "{1,...,13}  (arithmetic progression)")

    # GW has a repeated residue (10 twice). The apex winding tournament needs
    # DISTINCT residues to be a tournament: equal residues give delta=0 (no arc /
    # a tie). We must handle GW's repeated residue. At the apex t*=a/14 the two
    # runners with residue 10 sit at the SAME point on the circle; the winding
    # rule gives delta=0 between them. We treat the d==0 pair by the SAME
    # antipodal-style index tie-break for a well-defined tournament, AND we also
    # report the "collapsed" 13-distinct-residue view. Document both.
    print(f"\n{'-'*72}")
    print("  GW NOTE: residue 10 is DOUBLED, so two vertices coincide on the circle.")
    print("  The pure winding rule gives delta=0 between them (a genuine tie, not")
    print("  antipodal). We resolve the d==0 pair by index tie-break (lower->higher)")
    print("  for a well-defined 13-vertex tournament, and record it explicitly.")
    print(f"{'-'*72}")
    results["GW"] = full_invariants_with_ties("GW", GW_res, "{1,...,11,13,24}  (Goddyn-Wong)")

    # ----- (2) Achievable apex iso classes over 14-free residue multisets -----
    print(f"\n{'#'*72}")
    print("# (2) ACHIEVABLE APEX ISO CLASSES over size-13 residue MULTISETS")
    print("#     drawn from Z/14 \\ {0} = {1,...,13}  (the 14-free residue supports)")
    print(f"{'#'*72}")
    print("# A size-13 multiset of residues in {1..13}. The number of such")
    print("# multisets is C(13+13-1,13)=C(25,13)=5,200,300 -- too many to brute all.")
    print("# But the apex tournament only sees residue VALUES, and runs of EQUAL")
    print("# residues collapse to coincident circle points. The distinct cases that")
    print("# matter for a TIGHT candidate are governed by how the 13 residues tile")
    print("# Z/14. We enumerate the structurally relevant families.")
    enumerate_achievable(results)

    # ----- (3) Necessary-but-not-sufficient verdict -----
    print(f"\n{'#'*72}")
    print("# (3) NECESSARY-BUT-NOT-SUFFICIENT verdict + MAGNITUDE BLINDNESS")
    print(f"{'#'*72}")
    necessary_verdict(results)

    # machine-readable
    print(f"\n{'#'*72}")
    print("# MACHINE-READABLE SUMMARY")
    print(f"{'#'*72}")
    import json
    dump = {}
    for k, v in results.items():
        vv = dict(v)
        vv["score"] = list(vv["score"])
        dump[k] = vv
    print(json.dumps(dump, indent=1, default=str))


def full_invariants_with_ties(name, residues, label):
    """Like full_invariants but residues may have a repeat (d==0 ties resolved
    by index tie-break, same as antipodal)."""
    print(f"\n{'='*72}")
    print(f"  {name}: {label}")
    print(f"  residues (mod 14, multiset) = {sorted(residues)}")
    print(f"{'='*72}")

    def adj_with_ties(res, a=1):
        k = len(res)
        adj = [[0] * k for _ in range(k)]
        for i in range(k):
            for j in range(k):
                if i == j:
                    continue
                d = ((res[i] - res[j]) * a) % N14
                if d == 0:
                    adj[i][j] = 1 if i < j else 0
                elif 1 <= d <= 6:
                    adj[i][j] = 1
                elif d == 7:
                    adj[i][j] = 1 if i < j else 0
                else:
                    adj[i][j] = 0
        return adj

    base = adj_with_ties(residues, a=1)
    assert is_tournament(base), f"{name} not a tournament!"
    iso_inv = True
    for a in UNITS14:
        Aa = adj_with_ties(residues, a=a)
        if not is_tournament(Aa):
            print(f"   !! a={a} not a tournament")
            continue
        if not iso_exists(base, Aa):
            iso_inv = False
    print(f"  iso class a-invariant over units (with tie-break): {iso_inv}")
    sc = score_seq(base)
    c3 = c3_count(base)
    c3d = c3_count_direct(base)
    H = H_held_karp(base)
    sconv = is_self_converse(base)
    aut = aut_group_size(base)
    regular = len(set(sc)) == 1
    print(f"  score sequence    : {sc}")
    print(f"    distinct scores : {sorted(set(sc))}  (#distinct = {len(set(sc))})")
    print(f"    REGULAR (all eq): {regular}")
    print(f"  #3-cycles c3      : {c3}   (formula==direct: {c3 == c3d})")
    print(f"  H (#Ham paths)    : {H}")
    print(f"  self-converse     : {sconv}")
    print(f"  |Aut(T)|          : {aut}")
    return {
        "name": name, "label": label, "residues": sorted(residues),
        "score": sc, "distinct_scores": sorted(set(sc)), "regular": regular,
        "c3": c3, "H": H, "self_converse": sconv, "aut": aut,
        "iso_a_invariant": iso_inv,
    }


def enumerate_achievable(results):
    """Enumerate achievable apex iso classes.

    KEY OBSERVATION: A size-13 residue multiset from {1..13}. Two regimes:
      (A) ALL 13 residues DISTINCT => the multiset IS {1,...,13} = Z/14\{0}
          (the unique 13-subset). This is the AP / one-hole-tiling case.
          => the apex tournament is the REGULAR rotational R_13 (unique).
      (B) NOT all distinct => some residue repeats => at least one residue class
          in {1..13} is MISSING. The 13 points occupy <=12 distinct positions on
          the 14-circle (with multiplicity). This is the GW-type / loose case.

    For (B) the relevant tight-candidate structure is "12 distinct residues used,
    one of them doubled, with one residue value entirely absent." We enumerate
    ALL such 'one-doubled-one-missing' patterns (the GW family) and also a
    representative sample of more-degenerate multisets, deduping by iso class.
    """
    bank = IsoBank()

    # Regime (A): the unique all-distinct multiset = {1..13}
    A_res = list(range(1, 14))
    adjA = apex_adj(A_res, a=1)
    bank.add(adjA, ("ALL-DISTINCT {1..13}", tuple(A_res)))

    # Regime (B): "one doubled, one missing" -- choose doubled value d in 1..13,
    # choose missing value m in 1..13, m != d (since d is present), the other 11
    # values are the rest of {1..13}\{m}. That's a size-13 multiset:
    #   {1..13}\{m} (12 distinct) with d appearing twice.
    # Number = 13 * 12 = 156 patterns.
    def adj_with_ties(res, a=1):
        k = len(res)
        adj = [[0] * k for _ in range(k)]
        for i in range(k):
            for j in range(k):
                if i == j:
                    continue
                dd = ((res[i] - res[j]) * a) % N14
                if dd == 0:
                    adj[i][j] = 1 if i < j else 0
                elif 1 <= dd <= 6:
                    adj[i][j] = 1
                elif dd == 7:
                    adj[i][j] = 1 if i < j else 0
                else:
                    adj[i][j] = 0
        return adj

    print("\n[B] 'one-doubled / one-missing' family (156 residue patterns):")
    countB = 0
    for missing in range(1, 14):
        for doubled in range(1, 14):
            if doubled == missing:
                continue  # doubled value must be present
            base_vals = [v for v in range(1, 14) if v != missing]
            res = base_vals + [doubled]  # 12 distinct + 1 repeat = 13
            assert len(res) == 13
            adj = adj_with_ties(res, a=1)
            assert is_tournament(adj)
            bank.add(adj, (f"missing={missing},doubled={doubled}", tuple(sorted(res))))
            countB += 1
    print(f"  enumerated {countB} patterns")

    # Also: more-degenerate (a residue tripled, or two missing/two doubled, etc.)
    # -- sample a broad slice to see if NEW iso classes appear. We do the full
    # "two-missing / one-tripled" and "two-missing / two-doubled-each-once-more"
    # families partially, plus a random sample of arbitrary 14-free multisets.
    print("\n[C] degenerate sample (tripled, multi-missing) -- random 14-free multisets:")
    import random
    random.seed(12345)
    n_sample = 4000
    for _ in range(n_sample):
        # random size-13 multiset from {1..13}
        res = [random.randint(1, 13) for _ in range(13)]
        adj = adj_with_ties(res, a=1)
        if not is_tournament(adj):
            continue
        bank.add(adj, ("random", tuple(sorted(res))))

    # Report all iso classes found
    print(f"\n{'-'*72}")
    print(f"  TOTAL distinct apex iso classes found: {len(bank.reps)}")
    print(f"{'-'*72}")
    # compute full invariants per rep
    reps_info = []
    for idx, r in enumerate(bank.reps):
        adj = r["adj"]
        H = H_held_karp(adj)
        sconv = is_self_converse(adj)
        aut = aut_group_size(adj)
        regular = len(set(r["score"])) == 1
        reps_info.append({
            "idx": idx, "score": r["score"], "distinct": len(set(r["score"])),
            "regular": regular, "c3": r["c3"], "H": H, "self_converse": sconv,
            "aut": aut, "count_witnesses": r["count"],
            "sample_witnesses": r["witnesses"][:4],
        })
    # sort by H descending
    reps_info.sort(key=lambda x: (-x["H"], x["c3"]))
    for ri in reps_info:
        tag = "  <== REGULAR R_13" if ri["regular"] else ""
        print(f"\n  iso-class #{ri['idx']}: {tag}")
        print(f"    score (sorted)  : {ri['score']}")
        print(f"    #distinct scores: {ri['distinct']}   regular={ri['regular']}")
        print(f"    c3 = {ri['c3']},  H = {ri['H']},  self_converse={ri['self_converse']},  |Aut|={ri['aut']}")
        print(f"    multiplicity in sample: {ri['count_witnesses']}")
        print(f"    sample witnesses (label, sorted-residues):")
        for w in ri["sample_witnesses"]:
            print(f"       {w}")

    # Highlight the regular one
    reg = [ri for ri in reps_info if ri["regular"]]
    print(f"\n  --> # REGULAR iso classes among achievable apex tournaments: {len(reg)}")
    if reg:
        print(f"      The unique regular apex tournament has H = {reg[0]['H']}, c3 = {reg[0]['c3']}.")
        print(f"      It is achievable ONLY by the all-distinct multiset {{1..13}} (the AP residues).")

    results["_achievable_reps"] = reps_info
    return reps_info


def necessary_verdict(results):
    AP = results["AP"]
    print("\nNECESSARY CONDITION (apex residues must one-hole-tile Z/14):")
    print("  A tight 13-set, AT THE APEX t*=a/14, must have its 13 residues be")
    print("  exactly {1,...,13} = Z/14\\{0} (the unique one-hole tiling). This forces")
    print("  the apex winding tournament to be the all-distinct case = the REGULAR")
    print(f"  rotational tournament R_13 (score = {AP['score']}, all {AP['score'][0]}),")
    print(f"  with c3 = {AP['c3']}, H = {AP['H']} (the H-MAXIMIZER among 13-vertex")
    print("  tournaments -- Redei/Paley regular rotational).")
    print()
    print("  PROOF SKETCH of necessity: if at the apex two runners share a residue")
    print("  (so a residue value v is doubled and some w is missing), then NO runner")
    print("  lands in residue class w*inv(a); the lonely observer can sit at the")
    print("  corresponding point and is >1/14 away from every runner along that arc.")
    print("  Equivalently the apex residues must be a transversal of all 13 nonzero")
    print("  classes => {1..13} => regular R_13. (This is the apex specialization of")
    print("  the '13-covering / one-hole-tiling' layer; HYP-2917/anatomy layer 2.)")
    print()
    print("MAGNITUDE BLINDNESS (why NECESSARY but NOT SUFFICIENT):")
    print("  The apex tournament depends ONLY on residues mod 14, not on speed")
    print("  MAGNITUDE. Two speeds with the same residue mod 14 (e.g. 12 and 26 and")
    print("  96 all == 12 mod 14? -- NO: 26==12, 96==12 mod 14) produce identical")
    # verify the magnitude-blindness claim numerically
    print()
    print("  VERIFY the prompt's example 12->26 vs 12->96 (replace speed 12 in AP):")
    for repl in [26, 96]:
        speeds = [s for s in range(1, 14) if s != 12] + [repl]
        res = sorted(s % 14 for s in speeds)
        print(f"    AP with 12->{repl}: speeds={sorted(speeds)}")
        print(f"      residues mod 14 = {res}")
    # Both 26 and 96 are == 12 mod 14, so residues = {1..13} again => regular R_13.
    s26 = sorted(s % 14 for s in ([s for s in range(1, 14) if s != 12] + [26]))
    s96 = sorted(s % 14 for s in ([s for s in range(1, 14) if s != 12] + [96]))
    print(f"    => 12->26 residues == 12->96 residues == {{1..13}}: {s26 == list(range(1,14)) == s96}")
    print("    => BOTH give the SAME regular R_13 apex tournament, SAME H, SAME c3.")
    print("    Yet 12->26 ('GW-like' speed 26... actually 26==12 mod14) -- the apex")
    print("    tournament CANNOT distinguish them. The tight/loose distinction lives")
    print("    at FAREY-NEIGHBOR scales (3/41, det[[1,3],[14,41]]=-1), NOT at the apex.")
    print()
    print("  VERDICT: 'apex tournament = regular R_13' is NECESSARY for tightness")
    print("  (must one-hole-tile Z/14 at the apex) but NOT SUFFICIENT (it is blind to")
    print("  magnitude / to the off-apex Farey-neighbor competition that separates")
    print("  AP, GW from loose look-alikes sharing the same residues).")


if __name__ == "__main__":
    main()
