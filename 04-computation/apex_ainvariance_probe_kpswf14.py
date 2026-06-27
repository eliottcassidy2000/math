#!/usr/bin/env python3
"""apex_ainvariance_probe_kpswf14.py

PROBE: is the apex winding tournament's ISO CLASS really invariant under unit
multiplication a in (Z/14)^* ?  The main script flagged that for the AP, a=3 etc
give a NON-isomorphic tournament under the *index* tie-break of antipodal pairs.

This probe isolates WHY, and determines the CORRECT statement:
  - With antipodal pairs (delta==7) tie-broken by INDEX, multiplication by a is
    NOT a relabeling that respects the tie-break (the rotation permutes vertices
    but the index tie-break is anchored to the ORIGINAL labels). So the iso class
    can change. We quantify.
  - The mathematically natural invariant: the tournament is determined up to iso
    once we fix HOW antipodal ties break. If we break ALL antipodal ties the same
    consistent way relative to the residue VALUE (not the index), is the iso class
    a-invariant?  We test a residue-anchored tie-break:
       for delta==7 use: i->j iff residue[i] < residue[j].
    Under a*: delta'(i,j) = a*delta(i,j); delta==7 <=> delta'==7 (a odd). The
    residue-anchored rule is NOT a-invariant either in general, because a* permutes
    residues. The HONEST fact: the *unoriented* antipodal pairs are an a-invariant
    set, but their orientation is a convention; different conventions give possibly
    different iso classes, all sharing score/c3 PARITY structure.

We report, for AP and GW, ALL iso classes obtained as a ranges over units, under
each tie-break convention, with full invariants. This pins the exact dependence.
"""
import math
from math import comb
from itertools import combinations
from collections import defaultdict, Counter

N14 = 14
UNITS14 = [a for a in range(1, N14) if math.gcd(a, N14) == 1]


def adj_index_tiebreak(residues, a=1):
    """antipodal & d==0 ties broken by INDEX (lower->higher)."""
    k = len(residues)
    adj = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            d = ((residues[i] - residues[j]) * a) % N14
            if d == 0:
                adj[i][j] = 1 if i < j else 0
            elif 1 <= d <= 6:
                adj[i][j] = 1
            elif d == 7:
                adj[i][j] = 1 if i < j else 0
            else:
                adj[i][j] = 0
    return adj


def adj_residue_tiebreak(residues, a=1):
    """antipodal & d==0 ties broken by RESIDUE VALUE (lower residue -> higher)."""
    k = len(residues)
    adj = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            d = ((residues[i] - residues[j]) * a) % N14
            if d == 0:
                adj[i][j] = 1 if residues[i] < residues[j] else 0
            elif 1 <= d <= 6:
                adj[i][j] = 1
            elif d == 7:
                adj[i][j] = 1 if residues[i] < residues[j] else 0
            else:
                adj[i][j] = 0
    return adj


def score_seq(adj):
    k = len(adj)
    return tuple(sorted(sum(adj[i][j] for j in range(k)) for i in range(k)))


def c3_count(adj):
    k = len(adj)
    s = [sum(adj[i][j] for j in range(k)) for i in range(k)]
    return comb(k, 3) - sum(x * (x - 1) // 2 for x in s)


def H_held_karp(adj):
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


def iso_exists(A, B):
    k = len(A)
    sA = sorted(sum(A[i][j] for j in range(k)) for i in range(k))
    sB = sorted(sum(B[i][j] for j in range(k)) for i in range(k))
    if sA != sB:
        return False
    degA = [sum(A[i][j] for j in range(k)) for i in range(k)]
    degB = [sum(B[i][j] for j in range(k)) for i in range(k)]

    def refine(adj, init):
        colors = init[:]
        while True:
            sig = {}
            new = [0] * k
            for v in range(k):
                oc = tuple(sorted(colors[w] for w in range(k) if adj[v][w]))
                ic = tuple(sorted(colors[w] for w in range(k) if adj[w][v]))
                key = (colors[v], oc, ic)
                sig.setdefault(key, len(sig))
                new[v] = sig[key]
            if new == colors:
                return colors
            colors = new

    cA = refine(A, degA)
    cB = refine(B, degB)
    if Counter(cA) != Counter(cB):
        return False
    candB = defaultdict(list)
    for w in range(k):
        candB[cB[w]].append(w)
    order = sorted(range(k), key=lambda v: (len(candB[cA[v]]), v))
    mapping = {}
    used = [False] * k

    def ok(v, img):
        for u, iu in mapping.items():
            if A[v][u] != B[img][iu] or A[u][v] != B[iu][img]:
                return False
        return True

    def bt(idx):
        if idx == k:
            return True
        v = order[idx]
        for img in candB[cA[v]]:
            if used[img] or not ok(v, img):
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
    k = len(adj)
    deg = [sum(adj[i][j] for j in range(k)) for i in range(k)]

    def refine(colors):
        while True:
            sig = {}
            new = [0] * k
            for v in range(k):
                oc = tuple(sorted(colors[w] for w in range(k) if adj[v][w]))
                ic = tuple(sorted(colors[w] for w in range(k) if adj[w][v]))
                key = (colors[v], oc, ic)
                sig.setdefault(key, len(sig))
                new[v] = sig[key]
            if new == colors:
                return colors
            colors = new

    col = refine(deg[:])
    cand = {v: [w for w in range(k) if col[w] == col[v]] for v in range(k)}
    order = sorted(range(k), key=lambda v: (len(cand[v]), v))
    mp = {}
    used = [False] * k
    cnt = 0

    def ok(v, img):
        for u, iu in mp.items():
            if adj[v][u] != adj[img][iu] or adj[u][v] != adj[iu][img]:
                return False
        return True

    def bt(idx):
        nonlocal cnt
        if idx == k:
            cnt += 1
            return
        v = order[idx]
        for img in cand[v]:
            if used[img] or not ok(v, img):
                continue
            mp[v] = img
            used[img] = True
            bt(idx + 1)
            del mp[v]
            used[img] = False

    bt(0)
    return cnt


def analyze(name, residues, builder, builder_name):
    print(f"\n--- {name} under {builder_name} tie-break ---")
    classes = []  # list of (adj, [a-values])
    for a in UNITS14:
        adj = builder(residues, a)
        placed = False
        for c in classes:
            if iso_exists(adj, c["adj"]):
                c["a"].append(a)
                placed = True
                break
        if not placed:
            classes.append({"adj": adj, "a": [a]})
    print(f"  # distinct iso classes over a in {UNITS14}: {len(classes)}")
    for ci, c in enumerate(classes):
        adj = c["adj"]
        sc = score_seq(adj)
        print(f"   class {ci}: a-values {c['a']}")
        print(f"     score={sc}")
        print(f"     c3={c3_count(adj)}, H={H_held_karp(adj)}, "
              f"self_conv={iso_exists(adj, [[adj[j][i] for j in range(len(adj))] for i in range(len(adj))])}, "
              f"|Aut|={aut_size(adj)}")
    return classes


def main():
    AP = list(range(1, 14))
    GW = sorted([s % 14 for s in (list(range(1, 12)) + [13, 24])])  # multiset, 10 doubled
    print("AP residues:", AP)
    print("GW residues (multiset):", GW)

    print("\n### Why a-invariance FAILS under index tie-break ###")
    print("Antipodal pairs (delta==7): a is odd => 7*a == 7 mod 14, so the SET of")
    print("antipodal pairs is a-invariant. But the index tie-break orients them by")
    print("a FIXED labeling, while a* permutes positions on the circle. So a* is a")
    print("relabeling of the NON-antipodal arcs but does NOT commute with the")
    print("index-anchored tie-break => possibly different iso class. We quantify:")

    # Count antipodal pairs in AP and GW
    for nm, res in [("AP", AP), ("GW", GW)]:
        anti = [(i, j) for i in range(len(res)) for j in range(i+1, len(res))
                if (res[i]-res[j]) % N14 == 7]
        zero = [(i, j) for i in range(len(res)) for j in range(i+1, len(res))
                if (res[i]-res[j]) % N14 == 0]
        print(f"  {nm}: #antipodal(delta==7) pairs = {len(anti)} {anti}; "
              f"#coincident(delta==0) pairs = {len(zero)} {zero}")

    analyze("AP", AP, adj_index_tiebreak, "INDEX")
    analyze("AP", AP, adj_residue_tiebreak, "RESIDUE")
    analyze("GW", GW, adj_index_tiebreak, "INDEX")
    analyze("GW", GW, adj_residue_tiebreak, "RESIDUE")

    print("\n### KEY: AP has 0 coincident pairs but DOES have antipodal pairs.")
    print("### The 'rotational regular R_13' is the SAME iso class regardless of")
    print("### tie-break ONLY IF the antipodal orientation is itself rotationally")
    print("### symmetric. Check: does AP's score stay (6,...,6) for ALL a, both")
    print("### tie-breaks?  If yes, all a give REGULAR tournaments (score-invariant)")
    print("### even if not all literally isomorphic.")
    for builder, bn in [(adj_index_tiebreak, "INDEX"), (adj_residue_tiebreak, "RESIDUE")]:
        regs = []
        for a in UNITS14:
            adj = builder(AP, a)
            sc = score_seq(adj)
            regs.append((a, len(set(sc)) == 1, sc[0] if len(set(sc)) == 1 else None,
                         c3_count(adj), H_held_karp(adj)))
        print(f"  AP under {bn}: (a, regular?, common-score, c3, H):")
        for r in regs:
            print(f"     {r}")


if __name__ == "__main__":
    main()
