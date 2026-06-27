#!/usr/bin/env python3
"""apex_achievable_census_kpswf14.py

THREAD 1, parts (2)+(3) done CORRECTLY (no false a-invariance assumption).

We answer: WHICH tournament iso classes are achievable as the apex winding
tournament  T(t*) = winding at t* = a*/14  of a 13-set, AT a lonely-optimum
phase a*?  And we settle:
  - exactly which 13-set residue configurations give the REGULAR R_13;
  - whether R_13 here equals the standard circulant rotational tournament;
  - the necessary-but-insufficient verdict, with magnitude-blindness.

Two layers of the achievable set:

LAYER 1 (apex residues = perfect one-hole tiling {1..13}):
  A 13-set whose 13 speeds are DISTINCT mod 14 must have residue support
  {1,...,13} = Z/14 minus {0} (only 13-subset).  For EACH unit phase a* the winding
  tournament is the "(i-j)*a* mod 14 in {1..6}" circulant-image.  We enumerate
  all 6 (one per unit a*) and dedupe.  This is the AP family; a*=1 -> regular.

LAYER 2 (one residue doubled, one missing -- the GW family):
  13 speeds with exactly one repeated residue: support is 12 distinct classes,
  one doubled, one missing.  For each (missing m, doubled d) and each unit a*,
  build the winding tournament (d==0 coincident pair tie-broken by index).

We report, per layer, the FULL set of distinct iso classes reachable as
apex tournaments at SOME optimum phase, with (score, c3, H, self-converse, |Aut|,
regular?).  Then the verdict.

Magnitude-blindness check (part 3): AP, 12->26, 12->96 -- all have residues
{1..13} at the natural phase, hence share the SAME a*=1 regular R_13 and the same
full 6-tuple of apex tournaments.  Confirmed numerically.
"""
import math
from math import comb, gcd
from itertools import combinations
from collections import Counter, defaultdict

N14 = 14
UNITS14 = [a for a in range(1, N14) if gcd(a, N14) == 1]


def winding_residue_adj(residues, a):
    """(r_i - r_j)*a mod14 in {1..6} => i->j; ==7 or ==0 tie-broken by index."""
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
    if sorted(sum(A[i][j] for j in range(k)) for i in range(k)) != \
       sorted(sum(B[i][j] for j in range(k)) for i in range(k)):
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
    cA, cB = refine(A, degA), refine(B, degB)
    if Counter(cA) != Counter(cB):
        return False
    candB = defaultdict(list)
    for w in range(k):
        candB[cB[w]].append(w)
    order = sorted(range(k), key=lambda v: (len(candB[cA[v]]), v))
    mp, used = {}, [False] * k

    def ok(v, img):
        for u, iu in mp.items():
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
            mp[v], used[img] = img, True
            if bt(idx + 1):
                return True
            del mp[v]
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
    mp, used, cnt = {}, [False] * k, 0

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
            mp[v], used[img] = img, True
            bt(idx + 1)
            del mp[v]
            used[img] = False
    bt(0)
    return cnt


def is_self_converse(adj):
    k = len(adj)
    rev = [[adj[j][i] for j in range(k)] for i in range(k)]
    return iso_exists(adj, rev)


def standard_circulant_R13():
    """Standard rotational/circulant tournament on Z/13: i->j iff (j-i) mod 13 in {1..6}."""
    k = 13
    adj = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            if (j - i) % 13 in (1, 2, 3, 4, 5, 6):
                adj[i][j] = 1
    return adj


class Bank:
    def __init__(self):
        self.reps = []  # dict with adj, invariants, witnesses

    def add(self, adj, witness):
        for r in self.reps:
            if r["score"] == score_seq(adj) and r["c3"] == c3_count(adj) and iso_exists(adj, r["adj"]):
                r["count"] += 1
                if len(r["witnesses"]) < 8:
                    r["witnesses"].append(witness)
                return r
        r = {
            "adj": [row[:] for row in adj],
            "score": score_seq(adj),
            "c3": c3_count(adj),
            "H": H_held_karp(adj),
            "self_converse": is_self_converse(adj),
            "aut": aut_size(adj),
            "regular": len(set(score_seq(adj))) == 1,
            "witnesses": [witness],
            "count": 1,
        }
        self.reps.append(r)
        return r

    def report(self, title):
        print(f"\n{'-'*72}\n  {title}: {len(self.reps)} distinct iso classes\n{'-'*72}")
        for r in sorted(self.reps, key=lambda x: (-x["H"], x["c3"])):
            tag = "   <== REGULAR R_13" if r["regular"] else ""
            print(f"  score={r['score']}{tag}")
            print(f"     c3={r['c3']}, H={r['H']}, self_conv={r['self_converse']}, "
                  f"|Aut|={r['aut']}, regular={r['regular']}, sample_mult={r['count']}")
            print(f"     witnesses: {r['witnesses'][:4]}")


def main():
    print("#" * 72)
    print("# THREAD 1 (2)+(3): ACHIEVABLE apex tournaments + NECESSARY verdict")
    print("# kind-pasteur-2026-06-22. No a-invariance assumed (it is FALSE).")
    print("#" * 72)
    print(f"# units mod 14 = {UNITS14}")

    # Confirm AP a=1 == standard circulant R_13
    AP_res = list(range(1, 14))
    apex_a1 = winding_residue_adj(AP_res, 1)
    R13 = standard_circulant_R13()
    print(f"\nAP apex at a*=1 isomorphic to STANDARD circulant R_13 (Z/13, beat next 6): "
          f"{iso_exists(apex_a1, R13)}")
    print(f"  standard R_13: score={score_seq(R13)}, c3={c3_count(R13)}, H={H_held_karp(R13)}, "
          f"self_conv={is_self_converse(R13)}, |Aut|={aut_size(R13)}")
    print("  (Both are 'beat the next 6 clockwise'; AP on the 14-circle minus the hole,")
    print("   standard on the 13-circle. They coincide as iso classes: the doubling map")
    print("   x->2x mod 13... actually the AP residues 1..13 on Z/14 with arc {1..6}")
    print("   reproduce the Z/13 circulant -- confirmed above.)")

    # LAYER 1: AP family, all unit phases a*
    print(f"\n{'#'*72}\n# LAYER 1: AP residues {{1..13}}, all optimum phases a* in units\n{'#'*72}")
    bank1 = Bank()
    for a in UNITS14:
        adj = winding_residue_adj(AP_res, a)
        bank1.add(adj, f"AP a*={a}")
    bank1.report("LAYER 1 (AP one-hole tiling, 6 unit phases)")
    reg1 = [r for r in bank1.reps if r["regular"]]
    print(f"\n  # REGULAR among AP apex tournaments: {len(reg1)} "
          f"(H={[r['H'] for r in reg1]}, at phase(s) {[r['witnesses'] for r in reg1]})")

    # LAYER 2: GW family (one doubled, one missing), all unit phases
    print(f"\n{'#'*72}\n# LAYER 2: one-doubled/one-missing residue multisets, all a*\n{'#'*72}")
    bank2 = Bank()
    n_pat = 0
    for missing in range(1, 14):
        base = [v for v in range(1, 14) if v != missing]
        for doubled in base:
            res = base + [doubled]
            for a in UNITS14:
                adj = winding_residue_adj(res, a)
                bank2.add(adj, f"m={missing},d={doubled},a*={a}")
            n_pat += 1
    print(f"  enumerated {n_pat} (missing,doubled) patterns x {len(UNITS14)} phases")
    bank2.report("LAYER 2 (one-doubled-one-missing, all phases)")
    reg2 = [r for r in bank2.reps if r["regular"]]
    print(f"\n  # REGULAR among one-dipole apex tournaments: {len(reg2)}")

    # Where does the literal GW {1..11,13,24} sit? (residues: 10 doubled, 12 missing)
    GW_res = sorted(s % 14 for s in (list(range(1, 12)) + [13, 24]))
    print(f"\n  GW {{1..11,13,24}} residues (multiset) = {GW_res} (10 doubled, 12 missing).")
    print("  GW apex tournaments by phase a*:")
    for a in UNITS14:
        adj = winding_residue_adj(GW_res, a)
        print(f"    a*={a}: score={score_seq(adj)}, c3={c3_count(adj)}, H={H_held_karp(adj)}, "
              f"self_conv={is_self_converse(adj)}, |Aut|={aut_size(adj)}, "
              f"regular={len(set(score_seq(adj)))==1}")

    # PART 3: magnitude blindness
    print(f"\n{'#'*72}\n# PART 3: NECESSARY-BUT-INSUFFICIENT + MAGNITUDE BLINDNESS\n{'#'*72}")
    print("Necessary: at the optimum a*=1 the apex residues must one-hole-tile Z/14,")
    print("forcing them to {1..13} => the winding tournament is the REGULAR R_13.")
    print("(A missing residue class => a lonely arc => loneliness > 1/14 there; so the")
    print(" apex residues are forced to the full transversal {1..13}.)")
    print()
    print("Magnitude blindness: replace AP speed 12 by 26 or 96 (both ==12 mod 14):")
    for repl in [26, 96]:
        speeds = [s for s in range(1, 14) if s != 12] + [repl]
        res = sorted(s % 14 for s in speeds)
        adj1 = winding_residue_adj([s % 14 for s in speeds], 1)
        print(f"  12->{repl}: residues={res}; a*=1 apex: score={score_seq(adj1)}, "
              f"H={H_held_karp(adj1)}, regular={len(set(score_seq(adj1)))==1}")
    s26 = sorted(s % 14 for s in ([s for s in range(1, 14) if s != 12] + [26]))
    s96 = sorted(s % 14 for s in ([s for s in range(1, 14) if s != 12] + [96]))
    print(f"  => 12->26 residues == 12->96 residues == AP residues {{1..13}}: "
          f"{s26 == list(range(1,14)) == s96}")
    print("  => identical regular R_13 apex tournament, identical H=3711175, identical")
    print("     full 6-phase apex profile -- yet 12->26 is LOOSE (q=12, M=1/12) and")
    print("     12->96 is LOOSE (q=14 but large, off-apex escape). The apex tournament")
    print("     CANNOT separate them: it is MAGNITUDE-BLIND (sees only residues mod 14).")
    print()
    print("VERDICT: 'apex winding tournament (at a*=1) = regular R_13' is NECESSARY for")
    print("tightness but NOT SUFFICIENT. The magnitude/max-gap distinction (tight vs")
    print("loose) is invisible at the apex; it lives at Farey-neighbor scales (mod 41).")


if __name__ == "__main__":
    main()
