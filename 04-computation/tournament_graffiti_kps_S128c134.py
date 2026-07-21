#!/usr/bin/env python3
"""tournament_graffiti_kps_S128c134.py -- kind-pasteur-2026-07-21-S128c134

A TxGraffiti / Written-on-the-Wall style AUTOMATED CONJECTURE GENERATOR for TOURNAMENT
invariants -- the "front end" klein-S395 named as the one WOWII ingredient the repo lacks.

WOWII-103 shape: an EASY (structural) invariant bounds a HARD (extremal) one, refuted by a
tuned witness (triangle + leaves).  Tournament analog of the witness = a 3-CYCLE atom + a
transitive skeleton (THM-1830, my unstable-non-transitive family).  This computes a battery of
iso-invariant tournament invariants over all tournaments n=3..6 (labeled; invariants are
iso-invariant so testing is complete), machine-generates candidate linear inequalities between
them, and reports the TIGHT survivors (candidate theorems) plus which obvious bounds FAIL and
with what witness.
"""
import sys
import numpy as np
from itertools import combinations, permutations
from fractions import Fraction as Fr

# ---------------- tournament invariants ----------------
def from_bits(bits, n):
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def is_transitive_subset(A, S):
    # T[S] acyclic <=> no 3-cycle within S and it is a linear order: check it's a total order
    # acyclic tournament <=> scores within S are all distinct (0..|S|-1)
    sub = A[np.ix_(S, S)]
    sc = sorted(sub.sum(axis=1).tolist())
    return sc == list(range(len(S)))


def largest_transitive(A, n):
    for k in range(n, 0, -1):
        for S in combinations(range(n), k):
            if is_transitive_subset(A, list(S)):
                return k
    return 1


def domination_number(A, n):
    # min D s.t. every v not in D has an in-arc from some d in D (D out-dominates)
    for k in range(1, n + 1):
        for D in combinations(range(n), k):
            Ds = set(D)
            ok = True
            for v in range(n):
                if v in Ds:
                    continue
                if not any(A[d, v] for d in D):
                    ok = False
                    break
            if ok:
                return k
    return n


def num_kings(A, n):
    R2 = (A + A @ A) > 0
    cnt = 0
    for v in range(n):
        if all(R2[v, w] for w in range(n) if w != v):
            cnt += 1
    return cnt


def num_scc(A, n):
    reach = np.linalg.matrix_power(A + np.eye(n, dtype=np.int64), n) > 0
    mut = reach & reach.T
    seen = [False] * n
    c = 0
    for v in range(n):
        if not seen[v]:
            c += 1
            for w in range(n):
                if mut[v, w]:
                    seen[w] = True
    return c


def ham_paths(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if c:
                for w in range(n):
                    if not (mask >> w & 1) and A[last][w]:
                        dp[mask | (1 << w)][w] += c
    return sum(dp[full][last] for last in range(n))


def arborescences(A, n):
    # spanning out-arborescences rooted at 0 via Matrix-Tree (det of grounded Laplacian)
    Aout = A.copy()
    L = np.diag(Aout.sum(axis=0)) - Aout   # in-degree Laplacian for out-trees to root 0
    M = np.delete(np.delete(L, 0, 0), 0, 1)
    return int(round(np.linalg.det(M.astype(float))))


def invariants(A, n):
    c3 = int(round(np.trace(np.linalg.matrix_power(A, 3)) / 3))
    scores = A.sum(axis=1)
    inv = {
        "n": n,
        "c3": c3,
        "H": ham_paths([list(r) for r in A], n),
        "beta": largest_transitive(A, n),         # largest transitive subtournament
        "dom": domination_number(A, n),
        "kings": num_kings(A, n),
        "scc": num_scc(A, n),
        "smax": int(scores.max()),
        "smin": int(scores.min()),
        "srange": int(scores.max() - scores.min()),
        "sumC2": int(sum(s * (s - 1) // 2 for s in scores)),   # Schur term (c3 = C(n,3)-sumC2)
        "arb0": arborescences(A, n),
    }
    return inv


# ---------------- gather over n=3..6 ----------------
print("computing tournament invariants over n=3..6 ...", flush=True)
DATA = []
for n in range(3, 7):
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        DATA.append(invariants(from_bits(bits, n), n))
print("  tournaments: %d" % len(DATA), flush=True)

KEYS = ["c3", "H", "beta", "dom", "kings", "scc", "smax", "smin", "srange", "sumC2", "arb0", "n"]

# ---------------- TxGraffiti-lite: tightest linear pair bounds ----------------
print()
print("=" * 88)
print("SURVIVING PAIR INEQUALITIES  target <= a*source + b  (a in {1/2,1,2,3}, tight)")
print("=" * 88)
found = []
for t in KEYS:
    for x in KEYS:
        if t == x:
            continue
        for a in (Fr(1, 2), Fr(1), Fr(2), Fr(3)):
            # minimal integer b so that t <= a*x + b for all
            need = max(Fr(d[t]) - a * Fr(d[x]) for d in DATA)
            b = int(np.ceil(float(need)))
            if b > 6 or b < -6:
                continue
            # tightness = number of equality cases
            eq = sum(1 for d in DATA if Fr(d[t]) == a * Fr(d[x]) + b)
            # is it interesting? require many equalities and a 'clean' form
            if eq >= 20 and (a != Fr(1) or b != 0 or t != x):
                found.append((eq, t, a, x, b))
# dedup: keep the tightest per (t,x) direction, and only non-trivial
found.sort(reverse=True)
seen_tx = set()
shown = 0
for eq, t, a, x, b in found:
    key = (t, x)
    if key in seen_tx:
        continue
    seen_tx.add(key)
    astr = "%s" % a if a != 1 else ""
    astr = (astr + "*") if astr else ""
    bstr = (" + %d" % b) if b > 0 else ((" - %d" % (-b)) if b < 0 else "")
    print("  %-6s <= %s%-6s%-6s   [equality on %d tournaments]" % (t, astr, x, bstr, eq))
    shown += 1
    if shown >= 22:
        break

print()
print("=" * 88)
print("REFUTED 'obvious' bounds, with the tournament WITNESS (WOWII-style)")
print("=" * 88)
# test a few natural bounds and, if false, exhibit the extremal witness
tests = [
    ("H <= arb0", lambda d: d["H"] <= d["arb0"]),
    ("H >= arb0", lambda d: d["H"] >= d["arb0"]),
    ("c3 <= H", lambda d: d["c3"] <= d["H"]),
    ("beta >= n - c3", lambda d: d["beta"] >= d["n"] - d["c3"]),
    ("kings >= scc", lambda d: d["kings"] >= d["scc"]),
    ("H <= 2^(n-2)*c3 + 1", lambda d: d["H"] <= (1 << (d["n"] - 2)) * d["c3"] + 1),
]
for name, pred in tests:
    bad = [d for d in DATA if not pred(d)]
    if not bad:
        print("  HOLDS on all n<=6 : %s   (candidate theorem)" % name)
    else:
        w = min(bad, key=lambda d: (d["n"], d["c3"]))
        print("  FAILS: %-24s first witness n=%d c3=%d H=%d beta=%d arb0=%d scc=%d"
              % (name, w["n"], w["c3"], w["H"], w["beta"], w["arb0"], w["scc"]))
print()
print("  The 3-cycle-atom family (THM-1830: transitive skeleton + one 3-cycle) is the")
print("  tournament analog of WOWII-103's triangle+leaves witness -- a tunable non-transitive")
print("  obstruction for breaking tight tournament-invariant bounds.")
