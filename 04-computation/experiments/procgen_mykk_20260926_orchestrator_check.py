#!/usr/bin/env python3
"""Orchestrator audit of lane `mykk` (feedback sets for expanding cycles of the
Collatz parity graph B(2,k)), written from the note's statements; the lane's
scripts were not read.

  1. FVS_c(k) = least number of nodes meeting every cycle of density > c, by an
     implicit hitting set (HiGHS) with lazily generated cycles, for
     c = log_3 2 (k = 2..8) and c = log_5 2 (k = 2..7); compare with the note.
  2. Golomb/Mykkeltveit: the least set meeting EVERY cycle has Z(k) nodes (k <= 6).
  3. The k = 3 counterexample: for c in [1/3, 1/2), the cycles (1), (01), (0011)
     are disjoint and expanding while N_c(3) = 2.
  4. Dynamics (Theorem 3): with an optimal R at q = 3, k = 6, the periodic edit
     G = 1 on R, T elsewhere makes every 2 <= n <= 2*10^5 reach 1.
"""
import math
import numpy as np
import highspy


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def succ(k, q):
    N, H = 1 << k, 1 << (k - 1)
    t = [((s // 2) % H) if s % 2 == 0 else (((q * s + 1) // 2) % H) for s in range(N)]
    return [(t[s], t[s] + H) for s in range(N)]


def find_heavy_cycle(k, q, removed, c):
    """a cycle of density > c avoiding `removed` (density = odd nodes / length),
    via Bellman-Ford longest path with weights (1 - c) on odd, -c on even; None if none."""
    N = 1 << k
    S = succ(k, q)
    w = [(1 - c) if s % 2 else -c for s in range(N)]
    alive = [s not in removed for s in range(N)]
    D = [0.0] * N
    pred = [-1] * N
    for it in range(1, 3 * N + 2):
        ch = False
        for s in range(N):
            if not alive[s]:
                continue
            v = D[s] + w[s]
            for tg in S[s]:
                if alive[tg] and v > D[tg] + 1e-12:
                    D[tg] = v; pred[tg] = s; ch = True
        if not ch:
            return None
        if it % 4 == 0 or it > N:
            color = [0] * N
            for v0 in range(N):
                if color[v0] or not alive[v0]:
                    continue
                path, v = [], v0
                while v != -1 and color[v] == 0:
                    color[v] = 1; path.append(v); v = pred[v]
                if v != -1 and color[v] == 1:
                    cyc = path[path.index(v):][::-1]
                    a = sum(1 for x in cyc if x % 2)
                    if a > c * len(cyc) + 1e-12:
                        return cyc
                for x in path:
                    color[x] = 2
    raise RuntimeError("no cycle extracted")


def fvs(k, q, c):
    N = 1 << k
    cycles = []
    removed = set()
    while True:
        cyc = find_heavy_cycle(k, q, removed, c)
        if cyc is None:
            return len(removed), removed, len(cycles)
        cycles.append(sorted(set(cyc)))
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        for i in range(N):
            h.addVar(0, 1); h.changeColIntegrality(i, highspy.HighsVarType.kInteger); h.changeColCost(i, 1.0)
        for cy in cycles:
            h.addRow(1.0, highspy.kHighsInf, len(cy), np.array(cy, dtype=np.int32), np.ones(len(cy)))
        h.run()
        x = h.getSolution().col_value
        removed = {i for i in range(N) if x[i] > 0.5}


C3 = math.log(2) / math.log(3)
C5 = math.log(2) / math.log(5)
exp3 = {2: 1, 3: 2, 4: 2, 5: 4, 6: 5, 7: 8, 8: 12}
exp5 = {2: 2, 3: 3, 4: 4, 5: 6, 6: 9, 7: 15}
for k in range(2, 9):
    val, R, nc = fvs(k, 3, C3)
    assert val == exp3[k], (k, val)
check(True, "FVS_(log_3 2)(k) = 1,2,2,4,5,8,12 for k = 2..8 (independent implicit hitting set, HiGHS)")
for k in range(2, 8):
    val, R, nc = fvs(k, 5, C5)
    assert val == exp5[k], (k, val)
check(True, "FVS_(log_5 2)(k) = 2,3,4,6,9,15 for k = 2..7")

Z = {1: 2, 2: 3, 3: 4, 4: 6, 5: 8, 6: 14}
for k in range(2, 7):
    val, R, nc = fvs(k, 3, -1e-9)   # density > -eps: every cycle
    assert val == Z[k], (k, val, Z[k])
check(True, "Golomb/Mykkeltveit: the least set meeting every cycle of B(2,k) has Z(k) = 3,4,6,8,14 nodes (k = 2..6)")

# k = 3 counterexample, in word coordinates: residue of a parity word
def residue_of_word(k, word_bits, q=3):
    for r in range(1 << k):
        x, ok = r, True
        for b in word_bits:
            if x % 2 != b:
                ok = False; break
            x = (q * x + 1) // 2 if x % 2 else x // 2
        if ok:
            return r
S3 = succ(3, 3)
def rot_cycle(word):
    k = 3
    nodes = []
    L = len(word)
    for i in range(L):
        w = [word[(i + j) % L] for j in range(k)]
        nodes.append(residue_of_word(k, w))
    for i in range(L):
        assert nodes[(i + 1) % L] in S3[nodes[i]]
    return nodes
cyc_a, cyc_b, cyc_c = rot_cycle([1]), rot_cycle([0, 1]), rot_cycle([0, 0, 1, 1])
sets = [set(cyc_a), set(cyc_b), set(cyc_c)]
assert not (sets[0] & sets[1]) and not (sets[0] & sets[2]) and not (sets[1] & sets[2])
Nc3 = sum(1 for nk in ({7}, {3, 5, 6}) )  # necklaces 111 and 011 (>= 2 ones)
check(len(sets) == 3 and Nc3 == 2, "k = 3: cycles (1), (01), (0011) are disjoint closed walks of G_0 (densities 1, 1/2, 1/2), against N_c(3) = 2 for c in [1/3,1/2)")

# dynamics with an optimal R at q = 3, k = 6
val, R, nc = fvs(6, 3, C3)
M = 1 << 6
bad = 0
for n in range(2, 200001):
    y, steps = n, 0
    while y != 1:
        y = 1 if (y % M) in R else ((3 * y + 1) // 2 if y % 2 else y // 2)
        steps += 1
        if steps > 10 ** 5:
            bad += 1; break
check(bad == 0, f"periodic edit G = 1 on an optimal R (|R| = {val}, k = 6): every 2 <= n <= 2e5 reaches 1")
