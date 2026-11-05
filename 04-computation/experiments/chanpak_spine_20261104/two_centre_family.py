#!/usr/bin/env python3
"""two_centre_family.py -- the explicit two-centre family (OCF / Claim A):
T(w, s, dir) on {1..m} u {u, v}: TT_m on 1..m; u -> i iff w_i = 1 (else i -> u); v beats the suffix [s..m] and is beaten
by the prefix [1..s-1]; dir = 0: u -> v, dir = 1: v -> u.  Every odd cycle through v passes through u and has a
transitive complement, so H(T) = 1 + 2 N(w) + 2 r, with
  N(w) = sum_{i<j, w_i=1, w_j=0} c(j-i), c(1)=1, c(d)=2^(d-2)           (one-centre count, H(TT_m + u) = 1 + 2N)
  P_j  = #{increasing j-vertex paths in [s..m] ending at a vertex with w=0}  = sum_{t in [s..m], w_t=0} C(t-s, j-1)
  Q_l  = #{increasing l-vertex paths in [1..s-1] starting at a vertex with w=1} = sum_{t<s, w_t=1} C(s-1-t, l-1)
  r    = sum_{j+l odd} P_j Q_l with Q_0 = 1 if u -> v (dir 0), P_0 = 1 if v -> u (dir 1), (the empty path = the uv arc).
Part 1 checks the formula against brute-force Hamiltonian-path counts; part 2 measures coverage of the odd numbers."""
import sys, itertools
from math import comb
sys.path.insert(0, '.')
from spine_numeration import ham_paths

def c(d): return 1 if d == 1 else 2 ** (d - 2)

def N_of(w):
    m = len(w)
    return sum(c(j - i) for i in range(m) for j in range(i + 1, m) if w[i] == 1 and w[j] == 0)

def r_of(w, s, d):
    """w indexed 0..m-1 (vertex i+1); suffix = positions s-1..m-1 (vertices s..m); prefix = 0..s-2."""
    m = len(w)
    P = {}; Q = {}
    for t in range(s - 1, m):            # suffix positions
        if w[t] == 0:
            for j in range(1, t - (s - 1) + 2):
                P[j] = P.get(j, 0) + comb(t - (s - 1), j - 1)
    for t in range(0, s - 1):            # prefix positions
        if w[t] == 1:
            for l in range(1, (s - 2 - t) + 2):
                Q[l] = Q.get(l, 0) + comb(s - 2 - t, l - 1)
    if d == 0: Q[0] = Q.get(0, 0) + 1
    else: P[0] = P.get(0, 0) + 1
    return sum(pj * ql for j, pj in P.items() for l, ql in Q.items() if (j + l) % 2 == 1)

def build(w, s, d):
    m = len(w); n = m + 2; u = m; v = m + 1
    adj = [0] * n
    for i in range(m):
        for j in range(i + 1, m): adj[i] |= 1 << j
        if w[i] == 1: adj[u] |= 1 << i
        else: adj[i] |= 1 << u
        if i + 1 >= s: adj[v] |= 1 << i
        else: adj[i] |= 1 << v
    if d == 0: adj[u] |= 1 << v
    else: adj[v] |= 1 << u
    return n, adj

if __name__ == "__main__":
    bad = 0; tot = 0
    for m in range(1, 8):
        for w in itertools.product([0, 1], repeat=m):
            for s in range(1, m + 2):
                for d in (0, 1):
                    n, adj = build(w, s, d); H = ham_paths(n, adj)
                    Hf = 1 + 2 * N_of(w) + 2 * r_of(w, s, d); tot += 1
                    if H != Hf: bad += 1
    print("part 1: formula H = 1 + 2N + 2r checked on %d (w, s, dir), m <= 7: mismatches %d" % (tot, bad), flush=True)
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
    vals = set(); one_centre = set()
    for m in range(1, MMAX + 1):
        for w in itertools.product([0, 1], repeat=m):
            N = N_of(w); one_centre.add(1 + 2 * N)
            if 1 + 2 * N > LIMIT: continue
            for s in range(1, m + 2):
                for d in (0, 1):
                    H = 1 + 2 * N + 2 * r_of(w, s, d)
                    if H <= LIMIT: vals.add(H)
        miss = [x for x in range(1, LIMIT + 1, 2) if x not in vals]
        print("m <= %d: two-centre family covers %d odd values <= %d; missed (first 25): %s%s" % (m, len(vals), LIMIT, miss[:25], " ..." if len(miss) > 25 else ""), flush=True)
    miss1 = [x for x in range(1, LIMIT + 1, 2) if x not in one_centre]
    print("one-centre alone misses %d of %d odd values <= %d" % (len(miss1), (LIMIT + 1) // 2, LIMIT))
