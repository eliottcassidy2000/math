#!/usr/bin/env python3
"""self_lines_twisted_burnside_20261020.py -- opus-2026-10-20-S46

Self-lines of the tiling line metagraph WITHOUT enumerating the 2^m tilings.

Frame (THM-643/849): base path 0 -> 1 -> ... -> n-1; tiles = non-path pairs; kappa = flip all
tiles; sigma = staircase reflection (converse along the reversed path).  A tiling t is an
ordinary self-line endpoint iff T_{kappa t} is isomorphic to T_t (THM-849's set Q_n), and a
converse-merged self-loop endpoint iff T_{kappa t} is isomorphic to T_t or to T_t^op.  Since
(kappa T)^op = T with ONLY the path arcs reversed, the second case is "the Hamiltonian path is
reversal-neutral" (Q'_n).

For a fixed permutation g the condition  g . T_t = T_{kappa t}  is a linear system over F_2 in
the tile bits (one equation per pair, following g's cycles on unordered pairs), and likewise
g . T_t = (T_{kappa t})^op and g . T_t = T_t^op (the classical self-converse twist).  Each
system is solved by a parity union-find: it is consistent or not, and its solution set is an
affine space of dimension = number of pair-components not pinned by a path pair.  Summing
2^dim over g is a twisted Burnside sum; enumerating the solutions over all g gives the sets
Q_n, Q'_n exactly, and the number of g producing a given t is |Aut(T_t)| (R_kappa(t) is a coset
of Aut, THM-849 (11)).  Cost: n! systems of C(n,2) pairs -- no 2^m enumeration.

Output per n: |Q_n| (blue/black by sigma-fixedness), the |Aut| histogram of the endpoints,
self-lines = |Q|/2 by colour, the Klein-orbit counts, P_kappa = sum |Aut| (the pure Burnside
sum), the path-reversal set Q'_n, the converse-merged self-loop count, and the classical
SC(n) by the all-pairs twist as a control.

Usage: python3 self_lines_twisted_burnside_20261020.py n [n2 ...]
"""
import sys, itertools, time
from collections import Counter

def run(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    P = len(pairs)
    ispath = [j == i + 1 for (i, j) in pairs]
    tiles = [k for k in range(P) if not ispath[k]]
    tilepos = {k: a for a, k in enumerate(tiles)}
    m = len(tiles)
    ONE = P  # the constant-1 node
    # sigma on tiles: (i,j) -> (n-1-j, n-1-i)
    def sigma_bits(bits):
        out = 0
        for a, k in enumerate(tiles):
            if (bits >> a) & 1:
                i, j = pairs[k]; out |= 1 << tilepos[pid[(n - 1 - j, n - 1 - i)]]
        return out

    def solve(g, twist, pin=True):
        """twist(e, a, b) -> parity constant c for the equation y_e + y_{g^-1 e} = c.
        pin: path pairs are pinned to the base-path orientation (tilings); pin=False ranges over all
        labelled tournaments (the classical self-converse Burnside control).
        Returns None if inconsistent, else (components, root of the constant node ONE)."""
        parent = list(range(P + 1)); par = [0] * (P + 1)  # parity relative to parent
        def find(x):
            if parent[x] == x: return x, 0
            r, p = find(parent[x]); parent[x] = r; par[x] ^= p; return r, par[x]
        def union(x, y, c):  # enforce y_x + y_y = c
            rx, px = find(x); ry, py = find(y)
            if rx == ry: return (px ^ py) == c
            parent[rx] = ry; par[rx] = px ^ py ^ c; return True
        if pin:
            for k in range(P):
                if ispath[k]:
                    if not union(k, ONE, 0): return None   # y_e = 1 for path pairs (value of ONE is 1)
        ginv = [0] * n
        for v in range(n): ginv[g[v]] = v
        for k, (u, v) in enumerate(pairs):
            a, b = ginv[u], ginv[v]
            e2 = pid[(a, b) if a < b else (b, a)]
            c = twist(k, a, b)
            if not union(k, e2, c): return None
        comps = {}
        for k in range(P + 1):
            r, p = find(k); comps.setdefault(r, []).append((k, p))
        rootONE, pONE = find(ONE)
        return comps, rootONE, pONE

    def enumerate_solutions(res):
        comps, rootONE, pONE = res
        free = [r for r in comps if r != rootONE]
        sols = []
        for choice in itertools.product((0, 1), repeat=len(free)):
            val = {rootONE: 1 ^ pONE}      # the constant node ONE must take the value 1
            for r, c in zip(free, choice): val[r] = c
            bits = 0
            for r, members in comps.items():
                base = val[r]
                for k, p in members:
                    if k == ONE: continue
                    y = base ^ p
                    if not ispath[k] and y: bits |= 1 << tilepos[k]
            sols.append(bits)
        return sols

    # twists: ordinary kappa:  y_e + y_{g^-1 e} = [e tile] + [a > b]
    #         path reversal :  y_e + y_{g^-1 e} = [e path] + [a > b]
    #         self-converse :  y_e + y_{g^-1 e} = 1 + [a > b]
    tw_kappa = lambda k, a, b: ((0 if ispath[k] else 1) + (1 if a > b else 0)) & 1
    tw_rev = lambda k, a, b: ((1 if ispath[k] else 0) + (1 if a > b else 0)) & 1
    tw_sc = lambda k, a, b: (1 + (1 if a > b else 0)) & 1
    t0 = time.time()
    Q = Counter(); Qrev = Counter(); Pk = 0; Prev = 0; sc_sum = 0
    for g in itertools.permutations(range(n)):
        r = solve(g, tw_kappa)
        if r is not None:
            sols = enumerate_solutions(r); Pk += len(sols)
            for s in sols: Q[s] += 1
        r = solve(g, tw_rev)
        if r is not None:
            sols = enumerate_solutions(r); Prev += len(sols)
            for s in sols: Qrev[s] += 1
        r = solve(g, tw_sc, pin=False)
        if r is not None:
            comps, rootONE, pONE = r; sc_sum += 2 ** (len(comps) - 1)   # ONE is isolated when unpinned
    fact = 1
    for i in range(2, n + 1): fact *= i
    SC = sc_sum // fact
    assert sc_sum % fact == 0
    # colour and checks
    blueQ = sum(1 for t in Q if sigma_bits(t) == t); blackQ = len(Q) - blueQ
    autHist = Counter(Q.values()); autHistB = Counter(Q[t] for t in Q if sigma_bits(t) == t)
    # Klein orbits on black endpoints: {t, kappa t, sigma t, kappa sigma t} all in Q
    ALL = (1 << m) - 1
    kleinOK = all(((t ^ ALL) in Q) and (sigma_bits(t) in Q) for t in Q)
    union = set(Q) | set(Qrev)
    blueU = sum(1 for t in union if sigma_bits(t) == t)
    both = set(Q) & set(Qrev)
    print(f"n={n}: perms {fact}, pairs {P}, tiles {m}; time {time.time()-t0:.1f}s")
    print(f"  SC(n) by all-pairs twisted Burnside = {SC}")
    print(f"  Q_n (kappa t ~ t): |Q| = {len(Q)} (blue {blueQ}, black {blackQ}); self-lines blue {blueQ//2}, black {blackQ//2}; "
          f"P_kappa = sum|Aut| = {Pk}; |Aut| histogram all {dict(sorted(autHist.items()))}, blue {dict(sorted(autHistB.items()))}; "
          f"Klein closure {kleinOK}; |Q^K| mod 4 = {blackQ % 4}")
    print(f"  Q'_n (kappa t ~ t^op, i.e. path reversal neutral): |Q'| = {len(Qrev)}, P_rev = {Prev}, |Aut| histogram {dict(sorted(Counter(Qrev.values()).items()))}; "
          f"|Q cap Q'| = {len(both)}")
    print(f"  converse-merged self-loop endpoints |Q cup Q'| = {len(union)} (blue {blueU}, black {len(union)-blueU}); merged self-loops blue {blueU//2}, black {(len(union)-blueU)//2}")
    # score leg condition s_0 = s_{n-1} + 1 on Q' (proved) and e_0 = e_{n-1} on Q (THM-793)
    def scores(bits):
        s = [0] * n
        for i in range(n - 1): s[i] += 1
        for a, k in enumerate(tiles):
            i, j = pairs[k]
            if (bits >> a) & 1: s[i] += 1
            else: s[j] += 1
        return s
    bad_rev = sum(1 for t in Qrev if scores(t)[0] != scores(t)[n - 1] + 1)
    bad_k1 = sum(1 for t in Q if scores(t)[0] != scores(t)[n - 1] + 1)
    bad_k0 = sum(1 for t in Q if scores(t)[0] != scores(t)[n - 1])
    print(f"  leg conditions: Q' with s_0 != s_(n-1) + 1: {bad_rev}; Q with s_0 != s_(n-1) + 1: {bad_k1}; Q with s_0 != s_(n-1): {bad_k0}")
    sys.stdout.flush()
    return dict(n=n, SC=SC, Q=len(Q), blueQ=blueQ, blackQ=blackQ, Pk=Pk, Qrev=len(Qrev), union=len(union), blueU=blueU)

if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [4, 5, 6, 7]
    for n in ns: run(n)
