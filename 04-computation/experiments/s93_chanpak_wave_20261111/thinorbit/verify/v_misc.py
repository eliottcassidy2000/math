#!/usr/bin/env python3
"""Verifier (lane thinorbit, S93): three small independent checks.
(P) Lemma P with the path counts of the proof, by plain DFS, for EVERY marked block of order <= 5 (26,705 blocks),
    against my formula-based matrices (v_gen.all_blocks).
(D) Lemma D on all labelled tournaments with n <= 6 (odd cycles as directed cycles' vertex sets), and the number of
    tournaments at n = 6 with a component of Omega that is a chain of exactly three classes.
(S) two-apex letters: characteristic polynomial of words must be t^8 (t-1)^16 (t-2^n)^10 (t-3^n)^2 if the letters are
    simultaneously triangular with the common diagonal {0^8, 1^16, 2^10, 3^2}; commutators must be nilpotent.
"""
import sys
import itertools
import random
import numpy as np

sys.path.insert(0, "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit/verify")
from v_gen import all_blocks, tournament_from_bits


def check_P(kmax=5):
    bad = 0
    tot = 0
    for k in range(1, kmax + 1):
        arr = all_blocks(k)
        npairs = k * (k - 1) // 2
        for bits in range(1 << npairs):
            adj = tournament_from_bits(k, bits)
            # all directed paths of B (as tuples), by DFS
            paths = []

            def rec(p):
                paths.append(tuple(p))
                for t in range(k):
                    if t not in p and adj[p[-1]][t]:
                        p.append(t)
                        rec(p)
                        p.pop()

            for s in range(k):
                rec([s])
            for w in range(k):
                for u2 in range(k):
                    e = sum(1 for p in paths if p[0] == u2 and len(p) % 2 == 0)
                    o = sum(1 for p in paths if p[0] == u2 and len(p) % 2 == 1)
                    A = [0, 0]
                    for q in paths:
                        if q[0] == u2 and q[-1] == w:
                            sq = set(q)
                            A[len(q) % 2] += 1          # R empty
                            for r in paths:
                                if not (sq & set(r)):
                                    A[(len(q) + len(r)) % 2] += 1
                    N = arr[(bits * k + w) * k + u2]
                    pred = (1, 0, 0, e % 2, A[0] % 2, A[1] % 2, o % 2, A[1] % 2, A[0] % 2)
                    tot += 1
                    if tuple(int(x) % 2 for x in N) != pred:
                        bad += 1
    print(f"(P) Lemma P formula (e, o, A0, A1 by plain DFS) vs matrices mod 2: {tot} marked blocks of order <= {kmax}, mismatches {bad}")


def check_D(nmax=6):
    for n in range(3, nmax + 1):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        subsets = [c for r in range(3, n + 1, 2) for c in itertools.combinations(range(n), r)]
        fails_comp = fails_dom = fails_long = diam3 = chain3 = 0
        for bits in range(1 << len(pairs)):
            adj = [[0] * n for _ in range(n)]
            for t, (i, j) in enumerate(pairs):
                if (bits >> t) & 1:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
            # strong components by reachability
            reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
            for m in range(n):
                for i in range(n):
                    if reach[i][m]:
                        for j in range(n):
                            if reach[m][j]:
                                reach[i][j] = 1
            comp = [frozenset(j for j in range(n) if reach[i][j] and reach[j][i]) for i in range(n)]
            # odd cycle vertex sets: S carries a Hamiltonian directed cycle iff T[S] is strong (Camion) and |S| >= 3
            cyc = []
            for S in subsets:
                if all(comp[v] >= frozenset(S) for v in S):
                    # T[S] strong?  (S inside one strong component is not enough) -> test directly
                    sub = list(S)
                    r2 = {a: {a} for a in sub}
                    changed = True
                    while changed:
                        changed = False
                        for a in sub:
                            for b in list(r2[a]):
                                for c in sub:
                                    if adj[b][c] and c not in r2[a]:
                                        r2[a].add(c)
                                        changed = True
                    if all(len(r2[a]) == len(sub) for a in sub):
                        cyc.append(frozenset(S))
            if not cyc:
                continue
            # components of Omega
            unseen = set(range(len(cyc)))
            comps = []
            while unseen:
                s = unseen.pop()
                cur = {s}
                stack = [s]
                while stack:
                    x = stack.pop()
                    for y in list(unseen):
                        if cyc[x] & cyc[y]:
                            unseen.discard(y)
                            cur.add(y)
                            stack.append(y)
                comps.append(sorted(cur))
            strong_big = set(c for c in comp if len(c) >= 3)
            if len(comps) != len(strong_big):
                fails_comp += 1
            has_chain = False
            for cc in comps:
                W = comp[next(iter(cyc[cc[0]]))]
                if any(not (cyc[x] <= W) for x in cc):
                    fails_comp += 1
                dom = [x for x in cc if all(cyc[x] & cyc[y] for y in cc)]
                if not dom:
                    fails_dom += 1
                if not any(len(cyc[x]) >= len(W) - 1 for x in dom):
                    fails_long += 1
                rest = [x for x in cc if x not in dom]
                # diameter >= 3 ?  (two non-adjacent cycles without a common neighbour)
                for a in cc:
                    for b in cc:
                        if a < b and not (cyc[a] & cyc[b]) and not any(cyc[c] & cyc[a] and cyc[c] & cyc[b] for c in cc):
                            diam3 += 1
                if rest:
                    # rest = disjoint union of exactly two cliques with no edges between them?
                    a0 = rest[0]
                    c1 = [x for x in rest if x == a0 or (cyc[x] & cyc[a0])]
                    c2 = [x for x in rest if x not in c1]
                    ok = bool(c2)
                    ok = ok and all(cyc[x] & cyc[y] for x in c1 for y in c1) and all(cyc[x] & cyc[y] for x in c2 for y in c2)
                    ok = ok and not any(cyc[x] & cyc[y] for x in c1 for y in c2)
                    if ok:
                        has_chain = True
            if has_chain:
                chain3 += 1
        print(f"(D) n = {n}: {1 << len(pairs)} labelled tournaments; component/strong-component mismatches {fails_comp}; "
              f"components without a dominating cycle {fails_dom}; without a dominating cycle on >= |W|-1 vertices {fails_long}; "
              f"pairs at distance >= 3: {diam3}; tournaments with a three-class chain component: {chain3}", flush=True)


def check_S():
    sys.path.insert(0, "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit")
    import sympy as sp
    from two_apex_algebra import step_matrix
    letters = [(0, 0), (1, 0), (0, 1), (1, 1)]
    Ls = [sp.Matrix(36, 36, lambda i, j, L=step_matrix(a, b): int(L[i, j])) for (a, b) in letters]
    t = sp.symbols("t")
    rng = random.Random(7)
    words = [(i,) for i in range(4)] + [(i, j) for i in range(4) for j in range(4)]
    words += [tuple(rng.randrange(4) for _ in range(n)) for n in (3, 3, 3, 4, 4, 4, 5, 5, 6, 6)]
    bad = 0
    for wd in words:
        W = sp.eye(36)
        for l in wd:
            W = Ls[l] * W
        n = len(wd)
        cp = W.charpoly(t).as_expr()
        target = sp.expand(t ** 8 * (t - 1) ** 16 * (t - 2 ** n) ** 10 * (t - 3 ** n) ** 2)
        if sp.expand(cp - target) != 0:
            bad += 1
    nil = 0
    for i in range(4):
        for j in range(i + 1, 4):
            C = Ls[i] * Ls[j] - Ls[j] * Ls[i]
            if (C ** 36).is_zero_matrix:
                nil += 1
    print(f"(S) two-apex letters: charpoly of {len(words)} words = t^8 (t-1)^16 (t-2^n)^10 (t-3^n)^2: mismatches {bad}; "
          f"nilpotent commutators [L_i, L_j]: {nil} / 6")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "PDS"
    if "P" in which:
        check_P(5)
    if "S" in which:
        check_S()
    if "D" in which:
        check_D(6)
