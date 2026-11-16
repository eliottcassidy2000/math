#!/usr/bin/env python3
"""shuffle_check.py n: for every tournament of order n (gentourng), every pair of vertex subsets (S_A, S_B) with
|S_A|, |S_B| >= 2 and |S_A ∩ S_B| <= 2, check the path-splice bound
    H(U) >= sum over linear orders sigma of W = S_A ∩ S_B of h_A(sigma) h_B(sigma),
where U = T[S_A ∪ S_B], A = T[S_A], B = T[S_B], h_A(sigma) = number of Hamiltonian paths of A inducing sigma on W,
and also count the shuffles: Hamiltonian paths Z of U whose restrictions Z|A and Z|B are Hamiltonian paths of A and B.
Reports the number of instances, the number of violations, the number of equalities (|W| <= 1) and the min excess."""
import sys, subprocess, itertools
n = int(sys.argv[1])
out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
def parse(line):
    # gentourng default output: upper triangle bits, row by row: bit for (i,j), i<j: 1 means i->j
    bits = line.strip()
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            b = bits[k]; k += 1
            if b == '1': adj[i][j] = 1
            else: adj[j][i] = 1
    return adj
def ham_paths(adj, S):
    S = list(S)
    res = []
    def rec(path, used):
        if len(path) == len(S):
            res.append(tuple(path)); return
        for x in S:
            if x not in used and (not path or adj[path[-1]][x]):
                used.add(x); path.append(x); rec(path, used); path.pop(); used.discard(x)
    rec([], set())
    return res
inst = viol = eq1 = 0; minexc = None; shuffles_ok = True
for line in out:
    adj = parse(line)
    V = range(n)
    cache = {}
    def HP(S):
        key = frozenset(S)
        if key not in cache: cache[key] = ham_paths(adj, key)
        return cache[key]
    subsets = [frozenset(c) for r in range(2, n+1) for c in itertools.combinations(V, r)]
    for SA in subsets:
        for SB in subsets:
            W = SA & SB
            if len(W) > 2 or SA == SB: continue
            if SB < SA or SA < SB: continue
            U = SA | SB
            PA = HP(SA); PB = HP(SB); PU = HP(U)
            # h_A(sigma): order induced on W
            def ords(P):
                d = {}
                for p in P:
                    s = tuple(x for x in p if x in W); d[s] = d.get(s, 0) + 1
                return d
            dA, dB = ords(PA), ords(PB)
            bound = sum(dA[s]*dB.get(s, 0) for s in dA)
            HU = len(PU)
            # shuffles
            setA, setB = set(PA), set(PB)
            sh = sum(1 for z in PU if tuple(x for x in z if x in SA) in setA and tuple(x for x in z if x in SB) in setB)
            inst += 1
            if HU < bound: viol += 1; print("VIOLATION", line, sorted(SA), sorted(SB), HU, bound)
            if sh < bound: shuffles_ok = False; print("SHUFFLE COUNT BELOW BOUND", line, sorted(SA), sorted(SB), sh, bound)
            if len(W) <= 1:
                if HU == len(PA)*len(PB): eq1 += 1
                exc = HU - len(PA)*len(PB)
                if minexc is None or exc < minexc: minexc = exc
print(f"n={n}: instances {inst}, violations {viol}, |W|<=1 equalities {eq1}, min excess (|W|<=1) {minexc}, shuffle count >= bound everywhere: {shuffles_ok}")
