#!/usr/bin/env python3
"""additive.py: additive gluings on a shared arc.  For tournaments A, B of order <= 4 (all classes), every arc u -> v of A and
every arc u' -> v' of B (identified u = u', v = v'), and every orientation of the cross arcs between A - {u,v} and
B - {u,v}: U on (|A| + |B| - 2) vertices.  Record the cases with H(U) = H(A) + H(B) - 1 (and -3, +1, +3), whether
the path-splice bound sum_sigma h_A(sigma) h_B(sigma) is attained with equality, and whether every Hamiltonian path of
U restricts to Hamiltonian paths of A and B."""
import itertools, subprocess
from collections import Counter
def parse(s, n):
    adj = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj
def hampaths(adj, S):
    S = list(S); res = []
    def rec(path, used):
        if len(path) == len(S): res.append(tuple(path)); return
        for x in S:
            if x not in used and (not path or adj[path[-1]][x]):
                used.add(x); path.append(x); rec(path, used); path.pop(); used.discard(x)
    rec([], set()); return res
classes = {}
for n in (3, 4):
    classes[n] = [parse(s, n) for s in subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()]
results = Counter(); examples = {}; total = 0; eq_bound_all = True
for nA, nB in [(3,3),(3,4),(4,3),(4,4)]:
    for A in classes[nA]:
        for B in classes[nB]:
            HA = len(hampaths(A, range(nA))); HB = len(hampaths(B, range(nB)))
            for (u, v) in [(i, j) for i in range(nA) for j in range(nA) if A[i][j]]:
                for (u2, v2) in [(i, j) for i in range(nB) for j in range(nB) if B[i][j]]:
                    # build U: vertices: A's vertices 0..nA-1 (u, v among them); B's other vertices mapped to nA.. ; u2->u, v2->v
                    restB = [x for x in range(nB) if x not in (u2, v2)]; mapB = {u2: u, v2: v}
                    for i, x in enumerate(restB): mapB[x] = nA + i
                    n = nA + len(restB)
                    base = [[0]*n for _ in range(n)]
                    for i in range(nA):
                        for j in range(nA): base[i][j] = A[i][j]
                    for x in range(nB):
                        for y in range(nB):
                            if B[x][y] and (x in (u2, v2) and y in (u2, v2)): continue
                            if B[x][y]: base[mapB[x]][mapB[y]] = 1
                    crossA = [i for i in range(nA) if i not in (u, v)]; crossB = [nA + i for i in range(len(restB))]
                    pairs = [(a, b) for a in crossA for b in crossB]
                    for bits in range(1 << len(pairs)):
                        adj = [row[:] for row in base]
                        for t, (a, b) in enumerate(pairs):
                            if (bits >> t) & 1: adj[a][b] = 1
                            else: adj[b][a] = 1
                        PU = hampaths(adj, range(n)); HU = len(PU)
                        total += 1
                        diff = HU - (HA + HB)
                        if diff in (-3, -1, 1, 3):
                            SA = set(range(nA)); SB = set([u, v] + crossB)
                            PA = set(hampaths(adj, SA)); PB = set(hampaths(adj, SB))
                            restrict_all = all(tuple(x for x in p if x in SA) in PA and tuple(x for x in p if x in SB) in PB for p in PU)
                            def ords(P, W):
                                d = Counter()
                                for p in P: d[tuple(x for x in p if x in W)] += 1
                                return d
                            W = {u, v}; dA = ords(PA, W); dB = ords(PB, W)
                            bound = sum(dA[s]*dB.get(s, 0) for s in dA)
                            key = (diff, HA, HB); results[key] += 1
                            if key not in examples: examples[key] = (HU, bound, restrict_all, nA, nB)
                            if diff == -1 and HU != bound: eq_bound_all = False
print(f"instances {total}; equality cases by (H(U) - H(A) - H(B), H(A), H(B)):")
for key in sorted(results): print("  ", key, "count", results[key], "example (H(U), bound, all paths restrict, |A|, |B|):", examples[key])
print("every H(U) = H(A)+H(B)-1 case attains the path-splice bound:", eq_bound_all)
