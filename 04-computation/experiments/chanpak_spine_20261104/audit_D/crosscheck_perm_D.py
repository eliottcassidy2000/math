#!/usr/bin/env python3
"""crosscheck_perm_D.py -- audit D, independent brute force by permutations (no DP, no session imports).
Reads gentourng default output (upper triangle row by row, '1' at (i,j) = i->j) for n = 3..6 from stdin? No: calls gentourng itself.
Checks per class X and vertex a: e_a (paths ending at a), s_a, h; 2e_a <= h-1 for non-sinks, 2s_a <= h-1 for non-sources;
sole-exit law H(X+v) = 2h - e_a + sigma with sigma = [P A_a] - h computed from the definition by brute force over subsets
and permutations; max e_a/h; equality cases. Orders 3..6 complete (extension order <= 7); order 7 classes: bound + ratio only.
Also the labelled count at n <= 5 (all 2^C(n,2)) for the session's labelled numbers restricted to n <= 5."""
import itertools, subprocess, sys
from fractions import Fraction

def parse(line, n):
    adj = [[False]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            b = line[k] == '1'; k += 1
            adj[i][j] = b; adj[j][i] = not b
    return adj

def ham_paths(adj, verts):
    """list of all Hamiltonian paths (tuples) of the subtournament on verts, by permutations"""
    out = []
    for p in itertools.permutations(verts):
        if all(adj[p[i]][p[i+1]] for i in range(len(p)-1)):
            out.append(p)
    return out

def PA(adj, n, a):
    """[P A_a] = sum_{S subset V-a} h(X[S]) s_a(X[V-S]), h(empty) = 1, by brute force"""
    rest = [v for v in range(n) if v != a]
    tot = 0
    for r in range(len(rest)+1):
        for S in itertools.combinations(rest, r):
            hS = 1 if r == 0 else len(ham_paths(adj, S))
            comp = [v for v in range(n) if v not in S]
            sa = sum(1 for p in ham_paths(adj, comp) if p[0] == a)
            tot += hS * sa
    return tot

def run_classes(n, do_se):
    lines = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    best = Fraction(0); best_at = None; eq = []; viol = 0; viold = 0; se = [0, 0, 0, 0, 0]  # cases, lawfail, sigma<0, sigma=0, floorfail
    for line in lines:
        adj = parse(line, n)
        paths = ham_paths(adj, range(n)); h = len(paths)
        assert h % 2 == 1
        for a in range(n):
            ea = sum(1 for p in paths if p[-1] == a); sa = sum(1 for p in paths if p[0] == a)
            sink = not any(adj[a]); source = all(adj[a][j] for j in range(n) if j != a)
            if not sink:
                if 2*ea > h-1: viol += 1; print("VIOL", n, line, a)
                if Fraction(ea, h) > best: best = Fraction(ea, h); best_at = (line, a, h, ea)
                if 2*ea == h-1: eq.append((line, a, h))
                if do_se:
                    adj2 = [row + [False] for row in adj] + [[False]*(n+1)]
                    adj2[n][a] = True
                    for x in range(n):
                        if x != a: adj2[x][n] = True
                    H2 = len(ham_paths(adj2, range(n+1)))
                    sigma = PA(adj, n, a) - h
                    se[0] += 1
                    if H2 != 2*h - ea + sigma: se[1] += 1; print("LAWFAIL", n, line, a, h, ea, sigma, H2)
                    if sigma < 0: se[2] += 1
                    if sigma == 0: se[3] += 1
                    if 2*H2 < 3*h + 1: se[4] += 1; print("FLOORFAIL", n, line, a)
            if not source and 2*sa > h-1: viold += 1; print("VIOLDUAL", n, line, a)
    print(f"classes n={n}: {len(lines)} tournaments; violations={viol} dual={viold}; max e_a/h={best} at {best_at}; equality cases={len(eq)}: {eq}")
    if do_se: print(f"   sole-exit: cases={se[0]} lawfail={se[1]} sigma<0={se[2]} sigma=0:{se[3]} floorfail={se[4]}")

def run_labelled(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    cases = 0; sig0 = 0; maxsig = 0; viol = 0
    for bits in range(1 << len(pairs)):
        adj = [[False]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: adj[i][j] = True
            else: adj[j][i] = True
        paths = ham_paths(adj, range(n)); h = len(paths)
        for a in range(n):
            if not any(adj[a]): continue
            ea = sum(1 for p in paths if p[-1] == a)
            if 2*ea > h-1: viol += 1
            adj2 = [row + [False] for row in adj] + [[False]*(n+1)]
            adj2[n][a] = True
            for x in range(n):
                if x != a: adj2[x][n] = True
            H2 = len(ham_paths(adj2, range(n+1)))
            sigma = H2 - (2*h - ea)
            cases += 1; sig0 += (sigma == 0); maxsig = max(maxsig, sigma)
            assert sigma >= 0
    print(f"labelled n={n}: cases={cases} sigma=0:{sig0} maxsigma={maxsig} violations={viol}")

if __name__ == '__main__':
    for n in (3, 4, 5, 6): run_classes(n, True)
    run_classes(7, False)
    for n in (3, 4, 5): run_labelled(n)
