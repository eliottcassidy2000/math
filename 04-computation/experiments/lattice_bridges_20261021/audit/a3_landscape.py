#!/usr/bin/env python3
"""AUDIT of THM-4751 (3), HYP-9362 (a)/(b) and the proposed THM-4751 (6) (energy extremes).
Independent re-computation of the complete circulant class tables for odd n, using a3_hcount (own C code):
  * classes of connection sets S of Cay(Z_n, S) up to multipliers a in Z_n^x;
  * H by inclusion-exclusion over necklaces (mode ie) for every class, and by Held-Karp (mode dp) for n <= DPMAX;
  * class graph (one flip = reverse one difference class {k, n-k}); BFS distances to R_n and (prime n = 3 mod 4) P_n,
    compared with min over multipliers of the set distance;
  * strict local maxima (two conventions: against distinct neighbour classes; and counting self-loops);
  * additive energy E(S), directed 4- and 5-cycle counts c4, c5 by explicit DFS (not by zero-sum counting);
  * comparison with the session's data files ../../tournament_pseudoprimes_20261021/circ_*.out.
Usage: python3 a3_landscape.py NMAX [DPMAX]      Output: a3_landscape_NMAX.out (+ a3_classes_n.txt per n)
"""
import sys, os, re, glob, subprocess, json
from math import gcd
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'a3_hcount')
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 23
DPMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 21
NMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 3

def units(n): return [a for a in range(1, n) if gcd(a, n) == 1]
def canon(n, S, U):
    return min(tuple(sorted((a * s) % n for s in S)) for a in U)
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))

def run_bin(mode, n, sets):
    inp = ''.join(' '.join(map(str, (n,) + tuple(S))) + '\n' for S in sets)
    out = subprocess.run([BIN, mode], input=inp, capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        left, right = line.split(' : ')
        vals = tuple(int(x) for x in left.split())
        res[vals[1:]] = int(right.split()[0])
    return res

def energy(n, S):
    r = [0] * n
    for a in S:
        for b in S: r[(a + b) % n] += 1
    return sum(x * x for x in r)

def cycles_k(n, S, k):
    """number of directed k-cycles, by DFS over distinct vertices from 0 (cycles through 0), times n/k."""
    Sset = set(S); cnt = 0
    def dfs(v, depth, used):
        nonlocal cnt
        if depth == k:
            if (0 - v) % n in Sset: cnt += 1
            return
        for s in S:
            w = (v + s) % n
            if w not in used:
                used.add(w); dfs(w, depth + 1, used); used.discard(w)
    dfs(0, 1, {0})
    assert (cnt * n) % k == 0
    return cnt * n // k

def theirs(n):
    d = {}
    for f in glob.glob(os.path.join(HERE, '..', '..', 'tournament_pseudoprimes_20261021', 'circ_*.out')):
        for line in open(f):
            m_ = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
            if m_ and int(m_.group(1)) == n:
                d[tuple(sorted(int(x) for x in m_.group(2).split(',')))] = int(m_.group(4))
    return d

out = open(os.path.join(HERE, f'a3_landscape_{NMAX}.out' if NMIN == 3 else f'a3_landscape_{NMIN}_{NMAX}.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

for n in range(NMIN, NMAX + 1, 2):
    m = (n - 1) // 2; U = units(n)
    allsets = []
    for b in range(1 << m):
        S = tuple(sorted(((n - k) if (b >> (k - 1)) & 1 else k) for k in range(1, m + 1)))
        allsets.append(S)
    cls = {}
    for S in allsets: cls.setdefault(canon(n, S, U), []).append(S)
    reps = sorted(cls)
    R = canon(n, tuple(range(1, m + 1)), U)
    P = canon(n, tuple(sorted({x * x % n for x in range(1, n)})), U) if (isprime(n) and n % 4 == 3) else None
    H = run_bin('ie', n, reps)
    dpmatch = None
    if n <= DPMAX:
        Hd = run_bin('dp', n, reps); dpmatch = all(Hd[S] == H[S] for S in reps)
    # class graph
    nbr = {}
    for S in reps:
        Ss = set(S); lst = []
        for k in range(1, m + 1):
            T = set(Ss)
            if k in T: T.discard(k); T.add(n - k)
            else: T.discard(n - k); T.add(k)
            lst.append(canon(n, tuple(sorted(T)), U))
        nbr[S] = lst
    def bfs(src):
        dist = {src: 0}; dq = deque([src])
        while dq:
            x = dq.popleft()
            for y in nbr[x]:
                if y not in dist: dist[y] = dist[x] + 1; dq.append(y)
        return dist
    dR = bfs(R); dP = bfs(P) if P else None
    def setdist(A, B):
        A = set(A); return min(len(A - {(a * s) % n for s in B}) for a in U)
    bfs_eq_min = all(dR[S] == setdist(S, R) for S in reps) and (P is None or all(dP[S] == setdist(S, P) for S in reps))
    # local maxima
    strict_distinct, strict_loops, selfloops = [], [], 0
    for S in reps:
        others = [y for y in nbr[S] if y != S]
        if len(others) < len(nbr[S]): selfloops += 1
        if all(H[S] > H[y] for y in others): strict_distinct.append(S)
        if all(H[S] > H[y] for y in nbr[S]): strict_loops.append(S)
    gmax = max(reps, key=lambda S: H[S]); gmaxes = [S for S in reps if H[S] == H[gmax]]
    # energy / cycles
    Ev = {S: energy(n, S) for S in reps}
    c4 = {S: cycles_k(n, S, 4) for S in reps} if n >= 5 else {S: 0 for S in reps}
    c5 = {S: cycles_k(n, S, 5) for S in reps} if n >= 5 else {S: 0 for S in reps}
    rel4 = {4 * c4[S] - n * Ev[S] for S in reps}; rel5 = {2 * c5[S] + n * Ev[S] for S in reps}
    Emax = max(Ev.values()); Emin = min(Ev.values())
    argEmax = [S for S in reps if Ev[S] == Emax]; argEmin = [S for S in reps if Ev[S] == Emin]
    # comparison with session data
    th = theirs(n); cmp = None
    if th:
        thc = {}
        for S, h in th.items(): thc.setdefault(canon(n, S, U), set()).add(h)
        cmp = (len(thc), sum(1 for S in reps if S in thc and thc[S] == {H[S]}), sorted(set(reps) - set(thc)), sorted(set(thc) - set(reps)))
    def tag(S):
        t = []
        if S == R: t.append('R')
        if P is not None and S == P: t.append('P')
        return ('[' + ','.join(t) + ']') if t else str(list(S))
    say(f"n={n}: sets={len(allsets)} classes={len(reps)} self-loop classes={selfloops}  IE==DP: {dpmatch}  "
        f"BFS distance == min-over-multipliers set distance (to R{', P' if P else ''}): {bfs_eq_min}")
    if cmp: say(f"   session data: {cmp[0]} classes, H agrees on {cmp[1]}; classes missing in data: {cmp[2]}; extra in data: {cmp[3]}")
    say(f"   global max H={H[gmax]} at {[tag(S) for S in gmaxes]} dR={[dR[S] for S in gmaxes]}" + (f" dP={[dP[S] for S in gmaxes]}" if P else '')
        + f"   H(R)={H[R]}" + (f" H(P)={H[P]} d(R,P)={dR[P]}" if P else ''))
    say(f"   strict local maxima (vs distinct neighbour classes): {len(strict_distinct)}: " +
        '; '.join(f"{tag(S)} H={H[S]} dR={dR[S]}" + (f" dP={dP[S]}" if P else '') for S in sorted(strict_distinct, key=lambda S: -H[S])))
    if strict_loops != strict_distinct:
        say(f"   (counting self-loops as neighbours: {len(strict_loops)}: {[tag(S) for S in strict_loops]})")
    if P is not None:
        nb = sorted(set(nbr[P]) - {P}); say(f"   P neighbour classes: {len(nb)} -> {[tag(S) for S in nb]}, H = {[H[S] for S in nb]}, margin = {[f'{H[S]/H[P]:.6f}' for S in nb]}")
    if n >= 5:
        say(f"   energy: E in [{Emin},{Emax}]; argmax {[tag(S) for S in argEmax]}; argmin {[tag(S) for S in argEmin]}; "
            f"4c4-nE constant: {len(rel4)==1} {sorted(rel4)[:3]}; 2c5+nE constant: {len(rel5)==1} {sorted(rel5)[:3]}")
        hm = [S for S in gmaxes]
        say(f"   H-max is E-max: {all(S in argEmax for S in hm)}; H-max is E-min: {all(S in argEmin for S in hm)}; "
            f"E-min classes strict local max of H: {[S in strict_distinct for S in argEmin]}; H-rank of E-min: "
            f"{[1 + sum(1 for T in reps if H[T] > H[S]) for S in argEmin]}")
    with open(os.path.join(HERE, f'a3_classes_{n}.txt'), 'w') as fh:
        for S in sorted(reps, key=lambda S: -H[S]):
            fh.write(json.dumps({'S': list(S), 'H': H[S], 'dR': dR[S], 'dP': (dP[S] if P else None), 'E': Ev[S], 'c4': c4[S], 'c5': c5[S],
                                 'nbrs': [list(y) for y in nbr[S]], 'orbit': len(cls[S])}) + '\n')
out.close()
