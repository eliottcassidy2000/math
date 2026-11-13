"""Tail-law route to the C2 line.  For a digit-1 ternary clock B = S[u<-TT_j] (H(B) = A 3^j - 1) and a vertex y,
the one-out-neighbour tail (root r with N^+(r) = {y}, chain V of size k beaten by B and beating only r) has
H = H(B) + 2^k U_B(y)  (THM-5126 (A) / THM-5129).  If U_B(y)(j) = a' 3^j exactly, H = 3^j(a' 2^k + A) - 1.
Screen: both conventions (U = B_y, and the converse response E_y), j = 1..4."""
import sys, gzip, re, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digits import digits
from tern_check import parse, blowup
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def clocks():
    out = []
    L = gzip.open(os.path.join(ROOT, 'data', 'tern_all_orders.txt.gz'), 'rt').read().splitlines()
    for l in L:
        if not l.startswith('TERN'): continue
        p = l.split()
        if int(re.search(r' C=(-?\d+)', l).group(1)) != -1: continue
        n = int(p[1]); out.append((n, parse(p[2], n), int(p[3].split('=')[1])))
    seen = set()
    for l in gzip.open(os.path.join(ROOT, 'data', 'backscan_n10_b5.txt.gz'), 'rt').read().splitlines():
        if not l.startswith('TERN') or ' sigma=1 ' not in l: continue
        u = int(re.search(r' u=(\d+)', l).group(1))
        bk = tuple(sorted(tuple(map(int, x.split('>'))) for x in re.search(r'back=([\d>,]+)', l).group(1).strip(',').split(',')))
        if (u, bk) in seen: continue
        seen.add((u, bk)); n = 10; adj = [0]*n
        for i in range(n):
            for j in range(i+1, n):
                if (j, i) in bk: adj[j] |= 1 << i
                else: adj[i] |= 1 << j
        out.append((n, adj, u))
    return out
if __name__ == '__main__':
    C = clocks(); print('digit-1 clocks:', len(C), flush=True)
    hits = 0; near = []
    for (n, adj, u) in C:
        resp = {}
        for j in (1, 2, 3):
            A2, n2, block, newid = blowup(adj, n, u, j)
            H, sig, tau = digits(A2, n2)
            for y in range(n):
                if y == u: continue
                yy = newid[y]
                resp.setdefault(y, []).append((sig[yy] + H, tau[yy] + H))
        for y, vals in resp.items():
            for conv in (0, 1):
                r = [v[conv] for v in vals]
                if r[1] == 3*r[0] and r[2] == 9*r[0]:
                    hits += 1; print('PURE', n, u, y, 'conv', conv, r, flush=True)
                else:
                    d = abs(r[1] - 3*r[0]) + abs(r[2] - 9*r[0])
                    near.append((d, n, u, y, conv, r))
    near.sort()
    print('pure responses found:', hits, ' closest:', near[:5])
