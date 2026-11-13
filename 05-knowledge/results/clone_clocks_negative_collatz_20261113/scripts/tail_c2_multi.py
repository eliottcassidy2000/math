"""Tail-law route to the C2 line with a general out-set O of the new root (THM-5126 (A)):
   H(L o 0^k) = c0 + M 2^k,  c0 = H(B) - sum_{pairs {e,f} in O} W_ef,  M = sum_{y in O} U_B(y),
where U_B(y) = #(P1,P2) covering B, P2 starting at y (P1 possibly empty) and W_ef = # ordered two-path covers
(both nonempty) with the first ending at one of e,f and the second starting at the other.
With B = S[u<-TT_j] a digit-1 clock, the C2 line needs c0(j) = A' 3^j - 1 and M(j) = a' 3^j.  Screen |O| <= 4."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digits import subset_tables
from tern_check import blowup
from tail_c2_search import clocks
def responses(adj, n):
    end, start, h = subset_tables(adj, n)
    full = (1 << n) - 1; H = h[full]
    U = []
    for y in range(n):
        rest = full & ~(1 << y); tot = 0; A = rest
        while True:
            tot += h[A] * start[full & ~A][y]
            if A == 0: break
            A = (A - 1) & rest
        U.append(tot)
    W = {}
    for e in range(n):
        for f in range(e+1, n):
            tot = 0
            S = full
            sub = (full - 1) & full
            while sub:
                comp = full & ~sub
                if comp:
                    tot += end[sub][e] * start[comp][f] + end[sub][f] * start[comp][e]
                sub = (sub - 1) & full
            W[(e, f)] = tot
    return H, U, W
if __name__ == '__main__':
    C = clocks(); maxO = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    hits = []; best = []
    for ci, (n, adj, u) in enumerate(C):
        data = []
        for j in (1, 2, 3):
            A2, n2, block, newid = blowup(adj, n, u, j)
            H, U, W = responses(A2, n2)
            data.append((H, U, W, newid))
        verts = [y for y in range(n) if y != u]
        for size in range(1, maxO + 1):
            for O in itertools.combinations(verts, size):
                c0 = []; M = []
                for (H, U, W, nid) in data:
                    o = [nid[y] for y in O]
                    c0.append(H - sum(W[(min(a, b), max(a, b))] for a, b in itertools.combinations(o, 2)))
                    M.append(sum(U[y] for y in o))
                pureM = M[1] == 3*M[0] and M[2] == 9*M[0] and M[0] > 0
                c2off = c0[1] == 3*c0[0] + 2 and c0[2] == 9*c0[0] + 8
                if pureM and c2off:
                    hits.append((n, u, O, c0, M)); print('C2-LINE TAIL', n, u, O, c0, M, flush=True)
                if pureM:
                    # record the offset digit:  c0 = A' 3^j - d  if ternary-shaped
                    if (c0[2] - 3*c0[1]) == (c0[1] - 3*c0[0]):      # c0 = A' 3^j + C ; digit d = -C
                        d = (c0[1] - 3*c0[0]) // 2
                        best.append((d, n, u, O, c0, M))
    print('clocks', len(C), ' C2-line tails:', len(hits))
    import collections
    print('pure-response tails by offset digit:', sorted(collections.Counter(b[0] for b in best).items())[:20])
