# Stockmeyer's tournament family (D_{m,n}, D'_{m,n}) of "A census of non-reconstructable digraphs I" (1981), as described
# in Egrot-Hirsch arXiv:2008.01327 Def. 6.1 / Fig. 9 (tournament family D: (w,x,y,z) = (1,0,0,1) versus (0,1,1,0)).
# Checks: same deck, non-isomorphic, converse relation, self-converse, and H/hc (via hpaths.exe) -> Theorem H at odd order.
import sys, subprocess, itertools
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def odd(z):
    while z % 2 == 0: z //= 2
    return z

def build(m, n, wxyz):
    w, x, y, z = wxyz
    p = 2 ** m + 2 ** n if n >= 0 else 2 ** m
    M = [[0] * p for _ in range(p)]
    A = 2 ** m
    for i in range(1, p + 1):
        for j in range(1, p + 1):
            if i == j: continue
            if (i <= A and j <= A) or (i > A and j > A):
                M[i - 1][j - 1] = 1 if odd(j - i) % 4 == 1 else 0
            elif i <= A < j:
                M[i - 1][j - 1] = w if (i + j) % 2 == 0 else x
            else:
                M[i - 1][j - 1] = y if (i + j) % 2 == 0 else z
    for i in range(p):
        for j in range(p):
            if i != j: assert M[i][j] + M[j][i] == 1, (m, n, i, j)
    return M

def G(M):
    g = nx.DiGraph(); g.add_nodes_from(range(len(M)))
    g.add_edges_from((i, j) for i in range(len(M)) for j in range(len(M)) if M[i][j])
    return g

def iso(a, b):
    if sorted(d for _, d in a.out_degree()) != sorted(d for _, d in b.out_degree()): return False
    return DiGraphMatcher(a, b).is_isomorphic()

def deck_equal(a, b):
    ca = [a.subgraph([v for v in a if v != x]).copy() for x in a]
    cb = [b.subgraph([v for v in b if v != x]).copy() for x in b]
    used = [False] * len(cb)
    for c in ca:
        for k, d in enumerate(cb):
            if not used[k] and iso(c, d):
                used[k] = True; break
        else:
            return False
    return True

cases = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1)]
inp = []
rows = []
for m, n in cases:
    D = build(m, n, (1, 0, 0, 1)); Dp = build(m, n, (0, 1, 1, 0))
    a, b = G(D), G(Dp)
    p = len(D)
    res = dict(order=p, m=m, n=n, iso=iso(a, b), conv=iso(b, a.reverse(copy=True)), scD=iso(a, a.reverse(copy=True)), scDp=iso(b, b.reverse(copy=True)), deck=deck_equal(a, b))
    rows.append(res)
    for M in (D, Dp):
        inp.append(f"{p}\n" + "\n".join(" ".join(map(str, r)) for r in M))
out = subprocess.run(["./hpaths.exe"], input="\n".join(inp) + "\n", capture_output=True, text=True).stdout.split("\n")
for k, res in enumerate(rows):
    _, H1, c1 = out[2 * k].split(); _, H2, c2 = out[2 * k + 1].split()
    H1, c1, H2, c2 = int(H1), int(c1), int(H2), int(c2)
    p = res['order']
    thmH = (H1 == H2) if p % 2 == 0 else (H1 - H2 == 2 * (c1 - c2))
    print(f"order {p} (m={res['m']}, n={res['n']}): same deck {res['deck']}, isomorphic {res['iso']}, D' = converse(D) {res['conv']}, D self-converse {res['scD']}, D' self-converse {res['scDp']}; "
          f"H {H1}/{H2}, hc {c1}/{c2}; Theorem H {'OK' if thmH else 'FAIL'}")
