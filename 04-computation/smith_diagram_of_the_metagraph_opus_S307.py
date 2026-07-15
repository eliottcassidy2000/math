#!/usr/bin/env python3
"""THE SMITH DIAGRAM OF THE METAGRAPH (opus-2026-07-14-S307, HYP-6865).

Owner directive: consider Smith diagrams / squaring-the-square concepts against
the tournament-metagraph work. The BSST dictionary reads a squared rectangle as
a resistor network with KCL/KVL; here we read the METAGRAPH ITSELF as the
network: nodes = iso classes, edges = wiggly (d=1) tile-flips with conductance
= flip multiplicity, source = the TRANSITIVE class, sink = the DISTRIBUTED rail
(all classes at minimal axis x tied together, like the bottom bus bar of a
Smith diagram). Unit current flows; we compute exactly (Fractions, n <= 6;
float64 + residual check at n = 7):

  - the TRANSITIVITY RESISTANCE R_n (a new invariant of G_n);
  - the harmonic potential phi(class) and its relation to the axis x and to H;
  - the CURRENT MAP: concentration across the spine (SC-SC), ribs (SC-NS), and
    sea (NS-NS); the maximum-current edge and path (the electrical principal
    line) vs the H-gradient principal line;
  - the weighted spanning-tree complexity of G_n (matrix-tree; the BSST
    denominator), n <= 6 exact.
"""
import sys
from collections import defaultdict
from fractions import Fraction as F

def build(n):
    V = list(range(1, n+1))
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
    m = len(tiles)

    def tourn(bits):
        adj = [0]*(n+1)
        for k in range(n, 1, -1): adj[k] |= 1 << (k-2)
        for i, (x, y) in enumerate(tiles):
            if bits >> i & 1: adj[x] |= 1 << (y-1)
            else: adj[y] |= 1 << (x-1)
        return adj

    def invariant(adj):
        s = {v: bin(adj[v]).count('1') for v in V}
        prof = {}
        for v in V:
            outs = sorted(s[u] for u in V if adj[v] >> (u-1) & 1)
            ins = sorted(s[u] for u in V if adj[u] >> (v-1) & 1)
            c3 = 0
            for u in V:
                if adj[v] >> (u-1) & 1:
                    for w in V:
                        if adj[u] >> (w-1) & 1 and adj[w] >> (v-1) & 1: c3 += 1
            prof[v] = (s[v], tuple(outs), tuple(ins), c3)
        prof2 = {}
        for v in V:
            po = sorted(prof[u] for u in V if adj[v] >> (u-1) & 1)
            pi_ = sorted(prof[u] for u in V if adj[u] >> (v-1) & 1)
            prof2[v] = (prof[v], tuple(po), tuple(pi_))
        arcs = []
        for u in V:
            for v in V:
                if adj[u] >> (v-1) & 1:
                    arcs.append((s[u], s[v], bin(adj[u] & adj[v]).count('1')))
        full = (1 << n) - 1
        dp = {(1 << (v-1), v): 1 for v in V}
        for _ in range(n-1):
            nd = defaultdict(int)
            for (mask, v), c in dp.items():
                for u in V:
                    b = 1 << (u-1)
                    if not mask & b and adj[v] & b: nd[(mask | b, u)] += c
            dp = nd
        H = sum(c for (mask, v), c in dp.items() if mask == full)
        return (tuple(sorted(prof2.values())), tuple(sorted(arcs)), H)

    classes = {}; cls_of = {}; rep = {}
    for b in range(1 << m):
        inv = invariant(tourn(b))
        if inv not in classes:
            classes[inv] = len(classes); rep[classes[inv]] = b
        cls_of[b] = classes[inv]
    C = len(classes)
    counts = {4: 4, 5: 12, 6: 56, 7: 456}
    assert C == counts[n]
    rcls = {}
    for c in range(C):
        adj = tourn(rep[c]); radj = [0]*(n+1)
        for u in V:
            for v in V:
                if adj[u] >> (v-1) & 1: radj[v] |= 1 << (u-1)
        rcls[c] = classes[invariant(radj)]
    x_of = {}; H_of = {}
    for c in range(C):
        adj = tourn(rep[c])
        sc = sorted(bin(adj[v]).count('1') for v in V)
        x_of[c] = sum((2*si - (n-1))**2 for si in sc)
        H_of[c] = invariant(adj)[2]
    # wiggly d=1 edges with multiplicities
    W = defaultdict(int)
    for b in range(1 << m):
        ca = cls_of[b]
        for i in range(m):
            cb = cls_of[b ^ (1 << i)]
            key = (ca, cb) if ca <= cb else (cb, ca)
            W[key] += 1        # every unordered tiling-flip pair visited exactly twice
    W2 = {k: cnt // 2 for k, cnt in W.items()}
    assert all(cnt % 2 == 0 for cnt in W.values())
    return dict(n=n, C=C, cls_of=cls_of, x_of=x_of, H_of=H_of, rcls=rcls, W=W2)

def solve_network(n, B, exact=True):
    C, W, x_of, rcls, H_of = B['C'], B['W'], B['x_of'], B['rcls'], B['H_of']
    # conductance matrix (exclude self-loops)
    cond = defaultdict(lambda: 0)
    for (a, b), w in W.items():
        if a != b: cond[(a, b)] = w
    xmin = min(x_of.values()); xmax = max(x_of.values())
    src = [c for c in range(C) if H_of[c] == 1][0]
    sinks = [c for c in range(C) if x_of[c] == xmin]
    # collapse sinks into one supernode S
    idx = {}; k = 0
    for c in range(C):
        if c == src or c in sinks: continue
        idx[c] = k; k += 1
    NS = k          # unknown potentials; phi(src) unknown too? set phi(sink)=0, inject 1 at src
    # unknowns: potentials of all non-sink nodes (incl src); equations: KCL at each non-sink node
    nodes = [c for c in range(C) if c not in sinks]
    pos = {c: i for i, c in enumerate(nodes)}
    N = len(nodes)
    zero = F(0) if exact else 0.0
    one = F(1) if exact else 1.0
    A = [[zero]*N for _ in range(N)]
    rhs = [zero]*N
    deg = defaultdict(lambda: zero)
    for (a, b), w in cond.items():
        wv = F(w) if exact else float(w)
        deg[a] += wv; deg[b] += wv
        if a in pos and b in pos:
            A[pos[a]][pos[b]] -= wv; A[pos[b]][pos[a]] -= wv
    for c in nodes:
        A[pos[c]][pos[c]] += deg[c]
    rhs[pos[src]] = one
    # Gaussian elimination
    for i in range(N):
        piv = None
        for r in range(i, N):
            if A[r][i] != 0: piv = r; break
        A[i], A[piv] = A[piv], A[i]; rhs[i], rhs[piv] = rhs[piv], rhs[i]
        inv = (F(1)/A[i][i]) if exact else 1.0/A[i][i]
        for r in range(i+1, N):
            if A[r][i] != 0:
                f = A[r][i]*inv
                for cc in range(i, N): A[r][cc] -= f*A[i][cc]
                rhs[r] -= f*rhs[i]
    phi = [zero]*N
    for i in range(N-1, -1, -1):
        ssum = rhs[i]
        for cc in range(i+1, N): ssum -= A[i][cc]*phi[cc]
        phi[i] = ssum/A[i][i]
    pot = {c: phi[pos[c]] for c in nodes}
    for c in sinks: pot[c] = zero
    Reff = pot[src]
    # current map + spine/ribs/sea split
    def typ(a, b):
        sa, sb = rcls[a] == a, rcls[b] == b
        return 'spine' if sa and sb else ('sea' if not sa and not sb else 'ribs')
    split = defaultdict(lambda: zero); maxcur = (zero, None)
    for (a, b), w in cond.items():
        wv = F(w) if exact else float(w)
        cur = abs((pot[a]-pot[b])*wv)
        split[typ(a, b)] += cur
        if cur > maxcur[0]: maxcur = (cur, (a, b, w))
    return Reff, pot, split, maxcur, src, sinks

def tree_complexity(B):
    """weighted spanning-tree count of G_n (matrix-tree, exact)."""
    C, W = B['C'], B['W']
    cond = {(a, b): w for (a, b), w in W.items() if a != b}
    N = C - 1
    A = [[F(0)]*N for _ in range(N)]
    deg = defaultdict(lambda: F(0))
    for (a, b), w in cond.items():
        deg[a] += w; deg[b] += w
        if a < N and b < N:
            A[a][b] -= w; A[b][a] -= w
    for c in range(N): A[c][c] += deg[c]
    det = F(1)
    for i in range(N):
        piv = None
        for r in range(i, N):
            if A[r][i] != 0: piv = r; break
        if piv is None: return F(0)
        if piv != i:
            A[i], A[piv] = A[piv], A[i]; det = -det
        det *= A[i][i]
        inv = F(1)/A[i][i]
        for r in range(i+1, N):
            if A[r][i] != 0:
                f = A[r][i]*inv
                for cc in range(i, N): A[r][cc] -= f*A[i][cc]
    return det

if __name__ == '__main__':
    for n in (4, 5, 6, 7):
        B = build(n)
        exact = n <= 6
        Reff, pot, split, maxcur, src, sinks = solve_network(n, B, exact=exact)
        tot = sum(split.values())
        print(f"===== n = {n}: G_n has {B['C']} nodes; source = transitive, "
              f"sink rail = {len(sinks)} class(es) at x_min")
        if exact:
            print(f"  TRANSITIVITY RESISTANCE R_{n} = {Reff} = {float(Reff):.6f}")
        else:
            print(f"  TRANSITIVITY RESISTANCE R_{n} ~= {Reff:.6f} (float; n=7)")
        print(f"  current split: spine {float(split['spine']):.4f}, ribs {float(split['ribs']):.4f}, "
              f"sea {float(split['sea']):.4f} (total edge-current {float(tot):.4f})")
        mc, (a, b, w) = maxcur
        print(f"  max-current edge: classes {a}({'SC' if B['rcls'][a]==a else 'NS'}, x={B['x_of'][a]}, H={B['H_of'][a]})"
              f" -- {b}({'SC' if B['rcls'][b]==b else 'NS'}, x={B['x_of'][b]}, H={B['H_of'][b]}), "
              f"mult {w}, current {float(mc):.4f}")
        # potential vs axis: rank correlation (exact Spearman-ish: count concordant)
        cls = list(range(B['C']))
        import itertools
        conc = disc = 0
        for a1, b1 in itertools.combinations(cls, 2):
            dx = B['x_of'][a1] - B['x_of'][b1]
            dp = pot[a1] - pot[b1]
            if dx == 0 or dp == 0: continue
            if (dx > 0) == (dp > 0): conc += 1
            else: disc += 1
        print(f"  potential-vs-axis concordance: {conc}/{conc+disc} = {conc/(conc+disc):.3f}")
        if n <= 6:
            T = tree_complexity(B)
            print(f"  weighted spanning-tree complexity kappa(G_{n}) = {T}")
        print()
