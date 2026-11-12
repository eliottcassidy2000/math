# Verifier (lane density, S94). Independent recount of D-10 (balanced cell of F1) and of the quantity in the hypothesis
# of (BC) (D-06): ordered quadruples (h1,h2,h1',h2') of n-bit words with {h1,h2} != {h1',h2'} and equal N.
# Independence: states from the g-recursion (eta(h) = g(2h)); counting by sort + run lengths; small n also by Counter.
import sys, numpy as np
from math import gcd
from collections import Counter
def states(n):
    G = np.zeros(1 << (n + 1), np.int64); G[0] = 1
    for v in range(1, 1 << (n + 1)):
        u = v >> 1
        G[v] = G[u] if v & 1 else G[u] + u + (u & 1)
    h = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
    return G[2 * h], 2 * h + 1
def runs(sorted_vals):
    # run lengths of a sorted array
    if sorted_vals.size == 0: return np.zeros(0, np.int64)
    brk = np.flatnonzero(np.diff(sorted_vals)) + 1
    edges = np.concatenate(([0], brk, [sorted_vals.size]))
    return np.diff(edges)
sigma = int(sys.argv[1])
for n in [int(a) for a in sys.argv[2:]]:
    x, y = states(n); W = x.size
    cop = np.array([gcd(int(a), int(b)) == 1 for a, b in zip(x, y)])
    full = x[:, None] * x[None, :] + (1 << sigma) * (y[:, None] * y[None, :])
    assert int(full.max()) < 2 ** 62
    # unordered pairs h1 <= h2
    iu = np.triu_indices(W)
    vals = np.sort(full[iu]); M = vals.size
    r = runs(vals)
    distinct = r.size; E = int((r * r).sum()); maxr = int(r.max())
    span_odd = (int(vals[-1]) - int(vals[0])) // 2 + 1
    assert np.all(vals % 2 == 1)
    # ordered energy and the (BC) quadruple count
    allv = np.sort(full.ravel()); R = runs(allv)
    Eord = int((R * R).sum())
    trivial = 4 * (W * (W - 1) // 2) + W
    Q = Eord - trivial
    # collisions sharing an element: row collisions (same h1, h2 != h2'), split by coprimality of the shared state
    rc = 0; rc_cop = 0
    for i in range(W):
        rr = runs(np.sort(full[i])); k = int((rr * (rr - 1)).sum())
        rc += k
        if cop[i]: rc_cop += k
    print("sigma=%d n=%d M=%d (=2^(n-2)(2^(n-1)+1): %s) distinct=%d (%.4f of M, %.4f of span) E=%d E/M=%.4f max r=%d | coprime states %d/%d, same-row ordered collisions %d, in coprime rows %d | (BC): quadruples Q=%d, Q/4^n=%.5f"
          % (sigma, n, M, M == (1 << (n - 2)) * ((1 << (n - 1)) + 1), distinct, distinct / M, distinct / span_odd, E, E / M, maxr,
             int(cop.sum()), W, rc, rc_cop, Q, Q / 4 ** n))
    lo, hi = (1 << sigma) * 4 ** n, ((1 << (sigma + 2)) + 4) * 4 ** n
    assert lo < int(vals[0]) and int(vals[-1]) < hi
    if n <= 9:
        C = Counter()
        xs = [int(a) for a in x]; ys = [int(b) for b in y]
        for i in range(W):
            for j in range(i, W):
                C[xs[i] * xs[j] + (1 << sigma) * ys[i] * ys[j]] += 1
        assert len(C) == distinct and sum(v * v for v in C.values()) == E and max(C.values()) == maxr
        print("    Counter cross-check (pure Python integers): agrees")
    sys.stdout.flush()
