"""Dominance-pruned slow-factor tree, propagated, with inclusion dominance of column spaces (exploratory:
column spaces are computed modulo the prime 2147483629 in the minimal 61-dimensional representation; state
vectors are exact small integers held in float64).

A factor u is SLOW at rate rho if its operator has rank >= 2 and H(T[u']) < rho^|u'| for every prefix u'.
A finite tree of slow factors (depth L) gives H(W) >= rho^(|W| - L - 1) for every bottleneck-free word W.

Node = (L, x): L = column space of the operator of u, x = raw state vector of the standalone tournament T[u].
Dominance at equal depth: (L, x) dominates (L', x') if L' is contained in L and x <= x' coordinatewise; then for
every continuation v, rank(u'v) <= rank(uv) and H(uv) <= H(u'v) (the raw transfer matrices are nonnegative), so
every slow continuation of the dominated node is a slow continuation of the dominating one.  Only undominated
nodes are expanded.  MODE = 'global' compares all nodes of a level, 'suffix3' only nodes with equal last three letters.
usage: python3 domtree2.py RHO MAXDEPTH NODECAP MODE"""
import sys, time, collections, numpy as np
src = open('domtree.py').read().split("Af = [M.astype")[0].replace("RHO = float(sys.argv[1]); MAXD = int(sys.argv[2]); CAP = int(sys.argv[3])", "RHO = float(sys.argv[1]); MAXD = int(sys.argv[2]); CAP = int(sys.argv[3]); MODE = sys.argv[4]")
exec(src)
Af = [M.astype(np.float64) for M in A]; hf = h.astype(np.float64); v0f = v0.astype(np.float64)
def prune(nodes):
    """nodes: list of (w, b, x).  Returns the undominated ones."""
    nodes.sort(key=lambda t: (-t[1].shape[1], float(hf @ t[2])))
    keep = []; KX = np.zeros((0, len(v0f))); KR = np.zeros(0, dtype=np.int64)
    for w, b, x in nodes:
        dom = False
        if len(keep):
            cand = np.nonzero(np.all(KX <= x, axis=1) & (KR >= b.shape[1]))[0]
            for i in cand:
                kb = keep[i][1]
                if col_basis(np.hstack([kb, b])).shape[1] == kb.shape[1]: dom = True; break
        if not dom:
            keep.append((w, b, x)); KX = np.vstack([KX, x[None, :]]); KR = np.append(KR, b.shape[1])
    return keep
level = [((l,), basis((l,)), v0f.copy()) for l in range(8)]
t0 = time.time(); print("rate %.5f mode %s" % (RHO, MODE), flush=True)
for depth in range(1, MAXD + 1):
    raw = len(level)
    if MODE == 'global': level = prune(level)
    else:
        groups = collections.defaultdict(list)
        for nd in level: groups[nd[0][-3:]].append(nd)
        level = [nd for g in groups.values() for nd in prune(g)]
    ranks = collections.Counter(b.shape[1] for w, b, x in level)
    lo = min(float(hf @ x) for w, b, x in level) if level else 0
    print("depth %2d: %7d candidates -> %6d undominated; ranks %s; least H %.0f (rho^m = %.1f)  [%.0f s]" % (depth, raw, len(level), dict(sorted(ranks.items())), lo, RHO ** depth, time.time() - t0), flush=True)
    if not level: print("TREE CLOSED at depth", depth); break
    if len(level) > CAP: print("node cap reached"); break
    if depth == MAXD: break
    nxt = {}; m = depth + 1; thr = RHO ** m * (1 + 1e-9)
    for w, b, x in level:
        for l in range(8):
            b2 = col_basis(modmul(Mm[l], b))
            if b2.shape[1] <= 1: continue
            y = Af[l & ((1 << min(3, m - 1)) - 1)] @ x
            if hf @ y < thr: nxt[(b2.tobytes(), y.tobytes())] = ((w + (l,))[-3:], b2, y)
    level = list(nxt.values())
