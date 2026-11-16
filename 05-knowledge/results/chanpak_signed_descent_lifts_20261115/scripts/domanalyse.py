"""What do the column spaces of slow factors look like?  Build the slow tree (no pruning except exact duplicates of
(space, vector)) to a given depth and report, per depth: number of nodes, distribution of operator ranks, number of
distinct column spaces, and how many nodes survive when dominance uses general inclusion of column spaces inside
groups with the same last three letters."""
import pickle, sys, time, collections, numpy as np
src = open('domtree.py').read().split("Af = [M.astype")[0].replace("RHO = float(sys.argv[1]); MAXD = int(sys.argv[2]); CAP = int(sys.argv[3])", "RHO = float(sys.argv[1]); MAXD = int(sys.argv[2])")
exec(src)
Af = [M.astype(np.float64) for M in A]; hf = h.astype(np.float64); v0f = v0.astype(np.float64)
def key(b): return b.tobytes()
level = [((l,), basis((l,)), v0f.copy()) for l in range(8)]
t0 = time.time()
for depth in range(1, MAXD + 1):
    ranks = collections.Counter(b.shape[1] for w, b, x in level)
    spaces = {key(b) for w, b, x in level}
    print("depth %2d: %7d slow nodes; ranks %s; distinct column spaces %d  [%.0f s]" % (depth, len(level), dict(sorted(ranks.items())), len(spaces), time.time() - t0), flush=True)
    if depth >= 6:
        # general inclusion dominance inside suffix-3 groups
        groups = collections.defaultdict(list)
        for w, b, x in level: groups[w[-3:]].append((b, x))
        kept = 0
        for g, nodes in groups.items():
            nodes.sort(key=lambda t: (-t[0].shape[1], float(hf @ t[1])))
            keep = []
            for b, x in nodes:
                dom = False
                for kb, kx in keep:
                    if kb.shape[1] >= b.shape[1] and np.all(kx <= x):
                        # is col(b) inside col(kb)?  rank of the stacked basis
                        if col_basis(np.hstack([kb, b])).shape[1] == kb.shape[1]: dom = True; break
                if not dom: keep.append((b, x))
            kept += len(keep)
        print("          after inclusion dominance inside suffix-3 groups: %d nodes (%d groups)" % (kept, len(groups)), flush=True)
    if depth == MAXD: break
    nxt = {}; m = depth + 1; thr = RHO ** m * (1 + 1e-9)
    for w, b, x in level:
        for l in range(8):
            b2 = col_basis(modmul(Mm[l], b))
            if b2.shape[1] <= 1: continue
            y = Af[l & ((1 << min(3, m - 1)) - 1)] @ x
            if hf @ y < thr: nxt[(key(b2), y.tobytes())] = (w + (l,), b2, y)
    level = list(nxt.values())
