"""Adversarial search for counterexamples to HYP-9279 (Z) [I1: a cut vertex z is balanced against every vertex of the
interior blocks D_2..D_(k-1) of T - z] and (B) at cut vertices [I3], using the exact transfer chain of cutvertex_chain.py.
Hill-climbing over: block sizes/structures (random strong blocks), z's arcs.  Objective: minimum over interior-block vertices
of min(p, 1-p) (I1), or best partner balance of z (I3)."""
import random, sys
from fractions import Fraction as Fr
from cutvertex_chain import probs_cut, random_strong_block, components

def build(blocks_adj, zout):
    sizes = [len(B) for B in blocks_adj]; n = 1 + sum(sizes)
    A = [[False] * n for _ in range(n)]; idx = []; c = 1
    for B in blocks_adj:
        s = len(B); ids = list(range(c, c + s)); c += s; idx.append(ids)
        for i in range(s):
            for j in range(s):
                if B[i][j]: A[ids[i]][ids[j]] = True
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            for x in idx[a]:
                for y in idx[b]: A[x][y] = True
    for x in range(1, n):
        A[0][x] = zout[x - 1]; A[x][0] = not zout[x - 1]
    return A, idx

def objective(blocks_adj, zout, mode):
    A, idx = build(blocks_adj, zout)
    if not any(A[0][x] for x in idx[0]) or not any(A[x][0] for x in idx[-1]): return None
    H, blocks, res = probs_cut(A, 0)
    bal = {v: min(p, 1 - p) for v, p in res.items()}
    if mode == 'I1':
        if len(blocks) < 3: return None
        return min(bal[v] for b in blocks[1:-1] for v in b)
    return max(bal.values())

def search(seed, mode, steps):
    rng = random.Random(seed)
    k = rng.randint(3, 7) if mode == 'I1' else 2
    blocks = [random_strong_block(rng.choice([1, 1, 1, 3, 3, 4, 5] if mode == 'I1' else [1, 3, 4, 5, 6]), rng) for _ in range(k)]
    zout = [rng.random() < 0.5 for _ in range(sum(len(b) for b in blocks))]
    cur = objective(blocks, zout, mode)
    while cur is None:
        zout = [rng.random() < 0.5 for _ in range(sum(len(b) for b in blocks))]; cur = objective(blocks, zout, mode)
    for t in range(steps):
        nb = [b for b in blocks]; nz = list(zout); r = rng.random()
        if r < 0.6:
            i = rng.randrange(len(nz)); nz[i] = not nz[i]
        elif r < 0.8:
            j = rng.randrange(len(nb)); s = rng.choice([1, 1, 3, 3, 4, 5] if mode == 'I1' else [3, 3, 4, 5, 6, 7])
            old = len(nb[j]); nb[j] = random_strong_block(s, rng)
            off = sum(len(b) for b in blocks[:j]); nz = nz[:off] + [rng.random() < 0.5 for _ in range(s)] + nz[off + old:]
        elif r < 0.9 and len(nb) < 9 and mode == 'I1':
            j = rng.randrange(len(nb) + 1); s = rng.choice([1, 1, 3]); nb.insert(j, random_strong_block(s, rng))
            off = sum(len(b) for b in nb[:j]); nz = nz[:off] + [rng.random() < 0.5 for _ in range(s)] + nz[off:]
        elif len(nb) > 3 and mode == 'I1':
            j = rng.randrange(len(nb)); off = sum(len(b) for b in nb[:j]); s = len(nb[j]); del nb[j]; nz = nz[:off] + nz[off + s:]
        if sum(len(b) for b in nb) > 30: continue
        val = objective(nb, nz, mode)
        if val is not None and val <= cur:
            blocks, zout, cur = nb, nz, val
    return cur, blocks, zout

if __name__ == '__main__':
    mode = sys.argv[1]; runs = int(sys.argv[2]); steps = int(sys.argv[3])
    best = None
    for s in range(runs):
        val, blocks, zout = search(1000 + s, mode, steps)
        print(f"run {s}: {mode} objective {val} = {float(val):.6f}; block sizes {[len(b) for b in blocks]}; zout {''.join('1' if x else '0' for x in zout)}", flush=True)
        if best is None or val < best[0]: best = (val, blocks, zout)
    print("BEST", best[0], float(best[0]), [len(b) for b in best[1]], ''.join('1' if x else '0' for x in best[2]))
