#!/usr/bin/env python3
"""
ramsey_tabu_general.py — general-graph tabu search for R(B_{n-1}, B_n) >= 4n-1
witnesses: G on N = 4n-2 vertices, every edge with <= n-2 common neighbours,
every non-edge with <= n-1 common non-neighbours.

No structural assumption (witnesses need not be vertex-transitive).
Incremental exact maintenance of C = A^2 (common neighbours) under edge flips;
common non-neighbours derived: for i<j nonadjacent,
   cn(i,j) = N - 2 - deg(i) - deg(j) + C[i,j].
Energy = sum over violating pairs of the excess.
Tabu on flipped pairs; restarts; optional seeding from a two-level circulant.

Usage: python3 ramsey_tabu_general.py n [iters] [seed] [restarts]
"""
import numpy as np
import sys, time

def energy_full(A, n, N):
    M = A.astype(np.int64)
    C = M @ M
    deg = M.sum(1)
    e = 0
    on = C - (n - 2)
    e_on = np.where((M == 1) & (on > 0), on, 0)
    e += int(np.triu(e_on, 1).sum())
    cn = (N - 2) - deg[:, None] - deg[None, :] + C
    off = cn - (n - 1)
    mask_off = (M == 0)
    np.fill_diagonal(mask_off, False)
    e_off = np.where(mask_off & (off > 0), off, 0)
    e += int(np.triu(e_off, 1).sum())
    return e, C, deg

def verify(A, n, N):
    e, C, deg = energy_full(A, n, N)
    return e == 0

def run(n, iters, seed, restarts):
    N = 4 * n - 2
    rng = np.random.default_rng(seed)
    bestE_global = 10 ** 9
    for r in range(restarts):
        # seed: random regular-ish graph with density ~ (N-2)/2 / (N-1) ~ 1/2
        A = (rng.random((N, N)) < 0.5).astype(np.int8)
        A = np.triu(A, 1)
        A = A + A.T
        e, C, deg = energy_full(A, n, N)
        best = e
        tabu = {}
        t0 = time.time()
        it = 0
        last_improve = 0
        while it < iters:
            it += 1
            # candidate moves: sample a batch of random pairs, pick best delta
            bsz = 24
            us = rng.integers(0, N, size=bsz)
            vs = rng.integers(0, N, size=bsz)
            bestdelta = None; bestuv = None
            for u, v in zip(us, vs):
                if u == v:
                    continue
                if tabu.get((min(u, v), max(u, v)), 0) > it:
                    continue
                # compute delta by local formula
                delta = flip_delta(A, C, deg, u, v, n, N)
                if bestdelta is None or delta < bestdelta:
                    bestdelta = delta; bestuv = (u, v)
            if bestuv is None:
                continue
            u, v = bestuv
            # accept if improving or with small probability (noise)
            if bestdelta <= 0 or rng.random() < 0.02:
                apply_flip(A, C, deg, u, v)
                e += bestdelta
                tabu[(min(u, v), max(u, v))] = it + rng.integers(8, 40)
                if e < best:
                    best = e; last_improve = it
                if e == 0:
                    ok = verify(A, n, N)
                    print(f"[r{r}] ZERO at it={it} verify={ok} ({time.time()-t0:.0f}s)", flush=True)
                    if ok:
                        fn = f"ramsey_general_witness_n{n}_seed{seed}_r{r}"
                        np.save(fn + ".npy", A)
                        with open(fn + ".txt", "w") as f:
                            f.write("".join(str(A[i, j]) for j in range(N) for i in range(N)) + "\n")
                        print("WITNESS SAVED " + fn, flush=True)
                        return True
            if it - last_improve > 30000:
                break  # stagnated
            if it % 20000 == 0:
                print(f"[r{r}] it={it} e={e} best={best} ({time.time()-t0:.0f}s)", flush=True)
        bestE_global = min(bestE_global, best)
        print(f"[r{r}] done best={best}", flush=True)
    print(f"NO WITNESS: best {bestE_global}", flush=True)
    return False

def flip_delta(A, C, deg, u, v, n, N):
    """Exact energy delta for flipping edge (u,v). Recompute affected pairs:
    pairs (u,*), (v,*), and (x,y) with C[x,y] changed: only pairs where one of
    u,v is a common neighbour candidate: C[x,y] changes iff {x,y} ~ u,v link:
    C[x,y] += dA*(A[x,u]... ) — cheaper: pairs (x, u-or-v) energy and pairs
    (x,y) with x~u? Actually C[x,y] = sum_w A[x,w]A[w,y]; flipping A[u,v]
    changes terms w=v for pairs (x=u? no: C[x,y] with w=u needs A[x,u]A[u,y]:
    changes iff (x,y) has u as midpoint and one leg is v... Specifically
    dC[x,y] = dA * (A[x,u]*[y==v is wrong]...).  Do it by direct small recompute:
    affected pairs P = {(x,y): x in {u,v} or y in {u,v}} ∪ {(x,y): A[x,u]&A[..]}
    Simpler exact: dC[x,y] != 0 only if (x in nbr(u)∪{u,v}) and similar — to
    stay safe recompute energy contribution of all pairs involving u or v as
    endpoints, plus all pairs (x,y) with x,y in nbr(u)∪nbr(v)... that is O(N^2)
    worst case; N<=400 -> 160k ops, fine at Python+numpy row level."""
    # numpy vectorized local recompute: energy over pairs with endpoint in {u,v}
    # plus pairs whose C changes: those are (x,y) with (x,y) both adjacent... only
    # pairs where u or v is a MIDPOINT: (x,y) with A[x,u]=A[u,y]=1 changes via w=u
    # only if the flipped edge is a leg: legs are (u,v): so w=u contributes A[x,u]A[u,y]:
    # flipped entry A[u,v]: affects C[x,y] when (x==v or y==v) and w==u? no:
    # C[x,y] terms containing A[u,v]: w=u with x==v? A[v,u]A[u,y] -> pairs (v,y);
    # w=v: A[x,v]A[v,y] with x or y == u -> pairs (u,y),(x,u).
    # => C changes ONLY for pairs with endpoint u or v.  So local energy = pairs
    # with an endpoint in {u,v}.
    old = pair_energy(A, C, deg, u, v, n, N)
    apply_flip(A, C, deg, u, v)
    new = pair_energy(A, C, deg, u, v, n, N)
    apply_flip(A, C, deg, u, v)  # revert
    return new - old

def pair_energy(A, C, deg, u, v, n, N):
    """energy of all pairs with an endpoint in {u,v} (C only changes there,
    but deg(u),deg(v) change also affects non-edge counts for pairs with
    endpoint u or v — same set).  Exact, vectorized."""
    e = 0
    for w in (u, v):
        Crow = C[w].astype(np.int64); Arow = A[w]
        mask = np.ones(N, dtype=bool); mask[w] = False
        on = Crow - (n - 2)
        off = (N - 2) - (deg[w] + deg) + Crow - (n - 1)
        eon = np.where((Arow == 1) & mask & (on > 0), on, 0).sum()
        eoff = np.where((Arow == 0) & mask & (off > 0), off, 0).sum()
        e += int(eon + eoff)
    # pair (u,v) counted twice
    if A[u, v]:
        d = C[u, v] - (n - 2)
        if d > 0: e -= int(d)
    else:
        d = (N - 2) - deg[u] - deg[v] + C[u, v] - (n - 1)
        if d > 0: e -= int(d)
    return e

def apply_flip(A, C, deg, u, v):
    old = A[u, v]
    newv = 1 - old
    d = newv - old
    A[u, v] = A[v, u] = newv
    # C updates: C[v,y] += d*A[u,y] (w=u), C[x,v]... symmetric; C[u,y] += d*A[v,y]
    C[v, :] += d * A[u, :].astype(np.int64); C[:, v] = C[v, :]
    C[u, :] += d * A[v, :].astype(np.int64); C[:, u] = C[u, :]
    # C[u,u], C[v,v] diagonal adjust (self walks): C[u,u] = deg(u) etc; keep consistent:
    C[u, u] = 0; C[v, v] = 0  # never used (diagonal ignored)
    # note: C[v,y] change via w=u used A[u,y] AFTER flip for y=u? A[u,u]=0 fine.
    deg[u] += d; deg[v] += d
    return

if __name__ == "__main__":
    n = int(sys.argv[1])
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 300000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    print(f"general tabu: n={n} N={4*n-2} iters={iters} seed={seed}")
    run(n, iters, seed, restarts)
