"""S40 audit: an independent attempt to beat the natural order of P_p at the 'unbeaten' primes of HYP-9287.

Heuristics (all different from the session's best-improvement insertion + random kicks):
  A. rotation orders O_alpha: sort x in [0, p) by frac(alpha x) for alpha at every midpoint of consecutive Farey
     fractions of order p - 1 (all distinct such orders; includes the natural and all multiplier orders), each polished
     by first-improvement insertion local search;
  B. simulated annealing on insertion moves (exact O(|k - i|) deltas, Metropolis, geometric cooling), many restarts
     from random orders, in parallel, each finished by insertion local search and window-DP polishing;
  C. window-DP polishing: exact reordering of every window of W consecutive positions (subset DP), alternated with
     insertion local search until neither improves.
Every reported excess is recomputed from scratch.
Usage: python audit_search.py P1 [P2 ...] [--restarts R] [--steps S] [--window W] [--topk K]
"""
import sys, os, math, json, time
import numpy as np
from numba import njit, prange
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import legendre, h_forms


@njit(cache=True)
def leg_table(p):
    t = np.zeros(p, dtype=np.int64)
    for a in range(1, p):
        t[a] = legendre(a, p)
    return t


@njit(cache=True)
def excess(order, chi, p):
    s = 0
    for i in range(p):
        u = order[i]
        for j in range(i + 1, p):
            s += chi[(order[j] - u) % p]
    return s


@njit(cache=True)
def move(order, i, k):
    v = order[i]
    if k > i:
        for t in range(i, k):
            order[t] = order[t + 1]
    else:
        for t in range(i, k, -1):
            order[t] = order[t - 1]
    order[k] = v


@njit(cache=True)
def insertion_ls(order, chi, p):
    """First-improvement over vertices in random order; for each vertex take its best position. Returns excess."""
    improved = True
    while improved:
        improved = False
        perm = np.random.permutation(p)
        for r in range(p):
            i = 0
            v = perm[r]
            for t in range(p):
                if order[t] == v:
                    i = t
                    break
            bg = 0
            bk = i
            g = 0
            for k in range(i + 1, p):
                g -= 2 * chi[(order[k] - v) % p]
                if g > bg:
                    bg = g; bk = k
            g = 0
            for k in range(i - 1, -1, -1):
                g -= 2 * chi[(v - order[k]) % p]
                if g > bg:
                    bg = g; bk = k
            if bg > 0:
                move(order, i, bk)
                improved = True
    return excess(order, chi, p)


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def window_dp(order, chi, p, W):
    """Exact reordering of windows [s, s+W) for s = 0, W/2, W, ...; returns True if anything improved."""
    changed = False
    s = 0
    dp = np.zeros(1 << W, dtype=np.int32)
    par = np.zeros(1 << W, dtype=np.int32)
    inm = np.zeros(W, dtype=np.int64)
    verts = np.zeros(W, dtype=np.int64)
    while s < p:
        w = min(W, p - s)
        for a in range(w):
            verts[a] = order[s + a]
        for a in range(w):
            m = 0
            for b in range(w):
                if b != a and chi[(verts[a] - verts[b]) % p] == 1:   # verts[b] -> verts[a]
                    m |= 1 << b
            inm[a] = m
        size = 1 << w
        dp[0] = 0
        for S in range(1, size):
            best = -1
            bi = -1
            for a in range(w):
                if (S >> a) & 1:
                    R = S ^ (1 << a)
                    val = dp[R] + popcount(inm[a] & R)
                    if val > best:
                        best = val; bi = a
            dp[S] = best
            par[S] = bi
        # current internal forward count
        cur = 0
        for a in range(w):
            for b in range(a + 1, w):
                if chi[(verts[b] - verts[a]) % p] == 1:
                    cur += 1
        if dp[size - 1] > cur:
            S = size - 1
            pos = w - 1
            while S:
                a = par[S]
                order[s + pos] = verts[a]
                pos -= 1
                S ^= 1 << a
            changed = True
        s += W // 2
    return changed


@njit(cache=True)
def polish(order, chi, p, W):
    e = insertion_ls(order, chi, p)
    while True:
        ch = window_dp(order, chi, p, W)
        e2 = insertion_ls(order, chi, p)
        if not ch and e2 <= e:
            return e2
        e = e2


@njit(cache=True)
def anneal(order, chi, p, steps, T0, T1):
    cur = excess(order, chi, p)
    best = cur
    best_order = order.copy()
    lr = math.log(T1 / T0)
    for st in range(steps):
        T = T0 * math.exp(lr * st / steps)
        i = np.random.randint(p)
        k = np.random.randint(p)
        if k == i:
            continue
        v = order[i]
        d = 0
        if k > i:
            for j in range(i + 1, k + 1):
                d -= 2 * chi[(order[j] - v) % p]
        else:
            for j in range(k, i):
                d -= 2 * chi[(v - order[j]) % p]
        if d >= 0 or np.random.random() < math.exp(d / T):
            move(order, i, k)
            cur += d
            if cur > best:
                best = cur
                best_order[:] = order
    order[:] = best_order
    return best


@njit(parallel=True, cache=True)
def sa_batch(chi, p, R, steps, T0, T1, W, seed):
    res = np.zeros(R, dtype=np.int64)
    orders = np.zeros((R, p), dtype=np.int64)
    for r in prange(R):
        np.random.seed(seed + r)
        o = np.random.permutation(p)
        anneal(o, chi, p, steps, T0, T1)
        e = polish(o, chi, p, W)
        res[r] = e
        orders[r, :] = o
    return res, orders


def farey_midpoints(n):
    """Midpoints of consecutive Farey fractions of order n in (0, 1)."""
    a, b, c, d = 0, 1, 1, n
    mids = []
    while c <= n:
        mids.append((a / b + c / d) / 2)
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
    return mids


@njit(cache=True)
def rotation_scan(chi, p, alphas):
    vals = np.zeros(alphas.shape[0], dtype=np.int64)
    xs = np.arange(p)
    for t in range(alphas.shape[0]):
        al = alphas[t]
        key = np.empty(p)
        for x in range(p):
            y = al * x
            key[x] = y - math.floor(y)
        order = np.argsort(key)
        vals[t] = excess(order, chi, p)
    return vals


@njit(cache=True)
def rotation_order(p, al):
    key = np.empty(p)
    for x in range(p):
        y = al * x
        key[x] = y - math.floor(y)
    return np.argsort(key)


def main():
    args = sys.argv[1:]
    R = 96; steps_per = 3000; W = 14; topk = 60
    primes = []
    i = 0
    while i < len(args):
        if args[i] == '--restarts':
            R = int(args[i + 1]); i += 2
        elif args[i] == '--steps':
            steps_per = int(args[i + 1]); i += 2
        elif args[i] == '--window':
            W = int(args[i + 1]); i += 2
        elif args[i] == '--topk':
            topk = int(args[i + 1]); i += 2
        else:
            primes.append(int(args[i])); i += 1
    out = []
    for p in primes:
        t0 = time.time()
        chi = leg_table(p)
        h = h_forms(-p)
        nat = int(excess(np.arange(p), chi, p))
        assert nat == p * h
        # A. rotation family
        mids = np.array(farey_midpoints(p - 1))
        vals = rotation_scan(chi, p, mids)
        rot_raw = int(vals.max())
        idx = np.argsort(-vals)[:topk]
        np.random.seed(p)
        rot_best = -10 ** 9; rot_order = None
        for t in idx:
            o = rotation_order(p, mids[t]).astype(np.int64)
            e = polish(o, chi, p, W)
            if e > rot_best:
                rot_best = e; rot_order = o.copy()
        # B. simulated annealing batch
        res, orders = sa_batch(chi, p, R, steps_per * p, 3.0, 0.15, W, 1000 * p)
        b = int(np.argmax(res))
        sa_best = int(res[b])
        cand = [(rot_best, rot_order), (sa_best, orders[b].copy())]
        best, border = max(cand, key=lambda z: z[0])
        e_check = int(excess(np.array(border, dtype=np.int64), chi, p))
        assert e_check == best and sorted(border.tolist()) == list(range(p))
        row = dict(p=p, h=h, L=round(math.pi * h / math.sqrt(p), 4), natural=nat, rotation_raw_max=rot_raw,
                   rotation_polished=rot_best, sa_best=sa_best, sa_median=int(np.median(res)), restarts=R,
                   steps=steps_per * p, best=best, beaten=bool(best > nat), seconds=round(time.time() - t0, 1),
                   best_order=[int(x) for x in border])
        out.append(row)
        print({k: v for k, v in row.items() if k != 'best_order'}, flush=True)
    tag = '_'.join(str(p) for p in primes)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_search_%s.json' % tag[:60]), 'w'), indent=1)


if __name__ == '__main__':
    main()
