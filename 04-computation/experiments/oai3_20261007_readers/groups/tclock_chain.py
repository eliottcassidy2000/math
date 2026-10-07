"""T-clock relation chain for two 2-adic Collatz orbits (lane 'groups', session mac-mini-2026-10-07-oaimath3).

Setting.  T(z) = z/2 (z even), (3z+1)/2 (z odd) on Z_2.  Two orbits related by x = 3^j y + c with c in Z[1/3]
(an element of the affine group Z[1/6] x| <2,3> with no power of 2 in the multiplier).  Advance BOTH orbits by one
T-step.  With p = y mod 2 and e = c mod 2 (c is a 2-adic integer), x mod 2 = p XOR e, and the new relation is
x' = 3^j' y' + c' with
    (a) p=0, e=0:  j' = j,     c' = c/2
    (b) p=0, e=1:  j' = j + 1, c' = (3c + 1)/2
    (c) p=1, e=0:  j' = j,     c' = (3c + 1 - 3^j)/2
    (d) p=1, e=1:  j' = j - 1, c' = (c - 3^(j-1))/2
Representation: N = c * 3^d with d = max(0, -j); N is an integer (invariant).
For Haar y the parities p_0, p_1, ... are i.i.d. fair coins (Terras / Lagarias), so (j, c) is an exact Markov chain
driven by fair coins.  j is a martingale; its nonzero increments are i.i.d. +-1 (optional skipping).

Part 1 (FINITE-EXACT): check the four transitions against direct 2-adic iteration of x and y.
Part 2 (NUMERICAL, exact sampling): Monte Carlo of the merge time for
    M1:      x = y + 1                      start (0, 1), no forced prefix
    S19 lag-1 Mersenne pair: p = 3 q + 2  start (1, 2), q in 5 + 16 Z_2 (forced parities 1,0,0,0), then fair coins
             (S19's lag-1 non-merge share q1(T), T = template total = T-steps of p)
Usage: python3 tclock_chain.py [NSAMP SMAX SEED]   defaults 20000 200000 7
"""
import sys, random, math, time
from collections import Counter

FAIL = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAIL.append(msg)


def step(j, N, p):
    """one T-step of the relation chain; N = c * 3^max(0,-j).  Returns (j', N')."""
    e = N & 1
    if p == 0 and e == 0:
        return j, N >> 1                       # (a)  N even, exact
    if p == 0 and e == 1:                      # (b)
        if j >= 0:
            return j + 1, (3 * N + 1) >> 1
        d = -j
        return j + 1, (N + 3 ** (d - 1)) >> 1
    if p == 1 and e == 0:                      # (c)
        if j >= 0:
            return j, (3 * N + 1 - 3 ** j) >> 1
        d = -j
        return j, (3 * N + 3 ** d - 1) >> 1
    # p == 1 and e == 1                         # (d)
    if j >= 1:
        return j - 1, (N - 3 ** (j - 1)) >> 1
    return j - 1, (3 * N - 1) >> 1


def T(z):
    return z >> 1 if z & 1 == 0 else (3 * z + 1) >> 1


# ---------------------------------------------------------------- Part 1: exact transition check
def part1(seed=1, trials=300, K=400, S=300):
    print('Part 1 (FINITE-EXACT): relation chain vs direct 2-adic iteration (mod 2^(K-s))')
    rng = random.Random(seed)
    nsteps = 0
    ok_all = True
    for tr in range(trials):
        j0 = rng.randint(-6, 6)
        d0 = max(0, -j0)
        N0 = rng.randint(-10 ** 6, 10 ** 6)
        mod = 1 << K
        inv3 = pow(3, -1, mod)
        c_mod = (N0 * pow(inv3, d0, mod)) % mod
        y = rng.getrandbits(K)
        three_j = pow(3, j0, mod) if j0 >= 0 else pow(inv3, -j0, mod)
        x = (three_j * y + c_mod) % mod
        j, N = j0, N0
        for s in range(S):
            prec = K - s
            m = 1 << prec
            inv3m = pow(3, -1, m)
            cm = (N * pow(inv3m, max(0, -j), m)) % m
            tj = pow(3, j, m) if j >= 0 else pow(inv3m, -j, m)
            if (x - (tj * y + cm)) % m != 0:
                ok_all = False
                break
            if (x & 1) != ((y & 1) ^ (N & 1)):
                ok_all = False
                break
            p = y & 1
            j, N = step(j, N, p)
            x, y = T(x), T(y)
            nsteps += 1
    check(ok_all, 'x_s = 3^(j_s) y_s + c_s and x_s mod 2 = (y_s mod 2) XOR (c_s mod 2) at all %d steps of %d random relations'
          % (nsteps, trials))
    # merge precursors
    check(step(1, 1, 1) == (0, 0), '(j,c) = (1,1) with y odd merges: (0,0)')
    check(step(-1, -1, 0) == (0, 0), '(j,c) = (-1,-1/3) [N = -1] with y even merges: (0,0)')
    # (0,0) absorbing
    check(step(0, 0, 0) == (0, 0) and step(0, 0, 1) == (0, 0), '(0,0) is absorbing')
    # M1 classical: y = 4 mod 8 merges with y+1 in 3 steps
    st = (0, 1)
    for p in (0, 0, 1):
        st = step(st[0], st[1], p)
    check(st == (0, 0), 'x = y+1 with y = 4 mod 8 (parities 0,0,1) merges at T-step 3')


# ---------------------------------------------------------------- Part 2: Monte Carlo (exact sampling)
def run_chain(rng, j, N, forced, smax):
    """returns (merge_time or None, number of disagreement steps, j-path statistics)"""
    s = 0
    ndis = 0
    for p in forced:
        if j == 0 and N == 0:
            return s, ndis
        ndis += N & 1
        j, N = step(j, N, p)
        s += 1
    getbit = rng.getrandbits
    while s < smax:
        if j == 0 and N == 0:
            return s, ndis
        ndis += N & 1
        j, N = step(j, N, getbit(1))
        s += 1
    if j == 0 and N == 0:
        return s, ndis
    return None, ndis


def part2(nsamp, smax, seed):
    print('\nPart 2 (NUMERICAL, exact sampling of the chain): merge-time laws, %d samples each, horizon %d T-steps'
          % (nsamp, smax))
    grid = [27, 100, 400, 1600, 3200, 6400, 12800, 19000, 50000, 100000, 200000]
    grid = [g for g in grid if g <= smax]
    s19 = {27: 0.859, 100: 0.722, 400: 0.548, 1600: 0.356, 3200: 0.276, 6400: 0.202, 12800: 0.143, 19000: 0.118}
    cases = [('M1  x = y + 1', 0, 1, ()), ('S19 lag-1  p = 3q + 2, q in 5+16Z_2', 1, 2, (1, 0, 0, 0))]
    out = {}
    for name, j0, N0, forced in cases:
        rng = random.Random(seed * 7919 + j0 * 31 + N0)
        t0 = time.time()
        times = []
        dis_rate = []
        for _ in range(nsamp):
            mt, nd = run_chain(rng, j0, N0, forced, smax)
            times.append(mt)
        alive = lambda T: sum(1 for t in times if t is None or t > T) / nsamp
        print('  %s   (%.0f s)' % (name, time.time() - t0))
        print('    T        q(T) = P(no merge by T)    sqrt(T) q(T)     S19 q1(T)')
        for g in grid:
            q = alive(g)
            se = math.sqrt(q * (1 - q) / nsamp)
            extra = ('   %.3f' % s19[g]) if (name.startswith('S19') and g in s19) else ''
            print('    %-7d  %.4f +- %.4f         %6.2f%s' % (g, q, se, math.sqrt(g) * q, extra))
        # tail exponent on [smax/64, smax]
        xs, ys = [], []
        for g in [smax // 64, smax // 32, smax // 16, smax // 8, smax // 4, smax // 2, smax]:
            q = alive(g)
            if q > 0:
                xs.append(math.log(g)); ys.append(math.log(q))
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
        print('    tail exponent (least squares, T in [%d, %d]): %.3f' % (smax // 64, smax, -slope))
        merged = [t for t in times if t is not None]
        merged.sort()
        if merged:
            print('    merged %d/%d; median merge time %d' % (len(merged), nsamp, merged[len(merged) // 2]))
        out[name] = times
    return out


def part3(seed, nsamp=4000, smax=20000):
    """the disagreement clock and the exact variance identity Var(j_s - j_0) = E[N_s] (on unmerged paths frozen)"""
    print('\nPart 3 (NUMERICAL): disagreement clock N_s = #{r < s : c_r odd} and Var(j_s) = E[N_s] (exact identity)')
    rng = random.Random(seed + 99)
    marks = [100, 1000, 5000, 20000]
    js = {m: [] for m in marks}
    ns = {m: [] for m in marks}
    for _ in range(nsamp):
        j, N = 0, 1
        nd = 0
        s = 0
        for m in marks:
            while s < m:
                if not (j == 0 and N == 0):
                    nd += N & 1
                    j, N = step(j, N, rng.getrandbits(1))
                s += 1
            js[m].append(j); ns[m].append(nd)
    print('    s        E[N_s]/s    Var(j_s)/E[N_s]   E[j_s]')
    for m in marks:
        en = sum(ns[m]) / nsamp
        mj = sum(js[m]) / nsamp
        vj = sum((v - mj) ** 2 for v in js[m]) / nsamp
        print('    %-7d  %.4f      %.4f            %+.3f' % (m, en / m, vj / en if en else float('nan'), mj))
    # disagreement density on unmerged paths only (the clock while alive)
    rng = random.Random(seed + 7)
    alive_steps = 0; dis = 0
    for _ in range(2000):
        j, N = 0, 1
        for s in range(20000):
            if j == 0 and N == 0:
                break
            alive_steps += 1; dis += N & 1
            j, N = step(j, N, rng.getrandbits(1))
    print('    disagreement density while unmerged (M1, 2000 paths x 20000 steps): %.4f' % (dis / alive_steps))


if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:4]]
    NSAMP, SMAX, SEED = (args + [20000, 200000, 7][len(args):])[:3]
    part1()
    part2(NSAMP, SMAX, SEED)
    part3(SEED)
    print('\nALL CHECKS PASSED' if not FAIL else '\nFAILURES: %d' % len(FAIL))
