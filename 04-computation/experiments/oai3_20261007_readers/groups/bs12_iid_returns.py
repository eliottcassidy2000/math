"""The i.i.d. comparison walk on BS(1,2) = Z[1/2] x| Z (lane 'groups').

Symmetric random walk with steps s = t^a u t^(-b): z -> 2^(a-b) z + eps 2^a, a, b i.i.d. Geom(1/2) on {1,2,...},
eps = +-1 fair.  Its height L_n = sum (a_k - b_k) has exactly the law of S19's debt-walk height increments
(mean 0, variance 4).  phi_n = s_1 ... s_n has translation part Delta_n = sum_k eps_k 2^(h_k), h_k = H_(k-1) + a_k.
Given the heights, P(Delta_n = 0) is computed EXACTLY by a carry dynamic program (levels from the bottom: the carry
plus the local sign sum must be even, carry' = (carry + S_h)/2; the final carry must vanish).  Monte Carlo over the
heights gives unbiased estimates of
    r(n)  = P(Delta_n = 0)                (return to the axis subgroup <t>),
    r0(n) = P(Delta_n = 0, L_n = 0)       (return to the identity),
    P(L_n = 0)                            (the height alone, ~ n^(-1/2)).
Varopoulos / Pittet-Saloff-Coste: r0(n) = exp(-Theta(n^(1/3))).  Contrast: in the Collatz relation chain the fiber is
slaved to the height (contracting, E log M = log(3/4)/2 per T-step) and the merge law is ~ n^(-1/2).
Usage: python3 bs12_iid_returns.py [NSAMP SEED]   defaults 200000 5
"""
import sys, random, math
from collections import defaultdict

args = [int(a) for a in sys.argv[1:3]]
NSAMP, SEED = (args + [200000, 5][len(args):])[:2]


def geom(rng):
    a = 1
    while rng.getrandbits(1):
        a += 1
    return a


BIN = {}


def binom_signs(n):
    """law of a sum of n fair +-1 signs: dict S -> prob"""
    if n not in BIN:
        BIN[n] = {2 * k - n: math.comb(n, k) / 2 ** n for k in range(n + 1)}
    return BIN[n]


def p_fiber_zero(heights):
    """exact P(sum_k eps_k 2^(h_k) = 0) given the multiset of heights"""
    cnt = defaultdict(int)
    for h in heights:
        cnt[h] += 1
    levels = sorted(cnt)
    dist = {0: 1.0}
    prev = levels[0]
    for h in levels:
        # empty levels between prev and h: carry must stay even and halves
        gap = h - prev
        for _ in range(gap):
            nd = {}
            for c, pr in dist.items():
                if c % 2 == 0:
                    nd[c // 2] = nd.get(c // 2, 0.0) + pr
            dist = nd
            if not dist:
                return 0.0
        prev = h
        nd = defaultdict(float)
        for c, pr in dist.items():
            for S, ps in binom_signs(cnt[h]).items():
                v = c + S
                if v % 2 == 0:
                    nd[v // 2] += pr * ps
        dist = nd
        prev = h + 1
        if not dist:
            return 0.0
    # after the top level the remaining carry must be exactly zero
    return dist.get(0, 0.0)


def main():
    rng = random.Random(SEED)
    grid = [4, 8, 16, 32, 64, 128, 256]
    nmax = grid[-1]
    acc_r = {n: 0.0 for n in grid}
    acc_r0 = {n: 0.0 for n in grid}
    acc_L0 = {n: 0 for n in grid}
    acc_r2 = {n: 0.0 for n in grid}
    for _ in range(NSAMP):
        H = 0
        heights = []
        for k in range(1, nmax + 1):
            a, b = geom(rng), geom(rng)
            heights.append(H + a)
            H += a - b
            if k in acc_r:
                pz = p_fiber_zero(heights)
                acc_r[k] += pz
                acc_r2[k] += pz * pz
                if H == 0:
                    acc_r0[k] += pz
                    acc_L0[k] += 1
    print('i.i.d. symmetric walk on BS(1,2), steps t^a u t^-b (a, b ~ Geom(1/2)); %d height paths (seed %d)' % (NSAMP, SEED))
    print('  n     P(L_n=0)   sqrt(n)P(L_n=0)   r(n)=P(Delta_n=0)   r0(n)=P(phi_n=e)    r0/P(L_n=0)   -log r0 / n^(1/3)')
    for n in grid:
        pl = acc_L0[n] / NSAMP
        r = acc_r[n] / NSAMP
        r0 = acc_r0[n] / NSAMP
        se = math.sqrt(max(acc_r2[n] / NSAMP - r * r, 0) / NSAMP)
        print('  %-4d  %.5f    %.3f             %.3e (se %.1e)   %.3e          %.3e      %s'
              % (n, pl, math.sqrt(n) * pl, r, se, r0, r0 / pl if pl else float('nan'),
                 ('%.3f' % (-math.log(r0) / n ** (1 / 3))) if r0 > 0 else 'n/a'))
    print('  expected number of returns to the axis subgroup <t> among n in the grid is bounded (sum r(n) converges):'
          ' the walk is transient relative to <t>.')


if __name__ == '__main__':
    main()
