#!/usr/bin/env python3
"""audit_F (A): coalescence of y and y+1 by DIRECT big-integer orbits (no pair-chain formula, no Fractions).

Map T(x) = (m_i x + r_i)/d on x = i mod d.  y is uniform on [d^N, 2 d^N) with N = T + MARGIN, so y mod d^N is uniform and
(Terras / Matthews-Watts) the branch sequences of y and y+1 up to time T are exactly those of a Haar pair.  Both orbits are
iterated as Python integers; the first n <= T with T^n(y) = T^n(y+1) is the equal-time meeting time.  The debt is computed
only from the two branch sequences: M_n = prod m_{i_k} / prod m_{j_k}, stored as an integer exponent vector over the primes
of the multipliers (M_n = 1 iff the vector is zero).

Speed: the orbits are advanced in blocks of B steps: the branches of the next B steps are read off X mod d^B (a small
integer), then X <- (P X + Q) / d^B with the composed affine map (exact).  A block in which the two orbits become equal is
replayed step by step to find the exact meeting time.  --naive runs the plain one-step-at-a-time iteration (validation).

Recorded per sample: meeting time tau (or None), debt vector at the meeting (a meeting with nonzero debt would be a
non-Haar 'null' coincidence: flagged), visits of the debt to 0 at times 1..min(T, tau), and the coupling state
(a, b) = (M mod d, e mod d) with e mod d = (i - a j) mod d.
Output: q(T) = P(no meeting by T) with a binomial standard error, sqrt(T) q, 1/q; mean debt-zero visits over all samples and
among samples still unmerged at T, at T = 16, 64, ..., TMAX.

Usage: python3 bigint_pairs.py NAME NSAMP TMAX SEED [--naive] [--B 256] [--margin 64] [--offset 1]
       NAME is a key of MAPS below or a spec 'd:m0,r0;m1,r1;...'
"""
import sys, math, random, time, argparse

MAPS = {
    'x+1':      (2, [(1, 0), (1, 1)]),
    '3x+1':     (2, [(1, 0), (3, 1)]),
    '3x+5':     (2, [(1, 0), (3, 5)]),
    '5x+1':     (2, [(1, 0), (5, 1)]),
    'Z3_124':   (3, [(1, 0), (2, 1), (4, 1)]),
    'Z3_125':   (3, [(1, 0), (2, 1), (5, 2)]),
    'Z3_1416':  (3, [(1, 0), (4, 2), (16, 1)]),
    'Z3_157':   (3, [(1, 0), (5, 1), (7, 1)]),
    'Z5_12311': (5, [(1, 0), (2, 3), (3, 4), (1, 2), (1, 1)]),
    'Z5_12471': (5, [(1, 0), (2, 3), (4, 2), (7, 4), (1, 1)]),
    'Z5_12371': (5, [(1, 0), (2, 3), (3, 4), (7, 4), (1, 1)]),
}


def parse_map(name):
    if name in MAPS:
        return MAPS[name]
    dpart, rest = name.split(':')
    d = int(dpart)
    br = [tuple(int(t) for t in b.split(',')) for b in rest.split(';')]
    assert len(br) == d
    return d, br


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


class Map:
    def __init__(self, d, br):
        self.d = d
        self.ms = [m for m, r in br]
        self.rs = [r for m, r in br]
        for i, (m, r) in enumerate(br):
            assert math.gcd(m, d) == 1 and (m * i + r) % d == 0, (i, m, r)
        self.primes = sorted({p for m in self.ms for p in factor(m)})
        self.L = [tuple(factor(m).get(p, 0) for p in self.primes) for m in self.ms]
        # exponent-vector difference table and m mod d (for the coupling state)
        self.diff = [[tuple(a - b for a, b in zip(self.L[i], self.L[j])) for j in range(d)] for i in range(d)]
        self.lam = sum(math.log(m / d) for m in self.ms) / d

    def step(self, x):
        b = x % self.d
        return (self.ms[b] * x + self.rs[b]) // self.d, b

    def block(self, X, B, dB):
        """advance X by B steps; return (new X, list of B branches)"""
        d, ms, rs = self.d, self.ms, self.rs
        x = X % dB
        P, Q, D = 1, 0, 1
        brs = []
        for _ in range(B):
            b = x % d
            brs.append(b)
            m, r = ms[b], rs[b]
            x = (m * x + r) // d
            Q = m * Q + r * D
            P *= m
            D *= d
        num = P * X + Q
        Xn = num // dB
        return Xn, brs


def run_sample(mp, T, N, rnd, B, naive, offset):
    d = mp.d
    y = d ** N + rnd.randrange(d ** N)
    v, u = y, y + offset
    zero = tuple(0 for _ in mp.primes)
    debt = list(zero)
    k = len(debt)
    visits = []          # times t (1..) with debt == 0, up to and including the meeting
    tau = None
    debt_at_meet = None
    t = 0
    dB = d ** B
    while t < T:
        if naive:
            u2, i = mp.step(u)
            v2, j = mp.step(v)
            bu, bv = [i], [j]
            nsteps = 1
            u_new, v_new = u2, v2
        else:
            nsteps = min(B, T - t)
            if nsteps == B:
                u_new, bu = mp.block(u, B, dB)
                v_new, bv = mp.block(v, B, dB)
            else:
                u_new, bu = mp.block(u, nsteps, d ** nsteps)
                v_new, bv = mp.block(v, nsteps, d ** nsteps)
        met_in_block = (u_new == v_new)
        if met_in_block and nsteps > 1:
            # replay step by step to find the first meeting time inside the block
            uu, vv = u, v
            for s in range(nsteps):
                uu, i = mp.step(uu)
                vv, j = mp.step(vv)
                assert i == bu[s] and j == bv[s]
                dv = mp.diff[i][j]
                for c in range(k):
                    debt[c] += dv[c]
                if all(x == 0 for x in debt):
                    visits.append(t + s + 1)
                if uu == vv:
                    tau = t + s + 1
                    debt_at_meet = tuple(debt)
                    break
            assert tau is not None
            return tau, debt_at_meet, visits
        for s in range(nsteps):
            dv = mp.diff[bu[s]][bv[s]]
            for c in range(k):
                debt[c] += dv[c]
            if all(x == 0 for x in debt):
                visits.append(t + s + 1)
        if met_in_block:   # naive mode, nsteps == 1
            tau = t + 1
            debt_at_meet = tuple(debt)
            return tau, debt_at_meet, visits
        u, v = u_new, v_new
        t += nsteps
    return None, None, visits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('nsamp', type=int); ap.add_argument('tmax', type=int)
    ap.add_argument('seed', type=int)
    ap.add_argument('--naive', action='store_true')
    ap.add_argument('--B', type=int, default=0)
    ap.add_argument('--margin', type=int, default=64)
    ap.add_argument('--offset', type=int, default=1)
    ap.add_argument('--dump', default='')
    a = ap.parse_args()
    d, br = parse_map(a.name)
    mp = Map(d, br)
    B = a.B if a.B else max(16, int(256 * math.log(5) / math.log(d)) // 4 * 2)
    N = a.tmax + a.margin
    rnd = random.Random(a.seed)
    t0 = time.time()
    res = []
    for s in range(a.nsamp):
        res.append(run_sample(mp, a.tmax, N, rnd, B, a.naive, a.offset))
    el = time.time() - t0
    nullmeet = sum(1 for tau, dm, vis in res if tau is not None and any(dm))
    print(f"{a.name}: d={d} m={mp.ms} r={mp.rs} Lambda={mp.lam:+.4f} primes={mp.primes}; {a.nsamp} samples, offset {a.offset}, "
          f"N=T+{a.margin}, B={B}{' naive' if a.naive else ''}, seed {a.seed}; {el:.1f}s")
    print(f"   meetings with nonzero debt (non-Haar coincidences): {nullmeet}")
    checkpoints = []
    T = 4
    while T <= a.tmax:
        checkpoints.append(T); T *= 4
    if checkpoints[-1] != a.tmax:
        checkpoints.append(a.tmax)
    for T in checkpoints:
        alive = [r for r in res if r[0] is None or r[0] > T]
        q = len(alive) / a.nsamp
        se = math.sqrt(max(q * (1 - q), 1e-12) / a.nsamp)
        vis_all = sum(sum(1 for t in r[2] if t <= T) for r in res) / a.nsamp
        vis_alive = (sum(sum(1 for t in r[2] if t <= T) for r in alive) / len(alive)) if alive else float('nan')
        inv = 1 / q if q > 0 else float('inf')
        # window diagnostic: visits in (T/4, T] per chain still unmerged at T/4 (rank 1 ~ sqrt T, rank 2 ~ const, rank 3 -> 0)
        T0 = T // 4
        alive0 = [r for r in res if r[0] is None or r[0] > T0]
        win = (sum(sum(1 for t in r[2] if T0 < t <= T) for r in alive0) / len(alive0)) if alive0 else float('nan')
        print(f"  T={T:6d}  q={q:.4f} +- {se:.4f}  sqrt(T)q={math.sqrt(T)*q:7.3f}  1/q={inv:6.3f}  "
              f"visits(all)={vis_all:7.3f}  visits(unmerged)={vis_alive:7.3f}  window(T/4,T]/alive={win:7.3f}  win/sqrtT={win/math.sqrt(T):.3f}",
              flush=True)
    if a.dump:
        with open(a.dump, 'w') as f:
            for tau, dm, vis in res:
                f.write(f"{tau} {'' if dm is None else ','.join(map(str, dm))} {len(vis)} {' '.join(map(str, vis[:50]))}\n")


if __name__ == '__main__':
    main()
