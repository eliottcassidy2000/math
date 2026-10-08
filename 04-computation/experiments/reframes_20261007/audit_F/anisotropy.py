#!/usr/bin/env python3
"""audit_F (C3): does the adaptive coupling rule 'conspire'?  Lamperti-type test of the debt walk along real orbits.

The debt walk M_n (exponent vector x_n in Z^rho) is a martingale whose step law is chosen by the coupling state
(a, b) = (M mod d, e mod d).  If the frequencies of the coupling states (hence the step covariance) depend on the DIRECTION
of x_n, a zero-drift walk can change type (Georgiou-Menshikov-Mijatovic-Wade 2016 in d >= 2; Peres-Popov-Sousi 2013 in
d = 3).  The decisive quantity for a walk with bounded mean-zero steps is the effective dimension

        d_eff = E|xi_w|^2 / E[(xi_w . u_w)^2]      (time average over steps taken far from the origin),

where xi_w = W xi is the whitened step (W = Sigma^(-1/2), Sigma the time-averaged step covariance far from the origin),
and u_w = W x / |W x| the whitened radial direction.  Lamperti: d_eff > 2 + eps transient, d_eff < 2 - eps recurrent
(along the walk's own occupation measure).  An unconspired walk has d_eff = rho.

We run direct big-integer orbits of y, y+1 (as in bigint_pairs.py), keep only steps of unmerged chains with whitened
radius >= R0, and report d_eff overall and split by sign of log M (M > 1 vs M < 1) and by the half-space of the first
coordinate; also the frequency of the identity coupling (a, b) = (1, 0) in each half.
Two passes: pass 1 estimates Sigma, pass 2 measures d_eff with a per-chain jackknife-free standard error (chains are
independent; the SE is computed from per-chain sums by the delta method).
Usage: python3 anisotropy.py NAME NSAMP TMAX SEED [R0]
"""
import sys, math, random
sys.path.insert(0, '.')
from bigint_pairs import Map, parse_map


def chol_inv_sqrt(S):
    """return W with W S W^T = I (inverse Cholesky factor), S symmetric positive definite (small dim)"""
    n = len(S)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = S[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    # invert lower-triangular L
    Li = [[0.0] * n for _ in range(n)]
    for i in range(n):
        Li[i][i] = 1 / L[i][i]
        for j in range(i):
            Li[i][j] = -sum(L[i][k] * Li[k][j] for k in range(j, i)) / L[i][i]
    return Li


def matvec(A, v):
    return [sum(a * b for a, b in zip(row, v)) for row in A]


def walk_paths(mp, nsamp, T, seed, B=128):
    """generator over chains; each chain is a generator of (x_before, step, coupling_state) up to (about) the meeting"""
    d = mp.d
    rnd = random.Random(seed)
    N = T + 64
    dB = d ** B
    k = len(mp.primes)
    ratio_mod = [[(mp.ms[i] * pow(mp.ms[j], -1, d)) % d for j in range(d)] for i in range(d)]
    def chain(y):
        u, v = y + 1, y
        x = [0] * k
        Mmod = 1
        t = 0
        done = False
        while t < T and not done:
            un, bu = mp.block(u, B, dB)
            vn, bv = mp.block(v, B, dB)
            done = (un == vn)          # after a meeting all steps are zero steps at a fixed debt
            for s in range(B):
                i, j = bu[s], bv[s]
                b = (i - Mmod * j) % d
                step = mp.diff[i][j]
                yield tuple(x), step, (Mmod, b)
                for c in range(k):
                    x[c] += step[c]
                Mmod = Mmod * ratio_mod[i][j] % d
            u, v = un, vn
            t += B
    for _ in range(nsamp):
        y = d ** N + rnd.randrange(d ** N)
        yield chain(y)


def main():
    name, nsamp, T, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    R0 = float(sys.argv[5]) if len(sys.argv) > 5 else 6.0
    d, br = parse_map(name)
    mp = Map(d, br)
    k = len(mp.primes)
    # pass 1: covariance of steps (all steps of unmerged chains, |x| >= 3 in sup norm)
    S = [[0.0] * k for _ in range(k)]; n1 = 0
    for rec in walk_paths(mp, nsamp, T, seed):
        for x, st, cs in rec:
            if max(abs(c) for c in x) >= 3:
                n1 += 1
                for a in range(k):
                    for b in range(k):
                        S[a][b] += st[a] * st[b]
    S = [[S[a][b] / n1 for b in range(k)] for a in range(k)]
    W = chol_inv_sqrt(S)
    # pass 2: effective dimension
    def acc():
        return {'tot': 0.0, 'rad': 0.0, 'n': 0, 'id': 0, 'per': []}
    groups = {'all': acc(), 'M>1': acc(), 'M<1': acc(), 'x0>0': acc(), 'x0<0': acc()}
    logp = [math.log(p) for p in mp.primes]
    for rec in walk_paths(mp, nsamp, T, seed):
        loc = {g: [0.0, 0.0] for g in groups}
        for x, st, cs in rec:
            xw = matvec(W, x)
            r = math.sqrt(sum(c * c for c in xw))
            if r < R0:
                continue
            sw = matvec(W, st)
            tot = sum(c * c for c in sw)
            rad = (sum(a * b for a, b in zip(sw, xw)) / r) ** 2
            lm = sum(a * b for a, b in zip(x, logp))
            gs = ['all', 'M>1' if lm > 0 else 'M<1', 'x0>0' if x[0] > 0 else 'x0<0']
            for g in gs:
                G = groups[g]
                G['tot'] += tot; G['rad'] += rad; G['n'] += 1
                G['id'] += (cs == (1, 0))
                loc[g][0] += tot; loc[g][1] += rad
        for g in groups:
            if loc[g][1] > 0 or loc[g][0] > 0:
                groups[g]['per'].append(tuple(loc[g]))
    print(f"{name}: d={d} m={mp.ms} r={mp.rs} rank={k if k else 0} Lambda={mp.lam:+.4f}; {nsamp} chains, T={T}, seed {seed}; "
          f"whitening from {n1} steps; R0={R0}")
    print(f"   step covariance (lattice coords) = {[[round(v, 4) for v in row] for row in S]}")
    for g, G in groups.items():
        if G['n'] == 0 or G['rad'] == 0:
            print(f"   {g:5s}: no steps"); continue
        de = G['tot'] / G['rad']
        # delta-method SE over chains: ratio of sums
        per = G['per']; m = len(per)
        A = sum(a for a, b in per) / m; Bm = sum(b for a, b in per) / m
        var = sum((a - de * b) ** 2 for a, b in per) / (m - 1) / m / Bm ** 2 if m > 1 else float('nan')
        print(f"   {g:5s}: steps {G['n']:9d}  d_eff = {de:6.3f} +- {math.sqrt(var):.3f}   identity-coupling share {G['id']/G['n']:.4f}"
              f"   E|xi_w|^2 = {G['tot']/G['n']:.4f}")


if __name__ == '__main__':
    main()
