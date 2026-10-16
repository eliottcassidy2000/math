#!/usr/bin/env python3
"""S38 audit, item 3: independent check of dense-grid values in [2^34, 2^40) that are NOT on the S37 64-point grid
(5,760 such points rest on a single run).  Anchors: the S37 census values at N = 2^k (class sums mod 243, H, exact
prefix sums of H and H^2; equal to the dense run there, a5), k = 34..39.  From each anchor the next three dense grid
points N_(j+1), N_(j+2), N_(j+3) are reached by summing t(n) exactly with own code (parity table to 2^31, uint64
descent with exact Python fallback), and S, the class sums mod 3, H and the prefix sums of H, H^2 are compared with
dense40.npz.
Usage: python a4c_dense_only_points.py SCRATCH_S38_DIR REPO_ROOT
"""
import sys, os, json, time
import numpy as np
import numba as nb

SD, ROOT = sys.argv[1:3]
HERE = os.path.dirname(os.path.abspath(__file__))
LOWB = 31
GUARD = (2 ** 64 - 2) * 2 // 3 - 1


@nb.njit
def build(P, nmax):
    one = np.uint64(1)
    for n in range(2, nmax):
        if (n & 1) == 0:
            h = n >> 1
            par = 1 - int((P[h >> 6] >> np.uint64(h & 63)) & one)
        else:
            x = np.uint64(n)
            un = np.uint64(n)
            s = 0
            while x >= un:
                if x & one:
                    x = x + (x + one) // np.uint64(2)
                else:
                    x = x >> one
                s += 1
            xi = int(x)
            par = (s + int((P[xi >> 6] >> np.uint64(xi & 63)) & one)) & 1
        if par:
            P[n >> 6] |= one << np.uint64(n & 63)


@nb.njit
def chunk(P, a, b, low, guard, h0, flags):
    # n in [a, b): returns dS, dC0, dC1, dC2, sum hp, sum hp^2, hp_end; flags[i] = 1 where uint64 would overflow
    one = np.uint64(1)
    ul = np.uint64(low)
    dS = 0
    c0 = 0
    c1 = 0
    c2 = 0
    hp = 0
    s1 = 0
    s2 = 0
    r3 = a % 3
    nf = 0
    for n in range(a, b):
        x = np.uint64(n)
        s = 0
        bad = False
        while x >= ul:
            if x & one:
                if x > guard:
                    bad = True
                    break
                x = x + (x + one) // np.uint64(2)
            else:
                x = x >> one
            s += 1
        if bad:
            flags[nf] = n
            nf += 1
            t = 0
        else:
            xi = int(x)
            t = 1 - 2 * ((s + int((P[xi >> 6] >> np.uint64(xi & 63)) & one)) & 1)
        dS += t
        if r3 == 0:
            c0 += t
            hp -= t
        elif r3 == 1:
            c1 += t
            hp -= t
        else:
            c2 += t
            hp += 2 * t
        s1 += hp
        s2 += hp * hp
        r3 += 1
        if r3 == 3:
            r3 = 0
    return dS, c0, c1, c2, s1, s2, hp, nf


def tbig(P, n, low):
    s = 0
    x = n
    while x >= low:
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        s += 1
    return 1 - 2 * ((s + ((int(P[x >> 6]) >> (x & 63)) & 1)) & 1)


t0 = time.time()
P = np.zeros((1 << LOWB) // 64 + 1, np.uint64)
build(P, 1 << LOWB)
print(f'table to 2^{LOWB}: {time.time() - t0:.1f} s', flush=True)
c37 = np.load(os.path.join(ROOT, '05-knowledge/results/collatz_census_extension_2p40_20261015.npz'))
g37 = np.array(c37['grid'])
d = np.load(os.path.join(SD, 'dense40.npz'))
g = np.array(d['grid'])
res = []
bad_total = 0
for k in range(34, 40):
    N0 = 1 << k
    i37 = int(np.nonzero(g37 == N0)[0][0])
    cs = np.array(c37['classsums'][i37]).astype(object)
    C = [int(sum(cs[a::3])) for a in range(3)]
    S = sum(C)
    H = 3 * C[2] - S
    assert H == int(c37['H'][i37])
    p1 = int(c37['preH1_str'][i37])
    p2 = int(c37['preH2_str'][i37])
    jd = int(np.nonzero(g == N0)[0][0])
    pos = N0 + 1
    for step in (1, 2, 3):
        target = int(g[jd + step])
        CH = 1 << 22
        SUBCH = 1 << 12
        while pos <= target:
            b = min(pos + CH, target + 1)
            flags = np.zeros(4096, np.int64)
            dS, c0, c1, c2, s1, s2, hp, nf = chunk(P, pos, b, 1 << LOWB, np.uint64(GUARD), H, flags)
            if nf:
                # split the chunk into sub-chunks of 2^12; sub-chunks containing a flagged n are summed exactly in Python
                pieces = []
                q = pos
                while q < b:
                    e = min(q + SUBCH, b)
                    fl2 = np.zeros(SUBCH, np.int64)
                    r_ = chunk(P, q, e, 1 << LOWB, np.uint64(GUARD), 0, fl2)
                    if r_[7]:
                        dS_ = c0_ = c1_ = c2_ = s1_ = s2_ = hp_ = 0
                        for n in range(q, e):
                            t = tbig(P, n, 1 << LOWB)
                            dS_ += t
                            rr = n % 3
                            if rr == 0:
                                c0_ += t; hp_ -= t
                            elif rr == 1:
                                c1_ += t; hp_ -= t
                            else:
                                c2_ += t; hp_ += 2 * t
                            s1_ += hp_
                            s2_ += hp_ * hp_
                        r_ = (dS_, c0_, c1_, c2_, s1_, s2_, hp_, 0)
                    pieces.append((e - q,) + tuple(int(x) for x in r_[:7]))
                    q = e
                # combine the pieces (prefix sums of hp restart at 0 in each piece)
                dS = c0 = c1 = c2 = s1 = s2 = 0
                hp = 0
                for (ln_, a0, a1, a2, a3, a4, a5, a6) in pieces:
                    s1 += ln_ * hp + a4
                    s2 += ln_ * hp * hp + 2 * hp * a4 + a5
                    dS += a0; c0 += a1; c1 += a2; c2 += a3
                    hp += a6
            ln = b - pos
            p1 += ln * H + int(s1)
            p2 += ln * H * H + 2 * H * int(s1) + int(s2)
            S += int(dS)
            C = [C[0] + int(c0), C[1] + int(c1), C[2] + int(c2)]
            H += int(hp)
            pos = b
        j = jd + step
        ok = (S == int(d['S'][j]) and C == [int(x) for x in d['C3'][j]] and H == int(d['H'][j])
              and p1 == int(d['preH1_str'][j]) and p2 == int(d['preH2_str'][j]))
        bad_total += (not ok)
        res.append(dict(N=target, S=S, H=H, agrees=ok))
        print(f'   2^{k} + {step} dense steps: N = {target}, S = {S}, H = {H}, agrees with dense40: {ok}  ({time.time() - t0:.0f} s)', flush=True)
out = dict(points=res, mismatches=bad_total, seconds=time.time() - t0)
json.dump(out, open(os.path.join(HERE, 'a4c_dense_only_points.json'), 'w'), indent=1)
print(f'{len(res)} dense-only grid points in [2^34, 2^40) recomputed from S37 anchors: {bad_total} mismatches', flush=True)
