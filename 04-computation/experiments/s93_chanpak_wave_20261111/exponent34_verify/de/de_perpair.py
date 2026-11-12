#!/usr/bin/env python3
"""Independent verifier (part D/E) - the injection of (D) and the lattice-line bound of (C), pair by pair.

For every (w, t', s) with 1 <= w <= WMAX, 1 <= t' <= TMAX, 0 <= s <= SMAX and every odd pair (x, x') of w-bit words:
  * enumerate ALL (P, c) with P in [2^(t'-1), 2^t'), c in (2^s, 2^(s+1)] and count the solutions of
        c X - P Xi = gamma(x, x', P mod 2, c mod 2)
    in each of the four parity classes;
  * check, per parity class:  #solutions <= 1 + floor(g min((2^(t'-1)-1)/|X|, (2^s-1)/|Xi|))   (two solutions differ by j (Xi, X)/g)
    and the sharper parity-aware form with j even;
  * check, per pair:          #solutions <= 4 (1 + 2^max(t'-1, s) g / max(|X|, |Xi|))            (statement (C));
  * check that every solution is a genuine collision x(P|x|Z) = x(P|x'|Z) with Z the suffix of c, that the common
    prefix/suffix of the two words are exactly t' and s, and conversely that non-solutions are not collisions
    (so the map (h, h') -> (x, x', P, c) of (D) is a bijection onto the solution set);
  * return the total number of solutions = E(t', s) at n = t' + w + s, to be compared with the census.
Then: |X|, |Xi| >= 2^(w-1-R) for all odd pairs of windows without a run of more than R equal bits.
"""
import sys, json, time
import numpy as np
from numba import njit


def q_table(nmax):
    N = 1 << nmax
    q = np.zeros(N, dtype=np.int64)
    for k in range(nmax):
        lo, hi = 1 << k, 1 << (k + 1)
        h = np.arange(lo, hi, dtype=np.int64)
        half = h >> 1
        q[lo:hi] = q[half] + (h & 1) * ((half + 1) >> 1)
    return q


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def rev_bits(v, w):
    r = 0
    for i in range(w):
        r = (r << 1) | ((v >> i) & 1)
    return r


@njit(cache=True)
def test_cell(w, tp, s, q):
    """returns (total solutions, violations per-class integer bound, violations parity-aware bound,
                violations statement (C) bound, bijection failures, max over pairs of (#solutions / bound (C)) * 1e6,
                number of pairs with >= 1 solution, max #solutions of one pair)"""
    N = 1 << w
    total = 0
    v_int = 0
    v_par = 0
    v_C = 0
    v_bij = 0
    best = 0.0
    npos = 0
    maxsol = 0
    cnt = np.zeros(4, dtype=np.int64)
    A = 1 << (tp - 1)
    B = 1 << s
    for x in range(N):
        rx = rev_bits(x, w)
        for y in range(N):
            X = x - y
            Xi = rx - rev_bits(y, w)
            if (X & 1) == 0 or (Xi & 1) == 0:
                continue
            aX = abs(X)
            aXi = abs(Xi)
            g = gcd(aX, aXi)
            dq = 2 * (q[x] - q[y])
            dtop = ((x >> (w - 1)) & 1) - ((y >> (w - 1)) & 1)
            dbot = (x & 1) - (y & 1)
            cnt[:] = 0
            for P in range(A, 2 * A):
                for Z in range(B):
                    c = 2 * B - rev_bits(Z, s)
                    gam = dq + (P & 1) * dtop + (c & 1) * dbot
                    sol = (c * X - P * Xi == gam)
                    h1 = (((P << w) | x) << s) | Z
                    h2 = (((P << w) | y) << s) | Z
                    coll = (2 * h1 + 1 - 2 * q[h1] == 2 * h2 + 1 - 2 * q[h2])
                    if sol != coll:
                        v_bij += 1
                    if sol:
                        cnt[2 * (P & 1) + (c & 1)] += 1
            tot = cnt[0] + cnt[1] + cnt[2] + cnt[3]
            total += tot
            if tot > 0:
                npos += 1
            if tot > maxsol:
                maxsol = tot
            f1 = ((A - 1) * g) // aX
            f2 = ((B - 1) * g) // aXi
            bint = 1 + (f1 if f1 < f2 else f2)
            e1 = ((A - 1) * g) // (2 * aX)
            e2 = ((B - 1) * g) // (2 * aXi)
            bpar = 1 + (e1 if e1 < e2 else e2)
            for k in range(4):
                if cnt[k] > bint:
                    v_int += 1
                if cnt[k] > bpar:
                    v_par += 1
            mx = aX if aX > aXi else aXi
            M = A if A > B else B
            # statement (C): tot <= 4 (1 + M g / mx)  <=>  tot * mx <= 4 (mx + M g)
            if tot * mx > 4 * (mx + M * g):
                v_C += 1
            r = tot * mx / (4.0 * (mx + M * g))
            if r > best:
                best = r
    return total, v_int, v_par, v_C, v_bij, best, npos, maxsol


@njit(cache=True)
def max_run_w(x, w):
    best = 1
    cur = 1
    prev = x & 1
    for i in range(1, w):
        b = (x >> i) & 1
        if b == prev:
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
            prev = b
    return best


@njit(cache=True)
def runfree_sizes(w, R):
    """min over odd pairs of windows with no run > R of min(|X|,|Xi|); number of violations of >= 2^(w-1-R)"""
    N = 1 << w
    ok = np.zeros(N, dtype=np.bool_)
    rv = np.zeros(N, dtype=np.int64)
    for x in range(N):
        ok[x] = max_run_w(x, w) <= R
        rv[x] = rev_bits(x, w)
    mn = 1 << 60
    viol = 0
    npairs = 0
    for x in range(N):
        if not ok[x]:
            continue
        for y in range(N):
            if not ok[y]:
                continue
            X = x - y
            Xi = rv[x] - rv[y]
            if (X & 1) == 0 or (Xi & 1) == 0:
                continue
            npairs += 1
            a = abs(X)
            b = abs(Xi)
            m = a if a < b else b
            if m < mn:
                mn = m
            # need m >= 2^(w-1-R); for w-1-R < 0 this is m >= 1
            if w - 1 - R >= 0 and m < (1 << (w - 1 - R)):
                viol += 1
    return mn, viol, npairs


def main():
    WMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    TMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    SMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    col = json.load(open('de_collisions.json'))
    q = q_table(WMAX + TMAX + SMAX + 1)
    t0 = time.time()
    ncell = 0
    bad = 0
    worst = (0.0, None)
    mism = 0
    compared = 0
    for w in range(1, WMAX + 1):
        tot_w = 0
        for tp in range(1, TMAX + 1):
            for s in range(0, SMAX + 1):
                total, v_int, v_par, v_C, v_bij, best, npos, maxsol = test_cell(w, tp, s, q)
                ncell += 1
                tot_w += total
                if v_int or v_par or v_C or v_bij:
                    bad += 1
                    print("VIOLATION w=%d t'=%d s=%d: per-class int %d, parity-aware %d, statement (C) %d, bijection %d"
                          % (w, tp, s, v_int, v_par, v_C, v_bij))
                if best > worst[0]:
                    worst = (best, (w, tp, s))
                n = w + tp + s
                if str(n) in col:
                    c = col[str(n)]['cnt'].get('%d,%d' % (tp, s))
                    E = 2 * sum(c) if c else 0
                    compared += 1
                    if E != total:
                        mism += 1
                        print("MISMATCH with census: n=%d t'=%d s=%d: solutions %d, census %d" % (n, tp, s, total, E))
        print("w=%d: cells t'=1..%d, s=0..%d done, total solutions %d  (%.1fs)" % (w, TMAX, SMAX, tot_w, time.time() - t0))
        sys.stdout.flush()
    print("cells tested: %d ; cells with any violation: %d ; cells compared with the census: %d, mismatches: %d"
          % (ncell, bad, compared, mism))
    print("largest (#solutions of a pair)/(4(1 + 2^max(t'-1,s) g/max(|X|,|Xi|))) = %.4f at (w,t',s) = %s" % worst)
    print()
    print("run-free windows: min over odd pairs of min(|X|,|Xi|) against 2^(w-1-R)")
    for R in range(1, 7):
        line = []
        for w in range(1, 13):
            mn, viol, npairs = runfree_sizes(w, R)
            thr = 2 ** max(w - 1 - R, 0)
            if npairs == 0:
                line.append("w=%d:none" % w)
            else:
                line.append("w=%d:min=%d(thr %d)%s" % (w, mn, thr, '' if viol == 0 else ' VIOLATIONS=%d' % viol))
                if viol:
                    bad += 1
        print("  R=%d: %s" % (R, ' '.join(line)))
    print("TOTAL VIOLATIONS:", bad + mism)


if __name__ == '__main__':
    main()
