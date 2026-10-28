"""Local solvability and local densities for n = a^4 + b^4 + c^4 + d^4 over Z[w] (audit of HYP-9429 / THM-4826 (6)).

Independent of the author's code: residue sets by a meet-in-the-middle on pair sums (not iterated set addition);
densities by exact integer pair-count convolution (np.bincount), not FFT; every prime p <= 211 checked with the
Hensel criterion; Weil's bound for the rest.
"""
import math
import sys
import time

import numpy as np

T0 = time.time()
FAIL = []


def check(c, m):
    if not c:
        FAIL.append(m)
        print('  FAIL:', m)


def fourth_counts(M):
    """D[x, y] = #{alpha mod M : alpha^4 = x + y w (mod M)} as an int64 M x M array."""
    a = np.arange(M, dtype=np.int64)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()
    s, t = (A * A - B * B) % M, (2 * A * B - B * B) % M
    x, y = (s * s - t * t) % M, (2 * s * t - t * t) % M
    D = np.bincount(x * M + y, minlength=M * M).astype(np.int64)
    return D.reshape(M, M), (A % 2 == 0) & (B % 2 == 0), x, y


def pair_conv(D, M, chunk=4_000_000):
    """C2 = D * D (cyclic 2-d convolution mod M), exact, via the support (np.bincount with integer weights)."""
    idx = np.flatnonzero(D.ravel())
    w = D.ravel()[idx]
    sx, sy = idx // M, idx % M
    C = np.zeros(M * M, dtype=np.int64)
    n = len(idx)
    step = max(1, chunk // n)
    for i0 in range(0, n, step):
        i1 = min(n, i0 + step)
        X = (sx[i0:i1, None] + sx[None, :]) % M
        Y = (sy[i0:i1, None] + sy[None, :]) % M
        W = w[i0:i1, None] * w[None, :]
        # float64 weights are exact here: every partial sum is an integer below M^4 <= 2^36 < 2^53
        C += np.bincount((X * M + Y).ravel(), weights=W.ravel().astype(np.float64),
                         minlength=M * M).round().astype(np.int64)
    assert int(C.sum()) == M ** 4, (M, int(C.sum()))
    return C.reshape(M, M), len(idx)


def rational_counts(C2, M, ns):
    """N(n) = #{4-tuples mod M with sum = (n, 0)} = sum_P C2(P) C2((n,0) - P), exact integers (Python ints)."""
    out = {}
    flat = C2.astype(object)
    for n in ns:
        # shift: Q = (n - Px, -Py) mod M
        Qx = (n - np.arange(M)) % M
        Qy = (-np.arange(M)) % M
        S = C2[np.ix_(Qx, Qy)]
        out[n] = int(np.sum(C2.astype(np.int64) * S.astype(np.int64), dtype=np.int64)) if M <= 256 else \
            int(sum(int(v) for v in (C2.astype(object) * S.astype(object)).ravel()))
    return out


def residue_obstructions(M, nmax=200):
    D, _, _, _ = fourth_counts(M)
    S1 = (D > 0)
    # pair-sum support (boolean convolution via the support)
    idx = np.flatnonzero(S1.ravel())
    sx, sy = idx // M, idx % M
    P = np.zeros(M * M, dtype=bool)
    step = max(1, 2_000_000 // len(idx))
    for i0 in range(0, len(idx), step):
        X = (sx[i0:i0 + step, None] + sx[None, :]) % M
        Y = (sy[i0:i0 + step, None] + sy[None, :]) % M
        P[(X * M + Y).ravel()] = True
    P = P.reshape(M, M)
    reach = []
    for n in range(M):
        Qx = (n - np.arange(M)) % M
        Qy = (-np.arange(M)) % M
        if np.any(P & P[np.ix_(Qx, Qy)]):
            reach.append(n)
    reach = set(reach)
    bad = [n for n in range(nmax + 1) if n % M not in reach]
    return int(S1.sum()), reach, bad


def section_residues():
    print('== 1. rational residues reachable by four fourth powers (meet-in-the-middle on pair sums)')
    for M in (16, 32, 64, 128, 256, 512, 9, 27, 81, 243, 25, 125, 49, 343, 13, 169, 11, 121):
        ns, reach, bad = residue_obstructions(M)
        b8 = [n for n in range(201) if n % 8 == 5]
        if M % 2 == 0:
            check(bad == b8, f'M={M}')
            tag = 'exactly n = 5 mod 8' if bad == b8 else str(bad)
        else:
            check(bad == [], f'M={M}')
            tag = 'none' if not bad else str(bad)
        print(f'  M = {M:3d}: {ns:6d} fourth-power residues; obstructed n in 0..200: {tag}')
    # negative control: a wrong residue set (pretend fourth powers mod 16 were {0, 1} only, as over Z)
    print('  negative control: with the integer residues {0,1} mod 16 the obstructed set in 0..20 would be',
          [n for n in range(21) if n % 16 not in {a + b + c + d for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)}])


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def section_hensel(PMAX=211):
    """Local solvability of rational n (0 < n <= 200) at every prime p <= PMAX with the Hensel criterion; Weil above."""
    print(f'== 2. local solvability at every prime p <= {PMAX} (Hensel), Weil above')
    ns = list(range(1, 201))
    # p = 2: a solution mod 2^5 with a unit coordinate lifts (v(4 a^3) = 2 for a unit a, so precision 2*2 + 1 = 5)
    M = 32
    a = np.arange(M)
    A, B = np.meshgrid(a, a, indexing='ij')
    A, B = A.ravel(), B.ravel()
    s, t = (A * A - B * B) % M, (2 * A * B - B * B) % M
    X, Y = (s * s - t * t) % M, (2 * s * t - t * t) % M
    unit = ~((A % 2 == 0) & (B % 2 == 0))
    allres = set(zip(X.tolist(), Y.tolist()))
    unitres = set(zip(X[unit].tolist(), Y[unit].tolist()))
    pairs = {((x1 + x2) % M, (y1 + y2) % M) for (x1, y1) in allres for (x2, y2) in allres}
    upairs = {((x1 + x2) % M, (y1 + y2) % M) for (x1, y1) in unitres for (x2, y2) in allres}
    sol2 = {n for n in ns if any(((n - px) % M, (-py) % M) in pairs for (px, py) in upairs)}
    bad2 = [n for n in ns if n not in sol2]
    print(f'  p = 2: n with a unit-coordinate solution mod 32 (hence 2-adically representable): all except '
          f'{"exactly n = 5 mod 8" if bad2 == [n for n in ns if n % 8 == 5] else bad2}')
    check(bad2 == [n for n in ns if n % 8 == 5], 'p=2 Hensel')
    print('  p = 3: x^4 = 1 mod (1 - w) for units, so any n has a unit-coordinate solution mod (1 - w) (take n mod 3 '
          'or 3 units), and 4x^3 is a unit there: no 3-adic obstruction')
    worst = []

    def sumset(Aset, Bpts, shape):
        """Boolean sumset A + {b : b in Bpts} in the additive group (Z/p)^k given by shape; early exit if full."""
        out = np.zeros(shape, dtype=bool)
        for b in Bpts:
            out |= np.roll(Aset, shift=b, axis=tuple(range(len(shape))))
            if out.all():
                break
        return out
    for p in primes_upto(PMAX):
        if p in (2, 3):
            continue
        if p % 3 == 1:
            # split: Z[w]_p = Z_p x Z_p; n rational -> need n = sum of 4 fourth powers in F_p with a nonzero term
            shape = (p,)
            q4 = sorted({pow(x, 4, p) for x in range(p)})
            Q = [(v,) for v in q4]
            NZ = [(v,) for v in q4 if v]
            zero = np.zeros(shape, dtype=bool)
            zero[0] = True
        else:
            # inert: F_{p^2} = F_p[w]; rational n -> (n mod p, 0)
            shape = (p, p)
            a = np.arange(p, dtype=np.int64)
            A, B = np.meshgrid(a, a, indexing='ij')
            s, t = (A * A - B * B) % p, (2 * A * B - B * B) % p
            X, Y = (s * s - t * t) % p, (2 * s * t - t * t) % p
            q4 = sorted(set(zip(X.ravel().tolist(), Y.ravel().tolist())))
            Q = q4
            NZ = [v for v in q4 if v != (0, 0)]
            zero = np.zeros(shape, dtype=bool)
            zero[0, 0] = True
        s3 = zero
        for _ in range(3):
            s3 = sumset(s3, Q, shape)
        ok = sumset(s3, NZ, shape)
        if p % 3 == 1:
            bad = [n for n in ns if not ok[n % p]]
        else:
            bad = [n for n in ns if not ok[n % p, 0]]
        check(bad == [], f'p={p} local obstruction {bad}')
        if bad:
            worst.append((p, bad))
    print(f'  every n in 1..200 has a nonsingular solution mod p (hence is p-adically representable) for all primes '
          f'5 <= p <= {PMAX}: {not worst}')
    print('  p > 211: p does not divide n <= 200, the residue field has q >= p > 81^(2/3) = 18.7 elements, and Weil '
          '(N >= q^3 - 81 q^(3/2) > 0 for c != 0) gives a solution, which lifts (4x^3 a unit). [For p | n one would use '
          'x1^4 + x2^4 + x3^4 = -1, N >= q^2 - 27 q > 0 for q > 27; not needed above 211.]')
    print('  archimedean place: Z[w] (x) R = C, no condition')


def density_table():
    print('== 3. local densities (exact integer counts; pair-support convolution, no FFT)')
    ns = list(range(1, 201))
    res = {}
    for M in (32, 64, 128, 256, 512):
        D, _, _, _ = fourth_counts(M)
        C2, sup = pair_conv(D, M)
        N = rational_counts(C2, M, ns)
        res[M] = {n: N[n] / float(M) ** 6 for n in ns}
        print(f'  2-adic modulus {M}: support {sup}; done [{time.time() - T0:.0f} s]')
        sys.stdout.flush()
    res3 = {}
    for M in (9, 27, 81, 243):
        D, _, _, _ = fourth_counts(M)
        C2, sup = pair_conv(D, M)
        N = rational_counts(C2, M, ns)
        res3[M] = {n: N[n] / float(M) ** 6 for n in ns}
        print(f'  3-adic modulus {M}: support {sup}; done [{time.time() - T0:.0f} s]')
        sys.stdout.flush()
    # stability
    unstable2 = [n for n in ns if abs(res[128][n] - res[512][n]) > 1e-9]
    unstable2b = [n for n in ns if abs(res[256][n] - res[512][n]) > 1e-9]
    unstable3 = [n for n in ns if abs(res3[81][n] - res3[243][n]) > 1e-9]
    print(f'  2-adic: n <= 200 whose density at 2^7 differs from 2^9: {unstable2}')
    print(f'          n <= 200 whose density at 2^8 differs from 2^9: {unstable2b}')
    print(f'  3-adic: n <= 200 whose density at 3^4 differs from 3^5: {unstable3}')
    for n in unstable2:
        print(f'    n = {n}: sigma_2 at 32/64/128/256/512 = ' + '/'.join(f'{res[M][n]:.4f}' for M in (32, 64, 128, 256, 512)))
    allowed = [n for n in ns if n % 8 != 5]
    for label, s2 in (('author modulus 2^7', res[128]), ('stable modulus 2^9', res[512])):
        prod = {n: s2[n] * res3[243][n] for n in allowed}
        order = sorted(allowed, key=lambda n: (prod[n], n))
        rank = {n: i + 1 for i, n in enumerate(order)}
        med = float(np.median([prod[n] for n in allowed]))
        print(f'  [{label}] median {med:.4f}, minimum {prod[order[0]]:.4f}; misses/fills:')
        for n in (124, 140, 60, 6, 78, 75, 54, 102, 186):
            print(f'     n = {n:3d}: sigma_2 = {s2[n]:.4f}, sigma_3 = {res3[243][n]:.4f}, product {prod[n]:.4f}, rank {rank[n]} of {len(allowed)}')
        low = [n for n in allowed if abs(prod[n] - prod[order[0]]) < 1e-9]
        print(f'     lowest class: {low}')
        below = [n for n in (6, 54, 60, 75, 78, 124, 140) if prod[n] < med]
        print(f'     persistent misses below the median: {below}')
    return res, res3


if __name__ == '__main__':
    section_residues()
    section_hensel()
    density_table()
    print('FAILURES:', len(FAIL), FAIL[:10])
    print(f'[{time.time() - T0:.0f} s]')
