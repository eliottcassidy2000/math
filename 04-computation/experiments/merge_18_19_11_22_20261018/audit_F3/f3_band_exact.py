#!/usr/bin/env python3
"""Audit F3, items (d) and (e): exact ideal (N -> infinity) band computations, written from the definitions of
THM-4690 (sources, seeds, transport, literal clock), without importing the session's scripts.

Ideal source coefficient.  A seed s = r (mod 3) at scale s/N in (2^m, 2^(m+1)] (representative 3*2^e, e = m - 1)
carries h(r, e) = sum over its B-chain elements 4^k s:  kappa(r + k) while 4^k s <= 2N  (e + 2k <= -1),
lambda(r + k) when 2N < 4^k s <= 4N (e + 2k = 0), with kappa = (0, -3, 3), lambda = (-1, -1, 2).  Sources need e <= 0.

Preimages of a seed t (THM-4690(2)): for j >= 0 with B^j(t) = 1 (mod 3) [a = 2] or 2 (mod 3) [a = 1],
s_j = (2^a B^j(t) - 1)/3, sign (-1)^a, literal step clock a + 2j, scale 2^(a + 2j) t/3 in the ideal limit.

PART A (depth one, nodes mod 9, bands 3/2^b < t/N <= 6/2^b, weight 1/(8 2^b)):
  * exact lag constants M_D = sum over preimage pairs at clock lag D of sigma sigma' h h' (per N), by the closed form
      M_D = (4/3)(1/8)(64/63) sum_rho sum_m 2^(1 - a_i) sigma_i sigma_k sum_(e=-D-5)^(-D) 2^e h(r_i, e) h(r_k, e + D)
    (sum over i collapses because B^i permutes the classes mod 9; the e-sum is 6-periodic);
  * the recurrence M_(D+18) = 2^-18 M_D (D >= 2), the total E_1 - E_0 = 2 sum_D M_D, E_0, E_1 and the l1 mass, all
    as exact Fractions; E_1 also from the period-18 band coefficients.
PART B (depth two, nodes mod 27, bands 9/2^(b+1) < t/N <= 9/2^b, weight 1/(16 2^b)):
  * exact E_2 from the period-54 band coefficients;
  * the depth-2 equal-depth meet sums (pairs of words with different children) by clock lag, exact partial sums over
    b < BMAX (integers times 2^-b), compared with (E_2 - E_1)/2; child-scale cutoff test ('children up to scale 6').
Usage: python3 f3_band_exact.py [BMAX]
"""
import sys, os
from fractions import Fraction as Fr
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
KAPPA = (0, -3, 3)
LAM = (-1, -1, 2)


def h(r, e):
    """ideal source coefficient, scale 3*2^e (e <= 0)."""
    assert e <= 0
    tot, k = 0, 0
    while e + 2 * k <= -1:
        tot += KAPPA[(r + k) % 3]
        k += 1
    if e + 2 * k == 0:
        tot += LAM[(r + k) % 3]
    return tot


def B(x, mod):
    return (4 * x + 1) % mod


def pre_class(rho, mod):
    """a and the residue of s = (2^a rho - 1)/3 modulo mod/3, for a node residue rho modulo mod (or None)."""
    if rho % 3 == 1:
        a = 2
    elif rho % 3 == 2:
        a = 1
    else:
        return None
    v = (2 ** a * rho - 1) % mod
    assert v % 3 == 0
    return a, (v // 3) % (mod // 3)


def part_a(out):
    pr = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    pr('== PART A: depth one (first transport), exact ==')
    # periodicity check of h in e
    assert all(h(r, e) == h(r, e - 6) for r in range(3) for e in range(-30, 1))
    # closed-form lag constants
    def M(D):
        tot = Fr(0)
        for rho in range(9):
            pc = pre_class(rho, 9)
            if pc is None:
                continue
            ai, ri = pc
            x = rho
            for m in range(1, D + 2):
                x = B(x, 9)
                pk = pre_class(x, 9)
                if pk is None:
                    continue
                ak, rk = pk
                if ak - ai + 2 * m != D:
                    continue
                sg = (-1) ** (ai + ak)
                inner = sum(Fr(2) ** e * h(ri, e) * h(rk, e + D) for e in range(-D - 5, -D + 1))
                tot += sg * Fr(2) ** (1 - ai) * inner
        return Fr(4, 3) * Fr(1, 8) * Fr(64, 63) * tot
    Ms = {D: M(D) for D in range(1, 80)}
    nz = {D: v for D, v in Ms.items() if v != 0}
    for D in sorted(nz):
        if D <= 31:
            pr('  M_%-2d = %-22s = %+.12e   (denominator %s)' % (D, nz[D], float(nz[D]), nz[D].denominator))
    rec = all(Ms[D + 18] == Ms[D] / 2 ** 18 for D in range(2, 60))
    pr('  recurrence M_(D+18) = 2^-18 M_D for 2 <= D < 60:', rec, '; M_19 - 2^-18 M_1 =', Ms[19] - Ms[1] / 2 ** 18)
    dens = sorted({v.denominator for v in nz.values()})
    pr('  denominators of all nonzero M_D, D < 80: of the form 7*2^k only:',
       all(dd % 7 == 0 and ((dd // 7) & (dd // 7 - 1)) == 0 for dd in dens), '; any divisible by 19 or 73:',
       any(dd % 19 == 0 or dd % 73 == 0 for dd in dens))
    # total: sum_D M_D = M_1 + sum_(D=2..19) M_D / (1 - 2^-18)
    S = Ms[1] + sum(Ms[D] for D in range(2, 20)) / (1 - Fr(1, 2 ** 18))
    # diagonal E_0 by the same closed form
    E0 = Fr(0)
    for rho in range(9):
        pc = pre_class(rho, 9)
        if pc is None:
            continue
        a, r = pc
        E0 += Fr(2) ** (1 - a) * sum(Fr(2) ** e * h(r, e) ** 2 for e in range(-5, 1))
    E0 *= Fr(4, 3) * Fr(1, 8) * Fr(64, 63)
    pr('  E_0 =', E0, '(34/7:', E0 == Fr(34, 7), ')')
    E1_pairs = E0 + 2 * S
    pr('  E_0 + 2 sum_D M_D =', E1_pairs, '=', float(E1_pairs), '; equals 32898/9709:', E1_pairs == Fr(32898, 9709))
    pr('  E_1 - E_0 =', E1_pairs - E0, '=', float(E1_pairs - E0))
    pr('  partial sums: lags 1-13 give 2*sum = %s = %.10f; lags >= 17 give 2*sum = %s = %.3e'
       % (2 * sum(Ms[D] for D in range(1, 14)), float(2 * sum(Ms[D] for D in range(1, 14))),
          2 * (S - sum(Ms[D] for D in range(1, 14))), float(2 * (S - sum(Ms[D] for D in range(1, 14))))))
    # band coefficients C_b(rho) for nodes rho mod 9, b >= 0: sum_j sigma_j h(r_j, e_j), e_j = a_j + 2j - b - 1
    def C(rho, b):
        tot, x, j = 0, rho, 0
        while True:
            pc = pre_class(x, 9)
            if pc is not None:
                a, r = pc
                e = a + 2 * j - b - 1
                if e <= 0:
                    tot += (-1) ** a * h(r, e)
            if 2 * j - b - 1 >= 1:          # every further preimage has e >= 2j - b > 0
                break
            x = B(x, 9)
            j += 1
        return tot
    per = all(C(rho, b) == C(rho, b + 18) for rho in range(9) for b in range(0, 60))
    l1 = sum(Fr(1, 8 * 2 ** b) * sum(abs(C(rho, b)) for rho in range(9)) for b in range(18)) / (1 - Fr(1, 2 ** 18))
    E1 = sum(Fr(1, 8 * 2 ** b) * sum(C(rho, b) ** 2 for rho in range(9)) for b in range(18)) / (1 - Fr(1, 2 ** 18))
    pr('  band coefficients have period 18 in b (b < 60):', per)
    pr('  l1 first-transport mass =', l1, '(155896/87381:', l1 == Fr(155896, 87381), ')')
    pr('  E_1 from band coefficients =', E1, '(32898/9709:', E1 == Fr(32898, 9709), ')')
    pr('  9709 = (2^18-1)/27:', (2 ** 18 - 1) // 27 == 9709, '; E_1 - E_0 = -99820/67963:', E1 - Fr(34, 7) == Fr(-99820, 67963))
    return Ms


def words2(r27, b, child_cut=6):
    """depth-2 backward words of a node t = r27 (mod 27) in band b (representative t/N = 27/2^(b+2)):
    list of (child index, clock, signed ideal value).  Child scale 9*2^(tau1 - b - 2) <= child_cut."""
    out = []
    x, j1, ci = r27 % 27, 0, 0
    while True:
        pc = pre_class(x, 27)
        e1 = 2 * j1 - b - 2                     # child scale 9*2^(tau1-b-2) with tau1 = a1 + 2 j1 >= 2 j1 + 1
        if 9 * 2.0 ** (2 * j1 + 1 - b - 2) > child_cut:
            break
        if pc is not None:
            a1, p9 = pc
            tau1 = a1 + 2 * j1
            if 9 * 2.0 ** (tau1 - b - 2) <= child_cut:
                y, j2 = p9, 0
                while True:
                    pc2 = pre_class(y, 9)
                    if 2 * j2 + 1 + tau1 - b - 2 > 0:   # grandchild e = tau - b - 2 > 0 for all further j2
                        break
                    if pc2 is not None:
                        a2, s3 = pc2
                        tau = tau1 + a2 + 2 * j2
                        e = tau - b - 2
                        if e <= 0:
                            v = h(s3, e)
                            if v:
                                out.append((ci, tau, (-1) ** (a1 + a2) * v))
                    y = B(y, 9)
                    j2 += 1
            ci += 1
        x = B(x, 27)
        j1 += 1
    return out


def part_b(out, BMAX):
    pr = lambda *a: (print(*a, flush=True), out.append(' '.join(str(x) for x in a)))
    pr('== PART B: depth two, exact ==')
    E1 = Fr(32898, 9709)
    # E_2 from band coefficients (period 54)
    def K2(r, b):
        return sum(v for _, _, v in words2(r, b))
    per = all(K2(r, b) == K2(r, b + 54) for r in range(27) for b in range(0, 30))
    E2 = sum(Fr(1, 16 * 2 ** b) * sum(K2(r, b) ** 2 for r in range(27)) for b in range(54)) / (1 - Fr(1, 2 ** 54))
    pr('  K_2 has period 54 in b (b < 30 checked against b + 54):', per)
    pr('  E_2 =', E2, '; equals 572764898734162/222399981598543:', E2 == Fr(572764898734162, 222399981598543))
    target = (E2 - E1) / 2
    pr('  (E_2 - E_1)/2 = %.15f' % float(target))
    # lag-resolved meet sums, exact partial sums (integer correlations per band, then exact weights 1/(16 2^b))
    import numpy as np
    res = {}
    for cut, bmax in ((6, BMAX), (4, 60), (12, 60)):
        lagS = defaultdict(int)                       # lag -> exact Fraction-free numerator: sum_b S_D(b) 2^(bmax-b)
        for b in range(bmax):
            Sb = defaultdict(int)
            for r in range(27):
                ws = words2(r, b, child_cut=cut)
                if not ws:
                    continue
                tmin = min(t for _, t, _ in ws); tmax = max(t for _, t, _ in ws)
                chl = sorted({ci for ci, _, _ in ws}); cidx = {c: i for i, c in enumerate(chl)}
                Vc = np.zeros((len(chl), tmax - tmin + 1), dtype=np.int64)
                for ci, t, v in ws:
                    Vc[cidx[ci], t - tmin] += v
                V = Vc.sum(axis=0)
                L = V.shape[0]
                allc = np.correlate(V, V, 'full')[L - 1:]          # allc[D] = sum_t V[t] V[t+D]
                samec = np.zeros(L, dtype=np.int64)
                for row in Vc:
                    samec += np.correlate(row, row, 'full')[L - 1:]
                cross = allc - samec
                for D in range(L):
                    c = int(cross[D])
                    if c:
                        Sb[D] += c
            for D, c in Sb.items():
                if D == 0:
                    assert c % 2 == 0
                    c //= 2
                lagS[D] += c * 2 ** (bmax - b)
        den = 16 * 2 ** bmax
        lagF = {D: Fr(v, den) for D, v in lagS.items()}
        res[cut] = lagF
        tot = sum(lagF.values())
        if cut == 6:
            pr('  child cutoff 6, b < %d: total meet sum = %.15f ; |total - (E_2-E_1)/2| = %.3e'
               % (bmax, float(tot), abs(float(tot - target))))
            for D in sorted(lagF):
                if D <= 12 or abs(float(lagF[D])) > 1e-6:
                    pr('    lag %2d: %+.10f' % (D, float(lagF[D])))
            rest = sum(v for D, v in lagF.items() if D >= 10)
            pr('    all lags >= 10 together: %+.3e' % float(rest))
        else:
            ref = sum(v for D, v in res[6].items())
            pr('  child cutoff %s (b < %d): total %.15f, lag 0 %+.10f' % (cut, bmax, float(tot), float(lagF.get(0, 0))))


def main():
    BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    out = []
    part_a(out)
    part_b(out, BMAX)
    open(os.path.join(HERE, 'f3_band_exact.out'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
