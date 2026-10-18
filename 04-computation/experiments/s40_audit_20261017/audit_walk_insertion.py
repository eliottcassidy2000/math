"""S40 audit: THM-4666 statement 5 (insertion changes = -2 S(t), symmetry, L+ counts, Fekete identity, Mellin formula).

 (a) every single-vertex insertion move of the natural order of P_p, excess recomputed from scratch (p <= 71), compared
     with -2 S(t); for p <= 131 'no improving insertion' compared with 'all S(t) >= 0';
 (b) S(p - 1 - t) = S(t), S(p - 1) = 0 for all p = 3 mod 4 below 2e5;
 (c) L+ count in [7, 2e5], proportions in the L(1, chi) bins of THM-4666 (7), dyadic fractions (independent h);
 (d) Fekete identity f_p(x) = (1 - x) sum_t S(t) x^t as an integer polynomial identity (p < 3000);
 (e) Gamma(s) L(s, chi) = int_0^oo u^(s-1) f_p(e^-u) / (1 - e^(-p u)) du (mpmath, p = 23, 43; s = 0.5, 1, 2).
Usage: python audit_walk_insertion.py
"""
import sys, os, math, json, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import sieve, legendre, h_forms, paley_chi

OUT = {}


@njit(cache=True)
def excess_order(order, chi, p):
    s = 0
    for i in range(p):
        u = order[i]
        for j in range(i + 1, p):
            s += chi[(order[j] - u) % p]
    return s


@njit(cache=True)
def insertion_check(chi, p):
    """Return (#moves, #mismatches, max gain) over all moves of one vertex of the natural order to another position."""
    nat = np.arange(p)
    base = excess_order(nat, chi, p)
    S = np.zeros(p, dtype=np.int64)
    s = 0
    for t in range(1, p):
        s += chi[t]
        S[t] = s
    moves = 0
    bad = 0
    maxgain = -10 ** 9
    order = np.empty(p, dtype=np.int64)
    for i in range(p):
        for k in range(p):
            if k == i:
                continue
            # remove element at i, insert at k
            idx = 0
            for j in range(p):
                if j == i:
                    continue
                order[idx] = nat[j]
                idx += 1
            # now order[0..p-2] holds the rest; insert nat[i] at position k
            for j in range(p - 1, k, -1):
                order[j] = order[j - 1]
            order[k] = nat[i]
            e = excess_order(order, chi, p)
            t = k - i if k > i else i - k
            moves += 1
            if e - base != -2 * S[t]:
                bad += 1
            if e - base > maxgain:
                maxgain = e - base
    return moves, bad, maxgain


@njit(cache=True)
def walk_stats(chi, p):
    s = 0
    mn = 10 ** 9
    S = np.zeros(p, dtype=np.int64)
    for t in range(1, p):
        s += chi[t]
        S[t] = s
        if s < mn:
            mn = s
    sym_bad = 0
    for t in range(0, p):
        if S[p - 1 - t] != S[t]:
            sym_bad += 1
    return mn, S[p - 1], sym_bad


def main():
    t0 = time.time()
    SV = sieve(200000)
    P = [p for p in range(7, 200001, 4) if SV[p]]
    # (a)
    rows = []
    for p in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71]:
        chi = paley_chi(p)
        moves, bad, mg = insertion_check(chi, p)
        rows.append((p, moves, bad, mg))
    print('(a) all single-vertex moves of the natural order (excess recomputed): p, moves, mismatches with -2S(t), max gain:', rows, flush=True)
    OUT['insertion_moves'] = rows
    agree = 0; disagree = []
    for p in [p for p in P if p <= 131]:
        chi = paley_chi(p)
        moves, bad, mg = insertion_check(chi, p)
        mn, last, sb = walk_stats(chi, p)
        if (mg <= 0) == (mn >= 0) and bad == 0:
            agree += 1
        else:
            disagree.append(p)
    print('    p <= 131: "no improving move" <=> "min S >= 0" agrees at %d primes, disagrees at %s' % (agree, disagree), flush=True)
    # (b), (c)
    symbad = 0; lastbad = 0
    rec = []
    for p in P:
        chi = paley_chi(p)
        mn, last, sb = walk_stats(chi, p)
        symbad += (sb > 0); lastbad += (last != 0)
        h = h_forms(-p)
        rec.append((p, h, math.pi * h / math.sqrt(p), int(mn)))
    nL = sum(1 for r in rec if r[3] >= 0)
    print('(b) primes 7..2e5: %d; symmetry failures %d; S(p-1) != 0 failures %d' % (len(rec), symbad, lastbad), flush=True)
    print('(c) L+ count: %d of %d (%.4f)' % (nL, len(rec), nL / len(rec)), flush=True)
    bins = []
    for lo, hi in [(0, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 99)]:
        sel = [r for r in rec if lo <= r[2] < hi]
        inn = [r for r in sel if r[3] >= 0]
        bins.append((lo, hi, len(sel), len(inn), round(len(inn) / max(1, len(sel)), 3)))
    print('    L-bins (lo, hi, primes, in L+, share):', bins, flush=True)
    dy = []
    for k in range(3, 18):
        sel = [r for r in rec if 2 ** (k - 1) < r[0] <= 2 ** k]
        if sel:
            dy.append((k, len(sel), round(sum(r[3] >= 0 for r in sel) / len(sel), 4)))
    print('    dyadic (k, primes in (2^(k-1), 2^k], L+ share):', dy, flush=True)
    Lp = [r[0] for r in rec if r[3] >= 0]
    print('    first L+ members:', Lp[:30], '; 463 in L+:', 463 in Lp, '; 43 in L+:', 43 in Lp, flush=True)
    OUT['Lplus'] = dict(count=nL, total=len(rec), bins=bins, dyadic=dy, first=Lp[:40], has463=(463 in Lp))
    # (d) Fekete identity as integer polynomials
    fb = 0
    for p in [p for p in P if p < 3000]:
        chi = [legendre(t, p) for t in range(p)]
        S = np.cumsum([0] + chi[1:])        # S[t] for t = 0..p-1 (S[0] = 0)
        # (1 - x) sum_t S(t) x^t : coefficient of x^t is S(t) - S(t-1) for 1 <= t <= p-1, and -S(p-1) at x^p
        coeff = [0] * (p + 1)
        for t in range(1, p):
            coeff[t] += int(S[t]); coeff[t + 1] -= int(S[t])
        if coeff[:p] != chi or coeff[p] != 0:
            fb += 1
    print('(d) Fekete identity failures (p < 3000):', fb, flush=True)
    OUT['fekete_failures'] = fb
    # (e) Mellin formula
    import mpmath as mp
    mp.mp.dps = 30
    res = []
    for p in (23, 43):
        chi = [legendre(t, p) for t in range(p)]
        f = lambda x: sum(chi[t] * x ** t for t in range(1, p))
        for s in (0.5, 1, 2):
            lhs = mp.gamma(s) * mp.dirichlet(s, chi)
            g = lambda u: u ** (s - 1) * f(mp.e ** (-u)) / (1 - mp.e ** (-p * u))
            rhs = mp.quad(g, [0, 0.01, 0.1, 1, 10, mp.inf])
            res.append((p, s, float(lhs), float(rhs)))
    print('(e) Mellin check (p, s, Gamma(s)L(s,chi), integral):', res, flush=True)
    OUT['mellin'] = res
    # positivity of f_p on (0,1) for p = 43 (not in L+): sign of f_p on a grid
    chi = [legendre(t, 43) for t in range(43)]
    xs = np.linspace(0.001, 0.999, 999)
    vals = [sum(chi[t] * x ** t for t in range(1, 43)) for x in xs]
    print('    p = 43: min f_p on grid (0,1) = %.4f' % min(vals))
    print('elapsed %.1f s' % (time.time() - t0))
    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'audit_walk_insertion.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
