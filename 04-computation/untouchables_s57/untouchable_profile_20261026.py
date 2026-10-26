#!/usr/bin/env python3
"""Profiles of the even untouchables <= N from the exact multiplicity array (opus-2026-10-26-S57).

Input: mult_even_N.npy from untouchable_sieve_20261026.py (index j <-> even k = 2j, value = |s^{-1}(k)| capped at 255).
Computes, for even k in [N/10, N] (a late window) and in [2, N]:
  (a) the multiplicity distribution, mean, variance (Poisson comparison);
  (b) the untouchable share by residue class of k mod M for M in {3,4,5,7,8,9,14,15,30,105,210};
  (c) per-k Pollack-Pomerance test: bin even k by h(k) = s(k)/k and compare the untouchable share with e^{-1/h};
  (d) shape split: k-1 prime (preimage (k-1)^2), k-3 prime (preimage 2(k-3)), generic multiplicity g = mult - those;
  (e) the Pollack-Pomerance weighted sum  sum_{a even <= y} e^{-a/s(a)}/a  /  sum_{a <= y} 1/a  at y = 10^e.
Usage: python untouchable_profile_20261026.py N bigdir
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto


@nb.njit(cache=False)
def s_even_segment(lo, cnt, primes, out):
    """out[i] = s(m) for even m = lo + 2i."""
    rem = np.empty(cnt, np.int64)
    sig = np.ones(cnt, np.int64)
    for i in range(cnt):
        rem[i] = lo + 2 * i
    for i in range(cnt):
        pk = 1
        s = 1
        while rem[i] % 2 == 0:
            rem[i] //= 2
            pk *= 2
            s += pk
        sig[i] = s
    hi = lo + 2 * cnt
    for t in range(1, primes.shape[0]):
        p = primes[t]
        if p * p >= hi:
            break
        step = 2 * p
        first = ((lo + step - 1) // step) * step
        for m in range(first, hi, step):
            i = (m - lo) // 2
            pk = 1
            s = 1
            while rem[i] % p == 0:
                rem[i] //= p
                pk *= p
                s += pk
            sig[i] *= s
    for i in range(cnt):
        r = rem[i]
        sg = sig[i]
        if r > 1:
            sg *= (r + 1)
        out[i] = sg - (lo + 2 * i)


@nb.njit(cache=False)
def prime_flags_segment(lo, hi, base, flags):
    for i in range(hi - lo):
        flags[i] = 1
    for t in range(base.shape[0]):
        p = base[t]
        if p * p >= hi:
            break
        st = max(p * p, ((lo + p - 1) // p) * p)
        for m in range(st, hi, p):
            flags[m - lo] = 0
    for i in range(hi - lo):
        if lo + i < 2:
            flags[i] = 0


@nb.njit(cache=False)
def accumulate(lo_k, cnt, svals, mult, pf_lo, pflags, hbins, nb_h, mods, res_tot, res_unt, shape_tab,
               late_start, pp_sum, pp_marks, pp_idx):
    """k = lo_k + 2i (even). Accumulate per-k statistics."""
    for i in range(cnt):
        k = lo_k + 2 * i
        if k < 8:
            continue
        j = k // 2
        m = np.int64(mult[j])
        unt = 1 if m == 0 else 0
        late = 1 if k >= late_start else 0
        s = svals[i]
        h = s / k
        # PP weighted sum over a = k (even a), e^{-a/s(a)}/a
        if s > 0:
            pp_sum[0] += math.exp(-k / s) / k
        while pp_idx[0] < pp_marks.shape[0] and k >= pp_marks[pp_idx[0]]:
            pp_marks_val = pp_sum[0]
            pp_idx[0] += 1
        if late == 1:
            # h bins on 1/h in [0, 4): bin width 4/nb_h
            x = 1.0 / h if h > 0 else 99.0
            b = int(x / (4.0 / nb_h))
            if b >= nb_h:
                b = nb_h - 1
            hbins[b, 0] += 1
            hbins[b, 1] += unt
            for t in range(mods.shape[0]):
                M = mods[t]
                r = k % M
                res_tot[t, r] += 1
                res_unt[t, r] += unt
            # shapes: k-1 prime, k-3 prime
            a = np.int64(pflags[k - 1 - pf_lo])
            b3 = np.int64(pflags[k - 3 - pf_lo])
            g = np.int64(m) - a - b3
            if g < 0:
                g = 0
            if g > 3:
                g = 3
            shape_tab[a, b3, g] += 1


def run(N, bigdir, seg=1 << 22):
    mult = np.load(os.path.join(bigdir, f'mult_even_{N}.npy'))
    primes = primes_upto(int(math.isqrt(2 * N)) + 2)
    base = primes_upto(int(math.isqrt(N)) + 2)
    nb_h = 40
    hbins = np.zeros((nb_h, 2), np.int64)
    mods = np.array([3, 4, 5, 7, 8, 9, 11, 13, 14, 15, 30, 105, 210, 1155], np.int64)
    res_tot = np.zeros((mods.shape[0], 1155), np.int64)
    res_unt = np.zeros((mods.shape[0], 1155), np.int64)
    shape_tab = np.zeros((2, 2, 4), np.int64)
    late_start = N // 10
    pp_sum = np.zeros(1)
    marks = []
    x = 10
    while x <= N:
        marks.append(x)
        x *= 10
    pp_vals = {}
    lo = 2
    t0 = time.time()
    next_mark = 0
    while lo <= N:
        cnt = min(seg, (N - lo) // 2 + 1)
        sv = np.empty(cnt, np.int64)
        s_even_segment(lo, cnt, primes, sv)
        hi_k = lo + 2 * cnt
        pf_lo = max(0, lo - 4)
        pflags = np.empty(hi_k - pf_lo, np.uint8)
        prime_flags_segment(pf_lo, hi_k, base, pflags)
        dummy_marks = np.zeros(0, np.int64)
        idx = np.zeros(1, np.int64)
        before = pp_sum[0]
        accumulate(lo, cnt, sv, mult, pf_lo, pflags, hbins, nb_h, mods, res_tot, res_unt, shape_tab, late_start,
                   pp_sum, dummy_marks, idx)
        # record PP partial sums at powers of ten (segment granularity: exact if marks fall on segment ends)
        while next_mark < len(marks) and marks[next_mark] < hi_k:
            # recompute partial sum exactly up to the mark within this segment
            part = before
            for i in range(cnt):
                k = lo + 2 * i
                if k > marks[next_mark]:
                    break
                if k >= 8 and sv[i] > 0:
                    part += math.exp(-k / sv[i]) / k
            pp_vals[marks[next_mark]] = part
            next_mark += 1
        lo = hi_k
    return mult, hbins, mods, res_tot, res_unt, shape_tab, pp_vals, time.time() - t0


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    bigdir = sys.argv[2]
    mult, hbins, mods, res_tot, res_unt, shape_tab, pp_vals, dt = run(N, bigdir)
    out = {'N': N, 'seconds': dt}
    ev = mult[4:N // 2 + 1].astype(np.int64)     # even k >= 8
    late = mult[N // 20:N // 2 + 1].astype(np.int64)
    for name, arr in (('all_even_k>=8', ev), ('late_even_k_in_[N/10,N]', late)):
        mu = float(arr.mean())
        var = float(arr.var())
        z = float((arr == 0).mean())
        hist = np.bincount(np.minimum(arr, 12), minlength=13)
        out[name] = {'mean': mu, 'var': var, 'zero_share': z, 'poisson_zero': math.exp(-mu),
                     'hist_0_to_12plus': [int(v) for v in hist]}
        print(f'{name}: mean multiplicity {mu:.4f}, variance {var:.4f}, untouchable share {z:.5f}, '
              f'Poisson e^-mean {math.exp(-mu):.5f}; histogram {list(hist[:8])}...', flush=True)
    print('per-k Pollack-Pomerance test, late window: bins of x = 1/h(k) = k/s(k)')
    rows = []
    for b in range(hbins.shape[0]):
        tot, unt = int(hbins[b, 0]), int(hbins[b, 1])
        if tot < 1000:
            continue
        xm = (b + 0.5) * 4.0 / hbins.shape[0]
        rows.append((round(xm, 3), tot, unt / tot, math.exp(-xm)))
        print(f'   x in [{b * 0.1:.1f},{(b + 1) * 0.1:.1f}): n={tot:>10}  untouchable share {unt / tot:.4f}   '
              f'e^-x {math.exp(-xm):.4f}', flush=True)
    out['pp_per_k_bins'] = rows
    print('untouchable share among even k in the late window, by residue class:')
    rc = {}
    for t, M in enumerate(mods):
        shares = {}
        for r in range(M):
            if res_tot[t, r] > 0:
                shares[r] = res_unt[t, r] / res_tot[t, r]
        rc[int(M)] = {int(r): [int(res_tot[t, r]), round(v, 6)] for r, v in shares.items()}
        if M <= 30:
            print(f'   mod {M}: ' + ', '.join(f'{r}:{v:.4f}' for r, v in sorted(shares.items())), flush=True)
    out['residue_classes'] = rc
    print('shape table, late window: rows (k-1 prime, k-3 prime), columns generic multiplicity 0,1,2,3+:')
    st = {}
    for a in (0, 1):
        for b3 in (0, 1):
            row = [int(v) for v in shape_tab[a, b3]]
            st[f'{a}{b3}'] = row
            tot = sum(row)
            print(f'   k-1 prime={a}, k-3 prime={b3}: n={tot}, generic-zero share {row[0] / tot:.4f}, '
                  f'generic histogram {row}', flush=True)
    out['shape_table'] = st
    H = {x: v for x, v in pp_vals.items()}
    ppr = {}
    for x, v in sorted(H.items()):
        harm = math.log(x) + 0.5772156649 + 1 / (2 * x)
        ppr[str(x)] = v / harm
    print('Pollack-Pomerance weighted sum / harmonic sum at 10^e:', {k: round(v, 6) for k, v in ppr.items()})
    out['pp_ratio'] = ppr
    json.dump(out, open(os.path.join(HERE, f'untouchable_profile_{N}.json'), 'w'), indent=1)
