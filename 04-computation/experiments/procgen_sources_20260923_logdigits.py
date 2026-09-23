#!/usr/bin/env python3
"""procgen_sources_20260923_logdigits.py

The pasted {log n} construction: interleave the base-b digits of {log n},
n >= 2, round by round, so that after round r the prefix contains exactly the
rectangle {(n, j) : 2 <= n <= b^r, 1 <= j <= r}; T_r = r (b^r - 1).

The source was not located (see the note), so this script analyses the object
under both natural readings of "log" and two within-round orders:

  reading 'logb' : D(n, j) = j-th base-b digit of {log_b n}  (a log table:
                   for b^(r-1) < n <= b^r the row is the r-digit mantissa log
                   of the r-digit number n);
  reading 'ln'   : D(n, j) = j-th base-b digit of {ln n};
  order 'rows'   : round r = new column r for old rows (n increasing), then the
                   complete new rows n = b^(r-1)+1 .. b^r (row-major);
  order 'cols'   : round r = for j = 1..r, all new cells of column j (n increasing).

Outputs:
  L1 the rectangle invariant and T_r for every round;
  L2 separation: are the r-digit prefixes of the new rows distinct?
  L3 digit and block statistics of the whole prefix (k = 1, 2, 3), i.e. how
     far the number is from normal at the last checkpoint;
  L4 the column-1 law against the Benford-type prediction;
  L5 the 2-versus-3 slice: base 3, rows n = 2^k, column 1 = floor(3 {k log_3 2}),
     a 3-arc coding of the rotation by log_3 2 (the critical Collatz slope),
     with its factor complexity and a repetition (Dio) estimate; and the
     identity "top ternary digits of 2^k are read off {k log_3 2}".

Digits are computed in float64 with an exact fallback (mpmath, 60 digits) for
every value within 1e-6 (relative to the last digit) of a digit boundary.
Deterministic; peak memory about 150 MB.
"""
import math
import sys
from collections import Counter

import numpy as np
import mpmath

mpmath.mp.dps = 60


def say(*a):
    print(*a)
    sys.stdout.flush()


def digit_table(b, R, reading):
    """D[n, j-1] for 2 <= n <= b^R, 1 <= j <= R (rows 0, 1 unused)."""
    N = b ** R
    n = np.arange(N + 1, dtype=np.float64)
    n[:2] = 1.0
    if reading == 'logb':
        if b == 2:
            x = np.log2(n)
        elif b == 10:
            x = np.log10(n)
        else:
            x = np.log(n) / math.log(b)
    else:
        x = np.log(n)
    frac = x - np.floor(x)
    scaled = frac * float(b) ** R
    near = np.abs(scaled - np.rint(scaled)) < 1e-6
    # exact powers of b have fractional part 0 in the logb reading
    ints = np.floor(scaled).astype(np.int64)
    ints[np.rint(scaled) == float(b) ** R] = 0
    fixed = 0
    for m in np.nonzero(near)[0]:
        m = int(m)
        if m < 2:
            continue
        if reading == 'logb':
            k = 0
            t = m
            while t % b == 0:
                t //= b
                k += 1
            if t == 1:
                ints[m] = 0
                fixed += 1
                continue
            v = mpmath.log(m) / mpmath.log(b)
        else:
            v = mpmath.log(m)
        f = v - mpmath.floor(v)
        ints[m] = int(mpmath.floor(f * mpmath.mpf(b) ** R))
        fixed += 1
    D = np.zeros((N + 1, R), dtype=np.uint8)
    rem = ints.copy()
    for j in range(R - 1, -1, -1):
        D[:, j] = (rem % b).astype(np.uint8)
        rem //= b
    return D, fixed


def build(b, R, D, order):
    seq = []
    T = []
    for r in range(1, R + 1):
        lo, hi = b ** (r - 1), b ** r
        if order == 'rows':
            if lo >= 2:
                seq.append(D[2:lo + 1, r - 1])
            block = D[lo + 1:hi + 1, :r]
            seq.append(block.reshape(-1))
        else:
            for j in range(1, r + 1):
                if j == r and lo >= 2:
                    seq.append(D[2:lo + 1, r - 1])
                seq.append(D[lo + 1:hi + 1, j - 1])
        T.append(sum(len(s) for s in seq))
    return np.concatenate(seq), T


def block_dev(seq, b, k, chunk=1 << 21):
    """max_w |freq(w) b^k - 1| over overlapping k-blocks (chunked, int32)."""
    cnt = np.zeros(b ** k, dtype=np.int64)
    n = len(seq)
    for start in range(0, n - k + 1, chunk):
        stop = min(n - k + 1, start + chunk)
        code = np.zeros(stop - start, dtype=np.int32)
        for i in range(k):
            code = code * b + seq[start + i:stop + i].astype(np.int32)
        cnt += np.bincount(code, minlength=b ** k)
    f = cnt / cnt.sum()
    return float(np.max(np.abs(f * b ** k - 1))), int(np.argmax(np.abs(f * b ** k - 1)))


def section_L(b, R):
    say(f'--- base b = {b}, rounds r <= {R} (T_R = {R * (b ** R - 1)})')
    for reading in ('logb', 'ln'):
        D, fixed = digit_table(b, R, reading)
        # L2 separation of r-digit prefixes of the new rows of round R
        lo, hi = b ** (R - 1), b ** R
        rows = D[lo + 1:hi + 1, :R].astype(np.int64)
        codes = np.zeros(len(rows), dtype=np.int64)
        for j in range(R):
            codes = codes * b + rows[:, j]
        distinct = len(np.unique(codes)) / len(codes)
        # L4 column 1 of the new rows vs prediction
        col1 = Counter(D[lo + 1:hi + 1, 0].tolist())
        tot = hi - lo
        if reading == 'logb':
            pred = [(b ** ((d + 1) / b) - b ** (d / b)) / (b - 1) for d in range(b)]
        else:
            # exact finite-R density: n uniform on (b^(R-1), b^R]  <=>  u = ln n with density e^u
            Lb = math.log(b)
            u0, u1 = (R - 1) * Lb, R * Lb
            pred = []
            for d in range(b):
                tot_d = 0.0
                for m in range(int(math.floor(u0)) - 1, int(math.ceil(u1)) + 1):
                    a0, a1 = max(u0, m + d / b), min(u1, m + (d + 1) / b)
                    if a1 > a0:
                        tot_d += math.exp(a1 - u1) - math.exp(a0 - u1)
                pred.append(tot_d / (1 - math.exp(u0 - u1)))
        emp = [col1[d] / tot for d in range(b)]
        for order in ('rows', 'cols'):
            seq, T = build(b, R, D, order)
            okT = all(T[r - 1] == r * (b ** r - 1) for r in range(1, R + 1))
            devs = []
            for k in (1, 2, 3):
                dv, arg = block_dev(seq, b, k)
                devs.append(f'k={k}: {dv:.4f}')
            say(f'  [{reading:4s}|{order:4s}] T_r = r(b^r-1) for r=1..{R}: {okT}; len {len(seq)}; '
                f'max|freq*b^k-1|: ' + ', '.join(devs))
        say(f'  [{reading:4s}] fallback digits recomputed exactly: {fixed}; '
            f'distinct {R}-digit prefixes among rows ({lo},{hi}]: {distinct:.4f}')
        say(f'  [{reading:4s}] column-1 law of the new rows of round {R}, empirical/predicted '
            f'(logb: Benford law, R-independent; ln: exact finite-R law, oscillates with R): '
            + ', '.join(f'{d}:{e:.4f}/{p:.4f}' for d, (e, p) in enumerate(zip(emp, pred))))
    # convergence profile of the digit bias (rows order, logb reading) across checkpoints
    D, _ = digit_table(b, R, 'logb')
    seq, T = build(b, R, D, 'rows')
    prof = []
    for r in range(max(1, R - 4), R + 1):
        dv, _ = block_dev(seq[:T[r - 1]], b, 1)
        prof.append(f'r={r}: {dv:.4f} (r*dev={r * dv:.3f})')
    say('  [logb|rows] 1-block bias at the checkpoints T_r: ' + '; '.join(prof))
    seqc, Tc = build(b, R, D, 'cols')
    profc = []
    for r in range(max(2, R - 4), R + 1):
        dv, _ = block_dev(seqc[:Tc[r - 1]], b, 2)
        profc.append(f'r={r}: {dv:.3f}')
    say('  [logb|cols] 2-block bias at the checkpoints T_r (does not decay): ' + '; '.join(profc))
    del seqc
    # ln reading: column-1 law at two consecutive rounds (oscillation)
    D2, _ = digit_table(b, R, 'ln')
    for rr in (R - 1, R):
        lo, hi = b ** (rr - 1), b ** rr
        c = Counter(D2[lo + 1:hi + 1, 0].tolist())
        say(f'  [ln  ] column-1 frequencies of the new rows of round {rr}: '
            + ', '.join(f'{d}:{c[d] / (hi - lo):.4f}' for d in range(b)))


def lpf_dio(w, jmin):
    """Dio estimate 1 + max_{j >= jmin} LPF(j)/j with a simple O(n^2)-free method:
    LPF via Z-function on each suffix is too slow; use a suffix-automaton-free
    approach: for each j, LPF(j) = max_{a<j} LCE(a, j), computed with hashing
    and binary search over a restricted candidate set (previous occurrences of
    the next 12 letters)."""
    n = len(w)
    B1, M1 = 1315423911, (1 << 61) - 1
    h = [0] * (n + 1)
    pw = [1] * (n + 1)
    for i, c in enumerate(w):
        h[i + 1] = (h[i] * B1 + c + 1) % M1
        pw[i + 1] = pw[i] * B1 % M1

    def hs(i, L):
        return (h[i + L] - h[i] * pw[L]) % M1
    K = 12
    occ = {}
    best = 0.0
    lpf_at = []
    for j in range(n - K):
        key = hs(j, K)
        cands = occ.get(key, [])
        L_best = 0
        for a in cands[-64:]:
            lo, hi = K, n - j
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if hs(a, mid) == hs(j, mid):
                    lo = mid
                else:
                    hi = mid - 1
            L_best = max(L_best, lo)
        if j >= jmin and j > 0:
            best = max(best, L_best / j)
        lpf_at.append(L_best)
        occ.setdefault(key, []).append(j)
    return 1 + best


def section_L5():
    say('--- L5: the 2-versus-3 slice (base 3, rows n = 2^k, reading log_3)')
    K = 30000
    alpha = mpmath.log(2) / mpmath.log(3)
    word = []
    top_ok = True
    for k in range(1, K + 1):
        f = k * alpha - mpmath.floor(k * alpha)
        d1 = int(mpmath.floor(3 * f))
        word.append(d1)
        if k <= 400:
            # top ternary digit of 2^k from the integer itself vs from 3^f
            t = 2 ** k
            s = []
            while t:
                s.append(t % 3)
                t //= 3
            top = s[-1]
            mant = mpmath.mpf(3) ** f  # in [1, 3)
            if int(mpmath.floor(mant)) != top:
                top_ok = False
    say(f'  top ternary digit of 2^k equals floor(3^frac(k log_3 2)) for k <= 400: {top_ok}')
    comp = []
    for m in (1, 2, 3, 5, 8, 13, 21, 34, 55):
        comp.append((m, len({tuple(word[i:i + m]) for i in range(len(word) - m)})))
    say('  factor complexity p(m) of the column-1 word (3-arc coding of the rotation by log_3 2): '
        + ', '.join(f'p({m})={c}' for m, c in comp))
    dio = lpf_dio(word[:6000], 200)
    say(f'  finite-window repetition estimate 1 + max LPF(j)/j (j >= 200) on 6000 letters: {dio:.4f} '
        f'(Lemma J\' class; driven by the partial quotients of log_2 3 = [1;1,1,2,2,3,1,5,2,23,2,...])')
    # the critical Collatz Sturmian word of slope log_3 2 is the 2-arc coding of the same rotation
    stw = [int(mpmath.floor((k + 1) * alpha) - mpmath.floor(k * alpha)) for k in range(1, 2001)]
    say(f'  Sturmian word of slope log_3 2 (critical 3x+1 parity word), first 40 letters: '
        + ''.join(map(str, stw[:40])))
    say(f'  its complexity p(10) = {len({tuple(stw[i:i+10]) for i in range(len(stw)-10)})} (Sturmian: m+1 = 11)')


def main():
    say('procgen_sources_20260923_logdigits.py')
    say('=' * 78)
    say('L1-L4: the rectangle construction')
    section_L(2, 20)
    section_L(3, 12)
    section_L(10, 6)
    say('=' * 78)
    section_L5()
    say('=' * 78)
    say('logdigits: done')


if __name__ == '__main__':
    main()
