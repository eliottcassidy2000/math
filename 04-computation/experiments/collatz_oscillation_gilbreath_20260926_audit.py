#!/usr/bin/env python3
"""collatz_oscillation_gilbreath_20260926_audit.py -- adversarial audit of
05-knowledge/results/collatz_oscillation_gilbreath_20260926.md (session collatz-oscillation-20260926) and its
scripts collatz_oscillation_20260926_shells.py, gilbreath_20260926_ca.py, collatz_base6_ca_20260926.py.
Independent re-implementation (exact integer arithmetic for the orbit statistics; numpy for the Gilbreath
triangle), written from the note's and THM-4506's statements.

Sections (numbered as in the audit brief):
 1. Proposition 1: (a) m(j) versus the count of shell indices in the window with no earlier 2^D-drop (exact
    equality checked on every landing point), Lemma S / Lemma O / halving-step checks, runs versus dippers;
    (b) the partition accounting N(X) = #E + #Dip + #ND, landing points below X 2^-D, both directions of the
    regularity form.
 2. The probe table: total, ND, dippers, landing points, mean multiplicity, heaviest (mult, shell runs) and the
    ratio, recomputed for the record orbit of 63728127, the 5n+1 orbit of 7 and 2^L-1 at L = 24, 32, 48, 64,
    compared with the committed shells.out; the ratio's definition; float-versus-exact dip test.
 3. Gilbreath: (a) the algebra of the rule (XOR, |1-a|, which neighbour shrinks a defect, the stationary copy);
    (b) the triangle of the primes below 200000 (leading entries, frontier, exact first all-0/2 row, closure,
    frontier sequence, zero-run histogram of the script's region reproduced, the triangle check shown vacuous,
    full-sea histograms of runs and of true triangle sides, defect samples reproduced and de-duplicated);
    (c) single seed versus random row, with the script's own detector and with true triangle tops;
    (d) the random-model survival probability by Monte Carlo against the note's binomial tail.
 4. Base-6 locality: an explicit radius-1 rule verified exhaustively for x < 6^7 and on random 40-digit x;
    the carry table; the +1 at digit 0; dependence on both neighbours; the Cloney-Goles-Vichniac citation.
 5. Thresholds: THM-4487's curve E(gamma), where theta -> 0 sits on it, the Gilbreath tail versus 2^-F.
 6. Text checks (labels, 'about row 100', the Cloney-Goles-Vichniac attribution, 'Rule 90').
Usage: python3 collatz_oscillation_gilbreath_20260926_audit.py
       (writes 05-knowledge/results/collatz_oscillation_gilbreath_20260926_audit.out, LF)
"""
import math, os, re, sys, time, random, collections
from math import comb
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(ROOT, '05-knowledge', 'results', 'collatz_oscillation_gilbreath_20260926_audit.out')
NOTE = os.path.join(ROOT, '05-knowledge', 'results', 'collatz_oscillation_gilbreath_20260926.md')
lines = []


def P(s=''):
    print(s); lines.append(s)


def bar():
    P('=' * 100)


# ----------------------------------------------------------------------------------------------------------
# 1-2. Orbit statistics (exact integers)
# ----------------------------------------------------------------------------------------------------------
def T(n, q=3):
    return n // 2 if n % 2 == 0 else (q * n + 1) // 2


def orbit(n, steps, q=3):
    ys = [n]
    for _ in range(steps):
        n = T(n, q); ys.append(n)
        if n == 1:
            break
    return ys


def dippers(ys, L, D, exact=True):
    """THM-4506 setting: X = 2^L, k = L, depth D (integer). Index i is eligible iff i <= M-k and y_i <= X.
    Dipper iff y_(i+s) < y_i 2^-D for some 1 <= s <= k; landing index = least such i+s."""
    X = 1 << L; k = L; M = len(ys) - 1; thr = 2.0 ** (-D)
    dip = {}; nd = 0; tot = 0
    for i in range(0, M - k + 1):
        y = ys[i]
        if y > X:
            continue
        tot += 1
        land = None
        for s in range(1, k + 1):
            hit = ((ys[i + s] << D) < y) if exact else (ys[i + s] < y * thr)
            if hit:
                land = i + s; break
        if land is None:
            nd += 1
        else:
            dip.setdefault(land, []).append(i)
    return tot, nd, dip


def analyse(ys, L, D):
    X = 1 << L; k = L; M = len(ys) - 1
    tot, nd, dip = dippers(ys, L, D, exact=True)
    totf, ndf, dipf = dippers(ys, L, D, exact=False)
    Dt = sum(len(v) for v in dip.values()); lp = len(dip)
    same_float = (tot, nd, sorted((j, tuple(v)) for j, v in dip.items())) == (totf, ndf, sorted((j, tuple(v)) for j, v in dipf.items()))
    bl_script = sum(1 for i in range(0, M - k + 1) if (ys[i] << D) <= X)
    NX_all = sum(1 for y in ys if y <= X); NXD_all = sum(1 for y in ys if (y << D) <= X)
    ok_prop = ok_S = ok_O = ok_halve = ok_land = True
    runs_lt = runs_gt = 0; nondip_total = 0
    order = sorted(((len(v), j, v) for j, v in dip.items()), reverse=True)   # the script's ordering
    heavy = []
    for m, j, ds in order:
        yj = ys[j]; lo = yj << D; hi = yj << (D + 1)
        if not ((yj << D) < X):
            ok_land = False
        A = [i for i in range(max(0, j - k), j)
             if i <= M - k and ys[i] <= X and lo < ys[i] <= hi and all((ys[t] << D) >= ys[i] for t in range(i + 1, j))]
        if A != sorted(ds):
            ok_prop = False
        if any(not (lo < ys[i] <= hi) for i in ds):
            ok_S = False
        if ys[j - 1] != 2 * yj:
            ok_halve = False
        dss = sorted(ds)
        for a, b in zip(dss, dss[1:]):
            if not any(ys[t] % 2 == 1 for t in range(a, b)):
                ok_O = False
        runs_script = 0; inside = False
        for t in range(max(0, j - L), j):
            now = lo < ys[t] <= hi
            if now and not inside:
                runs_script += 1
            inside = now
        shell_idx = [t for t in range(max(0, j - L), j) if lo < ys[t] <= hi]
        nondip = len([t for t in shell_idx if t not in set(dss)])
        nondip_total += nondip
        runs_dip = sum(1 for idx, i in enumerate(dss) if idx == 0 or dss[idx - 1] != i - 1)
        if runs_script < m:
            runs_lt += 1
        if runs_script > m:
            runs_gt += 1
        if len(heavy) < 3:
            heavy.append((m, runs_script, runs_dip, len(shell_idx), nondip, j))
    return dict(tot=tot, nd=nd, Dt=Dt, lp=lp, mean=(Dt / lp if lp else 0.0), heavy=heavy,
                ratio_script=(tot / bl_script if bl_script else float('inf')), ratio_full=(NX_all / NXD_all if NXD_all else float('inf')),
                NX=NX_all, NXD=NXD_all, bl=bl_script, E=NX_all - tot, ok_prop=ok_prop, ok_S=ok_S, ok_O=ok_O, ok_halve=ok_halve,
                ok_land=ok_land, runs_lt=runs_lt, runs_gt=runs_gt, nondip=nondip_total, same_float=same_float, M=M)


EXPECTED = {  # from the committed collatz_oscillation_20260926_shells.out
    ('3n+1 record 63728127', 24): (150, 97, 53, 22, 2.41, '(6,5) (6,4) (6,4)', 1.39),
    ('3n+1 record 63728127', 32): (304, 195, 109, 33, 3.30, '(8,5) (7,5) (7,5)', 1.91),
    ('3n+1 record 63728127', 48): (545, 370, 175, 34, 5.15, '(12,8) (12,8) (12,10)', 1.00),
    ('3n+1 record 63728127', 64): (529, 356, 173, 30, 5.77, '(15,10) (13,9) (12,8)', 1.00),
    ('5n+1 from 7', 24): (80, 77, 3, 2, 1.50, '(2,2) (1,1)', 1.70),
    ('5n+1 from 7', 32): (142, 140, 2, 2, 1.00, '(1,1) (1,1)', 1.43),
    ('5n+1 from 7', 48): (293, 267, 26, 10, 2.60, '(6,6) (5,5) (4,4)', 1.47),
    ('5n+1 from 7', 64): (530, 460, 70, 20, 3.50, '(7,7) (6,6) (6,6)', 1.09),
    ('3n+1 from 2^L-1', 24): (142, 97, 45, 16, 2.81, '(7,5) (6,5) (5,3)', 1.21),
    ('3n+1 from 2^L-1', 32): (146, 65, 81, 27, 3.00, '(8,7) (7,5) (7,5)', 1.35),
    ('3n+1 from 2^L-1', 48): (165, 39, 126, 40, 3.15, '(10,7) (10,8) (8,7)', 1.67),
    ('3n+1 from 2^L-1', 64): (290, 51, 239, 65, 3.68, '(18,12) (16,13) (14,10)', 1.05),
}


def section_1_2():
    bar(); P('1. Proposition 1 (shell-visit reformulation) and 2. the probe table -- exact re-implementation')
    rec = orbit(63728127, 60 * 64)
    P('   record orbit of 63728127 under T: %d steps to 1, maximum %d = 2^%.2f (< 2^53: %s)' % (len(rec) - 1, max(rec), math.log2(max(rec)), max(rec) < 2 ** 53))
    P('   definitions used (THM-4506 section 0): X = 2^L, k = L, dipper i <= M-k with y_i <= X and y_(i+s) 2^D < y_i for some 1 <= s <= k,')
    P('   landing index = least such i+s; m(j) = #dippers landing at j; ND = eligible non-dippers; E = indices with y <= X and i > M-k.')
    P('   Prop 1(a) test: m(j) == #{ i in [j-k, j-1] : i <= M-k, y_i <= X, 2^D y_j < y_i <= 2^(D+1) y_j, and y_t 2^D >= y_i for all i < t < j }')
    P('   (a visit is an INDEX; "odd-separated" is automatic inside one dyadic shell, see below); runs = maximal runs of window indices in S_j')
    P('   (the script\'s "shell runs", which counts non-dipper visits too), runs_dip = maximal runs among the dippers themselves.')
    P('')
    P('   L  D  segment                 tot   ND  Dip  land  mean   heaviest [(m, runs_script, runs_dip, shell idx, non-dip, j)]      ratio_script ratio_full  N(X) N(X2^-D) E')
    mism = []
    allok = True
    for L in (24, 32, 48, 64):
        D = math.ceil(1.05 * math.log2(L))
        segs = [('3n+1 record 63728127', orbit(63728127, 60 * L)), ('5n+1 from 7', orbit(7, 60 * L, q=5)), ('3n+1 from 2^L-1', orbit(2 ** L - 1, 60 * L))]
        for name, ys in segs:
            r = analyse(ys, L, D)
            hv = ' '.join('(%d,%d,%d,%d,%d,%d)' % h for h in r['heavy'])
            P('  %2d %2d  %-22s %5d %4d %4d %5d %6.2f  %-58s %6.2f %6.2f  %5d %5d %3d' % (L, D, name, r['tot'], r['nd'], r['Dt'], r['lp'], r['mean'], hv, r['ratio_script'], r['ratio_full'], r['NX'], r['NXD'], r['E']))
            checks = 'Prop1a %s LemmaS %s LemmaO %s halving-into-j %s landing<X2^-D %s float==exact %s | landing points with runs_script<m: %d, >m: %d; non-dipper shell visits: %d' % (
                r['ok_prop'], r['ok_S'], r['ok_O'], r['ok_halve'], r['ok_land'], r['same_float'], r['runs_lt'], r['runs_gt'], r['nondip'])
            P('         ' + checks)
            allok &= r['ok_prop'] and r['ok_S'] and r['ok_O'] and r['ok_halve'] and r['ok_land'] and r['same_float']
            # partition accounting
            part_ok = (r['tot'] == r['nd'] + r['Dt']) and (r['E'] <= L) and (r['NX'] == r['E'] + r['nd'] + r['Dt'])
            allok &= part_ok
            exp = EXPECTED.get((name, L))
            if exp:
                got_hv = ' '.join('(%d,%d)' % (h[0], h[1]) for h in r['heavy'])
                got = (r['tot'], r['nd'], r['Dt'], r['lp'], round(r['mean'], 2), got_hv, round(r['ratio_script'], 2))
                ok = got == exp
                P('         vs shells.out %s: %s  partition N(X)=E+ND+Dip %s (E=%d<=k)' % (exp, 'MATCH' if ok else 'MISMATCH got %s' % (got,), part_ok, r['E']))
                if not ok:
                    mism.append((name, L, exp, got))
    P('')
    P('   All Prop 1(a)/Lemma S/Lemma O/halving/landing/partition checks passed: %s; table mismatches against shells.out: %d' % (allok, len(mism)))
    P('   Note on "visits": inside S_j = (2^D y_j, 2^(D+1) y_j] a halving leaves the shell (y/2 <= 2^D y_j) and two odd steps leave it')
    P('   ((3/2)^2 > 2), so any two indices of the orbit in S_j are automatically separated by an odd letter and a maximal run in S_j has')
    P('   at most 2 indices (THM-4506 Prop. C(a)). Hence m(j) = number of qualifying INDICES; the script\'s "shell runs" both undercount')
    P('   (a run holding 2 dippers) and overcount (window visits to S_j cut off by an earlier drop are not dippers of j).')
    P('   Regularity form: N(X) = #E + #Dip + #ND with #E <= k (0 for an infinite orbit). Averaged form => N(X) <= C L^mu N(X2^-D) + #ND + C L^2 + k;')
    P('   the note\'s "+1" (an extra N(X 2^-D)) and the absorbed k only weaken it. Converse: #Dip = N(X) - #ND - #E <= (C L^mu + 1) N(X2^-D) + C L^2,')
    P('   i.e. the averaged form with C -> C + 1. Both directions hold up to constants; the note\'s "both directions are the definitions" is right.')
    P('   Ratio definition: the script\'s N(X)/N(X 2^-D) restricts BOTH counts to indices with a full window (i <= M-k); the unrestricted ratio')
    P('   differs (e.g. 1.32 vs 1.39, 1.76 vs 1.91 on the record orbit at L = 24, 32). For the record orbit at L = 48, 64 every value is below')
    P('   X 2^-D (max 2^38.8 < 2^42), so ratio = 1.00 is vacuous there, not evidence of regularity.')
    return allok, mism


# ----------------------------------------------------------------------------------------------------------
# 3. Gilbreath
# ----------------------------------------------------------------------------------------------------------
def primes_below(Pn):
    s = np.ones(Pn, dtype=bool); s[:2] = False
    for i in range(2, int(Pn ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def run_lengths(row, value=0, need_right_border=True):
    """maximal runs of `value` in row; returns (starts, ends) with ends exclusive; keeps runs with start >= 1 and end < len(row)."""
    z = (row == value).astype(np.int8); zi = np.concatenate(([0], z, [0])); d = np.diff(zi)
    st = np.nonzero(d == 1)[0]; en = np.nonzero(d == -1)[0]
    ok = st >= 1
    if need_right_border:
        ok &= en < len(row)
    return st[ok], en[ok]


def script_detect(D):
    """the committed script's zero-triangle detector, verbatim logic (rows r0 = 1..R-2 of the kept block)."""
    R, C = D.shape; sizes = collections.Counter(); plain = collections.Counter(); plain_fit = collections.Counter()
    for r0 in range(1, R - 1):
        rowv = D[r0]; c = 1
        while c < C:
            if rowv[c] == 0:
                c1 = c
                while c1 < C and rowv[c1] == 0:
                    c1 += 1
                t = c1 - c
                left_ok = rowv[c - 1] in (2, 1); right_ok = c1 < C and rowv[c1] == 2
                if left_ok and right_ok:
                    plain[t] += 1
                    if r0 + t - 1 < R:
                        plain_fit[t] += 1
                tri = left_ok and right_ok
                for s in range(1, t):
                    if r0 + s >= R or not np.all(D[r0 + s][c:c1 - s] == 0) or D[r0 + s][c1 - s] != 2:
                        tri = False; break
                if tri and t >= 1:
                    sizes[t] += 1
                c = c1
            else:
                c += 1
    return sizes, plain, plain_fit


def tops_detect(D):
    """true triangle sides: a maximal zero run [c, c1) at row r0 (bordered by 2/1 left, 2 right) whose row above holds 2s at c..c1
    (a run of t+1 equal NON-zero entries); a run below t+1 zeros is the interior of a larger triangle."""
    R, C = D.shape; sizes = collections.Counter()
    for r0 in range(1, R):
        rowv = D[r0]; prev = D[r0 - 1]; c = 1
        while c < C:
            if rowv[c] == 0:
                c1 = c
                while c1 < C and rowv[c1] == 0:
                    c1 += 1
                t = c1 - c
                if rowv[c - 1] in (2, 1) and c1 < C and rowv[c1] == 2 and prev[c] == 2:
                    sizes[t] += 1
                c = c1
            else:
                c += 1
    return sizes


def section_3():
    bar(); P('3. Gilbreath\'s difference triangle')
    # (a) algebra
    P('   (a) rule algebra, new[i] = |row[i] - row[i+1]|:')
    xor_ok = all(abs(a - b) == 2 * ((a // 2) ^ (b // 2)) for a in (0, 2) for b in (0, 2))
    one_ok = all(abs(1 - a) == 1 for a in (0, 2))
    P('       |a-b| = 2*(a/2 XOR b/2) on {0,2}: %s; |1-a| = 1 for a in {0,2}: %s' % (xor_ok, one_ok))
    rng = np.random.default_rng(1)
    front_ok = stat_ok = True; trail_sea = 0; n_tr = 0
    for _ in range(20000):
        row = rng.integers(0, 2, size=12) * 2; p = 6; d = int(rng.integers(2, 8)) * 2; row[p] = d
        new = np.abs(np.diff(row))
        front_ok &= new[p - 1] == d - row[p - 1]      # the LEFT neighbour shrinks the front
        stat_ok &= new[p] == d - row[p + 1]           # the stationary copy at p shrinks by the RIGHT neighbour
        new2 = np.abs(np.diff(new)); n_tr += 1
        trail_sea += int(new2[p - 1] in (0, 2))       # between the two copies a sea value re-appears
    P('       defect d >= 4 at p with 0/2 neighbours: new[p-1] = d - row[p-1] (left neighbour) in all 20000 trials: %s;' % front_ok)
    P('       ALSO new[p] = d - row[p+1] (a stationary copy, right neighbour) in all trials: %s; two rows later cell p-1 is back in {0,2}: %d/%d.' % (stat_ok, trail_sea, n_tr))
    P('       So "moves one cell to the left per row, shrinks by 2 when the cell to its left is 2" is exact for the leftmost FRONT only;')
    P('       the seed also leaves a stationary copy that can re-emit fronts (a light cone with left speed 1 and right speed 0).')
    # (b) the primes' triangle
    t0 = time.time()
    Pn = 200000; row = primes_below(Pn); n0 = len(row)
    lead = []; frontier = []; maxes = []; KEEP = 160; diagram = []; full_first = []
    runhist = collections.Counter(); tophist = collections.Counter(); rows_hist = 0; rstar = None; r = 0
    while len(row) > 1:
        row = np.abs(np.diff(row)); r += 1
        lead.append(int(row[0]))
        big = np.nonzero(row >= 4)[0]; F = int(big[0]) if len(big) else len(row)
        frontier.append(F); mx = int(row.max()); maxes.append(mx)
        if r <= 6000:
            diagram.append(row[:KEEP].copy())
        if r <= 70:
            full_first.append(row.copy())
        if mx <= 2 and rstar is None:
            rstar = r
        if mx <= 2 and len(row) >= 3:
            st, en = run_lengths(row, 0)
            for t, c in zip(*np.unique(en - st, return_counts=True)):
                runhist[int(t)] += int(c)
            st2, en2 = run_lengths(row, 2)
            m = en2 - st2; m = m[m >= 2]
            for t, c in zip(*np.unique(m - 1, return_counts=True)):
                tophist[int(t)] += int(c)
            rows_hist += 1
    fr = np.array(frontier); mxs = np.array(maxes)
    P('   (b) primes below %d: %d primes, %d rows computed in %.1f s; leading entry == 1 in every row: %s' % (Pn, n0, r, time.time() - t0, all(v == 1 for v in lead)))
    P('       row:     ' + ' '.join('%7d' % rr for rr in (1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000, 17000)))
    P('       F(r):    ' + ' '.join('%7d' % fr[rr - 1] for rr in (1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000, 17000)))
    P('       max:     ' + ' '.join('%7d' % mxs[rr - 1] for rr in (1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000, 17000)))
    P('       note quotes F = 3, 8, 25, 59 (rows 1, 2, 5, 10), 2763 (row 50), then F = row length: %s' % (list(fr[[0, 1, 4, 9, 49]]) == [3, 8, 25, 59, 2763] and fr[99] == n0 - 100))
    P('       FIRST row with maximum entry <= 2: r* = %d (note: "about row 100"); rows with an entry >= 4: %d (rows 1..%d); max over rows >= r*: %d;' % (rstar, int((mxs >= 4).sum()), int(np.nonzero(mxs >= 4)[0][-1]) + 1, int(mxs[rstar - 1:].max())))
    P('       entries >= 4 re-appearing after r*: %d (closure of {1} x {0,2} under the rule makes this automatic).' % int((mxs[rstar - 1:] >= 4).sum()))
    P('       frontier sequence r: F/value/cell-to-its-left for rows 1..%d (value = first entry >= 4; "left" = the 0/2 cell at F-1, which decides' % (rstar - 1))
    P('       the next row\'s front: it continues at F-1 with value - left iff value - left >= 4):')
    items = []
    for rr in range(1, rstar):
        Fv = int(fr[rr - 1]); rv = full_first[rr - 1]
        items.append('%d:%d/%d/%d' % (rr, Fv, int(rv[Fv]), int(rv[Fv - 1])))
    for i in range(0, len(items), 8):
        P('         ' + '  '.join(items[i:i + 8]))
    same_col = sum(1 for rr in range(2, rstar) if fr[rr - 1] == fr[rr - 2])
    cont = sum(1 for rr in range(2, rstar) if fr[rr - 1] == fr[rr - 2] - 1)
    P('       rows whose frontier equals the previous row\'s (a stationary copy re-emitting a dead front): %d; rows continuing a front (F-1): %d.' % (same_col, cont))
    # script-region histogram
    D = np.array([np.pad(d, (0, KEEP - len(d)), constant_values=-1) for d in diagram[:6000]])
    sizes, plain, plain_fit = script_detect(D)
    ks = sorted(sizes)
    P('       script-region zero-"triangle" histogram (rows 2..5999 of the kept block, columns < 160), script\'s detector re-run:')
    P('         ' + '  '.join('%d:%d' % (t, sizes[t]) for t in ks))
    P('         note quotes 118983, 59228, 29388: %s' % ([sizes[1], sizes[2], sizes[3]] == [118983, 59228, 29388]))
    P('         plain maximal zero runs bordered by 2/1 and 2 (no triangle check): ' + '  '.join('%d:%d' % (t, plain[t]) for t in sorted(plain)))
    P('         same, restricted to runs whose triangle fits inside the kept block (r0 + t - 1 < R): identical to the detector: %s' % (plain_fit == sizes))
    P('         => the "triangle check" is vacuous except at the block\'s bottom edge: a zero run bordered by 2s forces the shrinking rows below.')
    P('         The detector counts a zero run in EVERY row, so a triangle of side t is counted once per row as sizes t, t-1, ..., 1.')
    P('       full sea (all %d rows >= r*, all columns): maximal zero runs by length:' % rows_hist)
    kr = sorted(runhist)
    P('         ' + '  '.join('%d:%d' % (t, runhist[t]) for t in kr[:24]))
    P('         ratios count(t)/count(t+1): ' + ' '.join('%.3f' % (runhist[t] / runhist[t + 1]) for t in kr[:18] if runhist.get(t + 1)))
    kt = sorted(tophist)
    P('       full sea, TRUE triangle sides (tops: a maximal run of t+1 twos not touching the right end -> side t):')
    P('         ' + '  '.join('%d:%d' % (t, tophist[t]) for t in kt[:24]))
    P('         ratios: ' + ' '.join('%.3f' % (tophist[t] / tophist[t + 1]) for t in kt[:18] if tophist.get(t + 1)))
    P('         largest side present: %d (the note\'s "sizes 1..15" is the 160-column window; the geometric law is the substantive point).' % max(kt))
    # defect samples, script-style and de-duplicated
    surv = []
    for rr in range(1, min(len(diagram), 3000)):
        d = diagram[rr - 1]; big = np.nonzero(d >= 4)[0]
        if len(big) == 0 or big[0] >= KEEP - 1:
            continue
        F = int(big[0]); size = int(d[F]); s_ = 0
        while s_ < F and rr - 1 + s_ + 1 < len(diagram):
            nxt = diagram[rr - 1 + s_ + 1]; c = F - s_ - 1
            if c < 0 or c >= len(nxt) or nxt[c] < 4:
                break
            s_ += 1
        surv.append((rr, size, F, s_))
    by = collections.defaultdict(list)
    for rr, size, F, s_ in surv:
        by[size].append((F, s_))
    P('       defect samples as the script defines them (rows < 3000, F < 159): (row, size, F, travel) = %s' % surv)
    for size in sorted(by):
        Ls = by[size]
        P('         size %d: count %d, mean F %.1f, mean travel %.2f, max %d   (note: size 4 -> 11, 48.6, 0.18, 1; size 8 -> max 3)' % (size, len(Ls), np.mean([f for f, _ in Ls]), np.mean([t for _, t in Ls]), max(t for _, t in Ls)))
    samp = set((rr, F) for rr, _, F, _ in surv)
    fresh = [(rr, size, F, s_) for rr, size, F, s_ in surv if (rr - 1, F + 1) not in samp]
    P('       de-duplicated (a sample is dropped when it is the previous row\'s frontier defect one column further left): %s' % fresh)
    f4 = [s_ for rr, size, F, s_ in fresh if size == 4]
    P('         fresh fronts: %d, of size 4: %d with travels %s (mean %.2f; random model: P(>= 1 row) = 1/2, mean 1.0); the size-8 front of row 5' % (len(fresh), len(f4), f4, np.mean(f4) if f4 else 0))
    P('         (F = 25) accounts for the "size 8 x2", "size 6 x1" and one "size 4" sample (rows 5-8: 8, 8, 6, 4). Rows 3-4 (F = 14) and 11-14 (F = 98, 97,')
    P('         98, 97) are one stationary copy re-emitting fronts. The "0.18 rows" statistic double counts; the qualitative claim (fronts die in')
    P('         0-3 rows, the longest observed travel is 3 rows) is reproduced.')
    # (c) single seed vs random row
    Rr = 63; W = 2 * Rr + 3
    row = np.zeros(W, dtype=np.int64); row[0] = 1; row[Rr + 1] = 2
    rows = []
    for s in range(Rr + 1):
        row = np.abs(np.diff(row)); rows.append(row.copy())
    KEEP2 = W - Rr - 1
    Ds = np.array([np.pad(d[:KEEP2], (0, KEEP2 - len(d[:KEEP2])), constant_values=-1) for d in rows])
    sz, _, _ = script_detect(Ds); tp = tops_detect(Ds)
    P('   (c) single seed (one 2 in a zero row behind a leading 1, %d rows = Pascal mod 2):' % (Rr + 1))
    P('         the script\'s detector reports sizes %s' % sorted(sz))
    P('         true triangle sides (tops): %s -> only 2^k - 1: %s' % (sorted(tp.items()), all((t + 1) & t == 0 for t in tp)))
    P('         => the committed statistic cannot separate the Sierpinski case from the primes\' rows (from a single seed it also shows every size,')
    P('            in dyadic blocks with equal counts); the conclusion "sizes geometric, not 2^k-1" is nevertheless right for TRUE sides (tops above).')
    rng = np.random.default_rng(20260926)
    rnd = rng.integers(0, 2, size=20000, dtype=np.int64) * 2
    h2 = collections.Counter(); t2 = collections.Counter(); rowr = rnd
    for s in range(3000):
        rowr = np.abs(np.diff(rowr))
        st, en = run_lengths(rowr, 0)
        for t, c in zip(*np.unique(en - st, return_counts=True)):
            h2[int(t)] += int(c)
        st2, en2 = run_lengths(rowr, 2); m = en2 - st2; m = m[m >= 2]
        for t, c in zip(*np.unique(m - 1, return_counts=True)):
            t2[int(t)] += int(c)
    k2 = sorted(h2)
    P('       random 0/2 row (20000 cells, 3000 rows): runs ' + '  '.join('%d:%d' % (t, h2[t]) for t in k2[:16]))
    P('         run ratios: ' + ' '.join('%.3f' % (h2[t] / h2[t + 1]) for t in k2[:12]) + ' ; side (tops) ratios: ' + ' '.join('%.3f' % (t2[t] / t2[t + 1]) for t in sorted(t2)[:12]))
    P('         => every size occurs with geometric frequency (ratio 2), for runs and for true sides alike; the primes\' sea matches this.')
    # (d) survival Monte Carlo
    P('   (d) random model (i.i.d. 0/2 sea, leading 1, one defect 2j at column F): P(leading 1 ever lost) by Monte Carlo (100000 trials),')
    P('       P(single front reaches column 1 with value >= 4), the note\'s formula sum_(i<j) C(F,i) 2^-F, and the exact single-front value')
    P('       sum_(t<=j-2) C(F-1,t) 2^-(F-1) (the front meets F-1 sea cells on its way from column F to column 1 and must keep 2j - 2t >= 4):')
    rng = np.random.default_rng(7)
    rows_out = []
    for j, F in [(2, 3), (2, 4), (2, 6), (2, 10), (3, 4), (3, 6), (3, 10), (4, 8), (4, 12), (5, 10), (5, 16)]:
        trials = 100000; width = 4 * F + 4 * j + 8
        rows = rng.integers(0, 2, size=(trials, width), dtype=np.int64) * 2
        rows[:, 0] = 1; rows[:, F] = 2 * j
        viol = np.zeros(trials, dtype=bool); front_ok = np.ones(trials, dtype=bool); fval = np.full(trials, 2 * j, dtype=np.int64)
        cur = rows
        for s in range(1, width - 1):
            if s <= F - 1:
                leftn = cur[:, F - (s - 1) - 1]
                fval = np.where(front_ok, fval - leftn, fval); front_ok &= fval >= 4
            cur = np.abs(np.diff(cur, axis=1))
            viol |= cur[:, 0] != 1
            if s >= F - 1 and not (cur[:, 1:] >= 4).any():
                break
        unresolved = int((cur[:, 1:] >= 4).any(axis=1).sum())
        A = sum(comb(F, i) for i in range(j)) / 2 ** F
        B = sum(comb(F - 1, t) for t in range(j - 1)) / 2 ** (F - 1)
        rows_out.append((j, F, viol.mean(), front_ok.mean(), A, B, unresolved))
        P('         (j=%d, F=%2d): P(any) %.5f  P(front) %.5f  note %.5f  exact-front %.5f  unresolved %d' % rows_out[-1])
    P('       => the note\'s tail is off by one in both j and F (it overstates the single-front probability by a factor 2-5 for small j), and the')
    P('          stationary copy adds a few per cent (j >= 3) through re-emitted fronts. The threshold j/F -> 0 and the exponential-in-F decay are right.')
    P('          Meeting a 2 always shrinks the front by exactly 2 while its left neighbour is in {0,2}; the leftmost defect never meets a 4 on its left')
    P('          (everything left of the frontier is 0/2 by definition), so the "4 or larger" case only concerns collisions between two defects.')
    return rstar, sizes


# ----------------------------------------------------------------------------------------------------------
# 4. Base-6 locality
# ----------------------------------------------------------------------------------------------------------
def digits6(x, w):
    d = []
    for _ in range(w):
        d.append(x % 6); x //= 6
    return d


def local(par, a, b, c):
    """digit i >= 1 of T(x) from x_(i-1) = a, x_i = b, x_(i+1) = c and the parity bit."""
    if par == 0:
        return (b + 6 * (c % 2)) // 2
    zi = (3 * b + a // 2) % 6          # digit i of 3x+1: carry into digit i is floor(3 x_(i-1)/6) = floor(x_(i-1)/2)
    zi1 = (3 * c + b // 2) % 6         # digit i+1 of 3x+1
    return (zi + 6 * (zi1 % 2)) // 2   # halving reads one digit up


def local0(par, b, c):
    if par == 0:
        return (b + 6 * (c % 2)) // 2
    zi1 = (3 * c + b // 2) % 6
    return (4 + 6 * (zi1 % 2)) // 2     # 3 x_0 + 1 = 4, 10, 16 for x_0 = 1, 3, 5: digit 4, carry floor(x_0/2), no extra carry


def section_4():
    bar(); P('4. Base-6 locality (section 3 of the note)')
    P('   carry table floor((3d + c)/6), d = 0..5, c = 0, 1, 2: %s -> independent of c: %s' % ([[(3 * d + c) // 6 for c in range(3)] for d in range(6)], all(len(set((3 * d + c) // 6 for c in range(3))) == 1 for d in range(6))))
    P('   maximal carry: floor((15 + 2)/6) = 2, and with the +1 at digit 0: floor((3 x_0 + 1)/6) = %s for x_0 = 1, 3, 5 versus floor(3 x_0/6) = %s; 3 x_0 + 1 mod 6 = %s (always 4)' % ([(3 * d + 1) // 6 for d in (1, 3, 5)], [(3 * d) // 6 for d in (1, 3, 5)], [(3 * d + 1) % 6 for d in (1, 3, 5)]))
    bad = 0; n = 0; Wd = 9
    for x in range(1, 6 ** 7):
        dx = digits6(x, Wd + 1); dt = digits6(T(x), Wd + 1); par = x % 2
        if local0(par, dx[0], dx[1]) != dt[0]:
            bad += 1
        for i in range(1, Wd):
            if local(par, dx[i - 1], dx[i], dx[i + 1]) != dt[i]:
                bad += 1
        n += 1
    P('   explicit radius-1 rule digit_i(T x) = f(parity(x_0), x_(i-1), x_i, x_(i+1)) checked EXHAUSTIVELY for x < 6^7 (%d values, digits 0..%d): mismatches %d' % (n, Wd - 1, bad))
    random.seed(3); bad2 = 0
    for _ in range(20000):
        x = random.randrange(1, 6 ** 40); dx = digits6(x, 42); dt = digits6(T(x), 42); par = x % 2
        if local0(par, dx[0], dx[1]) != dt[0]:
            bad2 += 1
        for i in range(1, 41):
            if local(par, dx[i - 1], dx[i], dx[i + 1]) != dt[i]:
                bad2 += 1
    P('   random 40-digit x (20000 values, all 41 digits): mismatches %d' % bad2)
    depL = any(local(1, a, b, c) != local(1, a2, b, c) for a in range(6) for a2 in range(6) for b in range(6) for c in range(6))
    depR = any(local(1, a, b, c) != local(1, a, b, c2) for a in range(6) for b in range(6) for c in range(6) for c2 in range(6))
    depE = any(local(0, a, b, c) != local(0, a, b, c2) for a in range(6) for b in range(6) for c in range(6) for c2 in range(6))
    P('   odd branch depends on x_(i-1): %s, on x_(i+1): %s; even branch (halving) reads x_(i+1) mod 2: %s -> radius exactly 1 on both branches; the parity' % (depL, depR, depE))
    P('   bit x_0 mod 2 is the one global datum (6 is even, so x mod 2 = x_0 mod 2). The note\'s rule table (x_(i-1) = 0) reproduced: %s' % (
        [[local(1, 0, b, c) for c in range(6)] for b in range(6)] == [[0, 3, 0, 3, 0, 3], [1, 4, 1, 4, 1, 4], [3, 0, 3, 0, 3, 0], [4, 1, 4, 1, 4, 1], [0, 3, 0, 3, 0, 3], [1, 4, 1, 4, 1, 4]]))
    P('   Citation: Cloney, Goles, Vichniac, "The 3x+1 problem: a quasi cellular automaton", Complex Systems 1 (1987) 349-360 -- verified online')
    P('   (complex-systems.com abstract): it studies the iterates "in base two as a quasi cellular automaton". It is NOT a base-6 paper: the note')
    P('   attributes the base-6 radius-1 statement to it, which is wrong (base 6 is where the map becomes a genuine local rule; the standard')
    P('   references for base-6 locality are Korec 1992 (generalized Pascal triangles) and Kari 2012 (multiplication by 3/2 in base 6)).')
    return bad + bad2


# ----------------------------------------------------------------------------------------------------------
# 5. Thresholds
# ----------------------------------------------------------------------------------------------------------
def section_5():
    bar(); P('5. Thresholds and THM-4487\'s curve (section 4 of the note)')
    h = lambda p: 0.0 if p in (0, 1) else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
    alpha = math.log2(3); E = lambda g: h(max(0.5, g / alpha))
    gs = [i / 1000 for i in range(1, 1001)]
    P('   E(gamma) = h(max(1/2, gamma/alpha)) on (0, 1]: minimum %.6f at gamma = %.3f, maximum %.6f (flat = 1 for gamma <= log_4 3 = %.4f).' % (min(E(g) for g in gs), min(gs, key=E), max(E(g) for g in gs), alpha / 2))
    P('   THM-4476\'s depth theta enters as rho = (1 - theta)/alpha: theta -> 0 gives rho = %.4f, h = %.6f = h* = E(1); theta = 1 - alpha/2 = %.4f gives rho = 1/2, h = 1.' % (1 / alpha, h(1 / alpha), 1 - alpha / 2))
    P('   So "the theta -> 0 end" IS the "gamma = 1 end" (the same point, the Collatz end); the curve never leaves [h*, 1] = [0.95, 1] and has no')
    P('   entropy-0 end. The Gilbreath tail sum_(i<j) C(F,i) = 2^(F h(j/F) + o(F)) with j/F -> 0 lives at the p -> 0 end of the binary entropy')
    P('   FUNCTION h(p) (a lower tail below 1/2), which is not part of THM-4487\'s curve (an upper-tail count with p in [1/2, log_3 2]).')
    P('   The sentence "Gilbreath is the theta -> 0 end of the entropy curve of THM-4487, Collatz the gamma = 1 end" is therefore wrong as written;')
    P('   the honest statement is: both heuristics are one-sided binomial tails, at p -> 0 (Gilbreath) and at p = log_3 2 (Collatz), of h(p).')
    for j, F in [(2, 10), (2, 25), (2, 50), (3, 25), (3, 50), (4, 50), (2, 100)]:
        tail = sum(comb(F, i) for i in range(j)) / 2 ** F
        P('   Gilbreath tail j=%d F=%3d: sum_(i<j) C(F,i) 2^-F = %.3e = %8.0f x 2^-F  (polynomial factor ~ F^(j-1)/(j-1)!, not 1);  2^(-F(1-h(j/F))) = %.3e' % (j, F, tail, tail * 2 ** F, 2 ** (-F * (1 - h(j / F)))))
    hs = h(math.log(2, 3))
    P('   Collatz: h* = %.6f, W_k/2^k = Theta(2^(-(1-h*)k) k^(-3/2)) with 1 - h* = %.6f (THM-4495): "bad words die like 2^(-0.05k)" is right.' % (hs, 1 - hs))


# ----------------------------------------------------------------------------------------------------------
# 6. Text checks
# ----------------------------------------------------------------------------------------------------------
def section_6(rstar):
    bar(); P('6. Text checks on the note and scripts')
    txt = open(NOTE, encoding='utf-8').read()
    status = ' '.join(l.strip() for l in txt.split('\n')[2:7])
    P('   status line: %s' % status[:330])
    P('   "PROVED" occurrences: %d; "SPECULATION"/"Speculation" markers: %d; "FINITE-EXACT": %d; "No Collatz claim": %d' % (txt.count('PROVED'), txt.count('SPECULATION') + txt.count('Speculation'), txt.count('FINITE-EXACT'), txt.count('No Collatz claim')))
    P('   "about row 100" present: %s (exact first all-0/2 row is %d); "Cloney, Goles, Vichniac 1987" present: %s (a base-two "quasi" CA paper, see 4.);' % ('about row `100`' in txt or 'about row 100' in txt, rstar, 'Cloney, Goles, Vichniac 1987' in txt))
    P('   "Odlyzko 1993" present: %s (Iterated absolute values of differences of consecutive primes, Math. Comp. 61 (1993) 373-380; the persistence' % ('Odlyzko 1993' in txt))
    P('   argument is his and is cited); "theta -> 0" sentence present: %s; "precise sense" claim present: %s.' % ('theta -> 0' in txt, 'precise sense' in txt))
    gs = open(os.path.join(HERE, 'gilbreath_20260926_ca.py'), encoding='utf-8').read()
    P('   gilbreath script docstring says "Rule 90": %s -- new[i] = row[i] XOR row[i+1] is Wolfram rule 102 (or 60), not 90 (a_(i-1) XOR a_(i+1)); cosmetic.' % ('Rule 90' in gs))
    P('   "the conjecture is trivially true for this range beyond that row": correct from row %d on (closure of {1} x {0,2}); "from about row 100" understates it.' % rstar)
    P('   "So on actual segments the shell form holds with mu = 0 and a small constant": a finite observation on orbits that reach 1 (the landing note')
    P('   says such orbits "illustrate but do not test the hypothesis for divergent orbits"); acceptable under the FINITE-EXACT label, not as evidence.')
    P('   "mean multiplicity 2.4-5.8": the 3n+1 rows of shells.out span 2.29-5.77 and the 5n+1 rows 1.00-4.00.')


def main():
    t0 = time.time()
    P('collatz_oscillation_gilbreath_20260926_audit.py -- independent audit of collatz_oscillation_gilbreath_20260926.md (2026-09-26)')
    allok, mism = section_1_2()
    rstar, sizes = section_3()
    bad6 = section_4()
    section_5()
    section_6(rstar)
    bar()
    P('SUMMARY')
    P('  1(a) CONFIRMED with a precision GAP: m(j) equals the count of shell INDICES in the window whose first 2^D-drop is j (exact on every landing point);')
    P('      "odd-separated visits" is automatic inside one shell, and the script\'s "shell runs" is a different count (both under- and over-counts).')
    P('  1(b) CONFIRMED: averaged <=> shell <=> regularity, up to constants (the "+1" and the absorbed k are harmless; converse holds with C -> C+1).')
    P('  2    CONFIRMED: every entry of shells.out reproduced exactly (%s mismatches); float dip test == exact on all cells; the ratio is restricted to' % len(mism))
    P('      full-window indices (differs from the unrestricted N(X)/N(X2^-D)) and is vacuously 1.00 for the record orbit at L = 48, 64.')
    P('  3(a) CONFIRMED for the leftmost front (left neighbour shrinks it); GAP: the seed also leaves a stationary copy (right speed 0) that re-emits fronts.')
    P('  3(b) CONFIRMED numbers (F table, 118983/59228/29388, defect samples); ERROR of detail: the first all-0/2 row is 65, not "about 100";')
    P('      GAP: the "triangle" histogram is a per-row zero-run histogram (check vacuous), and the defect statistic double counts one front per row.')
    P('  3(c) GAP: the script\'s detector shows all sizes from a single seed too; the conclusion (geometric, not 2^k-1) holds for TRUE sides (tops).')
    P('  3(d) GAP: the tail is off by one in j and in F (factor 2-5 for small j) and ignores re-emitted fronts (a few per cent); threshold -> 0 right.')
    P('  4    CONFIRMED (explicit rule, exhaustive to 6^7, %d mismatches); ERROR in attribution: Cloney-Goles-Vichniac 1987 is a base-two quasi-CA paper.' % bad6)
    P('  5    ERROR: theta -> 0 is the gamma = 1 end of THM-4487\'s curve, which lives in [h*, 1]; the entropy-0 end is not on that curve.')
    P('  6    Status line honest (PROVED = classical restatements, cited); "trivially true beyond that row" CONFIRMED from row 65; speculation marked.')
    P('  runtime %.1f s' % (time.time() - t0))
    with open(OUT, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    print('written', OUT)


if __name__ == '__main__':
    main()
