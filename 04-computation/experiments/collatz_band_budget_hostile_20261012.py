"""The paid dyadic-band exit hypothesis: least counterexample, glide-record
witnesses, the random-model threshold coefficient, and source 187 rules.

Run from the repository root with python3 and again with python3 -O.
Dependencies: numpy and numba (exhaustive scan), the maintained exit compiler
paid_exit_boundaries_20261012.py (and through it the charged-family engine).

Definitions (inherited).  For odd n>1 iterate U(x)=(3x+1)/2^a with letters
a=v2(3x+1).  The first dyadic half-band exit is the first state h with
2(h+1)<=n+1.  charge = number of odd letters in that prefix (last letter
included).  E(n)=floor(log2(n+1)).  Hypothesis (3) of
three_problems_marked_defects_20261012a.md section 3:

    charge <= 3 [E(n)^2 - E(h)^2].                       (3)
"""
from pathlib import Path
import importlib.util
import json
import math
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def E(n):
    return (n+1).bit_length()-1


def tz(y):
    return (y & -y).bit_length()-1


def first_band_exit(n):
    """Exact (big-integer) first half-band exit: (odd steps, charge, h, word)."""
    check(n > 1 and n % 2 == 1, 'odd source above one')
    x, word = n, []
    while True:
        y = 3*x+1
        a = tz(y)
        x = y >> a
        word.append(a)
        if 2*(x+1) <= n+1:
            return len(word), sum(b % 2 for b in word), x, word


def glide(n):
    """Classical glide: Collatz steps (3x+1, x/2) until the first value below n."""
    x, steps = n, 0
    while x >= n:
        x = x//2 if x % 2 == 0 else 3*x+1
        steps += 1
    return steps


# Glide records, fetched 2026-10-12 from Eric Roosendaal's table
# (www.ericr.nl/wondrous/glidrecs.html; records found by Tomas Oliveira e Silva).
# Each glide value is recomputed below; the table is used only as a candidate list.
GLIDE_RECORDS = [(3, 6), (7, 11), (27, 96), (703, 132), (10087, 171), (35655, 220),
                 (270271, 267), (362343, 269), (381727, 282), (626331, 287),
                 (1027431, 298), (1126015, 365), (8088063, 401), (13421671, 468),
                 (20638335, 476), (26716671, 486), (56924955, 502), (63728127, 613),
                 (217740015, 644), (1200991791, 649), (1827397567, 706),
                 (2788008987, 729), (12235060455, 892), (898696369947, 897),
                 (2081751768559, 988), (13179928405231, 1122), (31835572457967, 1161),
                 (70665924117439, 1177), (739448869367967, 1187),
                 (1008932249296231, 1445), (118303688851791519, 1471),
                 (180352746940718527, 1575), (1236472189813512351, 1614),
                 (2602714556700227743, 1639)]


def record_witnesses():
    rows, fails = [], {3: [], 4: [], 5: []}
    for n, g in GLIDE_RECORDS:
        check(glide(n) == g, 'fetched glide value recomputed exactly')
        steps, charge, h, word = first_band_exit(n)
        en, eh = E(n), E(h)
        d = en*en-eh*eh
        row = dict(n=n, glide=g, odd_steps=steps, charge=charge, h=h, E_n=en, E_h=eh,
                   drop=en-eh, budget3=3*d, ratio=charge/d, charge_per_2E_minus_1=charge/(2*en-1),
                   fails3=charge > 3*d)
        rows.append(row)
        for c in fails:
            if charge > c*d:
                fails[c].append(n)
    return rows, fails


def least_counterexample_scan(limit, nchunks=480):
    """Exhaustive numba scan of odd n<limit (int64 with overflow guard)."""
    import numpy as np
    from numba import njit, prange

    @njit(cache=True)
    def ntz(y):
        c = 0
        while (y & 1) == 0:
            y >>= 1
            c += 1
        return c

    @njit(cache=True)
    def nE(n):
        m = n+1
        e = -1
        while m:
            m >>= 1
            e += 1
        return e

    @njit(parallel=True, cache=True)
    def scan(lo, hi, nchunks):
        best_num = np.zeros(nchunks, np.int64); best_den = np.ones(nchunks, np.int64); best_n = np.zeros(nchunks, np.int64)
        fails = np.zeros(nchunks, np.int64); first_fail = np.zeros(nchunks, np.int64); flagged = np.zeros(nchunks, np.int64)
        fails4 = np.zeros(nchunks, np.int64); flagged_n = np.zeros((nchunks, 64), np.int64)
        size = (hi-lo)//nchunks+1
        for c in prange(nchunks):
            a0 = lo+c*size
            a1 = min(hi, a0+size)
            if a0 % 2 == 0:
                a0 += 1
            for n in range(a0, a1, 2):
                x = n; charge = 0; lim = n+1; ok = True
                while True:
                    y = 3*x+1
                    if y > (1 << 61):
                        ok = False
                        break
                    a = ntz(y)
                    x = y >> a
                    charge += a & 1
                    if 2*(x+1) <= lim:
                        break
                if not ok:
                    if flagged[c] < 64:
                        flagged_n[c, flagged[c]] = n
                    flagged[c] += 1
                    continue
                en = nE(n); eh = nE(x)
                d = en*en-eh*eh
                if charge > 3*d:
                    fails[c] += 1
                    if first_fail[c] == 0:
                        first_fail[c] = n
                if charge > 4*d:
                    fails4[c] += 1
                if charge*best_den[c] > best_num[c]*d:
                    best_num[c] = charge; best_den[c] = d; best_n[c] = n
        return best_num, best_den, best_n, fails, first_fail, flagged, fails4, flagged_n

    bn, bd, bnn, f, ff, fl, f4, fln = scan(3, limit, nchunks)
    i = max(range(nchunks), key=lambda k: bn[k]/bd[k])
    firsts = [int(x) for x in ff if x]
    flagged = sorted(int(x) for x in fln.ravel() if x)
    check(len(flagged) == int(fl.sum()), 'every overflow-flagged source is stored for exact recheck')
    rechecked = []
    for n in flagged:  # exact big-integer replay of the orbits that left int64
        steps, charge, h, _ = first_band_exit(n)
        d = E(n)**2-E(h)**2
        rechecked.append(dict(n=n, odd_steps=steps, charge=charge, E_n=E(n), E_h=E(h), ratio=charge/d, fails3=charge > 3*d))
    return dict(limit=limit, best_ratio=[int(bn[i]), int(bd[i]), int(bnn[i])],
                failures_coefficient3=int(f.sum()), failures_coefficient4=int(f4.sum()),
                least_failure=min(firsts) if firsts else None, overflow_flagged=int(fl.sum()),
                overflow_rechecked=rechecked, overflow_failures3=[r['n'] for r in rechecked if r['fails3']])


def threshold_coefficient():
    """Random model: letters i.i.d. P(a)=2^-a; glide of length k needs mean letter <= log2 3.

    Expected number of odd n<=N with a prefix of length k and empirical letter law nu
    is about N exp(-k KL(nu||P)); the charge of such a prefix is f(nu) k with f the
    odd-letter fraction.  Maximizing f k subject to N exp(-k KL)=1 gives
        charge_max ~ (f/KL) ln N = (f/KL)(ln 2) E(n),
    and the s=1 budget 3(2E-1) ~ 6E; so the coefficient needed in (3) tends to
        c* = (ln 2 / 2) max_nu f(nu)/KL(nu||P)  over mean(nu) <= log2 3.
    Computed by a two-parameter exponential tilt (mean and parity multipliers).
    """
    A = list(range(1, 81))
    P = [2.0**-a for a in A]
    L3 = math.log2(3)

    def law(t, u):
        w = [p*math.exp(t*a+u*(a % 2)) for p, a in zip(P, A)]
        z = sum(w)
        return [x/z for x in w]

    def stats(nu):
        f = sum(x for x, a in zip(nu, A) if a % 2)
        m = sum(x*a for x, a in zip(nu, A))
        kl = sum(x*math.log(x/p) for x, p in zip(nu, P) if x > 0)
        return f, m, kl

    def tilt_for_mean(u):
        lo, hi = -6.0, 1.0
        for _ in range(100):
            mid = (lo+hi)/2
            if stats(law(mid, u))[1] > L3:
                hi = mid
            else:
                lo = mid
        return lo

    best = None
    u = -1.0
    while u <= 3.0:
        t = tilt_for_mean(u)
        f, m, kl = stats(law(t, u))
        if best is None or f/kl > best[0]:
            best = (f/kl, u, t, f, m, kl)
        u += 0.005
    ratio, u, t, f, m, kl = best
    t0 = tilt_for_mean(0.0)
    f0, m0, kl0 = stats(law(t0, 0.0))
    return dict(c_star=math.log(2)/2*ratio, optimal_odd_fraction=f, optimal_KL=kl,
                optimal_letter_law_1_to_6=law(t, u)[:6], odd_steps_per_bit=math.log(2)/kl,
                pure_tilt=dict(odd_fraction=f0, KL=kl0, rate_per_odd_step=math.exp(-kl0),
                               glide_odd_steps_per_bit=math.log(2)/kl0, coefficient=math.log(2)/2*f0/kl0),
                three_to_h_star_minus_one=3**(0.950-1))


def source_187():
    spec = importlib.util.spec_from_file_location('pe', ROOT/'04-computation/experiments/paid_exit_boundaries_20261012.py')
    pe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pe)
    steps, charge, h, word = first_band_exit(187)
    check((tuple(word), h, charge) == ((1, 2, 1, 3, 1, 1, 3, 4), 19, 6), 'actual 187 first-band word')
    R187 = pe.ExitRule((1, 2, 1, 3, 1, 1, 3), 4, 3, 1)
    Q187 = pe.ExitRule((1, 2, 1, 3, 1, 1, 3), 4, 1, 2)
    rules = {}
    for label, rule in (('linear drop 3', R187), ('quadratic drop 1', Q187)):
        d = pe.rule_data(rule)
        check(d['residue'] == 187 and d['modulus'] == 2**16 and d['first'] == 187, 'rule owns source 187')
        check(pe.apply_rule(rule, 187) == (19, (1, 2, 1, 3, 1, 1, 3, 4), 6), 'rule replays 187 to 19 with charge 6')
        for t in list(range(65))+[2**90+5]:
            n = 187+(2**16)*t
            hh, ww, delta = pe.apply_rule(rule, n)
            check(pe.C.replay(n, ww) == hh and (hh+1)*(1 << rule.drop) <= n+1, 'all-height control member')
        rules[label] = {k: d[k] for k in ('residue', 'modulus', 'first', 'cut', 'charge_cap', 'drop', 'power')}
    for bad in (pe.ExitRule((1, 2, 1, 3, 1, 1, 3), 4, 1, 1), pe.ExitRule((1, 2, 1, 3, 1, 1, 3), 4, 2, 1)):
        try:
            pe.rule_data(bad)
        except ValueError:
            check(True, 'linear invoice rejects drops one and two for 187')
        else:
            check(False, 'unpaid linear 187 rule accepted')
    # grammar reach with the new rule, and an auto-compiler that gives every pending source its own rule
    base_lin = tuple(pe.RUNS+[pe.R127, pe.R135])
    base_quad = tuple(pe.RUNS+[pe.R127, pe.Q127, pe.R135, pe.Q159])
    reach = {}
    for label, extra, power in (('folded+R187', base_lin+(R187,), 1), ('quadratic+Q187', base_quad+(Q187,), 2)):
        ledger = pe.grammar(8191, extra, power)
        pending = [n for n in range(1, 8192, 2) if n not in ledger]
        reach[label] = dict(accepted=len(ledger), first_pending=pending[:12])
    check(reach['folded+R187']['accepted'] == 1668 and reach['quadratic+Q187']['accepted'] == 1702, 'matched reach')

    def compile_for(n, power):
        x, w = n, []
        while len(w) < 60:
            a = pe.C.valuation(3*x+1)
            x = (3*x+1) >> a
            w.append(a)
            if 2*(x+1) <= n+1:
                for drop in range(1, E(n)-E(x)+1):
                    for fm in range(a, 0, -1):
                        try:
                            rule = pe.ExitRule(tuple(w[:-1]), fm, drop, power)
                            d = pe.rule_data(rule)
                        except ValueError:
                            continue
                        if d['first'] <= n and n % d['modulus'] == d['residue']:
                            return rule, d
        return None, None

    auto = {}
    for label, rules_, power in (('auto-linear', list(base_lin), 1), ('auto-quadratic', list(base_quad), 2)):
        added, stop = [], None
        while len(added) < 12:
            ledger = pe.grammar(8191, tuple(rules_), power)
            pending = [n for n in range(1, 8192, 2) if n not in ledger]
            if not pending:
                break
            n = pending[0]
            rule, d = compile_for(n, power)
            if rule is None:
                stop = n
                break
            rules_.append(rule)
            added.append(dict(source=n, head=list(rule.head), final_min=rule.final_min, drop=rule.drop, charge_cap=d['charge_cap']))
        auto[label] = dict(added=added, accepted=len(ledger), pending=pending[:8], uncompilable_source=stop)
    check(auto['auto-linear']['uncompilable_source'] == 235, 'source 235 has no linear own rule')
    check(auto['auto-quadratic']['uncompilable_source'] is None, 'quadratic auto-compiler ran twelve rules')
    return dict(rules=rules, reach=reach, auto=auto)


def main(scan_limit):
    out = {}
    rows, fails = record_witnesses()
    out['glide_record_witnesses'] = rows
    out['record_failures_by_coefficient'] = fails
    check(fails[3][0] == 56924955 and 1008932249296231 in fails[4] and not fails[5], 'record witness summary')
    steps, charge, h, word = first_band_exit(56924955)
    check((steps, charge, h, E(56924955), E(h)) == (198, 151, h, 25, 24) and 2*(h+1) <= 56924956, 'least witness replay')
    out['least_witness_word'] = dict(n=56924955, odd_steps=steps, charge=charge, h=h, word=word)
    out['ratio_record_below_2_27'] = dict(n=96883183, **{k: v for k, v in zip(('odd_steps', 'charge', 'h'), first_band_exit(96883183)[:3])})
    check(out['ratio_record_below_2_27']['charge'] == 176 and E(96883183) == 26 and E(out['ratio_record_below_2_27']['h']) == 25, 'ratio record 176/51')
    out['threshold'] = threshold_coefficient()
    out['source_187'] = source_187()
    if scan_limit:
        out['exhaustive_scan'] = least_counterexample_scan(scan_limit)
        check(out['exhaustive_scan']['least_failure'] == 56924955, 'least counterexample below the scan limit')
    out['checks'] = CHECKS
    return out


if __name__ == '__main__':
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 2**27+1
    print(json.dumps(main(limit), indent=1, sort_keys=True, default=str))
