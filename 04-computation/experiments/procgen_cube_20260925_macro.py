#!/usr/bin/env python3
"""
procgen_cube_20260925_macro.py -- section I of the strategy-cube run: microcosm -> macrocosm.

  I1  uniform random strategies at levels k = 6..12 (each with its nu-partner), same
      certificates and the same cycle-search budget as levels <= 5;
  I2  class fractions per level (levels 1-5 exact, 6-12 sampled) and their trend;
  I3  drift concentration: the stationary odd density of the main bottom SCC per level;
  I4  the k = infinity model: an independent random sign for every odd integer
      (smallest extra cycle, P(no extra cycle with minimum <= X)).

Called from procgen_cube_20260925_run.py (section_I).  Randomness is seeded (reproducible).
"""
import os
import math
import random
from collections import Counter

LOG2, LOG3 = math.log(2), math.log(3)
CRIT = LOG2 / LOG3
NB = 131072
SAMPLES = {6: 4000, 7: 3000, 8: 2000, 9: 1500, 10: 1000, 11: 600, 12: 400}
SEED = 20260925


def negate_mask_int(k, mask):
    H = 1 << (k - 1)
    rev = int(format(mask, f'0{H}b')[::-1], 2)
    return ((1 << H) - 1) ^ rev


def sample_level(run_engine, k, S, rng):
    H = 1 << (k - 1)
    masks = [rng.getrandbits(H) for _ in range(S)]
    partners = [negate_mask_int(k, m) for m in masks]
    text = '\n'.join(format(m, 'x') for m in masks + partners) + '\n'
    env = {'CUBE_DIM': '1' if k <= 8 else '0'}
    recs, _ = run_engine(['list', k, NB], stdin_text=text, env_extra=env)
    assert len(recs) == 2 * S
    return masks, recs[:S], recs[S:]


def label_of(e):
    lab = []
    if e['ci']:
        lab.append('I')
    if e['cii']:
        lab.append('II')
    if any(c[0] > 1 for c in e['cyc']):
        lab.append('III')
    return '+'.join(lab) if lab else 'IV'


def main_bottom(e):
    return max(e['bots'], key=lambda b: b['mass'])


def dimension_check(check, hdr, store):
    """The engine's D (Perron formula inf_beta log2 rho(A(beta))) against the growth of the exact counts
    |Bad_L| (Python big integers) up to L = 600, for three strategies."""
    import procgen_cube_20260925_core as core
    hdr("F3. Exceptional dimension: Perron formula against exact counts |Bad_L| to L = 600")
    for (k, mask) in ((2, 0), (3, 2), (3, 8)):
        c = core.bad_counts_dp(k, mask, 600)
        D = store[(k, mask)]['eng']['D']
        sl = [math.log2(c[b] / c[a]) / (b - a) for (a, b) in ((200, 300), (300, 400), (400, 600))]
        check(sl[0] < sl[1] + 1e-3 and abs(sl[2] - D) < 0.01, "count slopes approach D")
        print(f"  {core.sigma_string(k, mask):5s} (level {k}): D = {D:.6f};  slopes of log2|Bad_L| on [200,300], [300,400], [400,600]:"
              f" {sl[0]:.4f}, {sl[1]:.4f}, {sl[2]:.4f}")
    print("  (the slopes carry the usual O(log L / L) correction; Collatz's D is h(log_3 2) = 0.949956)")


def single_flips(run_engine, check, hdr):
    """The 2^(k-1) single-flip neighbours of Collatz at levels 2..9: classes and extra cycles."""
    import procgen_cube_20260925_core as core
    hdr("H4. The single-flip neighbours of Collatz, levels 2-9 (same budget: starts <= 131072, cap 2^60)")
    for k in range(2, 10):
        H = 1 << (k - 1)
        text = '\n'.join(format(1 << i, 'x') for i in range(H)) + '\n'
        recs, _ = run_engine(['list', k, NB], stdin_text=text, env_extra={'CUBE_DIM': '0'})
        check(len(recs) == H, "single flips count")
        n1 = sum(1 for e in recs if e['ci'])
        n2 = sum(1 for e in recs if e['cii'])
        w = []
        for i, e in enumerate(recs):
            extra = sorted(c[0] for c in e['cyc'] if c[0] > 1)
            for (mn, ln, od, bs) in e['cyc']:
                if mn > 1:
                    sig = core.sig_table(k, 1 << i)
                    x, m, odd = mn, mn, 0
                    for _ in range(ln):
                        odd += x & 1
                        x = core.T(x, k, sig)
                        m = min(m, x)
                    check(x == mn and m == mn and odd == od, "single-flip cycle verified")
            if extra:
                w.append((2 * i + 1, extra[0]))
        print(f"  level {k}: {H} single flips; in (i): {n1}, in (ii): {n2}, with an extra cycle (iii): {len(w)}"
              f"  [flipped residue: smallest extra cycle] {w}")


def section_I(run_engine, check, hdr, store):
    dimension_check(check, hdr, store)
    single_flips(run_engine, check, hdr)
    fast = os.environ.get('CUBE_FAST') == '1'
    hdr("I. Microcosm -> macrocosm: class fractions as k grows")
    rng = random.Random(SEED)
    rows = {}
    # exact levels 1..5 from the store
    for k in (1, 2, 3, 4, 5):
        keys = [key for key in store if key[0] == k]
        labs = Counter(store[key]['label'] for key in keys)
        n = len(keys)
        opn = [store[key] for key in keys if store[key]['label'] == 'IV']
        obstructed = sum(1 for s in opn if store[(k, s['nu'])]['label'] != 'IV')
        drifts = [main_bottom(store[key]['eng'])['podd'] for key in keys]
        Dneg = [store[key]['eng']['D'] for key in keys if not store[key]['I'] and not store[key]['pos_bottom']]
        rows[k] = {'n': n, 'labs': labs, 'obstructed': obstructed, 'open': len(opn), 'exact': True,
                   'podd': drifts, 'Dneg': Dneg,
                   'pos': sum(1 for key in keys if store[key]['pos_bottom']),
                   'rmax1': sum(1 for key in keys if abs(store[key]['eng']['mumax'] - (LOG3 - LOG2)) < 1e-9),
                   'fix': sum(1 for key in keys if store[key]['sig'][0] == '-' or store[key]['sig'][-1] == '+')}
    levels = [6, 7, 8, 9, 10, 11, 12]
    if fast:
        levels = [6, 7, 8]
    for k in levels:
        S = SAMPLES[k] if not fast else 200
        masks, recs, precs = sample_level(run_engine, k, S, rng)
        labs = Counter()
        opn = obstructed = pos = rmax1 = fix = 0
        podd = []
        Dneg = []
        for m, e, pe in zip(masks, recs, precs):
            if k <= 8 and not e['ci'] and e['npos'] == 0:
                Dneg.append(e['D'])
            lab = label_of(e)
            plab = label_of(pe)
            # nu-invariance of the G-level data on the sample (sheet-blindness)
            check(e['ci'] == pe['ci'] and e['cii'] == pe['cii'] and e['bad'] == pe['bad'], "nu invariance (sample)")
            check(e['unres'] == 0, "no unresolved start (sample)")
            labs[lab] += 1
            if lab == 'IV':
                opn += 1
                if plab != 'IV':
                    obstructed += 1
            pos += e['npos'] > 0
            podd.append(main_bottom(e)['podd'])
            if k <= 10:
                rmax1 += abs(e['mumax'] - (LOG3 - LOG2)) < 1e-9
            else:
                rmax1 = None
            s1 = (m & 1)                       # sigma(1) = -1
            sm1 = not ((m >> ((1 << (k - 1)) - 1)) & 1)   # sigma(-1) = +1
            fix += bool(s1 or sm1)
        rows[k] = {'n': S, 'labs': labs, 'obstructed': obstructed, 'open': opn, 'exact': False, 'podd': podd,
                   'pos': pos, 'rmax1': rmax1, 'fix': fix, 'Dneg': Dneg}
    print(f"  levels 1-5: all strategies (exact); levels 6-12: uniform samples (seed {SEED}), each with its nu-partner;")
    print(f"  the same certificates and the same cycle budget (starts <= {NB}, cap 2^60) at every level.")
    print()
    print("   k   strategies   (i)      (ii)     (iii)    (iv)=OPEN  | OPEN sheet-obstructed | pos-drift  rho_max=1  sigma(1)=- or sigma(-1)=+")
    for k in sorted(rows):
        r = rows[k]
        n = r['n']
        labs = r['labs']
        fI = sum(v for lab, v in labs.items() if lab.split('+')[0] == 'I') / n
        fII = sum(v for lab, v in labs.items() if 'II' in lab.split('+')) / n
        fIII = sum(v for lab, v in labs.items() if 'III' in lab.split('+')) / n
        fIV = labs.get('IV', 0) / n
        ob = r['obstructed'] / max(1, r['open'])
        rm = f"{r['rmax1'] / n:.4f}" if r['rmax1'] is not None else "  n/a "
        tag = 'all' if r['exact'] else 'sample'
        print(f"  {k:2d}   {n:6d} {tag:6s} {fI:.4f}   {fII:.4f}   {fIII:.4f}   {fIV:.4f}     | {ob:.3f} of OPEN         | {r['pos'] / n:.4f}    {rm}     {r['fix'] / n:.4f}")
    # standard errors for the sampled OPEN fraction
    print()
    for k in sorted(rows):
        r = rows[k]
        if r['exact']:
            continue
        p = r['labs'].get('IV', 0) / r['n']
        print(f"   k = {k}: OPEN fraction {p:.4f} +- {math.sqrt(p * (1 - p) / r['n']):.4f} (1 s.e.)")
    hdr("I3. Drift concentration: stationary odd density pi_odd of the main bottom SCC")
    print(f"   critical density log_3 2 = {CRIT:.6f};  Collatz pi_odd = 1/2 (drift log(3/4) per odd step)")
    print("   k    mean pi_odd   s.d.     s.d.*2^(k/2)   frac(pi_odd > log_3 2)   mean D over negative-drift non-(i)")
    for k in sorted(rows):
        r = rows[k]
        xs = r['podd']
        m = sum(xs) / len(xs)
        sd = math.sqrt(sum((x - m) ** 2 for x in xs) / len(xs))
        above = sum(1 for x in xs if x > CRIT) / len(xs)
        Dn = r.get('Dneg') or []
        Ds = f"{sum(Dn) / len(Dn):.4f} (n={len(Dn)})" if Dn else "n/a"
        print(f"  {k:2d}    {m:.5f}      {sd:.5f}   {sd * 2 ** (k / 2):.3f}          {above:.4f}                   {Ds}")
    print(f"   (Collatz: pi_odd = 1/2, D = h(log_3 2) = 0.94996)")
    # ---------------------------------------------------------------- k = infinity
    hdr("I4. The k = infinity model: independent random sign for every odd integer")
    Ns = 10 ** 6 if not fast else 10 ** 5
    R = 1000 if not fast else 100
    recs = []
    import subprocess
    engine = run_engine.__globals__['ENGINE']
    out = subprocess.run([engine, 'rsign', str(Ns), str(R), str(SEED)], capture_output=True, text=True, check=True).stdout
    smallest = []
    ncycles = []
    escs = []
    for line in out.splitlines():
        if not line.startswith('R '):
            continue
        d = dict(p.split('=', 1) for p in line.split()[1:])
        cyc = []
        if d['cyc']:
            for c in d['cyc'].split(';'):
                mn, ln, od, bs = c.split(':')
                cyc.append((int(mn), int(ln), int(od), int(bs)))
        check(int(d['unres']) == 0, "rsign unresolved")
        extra = [c[0] for c in cyc if c[0] > 1]
        smallest.append(min(extra) if extra else None)
        ncycles.append(len(extra))
        escs.append(int(d['esc']))
    check(len(smallest) == R, "rsign count")
    print(f"   {R} independent sign fields, starts 1..{Ns}, cap 2^60")
    print("   X        P(no extra cycle with minimum <= X)")
    for X in (10, 30, 100, 300, 1000, 3000, 10 ** 4, 3 * 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6):
        if X > Ns:
            break
        p = sum(1 for s in smallest if s is None or s > X) / R
        print(f"   {X:<8d} {p:.4f}")
    cnt = Counter(ncycles)
    print(f"   number of extra cycles met (starts <= {Ns}): {dict(sorted(cnt.items()))}; mean {sum(ncycles) / R:.3f}")
    print(f"   escapes past 2^60: max over fields {max(escs)}")
    return rows
