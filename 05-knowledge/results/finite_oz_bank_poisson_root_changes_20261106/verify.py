"""verify.py - kind-pasteur-2026-11-06-S13, THM-4888 / HYP-9480: decide the pre-registered predictions P1-P6.

    python3 verify.py --out VALIDATION_NORMAL.json
    python3 -O verify.py --out VALIDATION_OPTIMIZED.json      (receipts must agree byte for byte)

Inputs (produced by the C programs of this packet, see README.md): p1_density.txt, p2_root.txt, p3_bankA_w26.txt,
p3_bankA_w27.txt (density.c), root_changes_1e8.txt (root_changes.c), p6_root_w25.txt, p6_root_w26.txt,
p6_seeds_w24.txt (residue_law.c), seeds_summary.json (seeds.py). P1-P5: PREREGISTRATION.md; P6: PREREGISTRATION_P6.md.
Identities that the proofs assert raise RuntimeError when they fail; prediction outcomes are recorded, not enforced.
"""
import sys, json, math, hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def q_array(n):
    """Q[0..n-1] by the recurrence Q(2m) = Q(m), Q(2m+1) = Q(m) + m + (m mod 2), level by level."""
    Q = np.zeros(n, dtype=np.int64)
    lo = 1
    while lo < n:
        hi = min(2 * lo, n)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx // 2
        Q[idx] = Q[half] + (idx & 1) * (half + (half & 1))
        lo = hi
    return Q


def T(m):
    """Session 206's converse involution on parameters m >= 1 (T(0) = 0)."""
    if m == 0:
        return 0
    n = m.bit_length()
    a = m - (1 << (n - 1))
    r = int(format(a, '0%db' % (n - 1))[::-1], 2) if n > 1 else 0
    return (1 << n) - 1 - r


def predict(h, R, e, d):
    """Local-factor Poisson prediction of the hit density among odd integers (PREREGISTRATION.md)."""
    lam_even = 2.0 / d * math.log(R / (R - d / 2.0))
    lam_odd = 2.0 / d * math.log((R - d / 2.0) / (R - d))
    P = 1
    for p in (3, 5, 7, 11, 13):
        if d % p == 0 and R % p == 0:
            P *= p
    if P == 1:
        return 1.0 - math.exp(-(lam_even + lam_odd)), P
    f = 1.0 / P
    c_even, c_odd = h % P, (h - e) % P
    if c_even == c_odd:
        return f * (1.0 - math.exp(-(lam_even + lam_odd) / f)), P
    return f * (1.0 - math.exp(-lam_even / f)) + f * (1.0 - math.exp(-lam_odd / f)), P


def main():
    out = HERE / 'VALIDATION.json'
    if '--out' in sys.argv:
        out = HERE / sys.argv[sys.argv.index('--out') + 1]
    rec = {}

    # ---------------- P5: identity A on all m < 2^26
    N = 1 << 26
    Q = q_array(N)
    m = np.arange(N, dtype=np.int64)
    runs = np.zeros(N, dtype=np.int64)
    starts = m & ~(m << 1)                  # bit j set iff a run of ones starts at j
    for j in range(27):
        runs += (starts >> j) & 1
    lhs = (Q - (m % 3) * (m % 3) + runs) % 3
    need(bool(np.all(lhs == 0)), 'mod-3 law Q(m) = m^2 - r(m) (mod 3)')
    need(bool(np.all(Q % 2 == 0)) and bool(np.all(Q[1:] <= m[1:] - 1)), 'Q even and Q(m) <= m - 1')
    rec['P5_identity_A'] = 'Q(m) = m^2 - r(m) (mod 3) for all m < 2^26 (and Q even, Q(m) <= m - 1)'

    # ---------------- P5: identity B on all physical seeds with U <= 64, m < 2^16
    sys.path.insert(0, str(HERE))
    import seeds as seedmod
    st = seedmod.states(64)
    M2 = 1 << 16
    Tm = [T(x) for x in range(M2)]
    for x in range(1, M2):
        need(Tm[Tm[x]] == x and Tm[x].bit_length() == x.bit_length(), 'T is an involution of each dyadic block')
        need(int(Q[Tm[x]]) - int(Q[x]) == 2 * (Tm[x] - x), 'Session 206 identity Q(T(m)) - Q(m) = 2(T(m) - m)')
    nseed = 0
    for (D, E, U, d) in st:
        h, R, e = D, U, E
        nseed += 1
        for x in range(M2):
            t = Tm[x]
            lhs_b = (h + R * t - d * int(Q[t]) - e * (t & 1)) - (h + R * x - d * int(Q[x]) - e * (x & 1))
            rhs_b = (R - 2 * d) * (t - x) - e * ((t & 1) - (x & 1))
            if lhs_b != rhs_b:
                raise RuntimeError('identity B fails')
    rec['P5_identity_B'] = {'seeds_U_le_64': nseed, 'parameters_per_seed': M2, 'status': 'holds'}
    summ = json.loads((HERE / 'seeds_summary.json').read_text())
    rec['P5_classification'] = {'states_U_le_4096': summ['states_U_le_4096'],
                                'states_with_U_eq_2d_and_E_eq_0': summ['states_with_U_eq_2d_and_E_eq_0'],
                                'only_rooted_edge': summ['states_with_U_eq_2d_and_E_eq_0'] == [[1, 0, 2, 1]]}
    need(rec['P5_classification']['only_rooted_edge'], 'classification of symmetric seeds')

    # ---------------- P1
    rows = []
    for line in (HERE / 'p1_density.txt').read_text().split('\n'):
        if not line.strip():
            continue
        h, R, e, d, hit, nodd = map(int, line.split())
        obs = hit / nodd
        pred, P = predict(h, R, e, d)
        rows.append([h, R, e, d, P, round(obs, 6), round(pred, 6), round(obs - pred, 6)])
    devs = [abs(r[7]) for r in rows]
    n_le4 = sum(1 for x in devs if x <= 0.004)
    rec['P1'] = {'seeds': len(rows), 'within_0.004': n_le4, 'max_abs_dev': max(devs),
                 'mean_abs_dev': round(sum(devs) / len(devs), 6),
                 'worst': sorted(rows, key=lambda r: -abs(r[7]))[:8],
                 'with_local_factor': sum(1 for r in rows if r[4] > 1),
                 'prediction': n_le4 >= 0.95 * len(rows) and max(devs) <= 0.010}
    rec['P1_table'] = rows

    # ---------------- P2
    h, R, e, d, hit, nodd = map(int, (HERE / 'p2_root.txt').read_text().split())
    dens = hit / nodd
    rec['P2'] = {'window': [2 ** 27, 2 ** 28], 'hit_density': round(dens, 6), 'prediction': 0.4978 <= dens <= 0.5005}

    # ---------------- P3
    p3 = {}
    ok = True
    for name, w in (('p3_bankA_w26.txt', [2 ** 26, 2 ** 27]), ('p3_bankA_w27.txt', [2 ** 27, 2 ** 28])):
        c, n, c0, n0, c1, n1, c2, n2 = map(int, (HERE / name).read_text().split())
        frac = c / n
        p3[str(w)] = {'miss_fraction': round(frac, 6), 'miss_by_class_mod3': [round(c0 / n0, 6), round(c1 / n1, 6),
                                                                                round(c2 / n2, 6)]}
        ok = ok and 0.0205 <= frac <= 0.0235
    rec['P3'] = {'windows': p3, 'prediction': ok}

    # ---------------- P4
    text = (HERE / 'root_changes_1e8.txt').read_text().split('\n')
    misses_line = [l for l in text if l.startswith('misses (<=2) above 10^5:')][0]
    big = [int(x) for x in misses_line.split(':')[1].split()]
    a_count = sum(1 for v in big if 4 * 10 ** 6 < v <= 10 ** 8)
    not1 = [v for v in big if v % 3 != 1]
    win = {}
    for l in text:
        l = l.strip()
        if l.startswith('['):
            rng, rest = l.split(':', 1)
            lo, hi = map(int, rng.strip('[)').split(','))
            vals = dict(tok.split(':') for tok in rest.split())
            win[(lo, hi)] = {k: int(v) for k, v in vals.items()}
    full = sorted(k for k in win if k[1] == 2 * k[0] and k[0] >= 2 ** 16)
    ratios = []
    for a, b in zip(full, full[1:]):
        ratios.append(round(win[b]['k1'] / win[a]['k1'], 4))
    totals = [l for l in text if l.startswith('<=')]
    rec['P4'] = {'summary_lines': totals, 'misses_above_1e5_upto_1e8': big,
                 'count_in_(4e6,1e8]': a_count, 'not_1_mod_3': not1, 'k1_window_ratios_from_2^16': ratios,
                 'a': a_count <= 6, 'b': len(not1) <= 2, 'c': all(1.45 <= r <= 1.75 for r in ratios),
                 'k2_misses_per_window': {f'[{a},{b})': win[(a, b)]['k2'] for (a, b) in sorted(win)}}
    rec['P4']['prediction'] = rec['P4']['a'] and rec['P4']['b'] and rec['P4']['c']

    # ---------------- P6 (PREREGISTRATION_P6.md): residue-resolved Poisson law, miss_c = exp(-kappa lambda_c)
    def load63(name):
        rows = {}
        for line in (HERE / name).read_text().split('\n'):
            if not line.strip():
                continue
            h, R, e, d, M, c, n, p, v = map(int, line.split())
            need(M == 63, 'P6 files are mod 63')
            rows.setdefault((h, R, e, d), {})[c] = (n, p, v)
        return rows

    p6a = {}
    ok_a = True
    for name, w in (('p6_root_w25.txt', [2 ** 25, 2 ** 26]), ('p6_root_w26.txt', [2 ** 26, 2 ** 27])):
        cls = load63(name)[(1, 2, 0, 1)]
        info = {}
        for M in (7, 9, 21, 63):
            agg = {}
            for c, (n, p, v) in cls.items():
                a = agg.setdefault(c % M, [0, 0, 0])
                a[0] += n
                a[1] += p
                a[2] += v
            miss = {c: 1 - a[2] / a[0] for c, a in agg.items()}
            pred = {c: math.exp(-a[1] / a[0] / 2) for c, a in agg.items()}
            info[str(M)] = {'max_abs_dev': round(max(abs(miss[c] - pred[c]) for c in agg), 6),
                            'mean_dev': round(sum(miss[c] - pred[c] for c in agg) / len(agg), 6),
                            'span_of_miss_rates': round(max(miss.values()) - min(miss.values()), 6)}
            if M == 7:
                info[str(M)]['miss'] = [round(miss[c], 5) for c in range(7)]
                info[str(M)]['pred'] = [round(pred[c], 5) for c in range(7)]
            ok_a = ok_a and info[str(M)]['max_abs_dev'] <= 0.006
        ok_a = ok_a and info['7']['span_of_miss_rates'] >= 0.05
        p6a[str(w)] = info
    rec['P6a'] = {'windows': p6a, 'prediction': ok_a}

    seeds6 = load63('p6_seeds_w24.txt')
    devs6, blind, worst6 = [], [], []
    zero_ok, nzero = True, 0
    for s, cls in sorted(seeds6.items()):
        live = [(n, p, v) for (n, p, v) in cls.values() if p > 0]
        lam_live = sum(p for n, p, v in live) / sum(n for n, p, v in live)
        for c, (n, p, v) in sorted(cls.items()):
            if p == 0:
                nzero += 1
                zero_ok = zero_ok and v == 0
                continue
            obs = 1 - v / n
            dv = obs - math.exp(-p / n)
            devs6.append(abs(dv))
            blind.append(abs(obs - math.exp(-lam_live)))
            worst6.append([list(s), c, round(p / n, 5), round(obs, 6), round(dv, 6)])
    srt = sorted(devs6)
    rec['P6b'] = {'seeds': len(seeds6), 'live_pairs': len(devs6), 'zero_intensity_classes': nzero,
                  'zero_intensity_classes_empty': zero_ok,
                  'within_0.004': sum(1 for x in devs6 if x <= 0.004), 'max_abs_dev': round(srt[-1], 6),
                  'p95_abs_dev': round(srt[int(0.95 * len(srt))], 6),
                  'mean_abs_dev_resolved': round(sum(devs6) / len(devs6), 6),
                  'mean_abs_dev_class_blind': round(sum(blind) / len(blind), 6),
                  'worst': sorted(worst6, key=lambda r: -abs(r[4]))[:6]}
    rec['P6b']['prediction'] = (rec['P6b']['within_0.004'] >= 0.95 * len(devs6) and srt[-1] <= 0.010 and zero_ok)

    # hashes of the inputs with line endings normalised to LF (the committed bytes on every platform)
    rec['input_sha256'] = {p.name: hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                           for p in sorted(HERE.glob('*.txt'))}
    text = json.dumps(rec, sort_keys=True, indent=1) + '\n'
    with open(out, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
