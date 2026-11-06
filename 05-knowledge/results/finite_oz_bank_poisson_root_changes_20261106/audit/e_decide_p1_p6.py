"""e_decide_p1_p6.py - auditor's decisions of P1, P2, P6a, P6b from the auditor's own counts (d_image.c outputs),
with the registered criteria re-implemented from PREREGISTRATION.md / PREREGISTRATION_P6.md (not from verify.py).

Also reports two alternative models that the registered tolerances should be able to reject:
  * INJ (no collisions): a seed's parameters give distinct values, miss_c = 1 - kappa*lambda_c (hit = parameters);
  * BLIND: class-blind Poisson, miss_c = exp(-kappa * lambda_live) with lambda_live pooled over live classes.
If INJ passes the same registered criteria, the criteria cannot distinguish the Poisson collision law from
"no collisions at all" and the confirmation is evidence only for the parameter intensities (capacity theorem).
"""
import math, json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

def gcd(a, b):
    while b: a, b = b, a % b
    return a

def registered_local_factor(h, R, e, d):
    """PREREGISTRATION.md 'Local-factor Poisson prediction', implemented from the text.
    Pieces: m even and m odd, intensities = parity-split capacities per odd integer (2/d) ln(R/(R-d/2)) and
    (2/d) ln((R-d/2)/(R-d)).  P = product of odd primes p <= 13 dividing gcd(d,R).  A piece landing in the classes
    of relative size f covers f (1 - exp(-lambda/f)); pieces landing in the same class mod P are merged."""
    le = 2.0 / d * math.log(R / (R - d / 2.0))
    lo = 2.0 / d * math.log((R - d / 2.0) / (R - d))
    P = 1
    for p in (3, 5, 7, 11, 13):
        if gcd(d, R) % p == 0: P *= p
    if P == 1:
        return 1 - math.exp(-(le + lo)), (le + lo), P
    f = 1.0 / P
    ce, co = h % P, (h - e) % P
    if ce == co:
        return f * (1 - math.exp(-(le + lo) / f)), min(f, le + lo), P
    return f * (1 - math.exp(-le / f)) + f * (1 - math.exp(-lo / f)), min(f, le) + min(f, lo), P

def load(fn):
    d = {}
    for l in open(HERE / fn):
        t = l.split()
        if t and t[0] == 'S':
            h, R, e, dd, M, c, n, p, v = map(int, t[1:])
            d.setdefault((h, R, e, dd), {})[c] = (n, p, v)
    return d

out = {}
# ---------------- P1
import verify as author   # only to compare the author's predict() with the registered text
rows = []
for s, cls in sorted(load('out_p1_w24_M1.txt').items()):
    n, p, v = cls[0]
    obs = v / n
    pred, inj, P = registered_local_factor(*s)
    apred, aP = author.predict(*s)
    rows.append((s, P, obs, pred, inj, apred, aP))
dev = [abs(r[2] - r[3]) for r in rows]
devinj = [abs(r[2] - r[4]) for r in rows]
out['P1'] = {
    'seeds': len(rows),
    'within_0.004': sum(x <= 0.004 for x in dev), 'within_2e-4': sum(x <= 2e-4 for x in dev), 'max_abs_dev': max(dev),
    'decision': sum(x <= 0.004 for x in dev) >= 0.95 * len(rows) and max(dev) <= 0.010,
    'author_predict_equals_registered_text': all(abs(r[3] - r[5]) < 1e-12 and r[1] == r[6] for r in rows),
    'seeds_with_local_factor': sum(r[1] > 1 for r in rows),
    'worst': [[list(r[0]), r[1], round(r[2], 6), round(r[3], 6)] for r in sorted(rows, key=lambda r: -abs(r[2] - r[3]))[:4]],
    'hit_density_range': [round(min(r[2] for r in rows), 5), round(max(r[2] for r in rows), 5)],
    'INJ_model_within_0.004': sum(x <= 0.004 for x in devinj), 'INJ_model_max_abs_dev': max(devinj),
    'INJ_model_would_pass_P1': sum(x <= 0.004 for x in devinj) >= 0.95 * len(rows) and max(devinj) <= 0.010,
}
# ---------------- P2
r = load('out_root_w27_M63.txt')[(1, 2, 0, 1)]
n = sum(x[0] for x in r.values()); v = sum(x[2] for x in r.values())
out['P2'] = {'hit_density': v / n, 'decision': 0.4978 <= v / n <= 0.5005}

# ---------------- P6a (root, kappa = 1/2)
def aggregate(cls, M):
    agg = {}
    for c, (n, p, v) in cls.items():
        a = agg.setdefault(c % M, [0, 0, 0]); a[0] += n; a[1] += p; a[2] += v
    return agg
p6a = {}; ok = True
spans = {}
for w in (20, 22, 24, 25, 26, 27):
    cls = load(f'out_root_w{w}_M63.txt')[(1, 2, 0, 1)]
    info = {}
    for M in (7, 9, 21, 63):
        agg = aggregate(cls, M)
        miss = {c: 1 - a[2] / a[0] for c, a in agg.items()}
        lam = {c: a[1] / a[0] for c, a in agg.items()}
        dpois = max(abs(miss[c] - math.exp(-lam[c] / 2)) for c in agg)
        dinj = max(abs(miss[c] - (1 - lam[c] / 2)) for c in agg)
        lam_live = sum(a[1] for a in agg.values()) / sum(a[0] for a in agg.values())
        dblind = max(abs(miss[c] - math.exp(-lam_live / 2)) for c in agg)
        info[M] = {'max_dev_poisson': round(dpois, 6), 'max_dev_INJ': round(dinj, 4), 'max_dev_class_blind': round(dblind, 4),
                   'span_miss': round(max(miss.values()) - min(miss.values()), 4),
                   'span_lambda': round(max(lam.values()) - min(lam.values()), 4),
                   'mean_signed_dev_poisson': round(sum(miss[c] - math.exp(-lam[c] / 2) for c in agg) / len(agg), 6)}
        if M == 7:
            info[M]['miss'] = [round(miss[c], 4) for c in range(7)]
            info[M]['lambda'] = [round(lam[c], 4) for c in range(7)]
    if w in (25, 26):
        ok = ok and all(info[M]['max_dev_poisson'] <= 0.006 for M in info) and info[7]['span_miss'] >= 0.05
    p6a[f'[2^{w},2^{w+1})'] = info
out['P6a'] = {'windows': p6a, 'decision_registered_windows_2^25_2^26': ok,
              'limit_lambda_by_Session206_root_residue_capacity': round(2 * math.log(2), 4)}

# ---------------- P6b (102 seeds, kappa = 1)
seeds6 = [tuple(map(int, l.split())) for l in open(HERE.parent / 'p6_seeds.txt') if l.strip()]
allc = load('out_p1_w24_M63.txt')
devs, devinj, devblind, zero_ok, nzero = [], [], [], True, 0
worst = []
for s in seeds6:
    cls = allc[s]
    live = [(n, p, v) for (n, p, v) in cls.values() if p > 0]
    lam_live = sum(p for n, p, v in live) / sum(n for n, p, v in live)
    for c, (n, p, v) in cls.items():
        if p == 0:
            nzero += 1; zero_ok = zero_ok and v == 0; continue
        obs = 1 - v / n; lam = p / n
        devs.append(abs(obs - math.exp(-lam))); devinj.append(abs(obs - (1 - lam))); devblind.append(abs(obs - math.exp(-lam_live)))
        worst.append((abs(obs - math.exp(-lam)), list(s), c, round(lam, 4), round(obs - math.exp(-lam), 6)))
srt = sorted(devs)
def summary(x):
    xs = sorted(x)
    return {'n': len(xs), 'within_0.004': sum(v <= 0.004 for v in xs), 'max': round(xs[-1], 6), 'p95': round(xs[int(0.95 * len(xs))], 6), 'mean': round(sum(xs) / len(xs), 7),
            'passes_registered_criteria': sum(v <= 0.004 for v in xs) >= 0.95 * len(xs) and xs[-1] <= 0.010}
out['P6b'] = {'seeds': len(seeds6), 'live_pairs': len(devs), 'zero_intensity_classes': nzero, 'zero_classes_empty': zero_ok,
              'poisson': summary(devs), 'INJ_no_collision_model': summary(devinj), 'class_blind_poisson': summary(devblind),
              'decision': summary(devs)['passes_registered_criteria'] and zero_ok,
              'worst5': [w[1:] for w in sorted(worst, reverse=True)[:5]]}
print(json.dumps(out, indent=1))
json.dump(out, open(HERE / 'out_decisions_p1_p2_p6.json', 'w'), indent=1)
