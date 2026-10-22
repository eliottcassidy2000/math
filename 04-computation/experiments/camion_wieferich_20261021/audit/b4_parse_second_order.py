# Audit b4: independent parser of the session's second_order_*.out files (THM-4750 (4) table, HYP-9361 (a),
# HYP-9358 (b) joint sporadic histogram).  Reads only the per-prime lines "p=.. eps=.. q2=.. [WIEFERICH] zero-classes: c:r[tags] ...".
# Recomputes: per class type (1, 2, 3, sporadic) the counts of 0/1/2 roots for 7 <= p <= 10^4; Wieferich lines;
# linear (L) and double (D) classes; odd lambda_C primes; the (class1, class2) contingency and chi-square;
# split rate by eps; the sporadic-count histogram, mean and variance; and the expected numbers of linear / double classes
# under the 1/p model.   Usage: python3 b4_parse_second_order.py <dir with second_order_*.out>
import sys, re, glob, os
from math import exp
d = sys.argv[1] if len(sys.argv) > 1 else '..'
files = sorted(glob.glob(os.path.join(d, 'second_order_*.out')))
lines = []
for f in files:
    for ln in open(f):
        if ln.startswith('p='): lines.append(ln.strip())
pat = re.compile(r'^p=(\d+) eps=([+-]1) q2=(\d+)( WIEFERICH)? zero-classes:(.*?)\s+lambda_C=(\d+)$')
recs = {}
for ln in lines:
    mt = pat.match(ln)
    assert mt, ln
    p = int(mt.group(1)); eps = int(mt.group(2)); wief = bool(mt.group(4)); body = mt.group(5).split(); lam = int(mt.group(6))
    cls = []
    for tok in body:
        mm = re.match(r'^(\d+):(\d+)([A-Z]*)$', tok); assert mm, tok
        cls.append((int(mm.group(1)), int(mm.group(2)), mm.group(3)))
    assert p not in recs, p
    recs[p] = dict(eps=eps, wief=wief, cls=cls, lam=lam)
ps = sorted(recs)
def isprime(n): return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))
allp = [q for q in range(3, 10001) if isprime(q)]
print(f"files: {[os.path.basename(f) for f in files]}")
print(f"primes present: {len(ps)} (min {ps[0]}, max {ps[-1]}); missing odd primes in [3,10^4]: {[q for q in allp if q not in recs]}")
P = [p for p in ps if 7 <= p <= 10000]
print(f"primes 7..10^4: {len(P)}")
tab = {t: [0, 0, 0, 0] for t in ('1', '2', '3', 'sporadic')}   # counts of 0,1,2 roots, >2
lin = []; dbl = []; zc = 0; wlines = []; oddlam = []; lamsum_mismatch = []
cont = {(0, 0): 0, (0, 2): 0, (2, 0): 0, (2, 2): 0}
split_by_eps = {1: [0, 0], -1: [0, 0]}
sporhist = {}
exp_lin = 0.0; exp_dbl = 0.0
for p in P:
    r = recs[p]; rootsum = 0; spor = 0; r1 = r2 = None
    for c, k, tag in r['cls']:
        zc += 1; rootsum += k
        typ = '1' if c == 1 else '2' if c == 2 else '3' if c == 3 else 'sporadic'
        if typ == 'sporadic': spor += 1
        tab[typ][min(k, 3)] += 1
        if 'L' in tag: lin.append((p, c, k, r['wief']))
        if 'D' in tag: dbl.append((p, c, k))
        if typ == '1': r1 = k
        if typ == '2': r2 = k
        if typ in ('1', '2', 'sporadic') and k in (0, 2) and not r['wief']:
            split_by_eps[r['eps']][0 if k == 0 else 1] += 1
        if not r['wief']:
            if typ != '3': exp_lin += 1.0 / p       # X(m-2) = 0 mod p (class 3 has X(1) = 4, never 0)
            exp_dbl += 1.0 / p                      # discriminant 0 (class 3: c1 = c2/2)
    if rootsum != r['lam']: lamsum_mismatch.append(p)
    if r['lam'] % 2: oddlam.append((p, r['lam'], [x for x in r['cls'] if x[2]]))
    if r['wief']: wlines.append((p, r['cls'], r['lam']))
    if r1 in (0, 2) and r2 in (0, 2): cont[(r1, r2)] += 1
    sporhist[spor] = sporhist.get(spor, 0) + 1
print(f"zero classes (7 <= p <= 10^4): {zc}")
print("table  class: [0 roots, 1 root, 2 roots, >2]")
for t in ('1', '2', '3', 'sporadic'): print(f"   {t:9s}: {tab[t]}  (total {sum(tab[t])})")
print("Wieferich lines:", wlines)
print(f"linear classes (tag L): {len(lin)} -> {lin}")
print(f"double-root classes (tag D): {len(dbl)} -> {dbl}")
print(f"expected under the 1/p model (non-Wieferich): linear {exp_lin:.2f}, double {exp_dbl:.2f}; "
      f"Poisson P(<=1 | mean {exp_dbl:.2f}) = {exp(-exp_dbl) * (1 + exp_dbl):.4f}")
print(f"lambda_C = sum of roots consistent on every line: {not lamsum_mismatch} {lamsum_mismatch[:5]}")
print(f"odd lambda_C primes: {oddlam}")
print(f"even lambda_C: {len(P) - len(oddlam)}, odd: {len(oddlam)}")
n = sum(cont.values()); r0 = cont[(0, 0)] + cont[(0, 2)]; c0 = cont[(0, 0)] + cont[(2, 0)]
chi = 0.0
for (a, b), o in cont.items():
    ra = r0 if a == 0 else n - r0; cb = c0 if b == 0 else n - c0; e = ra * cb / n; chi += (o - e) ** 2 / e
print(f"(class1, class2) contingency over {n} primes: {cont}; chi-square (1 df) = {chi:.3f}")
print(f"split rate (non-triangle classes, 0/2 roots, non-Wieferich) by eps: "
      f"eps=+1 inert/split {split_by_eps[1]} ({split_by_eps[1][1] / sum(split_by_eps[1]):.3f} split); "
      f"eps=-1 {split_by_eps[-1]} ({split_by_eps[-1][1] / sum(split_by_eps[-1]):.3f} split)")
H = [sporhist.get(i, 0) for i in range(max(sporhist) + 1)]
mean = sum(i * h for i, h in enumerate(H)) / len(P); var = sum(i * i * h for i, h in enumerate(H)) / len(P) - mean ** 2
print(f"sporadic classes per prime, histogram 0..: {H}; mean {mean:.4f}, variance {var:.4f}")
# Poisson(mean) chi-square over bins 0,1,2,3,>=4
from math import factorial
lam_ = mean; probs = [exp(-lam_) * lam_ ** i / factorial(i) for i in range(4)]; probs.append(1 - sum(probs))
obs = H[:4] + [sum(H[4:])]
chi2 = sum((o - len(P) * pr) ** 2 / (len(P) * pr) for o, pr in zip(obs, probs))
print(f"Poisson fit (bins 0,1,2,3,>=4): observed {obs}, expected {[round(len(P) * pr, 1) for pr in probs]}, chi-square {chi2:.2f}")
# dispersion test of the sporadic-class counts against the binomial model N_p ~ Bin(p-4, 1/p) (independent uniform residues)
mu = {p: (p - 4) / p for p in P}; sg = {p: (p - 4) * (1 / p) * (1 - 1 / p) for p in P}
Np = {p: sum(1 for c, k, tag in recs[p]['cls'] if c not in (1, 2, 3)) for p in P}
S1 = sum(Np[p] - mu[p] for p in P); S2 = sum((Np[p] - mu[p]) ** 2 for p in P); V = sum(sg.values())
# under the model E[S2] = V and Var[S2] = sum (E(N-mu)^4 - sg^2); for Bin(n, q): kappa4 = n q(1-q)(1 - 6q(1-q)), E(N-mu)^4 = kappa4 + 3 sg^2
var_S2 = sum(((p - 4) * (1 / p) * (1 - 1 / p) * (1 - 6 * (1 / p) * (1 - 1 / p)) + 3 * sg[p] ** 2) - sg[p] ** 2 for p in P)
print(f"binomial model: total observed {sum(Np.values())} vs expected {sum(mu.values()):.1f} (z = {S1 / V ** 0.5:.2f}); "
      f"sum of squared deviations {S2:.1f} vs expected {V:.1f}, sd {var_S2 ** 0.5:.1f}: dispersion z = {(S2 - V) / var_S2 ** 0.5:.2f}")
