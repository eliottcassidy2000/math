"""Audit check 7b: inert-prime valuation law 'v_p(Delta_p) = w at non-ordinary primes, apart from ~1/p coincidences'
(HYP-9330 (F); title of HYP-9330).  Probe all weight-2 survivors for SYSTEMATIC excess at non-ordinary primes,
for 7 <= p <= PMAX; report per point the count of non-ordinary primes with v = 2 and v >= 3, and the CM field
(discriminant D whose inert set matches the non-ordinary set, among small D)."""
import sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in
from a7_slopes_valuations import analyse_point
from a6_polys_modp import kron

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
PTS = {
    'Apery': [-1, 1], 'Domb': [-32, -8, -2, 1, 4, 8, 16, 64],
    's7': [-10648, -64, -8, -1, 2, 8, 24, 27, 125],
    's10': [-324, -64, -20, -9, -4, -2, 1, 4, 12, 16, 36, 196, 1296, 5776],
    'A274786': [-15228, -828, -28, -3, 4, 18, 36, 72, 147],
}
DS = [-3, -4, -7, -8, -11, -12, -15, -16, -19, -20, -24, -27, -28, -32, -35, -36, -40, -43, -48, -51, -52, -60, -64,
      -67, -72, -75, -88, -91, -99, -100, -112, -115, -120, -123, -147, -148, -163, -187, -232, -235, -267, -403, -427]
primes = primes_in(7, PMAX)
for fam, ms in PTS.items():
    u = seq(fam, 2 * PMAX + 2)
    for m in ms:
        rows = analyse_point(u, m, primes)
        nonord = [r for r in rows if not r['ordinary']]
        c = Counter(min(r['v'], 3) for r in nonord)
        hi = [r['p'] for r in nonord if r['v'] >= 3]
        npr = [r['p'] for r in nonord]
        Dm = [D for D in DS if [p for p in primes if m % p and D % p and kron(D, p) == -1] == [p for p in npr if D % p]]
        exp = sum(1 / p for p in npr)
        flag = '  <== SYSTEMATIC' if len(hi) > max(3, 3 * exp) else ''
        print(f'{fam:8s} m={m:>7}: non-ordinary primes {len(npr)} (CM D in {Dm[:3]}); v=2 at {c.get(2, 0)}, v>=3 at {len(hi)} '
              f'{hi[:12]}{"..." if len(hi) > 12 else ""} (chance expectation sum 1/p = {exp:.2f}){flag}', flush=True)
