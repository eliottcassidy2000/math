"""AUDIT2: compare the exhaustive labelled census (a2_census.c, n <= 8) with
 - THM-5042: E_odd(n) = {odd k <= 2^n - 1} minus {1, 5, 37};
 - THM-5040: E(n) = O(n) (cross-check only);
 - HYP-9610 (A): E0(n) = {odd k <= M0(n)} minus {3, 7, 13, M0(n) - 2} for n >= 5;
 - base_sets.json (the packet's census extraction) E and E0;
 - the exclusion facts: realisers of 1, 5, 37 (connected realisers, classes, parity).
usage: python a2_census_analyse.py NMAX
"""
import json, sys, os
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'egcp_20261109', 'base_sets.json')
js = json.load(open(BASE))

def load(n):
    rows = {}
    for line in open(os.path.join(HERE, 'a2_census_n%d.out' % n)):
        f = line.split()
        if f and f[0] == 'V':
            k = int(f[1]); v = list(map(int, f[2:]))
            rows[k] = dict(lab_odd=v[0], iso_odd=v[1], lab_even=v[2], iso_even=v[3], lab_even0=v[4], iso_even0=v[5],
                           lab_codd=v[6], iso_codd=v[7], lab_ceven=v[8], iso_ceven=v[9])
    return rows

def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def M0(n): return 2**n - n - 1 if n % 2 == 0 else 2**n - n - 4
def O(n):
    s = set(range(1, M(n) + 1, 2)) - {3, 7, 13}
    if n % 2 and n >= 5: s.discard(2**n - n - 2)
    return s

ok = True
allodd_cum = set(); conn37 = 0; conn5 = 0
for n in range(1, NMAX + 1):
    R = load(n)
    odd_vals = {k for k, r in R.items() if r['lab_odd']}
    even_vals = {k for k, r in R.items() if r['lab_even']}
    even0_vals = {k for k, r in R.items() if r['lab_even0']}
    allodd_cum |= odd_vals
    tgt_odd = set(range(1, 2**n, 2)) - {1, 5, 37}
    r1 = (odd_vals == tgt_odd) and (allodd_cum == tgt_odd)
    r2 = even_vals == O(n)
    r3 = even_vals == set(js['E'][str(n)]) if str(n) in js['E'] else None
    r4 = even0_vals == set(js['E0'][str(n)]) if str(n) in js['E0'] else None
    if n >= 5:
        tgt0 = set(range(1, M0(n) + 1, 2)) - {3, 7, 13, M0(n) - 2}
        r5 = even0_vals == tgt0
        extra = sorted(even0_vals - tgt0); miss = sorted(tgt0 - even0_vals)
    else:
        r5 = None; extra = miss = []
    # exclusion facts
    c37 = R.get(37, {}).get('iso_ceven', 0); c37o = R.get(37, {}).get('iso_codd', 0)
    c5 = R.get(5, {}).get('iso_ceven', 0); c5o = R.get(5, {}).get('iso_codd', 0)
    conn37 += c37; conn5 += c5
    odd_at_exc = [k for k in (1, 5, 37) if R.get(k, {}).get('lab_odd', 0)]
    print("n=%d  E_odd == odd<=2^n-1 minus{1,5,37}: %s | E == O(n): %s | E == base_sets E: %s | E0 == base_sets E0: %s | "
          "E0 == HYP-9610 formula: %s %s %s | conn. realisers of 37: even %d odd %d; of 5: even %d odd %d | odd graphs at {1,5,37}: %s"
          % (n, r1, r2, r3, r4, r5, ('extra ' + str(extra)) if extra else '', ('missing ' + str(miss)) if miss else '',
             c37, c37o, c5, c5o, odd_at_exc))
    ok &= r1 and r2 and (r3 is not False) and (r4 is not False) and (r5 is not False) and not odd_at_exc and c37o == 0 and c5o == 0
    if n >= 5:
        print("     M0(n) = %d, M0-2 = %d in E0: %s; M0 in E0: %s; max E0 = %d" % (M0(n), M0(n) - 2, (M0(n) - 2) in even0_vals,
              M0(n) in even0_vals, max(even0_vals)))
print("connected realisers of 37 over n <= %d (iso classes): %d (all even); of 5: %d" % (NMAX, conn37, conn5))
# base_sets.json E0 at n = 8, 9 against the HYP-9610 formula (consistency of the packet's own data)
for n in (8, 9):
    if str(n) in js['E0']:
        tgt0 = set(range(1, M0(n) + 1, 2)) - {3, 7, 13, M0(n) - 2}
        b = set(js['E0'][str(n)])
        print("base_sets.json E0[%d] == HYP-9610 formula: %s (extra %s, missing %s)" % (n, b == tgt0, sorted(b - tgt0)[:10], sorted(tgt0 - b)[:10]))
print("ALL:", "PASS" if ok else "FAIL")
