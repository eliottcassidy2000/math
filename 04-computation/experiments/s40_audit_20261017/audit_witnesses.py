"""S40 audit: recompute every witness order of the HYP-9287 local-search JSONs with independent code, and re-derive the
pre-registered prime lists and verdicts.
Usage: python audit_witnesses.py
"""
import sys, os, math, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import sieve, legendre, h_forms

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results')


@njit(cache=True)
def leg_table(p):
    t = np.zeros(p, dtype=np.int64)
    for a in range(1, p):
        t[a] = legendre(a, p)
    return t


@njit(cache=True)
def excess(order, chi, p):
    s = 0
    for i in range(p):
        for j in range(i + 1, p):
            s += chi[(order[j] - order[i]) % p]
    return s


def Lval(p, h):
    w = 6 if p == 3 else 2
    return 2 * math.pi * h / (w * math.sqrt(p))


def main():
    out = {}
    issues = []
    for fn in ['paley_ranking_local_search_20261017.json', 'paley_ranking_local_search_p132_500_20261017.json']:
        rows = json.load(open(os.path.join(RES, fn)))
        recs = []
        for r in rows:
            p = r['p']
            chi = leg_table(p)
            order = np.array(r['order'], dtype=np.int64)
            perm_ok = sorted(r['order']) == list(range(p))
            ex = int(excess(order, chi, p))
            h = int(h_forms(-p))
            nat = int(excess(np.arange(p), chi, p))
            Ltrue = Lval(p, h)
            beaten = ex > nat
            ok = perm_ok and ex == r['best'] and h == r['h'] and nat == r['natural'] and beaten == r['beaten'] and abs(r['L'] - Ltrue) < 1e-9
            if not ok:
                issues.append(dict(file=fn, p=p, perm_ok=perm_ok, excess=ex, stored_best=r['best'], h=h, stored_h=r['h'], natural=nat,
                                   stored_natural=r['natural'], beaten=beaten, stored_beaten=r['beaten'], L_true=Ltrue, L_stored=r['L']))
            recs.append(dict(p=p, h=h, L=round(Ltrue, 4), natural=nat, witness_excess=ex, beaten=beaten,
                             MAS_lower_bound=(p * (p - 1) // 2 + ex) // 2))
        out[fn] = recs
        print(fn, ': %d rows recomputed' % len(recs))
        for x in recs:
            print('   ', x)
    print('ISSUES:', json.dumps(issues, indent=1))
    # pre-registered lists on 131 < p <= 500
    S = sieve(1000)
    pr = [p for p in range(132, 501) if S[p] and p % 4 == 3]
    Ls = {p: math.pi * h_forms(-p) / math.sqrt(p) for p in pr}
    hi = [p for p in pr if Ls[p] >= 2.2]
    lo = [p for p in pr if Ls[p] <= 1.7]
    mid = [(p, round(Ls[p], 3)) for p in pr if 1.7 < Ls[p] < 2.2]
    print('primes p = 3 mod 4 in (131, 500]:', len(pr))
    print('L >= 2.2:', hi)
    print('L <= 1.7:', lo, len(lo))
    print('1.7 < L < 2.2:', mid)
    pre_hi = [167, 191, 239, 263, 311, 359, 383, 431, 439, 479]
    pre_lo = [139, 163, 179, 211, 223, 227, 251, 283, 307, 331, 347, 367, 379, 419, 443, 463, 467, 487, 491, 499]
    print('pre-registered high list matches:', hi == pre_hi, '; low list matches:', lo == pre_lo)
    beat = {x['p']: x['beaten'] for x in out['paley_ranking_local_search_p132_500_20261017.json']}
    p1 = [p for p in pre_hi if beat[p]]
    p2 = [p for p in pre_lo if not beat[p]]
    print('P1: beaten among high-L primes:', p1, '-> PASS' if len(p1) <= 1 else '-> FAIL')
    print('P2: unbeaten among low-L primes:', p2, '-> PASS' if not p2 else '-> FAIL')
    # every prime in range: largest beaten L, smallest unbeaten L
    allr = out['paley_ranking_local_search_p132_500_20261017.json']
    bL = max(x['L'] for x in allr if x['beaten']); uL = min(x['L'] for x in allr if not x['beaten'])
    print('largest beaten L = %.4f, smallest unbeaten L = %.4f' % (bL, uL))
    # correct L(1, chi_-3)
    print('L(1, chi_-3) = pi/(3 sqrt 3) = %.4f (session value 1.8138 = pi/sqrt 3)' % (math.pi / (3 * math.sqrt(3))))
    json.dump(dict(rows=out, issues=issues, hi=hi, lo=lo, mid=mid, P1_beaten=p1, P2_unbeaten=p2), open(os.path.join(HERE, 'audit_witnesses.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
