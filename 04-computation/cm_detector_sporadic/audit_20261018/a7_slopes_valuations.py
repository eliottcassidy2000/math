"""Audit checks 5, 7, 8(valuations).

For (family, m) and every prime 7 <= p <= PMAX with p not | m:
  Delta_p(m) = sum_{k=p}^{2p-1} u_k m^{-k} - u_1 m^{-1} sum_{k<p} u_k m^{-k}; v = v_p(Delta) (computed mod p^CAP, capped);
  ordinary <=> F_p(1/m) != 0 mod p;  slope = V_p(1/m) / (u_1 q_p(m) F_p(1/m)) mod p at ordinary primes.
Check 5: (B1') [v >= 2]  <=>  V_p(1/m) == u_1 q_p(m) F_p(1/m) (mod p)  at every (m, p) evaluated.
"""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in, FAMILIES

CAP = 8
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 199

POINTS = {
    'A': {'CM': [-4, 2], 'cusp': [-1, 8], 'ctrl': [3, 5, 10, 100]},
    'B': {'CM': [3, 9], 'cusp': [], 'ctrl': [2, 5, 10, 100]},
    'C': {'CM': [-3, 3], 'cusp': [1, 9], 'ctrl': [2, 5, 10, 100]},
    'D': {'CM': [], 'cusp': [], 'ctrl': [2, 3, 5, 100]},
    'E': {'CM': [], 'cusp': [4, 8], 'ctrl': [2, 3, 5, 100]},
    'F': {'CM': [6, 12], 'cusp': [8, 9], 'ctrl': [2, 3, 5, 100]},
    'Apery': {'CM': [-1, 1], 'cusp': [], 'ctrl': [2, 3, 5, 100]},
    'Domb': {'CM': [-32, -8, -2, 1, 8, 64], 'sing': [4, 16], 'cusp': [], 'ctrl': [2, 3, 5, 100]},
    's7': {'CM': [-10648, -64, -8, 2, 8, 24, 125], 'sing': [-1, 27], 'cusp': [], 'ctrl': [3, 5, 11, 100]},
    's10': {'CM': [-324, -64, -20, -9, -2, 1, 4, 12, 36, 196, 1296, 5776], 'sing': [-4, 16], 'cusp': [], 'ctrl': [2, 3, 5, 100]},
    'A274786': {'CM': [-15228, -828, -28, -3, 4, 18, 36, 72, 147], 'cusp': [], 'ctrl': [3, 5, 11, 100]},
    'HV': {'attr': [-7], 'sing': [1, 9, 25], 'ctrl': [5, 2, 3, 100]},
    'H333': {'attr': [-5832], 'sing': [729], 'ctrl': [5, 2, 7, 100]},
    'B2k4': {'z=-1': [-256], 'sing': [256], 'ctrl': [5, 2, 3, 100]},
    'Quintic': {'sing': [3125], 'ctrl': [5, 2, 3, 7, 100]},
}


def vp(n, p):
    if n == 0:
        return CAP
    v = 0
    while n % p == 0 and v < CAP:
        n //= p
        v += 1
    return v


def analyse_point(u, m, primes):
    rows = []
    for p in primes:
        if m % p == 0:
            continue
        P = p ** CAP
        p2 = p * p
        # Delta * m^{2p-1} = sum_{k=p}^{2p-1} u_k m^{2p-1-k} - u_1 sum_{k<p} u_k m^{2p-2-k}
        Nn = 0
        for k in range(p, 2 * p):
            Nn = (Nn + u[k] * pow(m, 2 * p - 1 - k, P)) % P
        for k in range(p):
            Nn = (Nn - u[1] * u[k] * pow(m, 2 * p - 2 - k, P)) % P
        v = vp(Nn, p)
        xinv = pow(m, -1, p)
        F = sum(u[k] * pow(xinv, k, p) for k in range(p)) % p
        V = sum(((u[p + j] - u[1] * u[j]) // p) * pow(xinv, j, p) for j in range(p)) % p
        qp = ((pow(m, p - 1, p2) - 1) // p) % p
        lemmaA_rhs = (u[1] * qp * F) % p
        b1 = v >= 2
        lemma_ok = (b1 == (V == lemmaA_rhs))
        slope = None
        if F % p and qp % p and u[1] % p:
            slope = V * pow(u[1] * qp * F, -1, p) % p
        rows.append(dict(p=p, v=v, ordinary=bool(F % p), slope=slope, lemmaA=lemma_ok, V0=(V == 0)))
    return rows


if __name__ == '__main__':
    primes = primes_in(7, PMAX)
    out = {}
    lemma_total = lemma_bad = 0
    for fam, groups in POINTS.items():
        w = FAMILIES[fam][1]
        u = seq(fam, 2 * PMAX + 2)
        for g, ms in groups.items():
            for m in ms:
                rows = analyse_point(u, m, primes)
                lemma_total += len(rows)
                lemma_bad += sum(not r['lemmaA'] for r in rows)
                ordv = Counter(r['v'] for r in rows if r['ordinary'])
                nordv = Counter(r['v'] for r in rows if not r['ordinary'])
                slopes = Counter(r['slope'] for r in rows if r['ordinary'] and r['slope'] is not None)
                s1 = slopes.get(1, 0)
                ns = sum(slopes.values())
                V0_nonord = Counter(r['V0'] for r in rows if not r['ordinary'])
                nonord_primes = [r['p'] for r in rows if not r['ordinary']]
                print(f'{fam:8s} w={w} {g:5s} m={m:>7}: ordinary v_p {dict(sorted(ordv.items()))}; non-ordinary v_p '
                      f'{dict(sorted(nordv.items()))}; slope=1 at {s1}/{ns} ordinary primes; V(1/m)==0 at non-ord: '
                      f'{dict(V0_nonord)}', flush=True)
                if fam == 'B2k4' and m == -256:
                    print('          non-ordinary primes:', nonord_primes, ' all == 3 mod 4:', all(p % 4 == 3 for p in nonord_primes),
                          ' all p==3 mod 4 are non-ordinary:', all((p in nonord_primes) for p in primes if p % 4 == 3))
                if g in ('attr', 'z=-1') or (fam == 'HV'):
                    print('          v_p by prime:', {r['p']: r['v'] for r in rows})
                out[f'{fam}:{m}'] = rows
    print(f'\nLemma (A) check: (B1\') <=> V == u_1 q_p F (mod p) agreed at {lemma_total - lemma_bad}/{lemma_total} (m,p) pairs')
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'a7_slopes_valuations_p{PMAX}.json'), 'w') as f:
        json.dump(out, f)
