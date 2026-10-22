# AUDIT task 1 (independent of forest_numbers.c / drt23_forest_exact.py / drt23_classes.py):
#  (a) required spot check: own DP aud_fwdpoly (layered pull DP) for class 0 (D23_0) and class 12 (P_23) at x = 2, 5
#      mod p = 1000000007, against sum_k f_k x^k mod p from drt23_forest_exact.json;
#  (b) extra: full exact reconstruction of F(T, x) for ALL 13 classes from own evaluations at x = 0..22 modulo the two
#      primes 2^61 - 1 and 4611686018427387847 (own Lagrange interpolation + CRT), compared with the stored f_k;
#  (c) own conversion (n - L)! N_L = sum_k f_k [t^L] (1 + t)^k (1 - t)^(n-1-k) by explicit polynomial convolution;
#  (d) own TT5 (recursion TT_k(T) = sum_v TT_(k-1)(T[N+(v)]) on bitmasks), own 4-census check;
#  (e) Delta N_8 = 7168 Delta TT5 for the 12 non-Paley classes; D23_1 vs D23_2; n = 19 pair.
import json, subprocess, sys, time
from math import factorial
from fractions import Fraction

HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
cls = json.load(open(f'{HERE}/drt23_classes.json'))
fex = {r['class']: r for r in json.load(open(f'{HERE}/drt23_forest_exact.json'))}
EXE = f'{HERE}/audit_loop4/aud_fwdpoly'
P1, P2 = 2 ** 61 - 1, 4611686018427387847

def F_mod(code, p, xs):
    out = subprocess.run([EXE, str(p)] + [str(x) for x in xs], input=code + '\n', capture_output=True, text=True, check=True).stdout
    return {int(a): int(b) for a, b in (l.split() for l in out.split('\n') if l.strip())}

def lagrange_coeffs(vals, p):
    """coefficients of the unique poly of degree < m with poly(x) = vals[x], x = 0..m-1, mod p (own Lagrange)"""
    m = len(vals)
    coef = [0] * m
    for i in range(m):
        # basis_i(x) = prod_(j != i) (x - j) / (i - j)
        num = [1]
        den = 1
        for j in range(m):
            if j == i: continue
            num = [((num[k - 1] if k >= 1 else 0) - j * (num[k] if k < len(num) else 0)) % p for k in range(len(num) + 1)]
            den = den * (i - j) % p
        s = vals[i] * pow(den, -1, p) % p
        for k in range(m): coef[k] = (coef[k] + s * num[k]) % p
    return coef

def crt2(a, b):
    t = (b - a) * pow(P1, -1, P2) % P2
    return a + P1 * t

def code_adj(code):
    m = len(code); n = int(round((1 + (1 + 8 * m) ** 0.5) / 2))
    out = [0] * n; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    return n, out

def tt(k, mask, out):
    """number of transitive k-subsets of the subtournament on mask"""
    if k == 1: return bin(mask).count('1')
    s = 0; m = mask
    while m:
        v = (m & -m).bit_length() - 1; m &= m - 1
        sub = out[v] & mask
        if bin(sub).count('1') >= k - 1: s += tt(k - 1, sub, out)
    return s

def census4(n, out):
    """counts of the 4 isomorphism types of 4-tournaments by score sequence"""
    from itertools import combinations
    c = {}
    for X in combinations(range(n), 4):
        msk = sum(1 << v for v in X)
        sc = tuple(sorted(bin(out[v] & msk).count('1') for v in X))
        c[sc] = c.get(sc, 0) + 1
    return c

def poly_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return r

def forest_N(f, n):
    """own conversion: (n-L)! N_L = sum_k f_k [t^L] (1+t)^k (1-t)^(n-1-k)"""
    tot = [0] * n
    onep, onem = [1, 1], [1, -1]
    powp = [[1]]
    for _ in range(n): powp.append(poly_mul(powp[-1], onep))
    powm = [[1]]
    for _ in range(n): powm.append(poly_mul(powm[-1], onem))
    for k, fk in enumerate(f):
        pk = poly_mul(powp[k], powm[n - 1 - k])
        for L, c in enumerate(pk): tot[L] += fk * c
    N = []
    for L in range(n):
        q, r = divmod(tot[L], factorial(n - L))
        assert r == 0, ('not divisible', L)
        N.append(q)
    return N

if __name__ == '__main__':
    t0 = time.time()
    print("(a) required spot check: own DP, p = 1000000007, x = 2, 5")
    p = 1000000007
    for k in (0, 12):
        c = cls[k]; f = [int(x) for x in fex[k]['f']]
        got = F_mod(c['code'], p, [2, 5])
        for x in (2, 5):
            want = sum(fk * pow(x, e, p) for e, fk in enumerate(f)) % p
            print(f"   class {k} ({c['label']}): x = {x}: own F mod p = {got[x]}, stored sum f_k x^k mod p = {want}, equal: {got[x] == want}")
    print(f"   [{time.time() - t0:.0f}s]", flush=True)

    print("\n(b)-(e) full independent reconstruction for all 13 classes (x = 0..22, primes 2^61-1 and 4611686018427387847)")
    n = 23
    res = {}
    for c in cls:
        k = c['class']
        v1 = F_mod(c['code'], P1, range(n)); v2 = F_mod(c['code'], P2, range(n))
        c1 = lagrange_coeffs([v1[x] for x in range(n)], P1); c2 = lagrange_coeffs([v2[x] for x in range(n)], P2)
        f = [crt2(a, b) for a, b in zip(c1, c2)]
        assert all(0 <= x <= factorial(n) for x in f)
        stored = [int(x) for x in fex[k]['f']]
        N = forest_N(f, n)
        Nst = [int(x) for x in fex[k]['N']]
        nn, out = code_adj(c['code'])
        T5 = tt(5, (1 << n) - 1, out)
        kap = sum(Fraction(N[L], factorial(n) // factorial(n - L)) for L in range(n))
        H = int(c['H'])
        checks = dict(f_equal_stored=(f == stored), N_equal_stored=(N == Nst), f0_eq_H=(f[0] == H), palin=(f == f[::-1]),
                      sum_nfact=(sum(f) == factorial(n)), odd_zero=all(N[L] == 0 for L in range(1, n, 2)),
                      forest_identity=(kap == Fraction(H * 2 ** (n - 1), factorial(n))), TT5_equal_json=(T5 == c['TT5']))
        res[k] = dict(N=N, TT5=T5, label=c['label'], f=f)
        print(f"   {c['label']:7s}: TT5 own {T5} (json {c['TT5']}); N_2,4,6,8,10 = {N[2]}, {N[4]}, {N[6]}, {N[8]}, {N[10]}; "
              + ', '.join(f"{a}={b}" for a, b in checks.items()) + f"  [{time.time() - t0:.0f}s]", flush=True)
    iP = [k for k in res if res[k]['label'] == 'P_23'][0]
    print("\n(e) Delta N_8 against 7168 Delta TT5 (own N, own TT5), relative to P_23:")
    allok = True
    for k in sorted(res):
        if k == iP: continue
        dN8 = res[k]['N'][8] - res[iP]['N'][8]; dT = res[k]['TT5'] - res[iP]['TT5']
        ok = dN8 == 7168 * dT; allok &= ok
        print(f"   {res[k]['label']:7s}: Delta N_8 = {dN8:9d}, Delta TT5 = {dT:4d}, 7168 Delta TT5 = {7168 * dT:9d}, equal: {ok};  Delta N_2,4,6 = "
              f"{[res[k]['N'][L] - res[iP]['N'][L] for L in (2, 4, 6)]}")
    print(f"   all 12 equal: {allok}")
    print(f"   D23_1 vs D23_2: N_8 equal {res[1]['N'][8] == res[2]['N'][8]}, N_10 {res[1]['N'][10]} vs {res[2]['N'][10]} (differ: {res[1]['N'][10] != res[2]['N'][10]}), "
          f"first differing L: {min(L for L in range(n) if res[1]['N'][L] != res[2]['N'][L])}")
    print("   own 4-census (score sequences) of all 13 classes identical:", len({tuple(sorted(census4(*code_adj(c['code'])).items())) for c in cls}) == 1,
          census4(*code_adj(cls[iP]['code'])))
    json.dump({str(k): {'label': v['label'], 'TT5': v['TT5'], 'N': [str(x) for x in v['N']]} for k, v in res.items()},
              open(f'{HERE}/audit_loop4/aud_task1_N.json', 'w'), indent=0)
    print(f"done [{time.time() - t0:.0f}s]")
