# checks.py -- exact checks of the lift lemma, the bijection Z_k = #{zeroless k-digit multiples of 2^k},
# the recursion Z_{k+1} = (9 Z_k + Delta_k)/2, the doubling criterion, and the unit/trace identity for small m.
import json, cmath, math, itertools
ok = True
def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name); ok &= bool(cond)
# (1) lift lemma: for k=1..7, n in a full period with n >= k+1, the 5 lifts n + j T_k have (k+1)-th digits = one parity class
for k in range(1, 8):
    T = 4 * 5**(k-1); M = 10**(k+1); good = True
    base = k + 1 + T  # start well past k+1
    for r in range(T):
        n = base + r
        digs = [(pow(2, n + j*T, M) // 10**k) % 10 for j in range(5)]
        lows = {pow(2, n + j*T, 10**k) for j in range(5)}
        par = digs[0] % 2
        if sorted(digs) != list(range(par, 10, 2)) or len(lows) != 1: good = False; break
    check('lift lemma k=%d (all %d residues): lifts run through the 5 digits of one parity, lower k digits fixed' % (k, T), good)
# (2) bijection: residues vs zeroless multiples
Z = {int(a): int(b) for a, b in json.load(open('zk_values.json'))['Z'].items()}
for k in range(1, 8):
    T = 4 * 5**(k-1); M = 10**k
    res = sum(1 for r in range(T) if '0' not in str(pow(2, k + T + r, M)).zfill(k))
    mult = sum(1 for x in range(2**k, M, 2**k) if '0' not in str(x).zfill(k))
    check('k=%d: #zeroless residues = #zeroless k-digit multiples of 2^k = Z_k = %d' % (k, Z[k]), res == mult == Z[k])
# (3) OEIS A181610 b-file agreement (n=1..26)
bfile = [4,18,81,364,1638,7371,33170,149268,671701,3022653,13601945,61208743,275439346,1239477074,5577646830,25099410745,
 112947348510,508263067945,2287183805359,10292327123878,46315472056678,208419624257654,937888309161430,4220497391215744,
 18992238260465327,85465072172060901]
check('Z_1..Z_26 agree with OEIS A181610 b-file (Yamanouchi)', all(Z[i+1] == v for i, v in enumerate(bfile)))
# (4) recursion parity counts from DFS for k<=12 (independent small DFS)
def dfs_counts(K):
    lev = [0]*(K+1); odd = [0]*(K+1); st = [(0, 0)]
    while st:
        i, q = st.pop(); lev[i] += 1; odd[i] += q & 1
        if i < K:
            for a in range(1 if q & 1 else 2, 10, 2): st.append((i+1, (q + a*5**i)//2))
    return lev, odd
lev, odd = dfs_counts(11)
check('Z_{k+1} = 4E_k + 5O_k = (9Z_k + Delta_k)/2 for k<=10', all(lev[k+1] == 4*(lev[k]-odd[k]) + 5*odd[k] for k in range(11)))
# (5) trace identity: G_m = sum_{t odd} P_m(zeta^t) = 2^(m-1) Delta_{m-1}, and |N(P_m(zeta))| = 1 (unit), for m <= 12
for m in range(1, 13):
    zs = [cmath.exp(2j*math.pi*t/2**m) for t in range(1, 2**m, 2)]
    vals = []
    for z in zs:
        p = 1
        for i in range(m):
            w = z**(10**i % 2**m) if True else None
            p *= sum(w**a for a in range(1, 10))
        vals.append(p)
    G = sum(vals).real; N = 1.0
    for v in vals: N *= abs(v)
    D = 2*Z[m] - 9*Z[m-1]
    check('m=%2d: Tr P_m(zeta_2^m) = %d = 2^(m-1) Delta_(m-1); |Norm| = %.12f (unit)' % (m, round(G), N), abs(G - 2**(m-1)*D) < 1e-6 and abs(N - 1) < 1e-8)
# (6) doubling criterion: for zeroless x, 2x is zeroless iff x contains none of 51,52,53,54
good = True
for x in range(1, 200000):
    s = str(x)
    if '0' in s: continue
    crit = not any(t in s for t in ('51','52','53','54')) and not s.endswith('5')
    if ('0' not in str(2*x)) != crit: good = False; print('counterexample', x); break
check('doubling criterion on zeroless x < 2e5: 2x zeroless iff x has no 51,52,53,54 and does not end in 5', good)
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
