#!/usr/bin/env python3
"""audit_C item 4: does the partner criterion (equal odd count, Terras difference K - K') agree with a DIRECT equal-time
merge of the orbits of M_K = 2^K - 1 and M_K' (K' = K - D), before 1?

Direct test: the deletion chain aligns T^D(M_K) = 3^D 2^K' - 1 with M_K' (state (D, 3^D - 1), anchored at -1).
Run both Terras orbits in lockstep with the debt k (start D, increment par(a) - par(b)) until they coincide or one of them
reaches 1. Record: coincidence time, debt at coincidence, coincidence value, first common odd value Z.
  * 'merge before 1' (Terras sense) = coincidence at a value >= 2 with zero debt (pair-chain absorption);
  * 'U-sense merge' (THM-4556: U^j(x) = U^j(y) != 1) additionally needs Z != 1.
Samples:
  (P) 90 random non-orphan odd K in [3, 12800]: the nearest partner and one random further partner;
  (N) random (K, D) with D NOT a partner (D <= 9), to test the converse; (N2) 120 more with D uniform in [1, K-2];
  (O) exhaustive: every D in [1, K-2] for the orphans K = 1039, 1081, 1113.
Also: a proof that the criterion is EXACTLY 'absorption at a value >= 8' for orbits that reach 1 is in the audit report.
"""
import random, time
from collections import defaultdict

data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split())
    data[K] = (o, t)
by_odd = defaultdict(list)
partners = {}
for K in range(2, 12801):
    o, t = data[K]
    ps = [Kp for Kp in by_odd[o] if t - data[Kp][1] == K - Kp]
    if K % 2 == 1 and K >= 3:
        partners[K] = ps
    by_odd[o].append(K)


def lockstep(K, D):
    Kp = K - D
    a = 3 ** D * (1 << Kp) - 1
    b = (1 << Kp) - 1
    k = D
    s = 0
    while True:
        if a == b:
            z = a
            while z % 2 == 0:
                z //= 2
            return dict(coincide=True, time=s, debt=k, value=a, Z=z)
        if a == 1 or b == 1:
            return dict(coincide=False, time=s, debt=k, value=None, Z=None, which=('a' if a == 1 else 'b'))
        k += (a & 1) - (b & 1)
        a = (3 * a + 1) >> 1 if a & 1 else a >> 1
        b = (3 * b + 1) >> 1 if b & 1 else b >> 1
        s += 1


rnd = random.Random(77)
t0 = time.time()
# (P)
nonorph = [K for K in partners if partners[K]]
sampP = rnd.sample([K for K in nonorph if K < 3000], 40) + rnd.sample([K for K in nonorph if K >= 3000], 50)
okP = badP = Zone = 0
nonzero_debt_coinc = 0
mintime = []
for K in sampP:
    ps = partners[K]
    tests = {max(ps)}
    if len(ps) > 1:
        tests.add(rnd.choice(ps))
    for Kp in tests:
        r = lockstep(K, K - Kp)
        if r['coincide'] and r['debt'] == 0 and r['value'] >= 2:
            okP += 1
            Zone += (r['Z'] == 1)
            mintime.append(r['value'])
        else:
            badP += 1
            print("criterion says partner but no merge:", K, Kp, r)
print(f"(P) partner pairs tested: {okP + badP}; merge (zero debt, value >= 2): {okP}; failures {badP}; "
      f"first common odd value Z = 1 in {Zone} of them; min merge value {min(mintime)}")
# (N)
okN = badN = 0
for _ in range(120):
    K = rnd.choice(list(partners))
    if K < 20:
        continue
    cand = [D for D in range(1, 10) if K - D >= 2 and (K - D) not in partners[K]]
    if not cand:
        continue
    D = rnd.choice(cand)
    r = lockstep(K, D)
    if r['coincide'] and r['value'] >= 2:
        badN += 1
        print("criterion says NOT partner but coincidence:", K, D, r)
    else:
        okN += 1
print(f"(N) non-partner pairs tested: {okN + badN}; no coincidence before 1: {okN}; contradictions {badN}")
# (N2) non-partner D drawn from the whole range [1, K-2]
okN2 = badN2 = 0
while okN2 + badN2 < 120:
    K = rnd.choice([K for K in partners if 50 <= K <= 8000])
    D = rnd.randint(1, K - 2)
    if (K - D) in partners[K]:
        continue
    r = lockstep(K, D)
    if r['coincide'] and r['value'] >= 2:
        badN2 += 1
        print("criterion says NOT partner but coincidence:", K, D, r)
    else:
        okN2 += 1
print(f"(N2) non-partner pairs with D uniform in [1, K-2]: {okN2 + badN2}; no coincidence before 1: {okN2}; contradictions {badN2}")
# (O)
for K in (1039, 1081, 1113):
    assert not partners[K]
    hits = []
    for D in range(1, K - 1):
        r = lockstep(K, D)
        if r['coincide'] and r['value'] >= 2:
            hits.append((D, r['debt'], r['value']))
    print(f"(O) orphan K = {K}: coincidences before 1 over all D in [1, {K-2}]: {hits}")
print(f"time {time.time()-t0:.1f}s")
