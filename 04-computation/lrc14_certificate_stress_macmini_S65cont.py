#!/usr/bin/env python3
"""
lrc14_certificate_stress_macmini_S65cont.py -- HYP-5730 phase 2: STRESS C1 (ledger) + C2 (descent)

Phase-1 result: C1 covers 100% of covering [1,18]; C2 covers 100% of random covering cap 60;
adversarial min #certified = 3 (never 0).  Phase 2:
  (A) structured adversaries: covering near-APs (dense blocks), monad-S2's detuned harmonics,
      prime-sum-rich sets (starve C2), the phase-1 min-certified sets as seeds;
  (B) longer, targeted hill-climbs at caps 60 / 100 / 150 with objective = #certified,
      tie-broken by total certificate margin;
  (C) firing-profile statistics: WHICH ruler fires C1 (its r(q), gcd profile, q/Vmax), and which
      k fires C2 -- the raw material for the a-priori theorem.
All exact integer arithmetic.
"""
from itertools import combinations
from math import gcd
from functools import reduce
import random

random.seed(67)

def covering(S):
    return all(any(v % k == 0 for v in S) for k in range(2, 15))

def primitive(S):
    return reduce(gcd, S) == 1

def rulers_of(S):
    return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})

def c1_data(S, q):
    """(fires, ledger_sum, n_classes, n_reps, cheap_count)"""
    m = -(-q // 14) - 1
    classes = {}
    for v in S:
        r = v % q
        if r == 0:
            return (False, None, None, None, None)
        key = min(r, q - r)
        classes.setdefault(key, gcd(r, q))
    total = sum(g * (2 * (m // g) + 1) - 1 for g in classes.values())
    nreps = 13 - len(classes)
    cheap = sum(1 for g in classes.values() if g > m)
    return (total < q - 1, total, len(classes), nreps, cheap)

def c2_data(S, q):
    if any(v % q == 0 for v in S):
        return None
    for k in sorted(kk for kk in range(15, q) if q % kk == 0):
        for s in range(1, k):
            if all(14 * (v * s % k) >= k and 14 * (k - v * s % k) >= k for v in S):
                return k
    return None

def live_exact(S, q):
    for p in range(1, q):
        if all(14 * (v * p % q) >= q and 14 * (q - v * p % q) >= q for v in S):
            return p
    return None

def ncert(S):
    tot = 0
    for q in rulers_of(S):
        f1 = c1_data(S, q)[0]
        f2 = (c2_data(S, q) is not None) if not f1 else False   # short-circuit: count certified rulers
        tot += bool(f1 or f2)
    return tot

# ---------------------------------------------------------------- (A) structured adversaries
print("=" * 100)
print("(A) STRUCTURED ADVERSARIES vs C1+C2")
print("=" * 100)

def report(name, S):
    S = sorted(S)
    if len(set(S)) != 13 or not covering(S) or not primitive(S):
        print(f"  {name}: SKIP (not a primitive covering 13-set) {S}")
        return
    certified, live_any = [], []
    for q in rulers_of(S):
        f1, ledger, ncls, nreps, cheap = c1_data(S, q)
        k2 = c2_data(S, q)
        if f1 or k2:
            certified.append((q, 'C1' if f1 else f'C2:k={k2}'))
    if not certified:
        # residual! find live rulers exactly
        for q in rulers_of(S):
            p = live_exact(S, q)
            if p:
                live_any.append((q, p))
        print(f"  {name}: *** RESIDUAL *** live rulers {live_any[:5]}  S={S}")
    else:
        print(f"  {name}: certified x{len(certified)}, first {certified[:3]}  S={S}")

# monad-S2 detuned harmonics
for repl in (83, 85):
    S = [14 * k for k in range(1, 14)]; S[5] = repl
    report(f"detuned-harmonic 84->{repl}", S)
# covering near-APs: dense block + patches for missing divisors
report("near-AP block 7..19 patched", [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 22])
report("near-AP block 2..14", [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
report("min-supply set S65", [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14])
report("phase1 adversarial min", [7, 11, 14, 17, 22, 23, 24, 25, 26, 27, 29, 30, 39])
report("worst7Struct @91", sorted(91 - e for e in [0, 7, 14, 21, 26, 29, 37, 44, 51, 58, 67, 75, 82]))
# dilated-AP-adjacent: c*{1..13} not primitive; perturb two elements to force primitivity+covering
S = [2 * k for k in range(1, 14)]; S[0] = 1; S[6] = 15   # {1,4,6,8,10,12,15,16,18,20,22,24,26}
report("even-dilate 2-perturbed", S)

# ---------------------------------------------------------------- (B) hard hill-climbs
print()
print("=" * 100)
print("(B) HARD ADVERSARIAL: minimize #certified rulers (caps 60/100/150), seeded + long")
print("=" * 100)
seeds = [
    [7, 11, 14, 17, 22, 23, 24, 25, 26, 27, 29, 30, 39],
    [6, 11, 19, 23, 24, 25, 26, 27, 28, 31, 34, 37, 40],
    [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14],
]
overall = None
for cap in (60, 100, 150):
    best = None
    starts = [s for s in seeds if max(s) <= cap]
    while len(starts) < 5:
        S = sorted(random.sample(range(1, cap + 1), 13))
        if covering(S) and primitive(S):
            starts.append(S)
    for si, S0 in enumerate(starts):
        S, cur = list(S0), ncert(S0)
        stall = 0
        for step in range(300):
            T = list(S)
            if random.random() < 0.5:
                T[random.randrange(13)] = random.randrange(1, cap + 1)
            else:                                   # local move: +-1..3 on one speed
                idx = random.randrange(13)
                T[idx] = max(1, min(cap, T[idx] + random.choice([-3, -2, -1, 1, 2, 3])))
            T = sorted(set(T))
            if len(T) != 13 or not covering(T) or not primitive(T):
                continue
            nc = ncert(T)
            if nc <= cur:
                if nc < cur: stall = 0
                S, cur = T, nc
            stall += 1
            if cur == 0:
                break
        if best is None or cur < best[0]:
            best = (cur, S)
        if cur == 0:
            break
    print(f"cap {cap}: adversarial min #certified = {best[0]}  S = {best[1]}")
    if best[0] == 0:
        report(f"cap-{cap} DEFEATER", best[1])
    if overall is None or best[0] < overall[0]:
        overall = best
print(f"OVERALL adversarial min #certified rulers = {overall[0]}")

# ---------------------------------------------------------------- (C) firing profile
print()
print("=" * 100)
print("(C) FIRING PROFILE over the 966 covering [1,18]: the ruler that fires C1 -- why")
print("=" * 100)
prof = {"reps": [], "cheap": [], "qq": [], "firstq_over_V": []}
for S in combinations(range(1, 19), 13):
    if not (covering(S) and primitive(S)):
        continue
    S = list(S)
    for q in rulers_of(S):
        f1, ledger, ncls, nreps, cheap = c1_data(S, q)
        if f1:
            prof["reps"].append(nreps)
            prof["cheap"].append(cheap)
            prof["qq"].append(q)
            prof["firstq_over_V"].append(q / S[-1])
            break
from collections import Counter
print(f"first-firing ruler: reps distribution {Counter(prof['reps'])}")
print(f"                    cheap-count distribution {Counter(prof['cheap'])}")
print(f"                    q/Vmax: min {min(prof['firstq_over_V']):.3f} "
      f"median {sorted(prof['firstq_over_V'])[len(prof['firstq_over_V'])//2]:.3f} "
      f"max {max(prof['firstq_over_V']):.3f}")
print(f"                    q values {Counter(prof['qq']).most_common(8)}")
print()
print("Done.")
