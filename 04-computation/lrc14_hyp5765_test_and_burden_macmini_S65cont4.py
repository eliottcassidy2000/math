#!/usr/bin/env python3
"""
lrc14_hyp5765_test_and_burden_macmini_S65cont4.py -- test HYP-5765 AS STATED + THM-676 burden

PART 1 -- HYP-5765 STRESS TEST (before any proof attempt: is it even true as stated?).
  HYP-5765: no primitive covering 13-set can simultaneously
    (a) block/kill every window descent k in [15,28] dividing a pair sum, AND
    (b) block/kill every PRIME p >= 29 dividing a pair sum.
  SUSPECTED DEFECT (self-flagged): (b) omits COMPOSITE k >= 29 (the cont-3 dodger was caught at
  k = 49 = 7^2), and smooth-sum adversaries could starve (b).  Adversarial hunt decides.

PART 2 -- THM-676 (the descent-burden theorem) verification:
  (i)  some parity class has >= 7 members  =>  >= 21 even pair sums;
  (ii) their halves h = q/2 take >= 11 DISTINCT values (Freiman |A+A| >= 2|A|-3 for i<j sums);
       h > 14 whenever q > 28: >= 11 unavoidable proper descent moduli;
  (iv) equality (exactly 11) iff the majority-parity class is an AP;
  (v)  no-primality-escape: any COMPOSITE pair sum q > 196 has a proper divisor > 14
       (q = ab, b >= sqrt(q) > 14); even sums always have h = q/2 > 14 once q > 28.
All exact integers.
"""
from math import gcd, isqrt
from functools import reduce
import random

random.seed(76)

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))

def danger_j(k):
    return -(-k // 14) - 1

def blocked(R, k):
    j = danger_j(k)
    bad = set(range(0, j + 1)) | set(range(k - j, k))
    for s in range(1, k):
        if all((r * s) % k not in bad for r in R):
            return False
    return True

def covering(S):  return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1
def rulers_of(S): return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})

def hyp5765_open(S):
    """#open descents among [(a) window k in [15,28]] u [(b) primes p >= 29 dividing sums]."""
    n = 0
    checked = set()
    for q in rulers_of(S):
        # (a) window divisors
        for k in range(15, 29):
            if q % k == 0 and q > k and k not in checked:
                checked.add(k)
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    n += 1
        # (b) prime divisors >= 29
        qq = q
        for p in range(2, isqrt(q) + 1):
            while qq % p == 0:
                if p >= 29 and p < q and p not in checked:
                    checked.add(p)
                    R = [v % p for v in S]
                    if 0 not in R and not blocked(R, p):
                        n += 1
                qq //= p
        if qq >= 29 and qq < q and qq not in checked:   # remaining prime factor
            checked.add(qq)
            R = [v % qq for v in S]
            if 0 not in R and not blocked(R, qq):
                n += 1
    return n

def full_capture(S):
    """If (a)&(b) fully dodged: what still catches S? composite k >= 29, or exact-only."""
    for q in rulers_of(S):
        for k in range(15, q):
            if q % k == 0:
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    kind = 'window' if k <= 28 else ('prime>=29' if is_prime(k) else 'COMPOSITE>=29')
                    return f"C2 k={k} ({kind}, q={q})"
    for q in rulers_of(S):
        for p in range(1, q):
            if all(14 * (v * p % q) >= q and 14 * (q - v * p % q) >= q for v in S):
                return f"exact-only q={q} p={p}"
    return "NOT LONELY?!"

print("=" * 100)
print("PART 1 -- HYP-5765 STRESS TEST: hunt covering sets achieving (a) AND (b)")
print("=" * 100)
refuted = None
for cap in (150, 250, 400):
    best = None
    for restart in range(6):
        while True:
            S = sorted(random.sample(range(1, cap + 1), 13))
            if covering(S) and primitive(S):
                break
        cur = hyp5765_open(S)
        for step in range(220):
            T = list(S)
            T[random.randrange(13)] = random.randrange(1, cap + 1)
            T = sorted(set(T))
            if len(T) != 13 or not covering(T) or not primitive(T):
                continue
            sc = hyp5765_open(T)
            if sc <= cur:
                S, cur = T, sc
            if cur == 0:
                break
        if best is None or cur < best[0]:
            best = (cur, S)
        if cur == 0:
            break
    print(f"cap {cap}: min open[(a)+(b)] = {best[0]}  S = {best[1]}")
    if best[0] == 0:
        cap_diag = full_capture(best[1])
        print(f"   *** (a) AND (b) BOTH ACHIEVED -- HYP-5765 REFUTED AS STATED ***")
        print(f"   capture beyond the conjecture's scope: {cap_diag}")
        refuted = (best[1], cap_diag)
        break
print(f"HYP-5765 as stated: {'REFUTED -- needs composite moduli in (b)' if refuted else 'survived this hunt (no (a)&(b) achiever found)'}")

# ================================================================ PART 2: burden theorem
print()
print("=" * 100)
print("PART 2 -- THM-676 descent-burden verification")
print("=" * 100)
# (ii)+(iv): distinct half-sums >= 11, equality iff AP -- exhaustive-ish check on 7-subsets
viol_lb = viol_eq = 0
for trial in range(200000):
    A = sorted(random.sample(range(1, 300), 7))
    sums = {A[i] + A[j] for i in range(7) for j in range(i + 1, 7)}
    d = len(sums)
    isAP = len({A[i+1] - A[i] for i in range(6)}) == 1
    if d < 11:
        viol_lb += 1
    if (d == 11) != isAP:
        viol_eq += 1
print(f"200k random 7-sets: #distinct pairwise sums >= 11 violations: {viol_lb}; "
      f"(= 11 <=> AP) violations: {viol_eq}")

# (v): composite q > 196 has a proper divisor > 14 (q = p*(q/p), q/p >= sqrt(q) > 14)
viol_c = 0
for q in range(197, 20000):
    if not is_prime(q):
        p = next(d for d in range(2, isqrt(q) + 1) if q % d == 0)
        if q // p <= 14:
            viol_c += 1
print(f"composite q in (196, 20000]: proper-divisor>14 violations: {viol_c}")

# burden distribution on covering sets
print()
print("Burden on covering sets (cap 250): distinct half-sum moduli h > 14 from the majority parity")
from collections import Counter
dist = Counter()
apcount = 0
for _ in range(300):
    while True:
        S = sorted(random.sample(range(1, 251), 13))
        if covering(S) and primitive(S):
            break
    evens = [v for v in S if v % 2 == 0]
    odds = [v for v in S if v % 2 == 1]
    maj = evens if len(evens) >= len(odds) else odds
    hs = {(maj[i] + maj[j]) // 2 for i in range(len(maj)) for j in range(i + 1, len(maj))}
    hs = {h for h in hs if h > 14}
    dist[len(hs)] += 1
lo, hi = min(dist), max(dist)
import statistics
print(f"forced half-sum moduli: min {lo}, max {hi} over 300 covering sets "
      f"(theorem floor: 11 when majority class has 7 members and sums > 28)")
print()
print("Done.")
