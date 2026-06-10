#!/usr/bin/env python3
"""
taxicab_inert_audit_kpc1.py
THREAD A (kind-pasteur-2026-06-10-S1): HYP-2367 / THM-461 -- full independent audit.

Companion to taxicab_moser_bridge_kpc1.py.  That script established the
doubly-primitive 2-rep list <= 10^12 and tested inert primes THROUGH the
theorem (inert p | n  =>  p | gcd(d_i)).  This script removes the circularity:
it regenerates the list, then DIRECTLY factors every doubly-primitive n and
every cofactor, and audits, for every primitive representation n = x^3 + y^3
(gcd(x,y) = 1, d = x + y, m = n/d = x^2 - xy + y^2):

  AUDIT-1 (split lemma, THM-461.B): m = 3^eps * (split primes only), eps <= 1,
          and eps = 1 iff 3 | d.
  AUDIT-2 (Cor C exact form): for every inert prime p (p == 2 mod 3) of n,
          v_p(d) = v_p(n)  (the FULL inert content of n sits in EVERY
          primitive d; hence in gcd of the d's).
  AUDIT-3: v_3(n) - v_3(d) <= 1 for every primitive rep.
  AUDIT-4: recount of "n containing an inert prime" by direct factorization,
          to compare against the gcd-derived count (1952 expected).

Exact integer arithmetic; numpy int64 only for the collision enumeration
(max 2*10^12 << 2^63), all hits re-derived in pure ints.
"""

import time
from math import isqrt, gcd
from array import array

T0 = time.time()


def icbrt_floor(m):
    if m < 0:
        raise ValueError
    r = int(round(m ** (1.0 / 3.0)))
    while r * r * r > m:
        r -= 1
    while (r + 1) ** 3 <= m:
        r += 1
    return r


# ---- primes <= 10^6 (sieve) -------------------------------------------------
M = 10 ** 6
sieve = bytearray([1]) * (M + 1)
sieve[0] = sieve[1] = 0
i = 2
while i * i <= M:
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    i += 1
PRIMES = [p for p in range(2, M + 1) if sieve[p]]
print(f"primes <= 10^6: {len(PRIMES)}   [{time.time()-T0:.1f}s]")


def factorize(n):
    fac = []
    m = n
    for p in PRIMES:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            fac.append((p, e))
    if m > 1:
        fac.append((m, 1))  # prime: m <= 10^12 with no factor <= 10^6 = sqrt(10^12)
    return fac


# ---- regenerate doubly-primitive 2-rep list <= 10^12 ------------------------
import numpy as np

N_C = 10 ** 12
X_C = 10 ** 4
cb = np.arange(0, X_C + 1, dtype=np.int64) ** 3
chunks = []
for x in range(1, X_C + 1):
    rem = N_C - x * x * x
    if rem < 1:
        break
    ymax = min(x, icbrt_floor(rem))
    if ymax >= 1:
        chunks.append(x * x * x + cb[1:ymax + 1])
S = np.concatenate(chunks)
del chunks
S.sort()
dup_vals = np.unique(S[1:][S[1:] == S[:-1]])
del S
cube_set = set(i * i * i for i in range(1, X_C + 1))


def positive_reps(n):
    out = []
    xhi = icbrt_floor(n - 1)
    xlo = max(1, icbrt_floor(n // 2) - 1)
    for x in range(xlo, xhi + 1):
        r = n - x * x * x
        if 1 <= r <= x * x * x and r in cube_set:
            y = icbrt_floor(r)
            if y * y * y == r and 1 <= y <= x:
                out.append((x, y))
    return out


dp = []
for v in dup_vals:
    n = int(v)
    prr = [(x, y) for (x, y) in positive_reps(n) if gcd(x, y) == 1]
    if len(prr) >= 2:
        dp.append((n, prr))
dp.sort()
print(f"doubly-primitive 2-rep n <= 10^12: {len(dp)} (expect 5464)   [{time.time()-T0:.1f}s]")

# ---- the audits --------------------------------------------------------------
viol1 = viol2 = viol3 = 0
n_with_inert_direct = 0
inert_prime_hist = {}
total_prim_reps = 0
for n, prr in dp:
    nf = factorize(n)
    inert_of_n = [(p, e) for p, e in nf if p % 3 == 2]
    if inert_of_n:
        n_with_inert_direct += 1
        for p, e in inert_of_n:
            inert_prime_hist[p] = inert_prime_hist.get(p, 0) + 1
    v3n = next((e for p, e in nf if p == 3), 0)
    for (x, y) in prr:
        total_prim_reps += 1
        d = x + y
        m = n // d
        mf = factorize(m)
        # AUDIT-1: split lemma on the cofactor
        eps = next((e for p, e in mf if p == 3), 0)
        bad_inert_in_m = [p for p, e in mf if p % 3 == 2]
        if bad_inert_in_m or eps > 1 or ((eps == 1) != (d % 3 == 0)):
            viol1 += 1
            print(f"  AUDIT-1 VIOLATION n={n} rep=({x},{y}) m={m} fac={mf}")
        # AUDIT-2: full inert content of n in d
        for p, e in inert_of_n:
            vd = 0
            dd = d
            while dd % p == 0:
                dd //= p
                vd += 1
            if vd != e:
                viol2 += 1
                print(f"  AUDIT-2 VIOLATION n={n} rep=({x},{y}) inert p={p} "
                      f"v_p(n)={e} v_p(d)={vd}")
        # AUDIT-3: v3
        vd3 = 0
        dd = d
        while dd % 3 == 0:
            dd //= 3
            vd3 += 1
        if v3n - vd3 > 1:
            viol3 += 1
            print(f"  AUDIT-3 VIOLATION n={n} rep=({x},{y}) v3(n)={v3n} v3(d)={vd3}")

print()
print(f"primitive representations audited: {total_prim_reps} over {len(dp)} numbers")
print(f"AUDIT-1 (cofactor m = 3^(0|1) * split primes, eps=1 iff 3|d): {viol1} violations")
print(f"AUDIT-2 (v_p(d) = v_p(n) for every inert p | n, every primitive rep): {viol2} violations")
print(f"AUDIT-3 (v_3(n) - v_3(d) <= 1): {viol3} violations")
print(f"AUDIT-4 direct count of doubly-primitive n containing an inert prime: "
      f"{n_with_inert_direct} (gcd-derived count in bridge script: 1952)")
hist = sorted(inert_prime_hist.items())
print(f"inert primes occurring (p: #n): {hist[:15]}{' ...' if len(hist) > 15 else ''}")
print(f"ALL AUDITS {'PASS' if viol1 == viol2 == viol3 == 0 and n_with_inert_direct == 1952 else 'CHECK'}"
      f"   [{time.time()-T0:.1f}s total]")
