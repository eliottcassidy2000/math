# taxicab_inert_meet_kpc1_verify.py
# ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1) for THM-461 claims A4/A5 (and the
# counting half of A1's part (c) numbers).
# FRESH method: numpy packed-key sort (key = n<<14 | x) over all 1<=y<=x with x^3+y^3<=10^12,
# instead of the worker's dict-based meet-in-the-middle. Exact: all values < 2^54 in int64.
# Then exact pure-python post-processing: primitivity, Cor C inert confinement, v3 bound,
# Pollard-rho factorizations (exact integers).
import math
import random
import numpy as np
from math import gcd

random.seed(20260610)
LIM = 10**12
XMAX = 10**4

def icbrt(m):
    r = round(m ** (1.0 / 3.0))
    while r * r * r > m:
        r -= 1
    while (r + 1) ** 3 <= m:
        r += 1
    return r

# ---------- build all pairs
chunks = []
for x in range(1, XMAX + 1):
    x3 = x * x * x
    if x3 + 1 > LIM:
        break
    ymax = min(x, icbrt(LIM - x3))
    if ymax < 1:
        continue
    y = np.arange(1, ymax + 1, dtype=np.int64)
    n = y * y * y + x3
    chunks.append((n << np.int64(14)) | np.int64(x))
keys = np.concatenate(chunks)
chunks = None
total_pairs = int(keys.size)
print(f"total pairs 1<=y<=x, x^3+y^3 <= 10^12: {total_pairs}")
keys.sort()
nvals = keys >> np.int64(14)
dup = np.flatnonzero(nvals[1:] == nvals[:-1])

runs = []
if dup.size:
    br = np.flatnonzero(np.diff(dup) > 1)
    starts = np.concatenate(([0], br + 1))
    ends = np.concatenate((br, [dup.size - 1]))
    for s_, e_ in zip(starts, ends):
        runs.append((int(dup[s_]), int(dup[e_]) + 1))
print(f"n with >= 2 positive reps (taxicab-style): {len(runs)}")

# ---------- extract reps; doubly-primitive numbers
doubly = []  # (n, [(x,y),...] primitive reps)
for i0, i1 in runs:
    n = int(nvals[i0])
    prims = []
    for j in range(i0, i1 + 1):
        x = int(keys[j]) & 16383
        y = icbrt(n - x * x * x)
        assert x * x * x + y * y * y == n and 1 <= y <= x
        if gcd(x, y) == 1:
            prims.append((x, y))
    if len(prims) >= 2:
        doubly.append((n, prims))
doubly.sort()
nprims = sum(len(p) for _, p in doubly)
print(f"doubly-primitive numbers (>=2 primitive reps) <= 10^12: {len(doubly)}")
print(f"total primitive reps over those numbers: {nprims}")

# ---------- exact factorization (Miller-Rabin + Pollard rho)
def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def rho(n):
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d

def factorize(n):
    fac = {}
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        while n % p == 0:
            fac[p] = fac.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            fac[m] = fac.get(m, 0) + 1
        else:
            d = rho(m)
            stack += [d, m // d]
    return fac

def facstr(fac):
    return "*".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(fac.items()))

# ---------- audits
corC_viol = 0          # inert p: v_p(d) != v_p(n) for some primitive rep
v3_viol = 0            # v3(d) < v3(n) - 1 for some primitive rep
inert_carrying = 0     # n has a prime == 2 (mod 3)
gcd_gt1 = 0            # gcd of primitive-rep d's > 1
inert_via_gcd = 0      # gcd of d's contains an inert prime
first_inert = None
hist_all = {}          # inert prime p -> #n divisible by p
hist_min = {}          # smallest inert prime of n -> count
for n, prims in doubly:
    ds = [x + y for x, y in prims]
    G = 0
    for d in ds:
        G = gcd(G, d)
    if G > 1:
        gcd_gt1 += 1
    fac = factorize(n)
    assert math.prod(p ** e for p, e in fac.items()) == n
    inerts = sorted(p for p in fac if p % 3 == 2)
    if inerts:
        inert_carrying += 1
        if first_inert is None:
            first_inert = (n, facstr(fac))
        for p in inerts:
            hist_all[p] = hist_all.get(p, 0) + 1
        hist_min[inerts[0]] = hist_min.get(inerts[0], 0) + 1
        # Cor C: every inert prime's FULL valuation must sit in every primitive d
        for p in inerts:
            pe = p ** fac[p]
            for d in ds:
                if d % pe != 0 or (d // pe) % p == 0:  # v_p(d) == v_p(n) exactly (d|n forces <=)
                    if (d % pe) != 0:
                        corC_viol += 1
        if G > 1 and any(p % 3 == 2 for p in factorize(G)):
            inert_via_gcd += 1
    else:
        if G > 1 and any(p % 3 == 2 for p in factorize(G)):
            inert_via_gcd += 1
    v3n = fac.get(3, 0)
    for d in ds:
        v3d = 0
        dd = d
        while dd % 3 == 0:
            dd //= 3
            v3d += 1
        if v3d < v3n - 1:
            v3_viol += 1

print(f"n with gcd(primitive d's) > 1: {gcd_gt1}")
print(f"inert-carrying doubly-primitive n (some p == 2 mod 3 divides n): {inert_carrying}")
print(f"inert-carrying count derived via gcd(d's) factorization instead: {inert_via_gcd}")
print(f"first inert-carrying doubly-primitive n: {first_inert}")
print(f"[A4] Cor C violations (v_p(d) != v_p(n) at an inert p, any primitive rep): {corC_viol}")
print(f"[A4] v3 violations (v3(d) < v3(n)-1): {v3_viol}")
print("histogram (inert p -> #n divisible by p), p <= 50:",
      {p: hist_all[p] for p in sorted(hist_all) if p <= 50})
print("histogram (smallest inert prime -> count), p <= 50:",
      {p: hist_min[p] for p in sorted(hist_min) if p <= 50})

print("\nfirst 25 doubly-primitive numbers:")
for n, prims in doubly[:25]:
    fac = factorize(n)
    split = all(p % 3 == 1 for p in fac) and 3 not in fac
    splitor3 = all(p % 3 == 1 or p == 3 for p in fac)
    print(f"  {n} = {facstr(fac)}   reps {prims}   d's {[x+y for x,y in prims]}"
          f"   completely-split={split} split-or-3={splitor3}")
print("DONE taxicab_inert_meet_kpc1_verify")
