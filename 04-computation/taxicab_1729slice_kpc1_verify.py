# taxicab_1729slice_kpc1_verify.py — ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1)
# Independent check of worker claim A6 (the 1729 slice):
#   - full divisor scan of 1729: good divisors exactly {13,19}; mod-3 gate
#     (3d | d^3-n) open at ALL 8 divisors; cofactors 133=7*19, 91=7*13;
#     both reps primitive; gcd(13,19)=1
#   - r_E(1729) = 48 by direct norm-form count; B(1729) = 8 = tau(1729);
#     #units(L_1729) = 12 + r_E = 60 (THM-434 formula)
#   - B-records over all t: (1,1)(7,2)(49,3)(91,4)(637,6)(1729,8); no B=7
#     below 7^6; degenerate rungs (4t-1 = 3*square): 1, 7, 91 degenerate,
#     {3,13,49,133,637,1729} non-degenerate; records over non-degenerate t
#     = THM-434's published list 3, 13, 49, 133, 637, 1729
#   - Ta(2) = 1729 (brute force below 1729); 9^3+10^3-12^3 = 1 near-miss
#
# FRESH METHOD: explicit table from the criterion; B via character sieve
# (numpy) cross-checked against brute-force norm-form counts.
import math
import numpy as np

n = 1729
print(f"n = {n} = ", end="")
m, fs = n, []
for p in range(2, n + 1):
    while m % p == 0:
        fs.append(p)
        m //= p
    if m == 1:
        break
print(" * ".join(map(str, fs)), f"; all primes == 1 mod 3: {all(p % 3 == 1 for p in fs)}")

divs = [d for d in range(1, n + 1) if n % d == 0]
print(f"divisors ({len(divs)} = tau): {divs}")
print(f"{'d':>6} {'3d|(d^3-n)':>11} {'s=(d^3-n)/3d':>13} {'Delta':>10} {'square?':>8} roots")
good = []
for d in divs:
    num = d ** 3 - n
    gate = (num % (3 * d) == 0)
    if not gate:
        print(f"{d:>6} {'NO':>11}")
        continue
    s = num // (3 * d)
    Delta = d * d - 4 * s
    sq = Delta >= 0 and math.isqrt(Delta) ** 2 == Delta
    roots = ""
    if sq:
        e = math.isqrt(Delta)
        assert (e - d) % 2 == 0   # A2 automatic parity
        x, y = (d + e) // 2, (d - e) // 2
        assert x ** 3 + y ** 3 == n
        roots = f"{{{x},{y}}} gcd={math.gcd(x, y)} cofactor={n // d}"
        good.append(d)
    print(f"{d:>6} {'yes':>11} {s:>13} {Delta:>10} {('YES' if sq else 'no'):>8} {roots}")
print(f"good divisors: {good} (worker claims exactly [13, 19])")
print(f"cofactors: {n//13} = 7*19 -> {n//13 == 7*19}; {n//19} = 7*13 -> {n//19 == 7*13}")
print(f"gcd(13,19) = {math.gcd(13, 19)}")

# r_E(1729) by direct norm-form count (both sign conventions)
R = math.isqrt(4 * n // 3) + 2
c_plus = sum(1 for a in range(-R, R + 1) for b in range(-R, R + 1)
             if a * a + a * b + b * b == n)
c_minus = sum(1 for a in range(-R, R + 1) for b in range(-R, R + 1)
              if a * a - a * b + b * b == n)
print(f"r_E(1729): x^2+xy+y^2 count = {c_plus}, x^2-xy+y^2 count = {c_minus} "
      f"(worker: 48); 12 + r_E = {12 + c_plus} (worker/THM-434: 60 units)")

# B sieve up to 7^6 (character divisor sum), exact ints
T = 7 ** 6
B = np.zeros(T + 1, dtype=np.int64)
for d in range(1, T + 1):
    r = d % 3
    if r == 1:
        B[d::d] += 1
    elif r == 2:
        B[d::d] -= 1
print(f"B(1729) = {int(B[1729])} (worker: 8); tau(1729) = {len(divs)}")
# cross-check sieve against brute force for t <= 200
Rs = math.isqrt(4 * 200 // 3) + 2
bf = [0] * 201
for a in range(-Rs, Rs + 1):
    for b in range(-Rs, Rs + 1):
        v = a * a + a * b + b * b
        if 0 < v <= 200:
            bf[v] += 1
bad = [t for t in range(1, 201) if bf[t] != 6 * int(B[t])]
print(f"sieve vs brute force r=6B, t<=200: {len(bad)} mismatches")

# records of B over ALL t >= 1
recs = []
best = -1
for t in range(1, 1730):
    if int(B[t]) > best:
        best = int(B[t])
        recs.append((t, best))
print(f"B-records over all t in [1,1729]: {recs} "
      f"(worker: (1,1)(7,2)(49,3)(91,4)(637,6)(1729,8))")

# no B = 7 strictly below 7^6; B(7^6) = 7
where7 = np.flatnonzero(B[:T] == 7)
print(f"t < 7^6 with B(t) = 7: {len(where7)}; B(7^6) = {int(B[T])}")

# degeneracy: 4t-1 = 3 * square
def degenerate(t):
    v = 4 * t - 1
    if v % 3:
        return False
    s = math.isqrt(v // 3)
    return s * s == v // 3
print("degenerate (4t-1 = 3*sq): "
      + ", ".join(f"{t}:{degenerate(t)}" for t in [1, 7, 91, 3, 13, 49, 133, 637, 1729]))

# records of B over NON-degenerate t >= 2 (THM-434's ladder domain)
recs_nd = []
best = 0
for t in range(2, 1730):
    if not degenerate(t) and int(B[t]) > best:
        best = int(B[t])
        recs_nd.append((t, best))
print(f"B-records over non-degenerate t in [2,1729]: {recs_nd}")
print(f"matches THM-434 published rungs 3,13,49,133,637,1729: "
      f"{[t for t, _ in recs_nd] == [3, 13, 49, 133, 637, 1729]}")

# Ta(2): brute force smallest n with two positive reps
from collections import Counter as Ctr
c = Ctr()
for x in range(1, 14):
    for y in range(1, x + 1):
        c[x ** 3 + y ** 3] += 1
two = sorted(k for k, v in c.items() if v >= 2)
print(f"smallest n with two positive reps (brute force x,y<=13): {two[:3]}")
print(f"no n < 1729 with two positive reps: {two[0] == 1729}")

# Ramanujan near-miss
print(f"9^3 + 10^3 - 12^3 = {9**3 + 10**3 - 12**3} (Klein: 12^3 = {12**3})")
print(f"1729 = 9^3+10^3 = {9**3 + 10**3 == 1729}; = 12^3+1 = {12**3 + 1 == 1729}")

ok = (good == [13, 19] and c_plus == 48 and c_minus == 48 and int(B[1729]) == 8
      and len(bad) == 0 and recs == [(1, 1), (7, 2), (49, 3), (91, 4), (637, 6), (1729, 8)]
      and len(where7) == 0 and int(B[T]) == 7
      and degenerate(1) and degenerate(7) and degenerate(91)
      and not any(degenerate(t) for t in [3, 13, 49, 133, 637, 1729])
      and [t for t, _ in recs_nd] == [3, 13, 49, 133, 637, 1729]
      and two[0] == 1729)
print("VERDICT(A6 computational side):", "PASS" if ok else "FAIL")
