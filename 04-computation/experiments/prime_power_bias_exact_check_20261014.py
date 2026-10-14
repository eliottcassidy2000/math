"""Exact finite checks for the mod-4 prime-power correction; standard library only."""
from fractions import Fraction
from math import isqrt, lcm, log, sqrt
import json

N = 100_000
sieve = bytearray(b"\x01") * (N + 1)
sieve[:2] = b"\x00\x00"
for p in range(2, isqrt(N) + 1):
    if sieve[p]:
        sieve[p * p : N + 1 : p] = b"\x00" * ((N - p * p) // p + 1)
primes = [p for p in range(2, N + 1) if sieve[p]]
K = N.bit_length()
scale = lcm(*range(1, K + 1))
pi1 = [0] * (N + 1)
pi3 = [0] * (N + 1)
odd = [0] * (N + 1)
j_jump = [0] * (N + 1)
higher_jump = [0] * (N + 1)
for p in primes:
    if p == 2:
        continue
    pi1[p] = int(p % 4 == 1)
    pi3[p] = int(p % 4 == 3)
    odd[p] = 1
    pk, k = p, 1
    while pk <= N:
        sign = 1 if pk % 4 == 3 else -1
        j_jump[pk] += sign * (scale // k)
        if k >= 3:
            higher_jump[pk] += sign * (scale // k)
        pk *= p
        k += 1

first = None
j = higher = checks = 0
snapshots = {}
for x in range(1, N + 1):
    for a in (pi1, pi3, odd):
        a[x] += a[x - 1]
    j += j_jump[x]
    higher += higher_jump[x]
    delta = pi3[x] - pi1[x]
    rhs_scaled = j + (scale // 2) * odd[isqrt(x)] - higher
    if delta * scale != rhs_scaled:
        raise RuntimeError((x, delta, rhs_scaled, scale))
    checks += 1
    if pi1[x] > pi3[x] and first is None:
        first = {"x": x, "pi1": pi1[x], "pi3": pi3[x]}
    if x in (26833, 26849, 26861, 26863, 26879):
        snapshots[str(x)] = {
            "pi1": pi1[x], "pi3": pi3[x],
            "delta_J": str(Fraction(j, scale)),
            "square_correction": str(Fraction(odd[isqrt(x)], 2)),
            "higher_correction": str(Fraction(-higher, scale)),
        }
if first != {"x": 26861, "pi1": 1473, "pi3": 1472}:
    raise RuntimeError(first)
print(json.dumps({
    "scope": f"Exact identity at every integer 1 <= x <= {N}",
    "checks": checks,
    "first_strict_pi1_lead": first,
    "snapshots": snapshots,
    "note": "Finite corroboration; asymptotics and zero-distribution claims need separate proofs.",
}, indent=2))
