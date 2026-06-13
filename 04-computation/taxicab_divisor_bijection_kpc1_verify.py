# taxicab_divisor_bijection_kpc1_verify.py
# ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1) for THM-461 claims A1/A2.
# FRESH method, deliberately different from the worker's |x|,|y|<=1000 box scan:
#   reps are enumerated by the (d = x+y, x) parametrization n = d*(3x^2-3dx+d^2),
#   good divisors are counted by a per-d arithmetic-progression sweep.
# A small third-method box brute force cross-checks my own enumeration on n<=20000.
import math

N = 10**6
DMAX = 1
while (DMAX + 1) ** 3 <= 4 * N:
    DMAX += 1  # d^3 <= 4n <= 4N is necessary (n = d*q >= d^3/4)
print(f"N = {N}, DMAX (largest d with d^3 <= 4N) = {DMAX}")

# ---------- Side 1: all unordered integer reps {x,y}, x >= y, x^3+y^3 = n in [1,N]
# q = x^2-xy+y^2 > 0 for (x,y) != (0,0), so d = x+y = n/q > 0; x >= ceil(d/2);
# q = 3x^2-3dx+d^2 strictly increasing in x for x > d/2.
repcnt = bytearray(N + 1)
totreps = 0
maxcoord = 0
for d in range(1, DMAX + 2):
    x = (d + 1) // 2
    while True:
        n = d * (3 * x * x - 3 * d * x + d * d)
        if n > N:
            break
        if n >= 1:
            repcnt[n] += 1
            totreps += 1
            mc = max(abs(x), abs(d - x))
            if mc > maxcoord:
                maxcoord = mc
        x += 1
print(f"[reps]  total unordered integer reps over n in [1,{N}]: {totreps}")
print(f"[reps]  max |coordinate| over all reps: {maxcoord}  (claimed bound/observed: 577)")

# ---------- Side 2: good divisors per the THM-461 criterion
# d>0, d|n, 3d|(d^3-n), Delta=(4n-d^3)/(3d) a perfect square (with parity check counted separately)
gdcnt = bytearray(N + 1)
totgood = 0
parity_viol = 0
for d in range(1, DMAX + 2):
    d3 = d * d * d
    td = 3 * d
    n0 = -(-d3 // 4)            # ceil(d^3/4): below this Delta < 0
    n0 = -(-n0 // d) * d        # first multiple of d >= n0
    if n0 < d:
        n0 = d
    for n in range(n0, N + 1, d):
        if (d3 - n) % td:
            continue
        delta = (4 * n - d3) // td
        if delta < 0:
            continue
        e = math.isqrt(delta)
        if e * e == delta:
            gdcnt[n] += 1
            totgood += 1
            if (e - d) % 2:
                parity_viol += 1
print(f"[divs]  total good divisors over n in [1,{N}]: {totgood}")
print(f"[A2]    parity violations (sqrt(Delta) !== d mod 2 among square-Delta divisors): {parity_viol}")

# ---------- Compare
if repcnt == gdcnt:
    print(f"[A1]    BIJECTION HOLDS: repcnt == gdcnt at every n in [1,{N}]  (0 mismatches)")
else:
    bad = [n for n in range(1, N + 1) if repcnt[n] != gdcnt[n]]
    print(f"[A1]    MISMATCHES at {len(bad)} values of n; first 20: {bad[:20]}")

# ---------- Third-method sanity: raw box brute force for n <= 20000
NB = 20000
B = 100  # max coord for n<=NB is d/2+sqrt(NB/(3d)) <= 82.2 at d=1
brute = bytearray(NB + 1)
for x in range(-B, B + 1):
    for y in range(-B, x + 1):  # y <= x: unordered
        n = x * x * x + y * y * y
        if 1 <= n <= NB:
            brute[n] += 1
ok = all(brute[n] == repcnt[n] for n in range(1, NB + 1))
print(f"[sanity] box brute force |x|,|y|<={B} agrees with (d,x)-enumeration on n<=,{NB}: {ok}")

# edge cases explicitly
for n in (1, 2, 1729):
    reps = []
    for d in range(1, DMAX + 2):
        if d * d * d > 4 * n or n % d:
            continue
        if (d * d * d - n) % (3 * d):
            continue
        delta = (4 * n - d * d * d) // (3 * d)
        e = math.isqrt(delta)
        if e * e == delta:
            x, y = (d + e) // 2, (d - e) // 2
            reps.append((x, y, d))
    print(f"[edge]  n={n}: reps (x,y,d) = {reps}")
print("DONE taxicab_divisor_bijection_kpc1_verify")
