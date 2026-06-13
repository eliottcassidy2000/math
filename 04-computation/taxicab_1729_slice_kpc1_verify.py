# taxicab_1729_slice_kpc1_verify.py
# ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1) for THM-461 claim A6 (the 1729 slice)
# plus the B-record / THM-434 record-rung cross-checks.
import math

# ---------- full divisor scan of 1729
n = 1729
divs = sorted(d for d in range(1, n + 1) if n % d == 0)
print(f"divisors of {n}: {divs} (tau = {len(divs)})")
good = []
for d in divs:
    d3 = d ** 3
    gate = (d3 - n) % (3 * d) == 0
    row = f"  d={d:5d} gate(3d | d^3-n)={'open' if gate else 'CLOSED'}"
    if gate:
        s = (d3 - n) // (3 * d)
        delta = d * d - 4 * s
        row += f" s={s} Delta={delta}"
        if delta >= 0 and math.isqrt(delta) ** 2 == delta:
            e = math.isqrt(delta)
            x, y = (d + e) // 2, (d - e) // 2
            row += f" SQUARE e={e} -> (x,y)=({x},{y}), gcd={math.gcd(x,y)}, cofactor n/d={n//d}"
            good.append(d)
        else:
            row += " not a square"
    print(row)
print(f"good divisors of 1729: {good} (claimed: [13, 19]); gcd = {math.gcd(*good) if len(good)==2 else '?'}")
print(f"cofactors: 1729/13 = {1729//13} (= 7*19: {7*19}), 1729/19 = {1729//19} (= 7*13: {7*13})")

# ---------- r_E(1729) by direct norm-form count, BOTH forms
cnt_plus = sum(1 for a in range(-100, 101) for b in range(-100, 101)
               if a * a + a * b + b * b == 1729)
cnt_minus = sum(1 for a in range(-100, 101) for b in range(-100, 101)
                if a * a - a * b + b * b == 1729)
print(f"r_E(1729): x^2+xy+y^2 count = {cnt_plus}, x^2-xy+y^2 count = {cnt_minus} (claimed 48)")
print(f"THM-434 units(L_1729) = 12 + r_E = {12 + cnt_plus} (claimed 60)")

# ---------- B(t) sieve and records up to 7^6
T = 7 ** 6
B = [0] * (T + 1)
for d in range(1, T + 1):
    c = 1 if d % 3 == 1 else (-1 if d % 3 == 2 else 0)
    if c:
        for m in range(d, T + 1, d):
            B[m] += c
print(f"B(1729) = {B[1729]} (claimed 8 = tau(1729))")

# independent cross-check of the divisor-character formula: direct norm-form count for t<=300
for t in range(1, 301):
    r = sum(1 for a in range(-25, 26) for b in range(-25, 26) if a * a + a * b + b * b == t)
    assert r == 6 * B[t], (t, r, B[t])
print("cross-check r_E(t) = 6*B(t) by direct count: OK for all t <= 300")

def degenerate(t):  # 4t-1 = 3*square  <=>  omega_t in Q(sqrt-3) (THM-434 exclusion)
    m = 4 * t - 1
    if m % 3:
        return False
    s = math.isqrt(m // 3)
    return s * s == m // 3

rec_all, rec_nd = [], []
best_a = best_n = 0
for t in range(1, 1730):
    if B[t] > best_a:
        best_a = B[t]
        rec_all.append((t, B[t]))
    if not degenerate(t) and B[t] > best_n:
        best_n = B[t]
        rec_nd.append((t, B[t]))
print(f"B-records over ALL t <= 1729: {rec_all}")
print(f"B-records over NON-DEGENERATE t <= 1729 (4t-1 != 3*sq): {rec_nd}")
print(f"degenerate among record t's: {[t for t,_ in rec_all if degenerate(t)]}")
first7 = next(t for t in range(1, T + 1) if B[t] == 7)
print(f"first t with B(t) = 7: {first7} (7^6 = {7**6})")

# ---------- Ta(2) brute force and Ramanujan near-miss
pos = {}
for x in range(1, 13):
    for y in range(1, x + 1):
        m = x ** 3 + y ** 3
        if m <= 1729:
            pos.setdefault(m, []).append((x, y))
two = sorted(m for m, v in pos.items() if len(v) >= 2)
print(f"n <= 1729 with >= 2 positive reps: {two}, reps of 1729: {pos.get(1729)}")
print(f"9^3 + 10^3 - 12^3 = {9**3 + 10**3 - 12**3}  (12^3 = {12**3})")
print(f"1729 degenerate as a rung? {degenerate(1729)} (4t-1 = {4*1729-1} = 3*{(4*1729-1)//3}; "
      f"isqrt = {math.isqrt((4*1729-1)//3)})")
print("DONE taxicab_1729_slice_kpc1_verify")
