# taxicab_bijection_kpc1_verify.py — ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1)
# Independent check of worker claims A1 + A2 (THM-461):
#   A1: unordered integer reps {x,y} of n = x^3+y^3 biject with good divisors d:
#       d>0, d|n, 3d | (d^3-n), Delta = (4n-d^3)/(3d) a perfect square;
#       {x,y} = roots of z^2 - dz + s, s = (d^3-n)/(3d); only d <= (4n)^(1/3) can pass.
#   A2: parity condition sqrt(Delta) == d (mod 2) is AUTOMATIC.
#
# FRESH METHOD (different from worker):
#   * rep side: enumerate unordered pairs y <= x over a DIFFERENT box [-1100,1100]
#     (worker used |x|,|y| <= 1000); bucket by n.
#   * divisor side: DIVISOR-MAJOR sieve (d outer, multiples n inner) instead of
#     n-major divisor scan; reconstruct roots and re-assert x^3+y^3 = n.
#   * box completeness re-derived independently: for a rep of n >= 1 with larger
#     coordinate X (necessarily X >= 1) and d = x+y (necessarily d >= 1 since
#     n = d*q with q = x^2-xy+y^2 > 0), we have n = f(d) = d(3X^2 - 3dX + d^2)
#     with f'(d) = 3(X-d)^2 >= 0, so n >= f(1) = 3X^2 - 3X + 1.
#     X = 578 forces n >= 1002253 - 1734 = 1000519 > 10^6. So every coordinate
#     of every rep of n <= 10^6 satisfies |coord| <= 577; box 1100 is complete.
#     (smaller coordinate: |y| <= |x| by choice of X = max(|x|,|y|).)
import math

N = 10**6
B = 1100

# ---------- rep side ----------
repcnt = bytearray(N + 1)   # counts are tiny (max ~ a few); bytearray is fine
maxc_rep = 0
tot_rep = 0
for x in range(-B, B + 1):
    x3 = x * x * x
    for y in range(-B, x + 1):          # unordered: y <= x
        n = x3 + y * y * y
        if 1 <= n <= N:
            repcnt[n] += 1
            tot_rep += 1
            mc = max(abs(x), abs(y))
            if mc > maxc_rep:
                maxc_rep = mc

print(f"[rep side] box [-{B},{B}], unordered pairs y<=x")
print(f"[rep side] total reps of n in [1,{N}]: {tot_rep}")
print(f"[rep side] max |coordinate| seen: {maxc_rep} (proved bound: 577)")
assert maxc_rep <= 577, "box-completeness derivation violated!"

# ---------- divisor side (divisor-major sieve) ----------
divcnt = bytearray(N + 1)
tot_div = 0
parity_viol = 0
recon_fail = 0
maxc_div = 0
Dmax = round((4 * N) ** (1 / 3)) + 2
while Dmax**3 > 4 * N:
    Dmax -= 1
print(f"[div side] only d <= (4n)^(1/3) can pass (Delta >= 0); Dmax = {Dmax}")
for d in range(1, Dmax + 1):
    d3 = d * d * d
    td = 3 * d
    # start n at first multiple of d with 4n >= d^3
    nstart = ((d3 + 4 * d - 1) // (4 * d)) * d
    if nstart < d:
        nstart = d
    for n in range(nstart, N + 1, d):
        num = 4 * n - d3          # = 3d * Delta ; note 3d | (d^3-n) <=> 3d | (4n-d^3) given d|n
        if num < 0:
            continue
        if num % td:
            continue
        Delta = num // td
        e = math.isqrt(Delta)
        if e * e != Delta:
            continue
        # good divisor found
        if (e - d) % 2:                      # A2: should be impossible
            parity_viol += 1
            continue
        xx = (d + e) // 2
        yy = (d - e) // 2
        if xx * xx * xx + yy * yy * yy != n:  # re-assert the reconstruction
            recon_fail += 1
            continue
        divcnt[n] += 1
        tot_div += 1
        mc = max(abs(xx), abs(yy))
        if mc > maxc_div:
            maxc_div = mc

print(f"[div side] total good divisors over n in [1,{N}]: {tot_div}")
print(f"[div side] parity violations (A2 says 0): {parity_viol}")
print(f"[div side] reconstruction failures x^3+y^3 != n: {recon_fail}")
print(f"[div side] max |coordinate| from reconstructed roots: {maxc_div}")

# ---------- compare ----------
mismatch = [n for n in range(1, N + 1) if repcnt[n] != divcnt[n]]
print(f"[compare] n with repcnt != divcnt: {len(mismatch)}"
      + (f" first few: {mismatch[:10]}" if mismatch else ""))

# ---------- explicit edge cases ----------
print("[edges] n=1 (cube, {1,0}):", repcnt[1], divcnt[1])
print("[edges] n=2 (x=y=1, Delta=0, d^3=4n boundary):", repcnt[2], divcnt[2])
print("[edges] n=7 (negative member {2,-1}):", repcnt[7], divcnt[7])
print("[edges] n=8 (cube, {2,0}):", repcnt[8], divcnt[8])
print("[edges] n=1729:", repcnt[1729], divcnt[1729])

ok = (len(mismatch) == 0 and parity_viol == 0 and recon_fail == 0
      and tot_rep == tot_div)
print("VERDICT(A1,A2):", "PASS" if ok else "FAIL")
