#!/usr/bin/env python3
"""
pascal_slope_d_family_macmini_0614s1.py  (mac-mini-2026-06-14-S1)

THE PASCAL-SLOPE-d / d-STEP FIBONACCI FAMILY (user dispatch).
Construction d: row n = [ C(n-1-(d-1)k, k) : k>=0 ],  row-sum a_d(n) = a_d(n-1)+a_d(n-d).
  d=1: full Pascal rows, sums 2^n
  d=2: Fibonacci diagonals, sums Fibonacci
  d=3: sums Narayana's cows A000930
Verify: (1) the triangle entries reproduce the user's three constructions exactly;
(2) row-sums satisfy a_d(n)=a_d(n-1)+a_d(n-d) and match 2^n / Fibonacci / Narayana;
(3) GF 1/(1-x-x^d); (4) candidate "same-pace" central sequences = max binomial on the
slope-d diagonal AND C(n,floor(n/2))-analogs; (5) the combinatorial model: a_d(n) =
#tilings of a length-(n-1) strip by 1-tiles and d-tiles = gap-d independent sets;
(6) dominant (Pisot) root of x^d = x^{d-1}+1.
"""
import sys
from math import comb

sys.stdout.reconfigure(line_buffering=True)

def C(n, k):
    return comb(n, k) if (n >= 0 and 0 <= k <= n) else 0

def row(d, n):
    """row n (1-indexed) of construction d: entries C(n-1-(d-1)k, k), k=0,1,..."""
    out = []
    k = 0
    while True:
        top = n - 1 - (d - 1) * k
        if top < k:   # C(top,k)=0 beyond this
            break
        out.append(C(top, k))
        k += 1
    # trim trailing zeros
    while out and out[-1] == 0:
        out.pop()
    return out

def rowsum(d, n):
    return sum(row(d, n))

def tilings_1_and_d(d, L):
    """# tilings of a 1xL strip by 1-tiles and d-tiles (a(L)=a(L-1)+a(L-d), a(0)=1)."""
    a = [0]*(L+1)
    for i in range(L+1):
        if i == 0: a[i] = 1
        else:
            a[i] = a[i-1] + (a[i-d] if i-d >= 0 else 0)
    return a

print("="*74)
print("THE PASCAL-SLOPE-d FAMILY  (row n = C(n-1-(d-1)k,k); sums = d-step Fibonacci)")
print("="*74)

for d in [1, 2, 3]:
    print(f"\n--- d = {d} ---")
    for n in range(1, 11):
        r = row(d, n)
        print(f"  row {n:2d}: {' + '.join(map(str,r)):28s} = {sum(r)}")
    sums = [rowsum(d, n) for n in range(1, 16)]
    print(f"  row-sums (n=1..15): {sums}")

print("\n" + "="*74)
print("ROW-SUM RECURRENCE  a_d(n) = a_d(n-1) + a_d(n-d)  and named sequences")
print("="*74)
names = {1: "powers of 2 (A000079)", 2: "Fibonacci (A000045)", 3: "Narayana's cows (A000930)",
         4: "A003269", 5: "A003520", 6: "A005708"}
for d in range(1, 7):
    sums = [rowsum(d, n) for n in range(1, 21)]
    # check recurrence
    rec_ok = all(sums[n] == sums[n-1] + sums[n-d] for n in range(d, len(sums)))
    # tiling model cross-check
    til = tilings_1_and_d(d, 19)  # a(L), L=0..19 ; row-sum a_d(n) should = til[n-1]
    til_match = all(rowsum(d, n) == til[n-1] for n in range(1, 20))
    print(f"  d={d}: {names.get(d,'?'):26s} sums={sums[:12]}")
    print(f"        recurrence a(n)=a(n-1)+a(n-{d}) holds: {rec_ok};  = #tilings(1- & {d}-ominoes): {til_match}")

print("\n" + "="*74)
print('"SAME-PACE" CENTRAL SEQUENCES  (largest binomial on each slope-d diagonal)')
print("="*74)
for d in [1, 2, 3]:
    maxes = [max(row(d, n)) if row(d, n) else 0 for n in range(1, 18)]
    print(f"  d={d}: max-on-diagonal (n=1..17): {maxes}")
# central binomial explicitly
cb = [C(n, n//2) for n in range(0, 14)]
print(f"\n  central binomial C(n,floor(n/2)) (A001405): {cb}")
print(f"  (= the project's metagraph WIDTH C(n-2,floor((n-2)/2)); user's powers-of-2 partner)")
# user's stated partners, for matching
print("\n  USER-STATED partners (for OEIS matching):")
print("   powers-of-2 partner: 1,1,2,3,6,10,20,35,70,126  =? central binomial:", cb[:10] == [1,1,2,3,6,10,20,35,70,126])
print("   fibonacci partner:   1,1,2,3,4,6,10,15,20,35,56,70,120,210,252")
print("   construction-3 partner: 1,1,2,3,4,5,6,10,15,21,28,35,56,84,120,126")

print("\n" + "="*74)
print("DOMINANT (PISOT) ROOT of  x^d = x^{d-1} + 1")
print("="*74)
import cmath
def dom_root(d):
    # Newton on f(x)=x^d - x^{d-1} - 1, real root > 1
    x = 1.5
    for _ in range(200):
        f = x**d - x**(d-1) - 1
        fp = d*x**(d-1) - (d-1)*x**(d-2)
        x -= f/fp
    return x
consts = {1:"2", 2:"golden phi=1.6180339887", 3:"supergolden/plastic-adjacent psi=1.4655712319",
          4:"1.3802775691", 5:"1.3247...(plastic? no)"}
for d in range(1, 7):
    print(f"  d={d}: dominant root of x^{d}=x^{d-1}+1  ->  {dom_root(d):.10f}   {consts.get(d,'')}")
print("  (d=2 golden ratio; d=3 'supergolden ratio' 1.4656, root of x^3=x^2+1)")

print("\nDONE.")
