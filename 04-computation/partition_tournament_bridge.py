#!/usr/bin/env python3
"""
partition_tournament_bridge.py — opus-2026-03-14-S80

Exploring the partition function p(n) and Young tableaux
through the (2,3) lens. Key questions:
- p(n) mod 2 and mod 3 — what do the tournament keys say?
- Standard Young tableaux of shape (n) and the hook length formula
- The partition function evaluated at tournament-meaningful n
- Ramanujan congruences p(5n+4)=0 mod 5, etc. and (2,3,5)
- Robinson-Schensted correspondence and tournament structure
"""

from math import comb, factorial, gcd, sqrt
from fractions import Fraction
from functools import lru_cache

def section(title, n):
    print(f"\n{'='*70}")
    print(f"{n}. {title}")
    print(f"{'='*70}\n")

KEY1, KEY2 = 2, 3

# ============================================================
section("PARTITION NUMBERS AND (2,3)", 1)
# ============================================================

@lru_cache(maxsize=1000)
def p(n):
    """Number of partitions of n."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    result = 0
    for k in range(1, n+1):
        # Euler's pentagonal theorem
        for sign in [1, -1]:
            m = k * (3*k + sign) // 2  # generalized pentagonal numbers
            if m > n:
                break
            if k % 2 == 1:
                result += p(n - m)
            else:
                result -= p(n - m)
            if m == n:
                break
    return result

# Actually let me just use the recurrence properly
@lru_cache(maxsize=1000)
def partition(n):
    """Partition function via pentagonal number theorem."""
    if n < 0: return 0
    if n == 0: return 1
    result = 0
    k = 1
    while True:
        # Generalized pentagonal numbers: k(3k-1)/2 and k(3k+1)/2
        pent1 = k * (3*k - 1) // 2
        pent2 = k * (3*k + 1) // 2
        if pent1 > n and pent2 > n:
            break
        sign = (-1)**(k+1)
        if pent1 <= n:
            result += sign * partition(n - pent1)
        if pent2 <= n:
            result += sign * partition(n - pent2)
        k += 1
    return result

print("Partition numbers p(n) and their (2,3) content:")
print(f"{'n':>4s} {'p(n)':>10s} {'mod 2':>6s} {'mod 3':>6s} {'mod 5':>6s} {'mod 7':>6s} {'notes':>30s}")
print("-"*75)

for n in range(31):
    pn = partition(n)
    notes = ""
    if pn == 2: notes = "= KEY1"
    elif pn == 3: notes = "= KEY2"
    elif pn == 5: notes = "= KEY1+KEY2"
    elif pn == 7: notes = "= H_forb_1"
    elif pn == 11: notes = "prime"
    elif pn == 15: notes = "= C(6,2)"
    elif pn == 22: notes = "= 2*11"
    elif pn == 30: notes = "= h(E8)"
    elif pn == 42: notes = "= f(9)"
    elif pn == 56: notes = "= f(10) = dim(V_E7)"
    elif pn == 77: notes = "= 7*11"
    elif pn == 101: notes = "prime"
    elif pn == 135: notes = "= 5*27 = (KEY1+KEY2)*KEY2^3"
    elif pn == 176: notes = "= 16*11"
    elif pn == 231: notes = "= C(22,2)"
    elif pn == 297: notes = "= 27*11"
    elif pn == 385: notes = "= 5*7*11"
    elif pn == 490: notes = "= 2*5*7^2"
    elif pn == 627: notes = "= 3*11*19"
    elif pn == 792: notes = "= C(12,5) = C(h(E6),5)"
    elif pn == 1002: notes = "= 2*3*167"
    elif pn == 1255: notes = "= 5*251"
    elif n == 30: notes = f"p(h(E8)) = {pn}"

    print(f"{n:4d} {pn:10d} {pn%2:6d} {pn%3:6d} {pn%5:6d} {pn%7:6d}  {notes}")

# ============================================================
section("RAMANUJAN CONGRUENCES AND (2,3,5)", 2)
# ============================================================

print("Ramanujan's partition congruences:")
print(f"  p(5n+4) = 0 (mod 5)   for all n >= 0")
print(f"  p(7n+5) = 0 (mod 7)   for all n >= 0")
print(f"  p(11n+6) = 0 (mod 11) for all n >= 0")
print()
print("The moduli: 5, 7, 11")
print(f"  5 = KEY1+KEY2")
print(f"  7 = H_forb_1 = Phi_3(KEY1)")
print(f"  11 = 2*KEY1+KEY1+KEY2 = 2*KEY1+KEY2+KEY1+KEY2... = f(KEY1*KEY2) + KEY1+KEY2? No...")
print(f"  Actually: 5, 7, 11 are the 3rd, 4th, 5th primes")
print()
print("The residues: 4, 5, 6")
print(f"  These are delta(p) = (p^2-1)/24:")
print(f"  (5^2-1)/24 = 24/24 = 1... no, delta = (p-1)/? ")
print()

# Verify: the residue is (p²-1)/24 mod p, which is a standard result
# Actually the residue r in p(pn+r) = 0 mod p is r = 24^{-1}*(p^2-1)/24...
# Standard: r = (p-1)*24^{-1} mod p where 24^{-1} is modular inverse
# For p=5: 24 ≡ 4 (mod 5), 4^{-1} = 4 (mod 5), r = 4*4 = 16 ≡ 1 (mod 5)... no
# Actually r = -1/24 mod p = -(24^{-1}) mod p
# For p=5: 24^{-1} ≡ 4 (mod 5), so r = -4 ≡ 1 (mod 5). But r=4, not 1.
# Hmm, let me just verify directly

print("Verification:")
for n in range(10):
    print(f"  p(5*{n}+4) = p({5*n+4}) = {partition(5*n+4)} = 0 mod 5? {partition(5*n+4) % 5 == 0}")

print()
for n in range(8):
    print(f"  p(7*{n}+5) = p({7*n+5}) = {partition(7*n+5)} = 0 mod 7? {partition(7*n+5) % 7 == 0}")

print()
for n in range(6):
    print(f"  p(11*{n}+6) = p({11*n+6}) = {partition(11*n+6)} = 0 mod 11? {partition(11*n+6) % 11 == 0}")

print()
print("All verified! The (2,3,5) triple appears as the modular structure:")
print("  - The Ramanujan congruences use moduli KEY1+KEY2, Phi_3(KEY1), 11")
print("  - The generating function is 1/prod(1-q^n) = 1/Euler function")
print("  - The Dedekind eta function: eta(tau) = q^{1/24} * prod(1-q^n)")
print("  - The 1/24 relates to |BT| = 24")

# ============================================================
section("PARTITION AT TOURNAMENT-MEANINGFUL VALUES", 3)
# ============================================================

special_n = [
    (2, "KEY1"),
    (3, "KEY2"),
    (5, "KEY1+KEY2"),
    (6, "h(G2)"),
    (7, "H_forb_1"),
    (8, "rank(E8)"),
    (10, "V(Petersen)"),
    (12, "h(E6)"),
    (14, "dim(G2)"),
    (18, "h(E7)"),
    (21, "H_forb_2"),
    (24, "|BT|"),
    (30, "h(E8)"),
    (48, "|BO|"),
    (56, "dim(V_E7)"),
    (78, "dim(E6)"),
    (120, "|BI|"),
]

print("p(n) at tournament-meaningful n values:")
for n, meaning in special_n:
    pn = partition(n)
    # Factor pn
    factors = []
    temp = pn
    for prime in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
        while temp % prime == 0:
            factors.append(prime)
            temp //= prime
    if temp > 1:
        factors.append(temp)
    factor_str = " * ".join(str(f) for f in factors) if len(factors) > 1 else str(pn)

    print(f"  p({n:3d}) = {pn:>15d}  ({meaning:15s})  = {factor_str}")

# ============================================================
section("PENTAGONAL NUMBERS AND THE TOURNAMENT POLYNOMIAL", 4)
# ============================================================

print("The pentagonal number theorem connects to the tournament polynomial:")
print("  Euler: prod(1-q^n) = sum (-1)^k q^{k(3k-1)/2}")
print("  The generalized pentagonal numbers: omega(k) = k(3k-1)/2")
print()

print("Generalized pentagonal numbers (both signs of k):")
for k in range(-6, 7):
    if k == 0: continue
    pent = abs(k) * (3*abs(k) - 1) // 2
    pent_neg = abs(k) * (3*abs(k) + 1) // 2
    if k > 0:
        print(f"  k={k:+d}: omega = {pent:4d}, omega' = {pent_neg:4d}")

print()
print("Pentagonal numbers omega_k = k(3k-1)/2:")
pents = []
for k in range(1, 10):
    omega = k * (3*k - 1) // 2
    pents.append(omega)
    star = ""
    if omega == 1: star = " = 1"
    elif omega == 5: star = " = KEY1+KEY2"
    elif omega == 12: star = " = h(E6)"
    elif omega == 22: star = " = 2*11"
    elif omega == 35: star = " = 5*7 = (KEY1+KEY2)*H_forb_1"
    elif omega == 51: star = " = 3*17"
    elif omega == 70: star = " = C(8,4) = C(rank(E8),4)"
    elif omega == 92: star = " = 4*23"
    print(f"  k={k}: {omega}{star}")

print()
print("The coefficient 3 in k(3k-1)/2 is KEY2!")
print("Pentagonal numbers = k(KEY2*k - 1)/2")
print()
print("The first few: 1, 5, 12, 22, 35, 51, 70, 92, ...")
print(f"  Differences: 4, 7, 10, 13, 16, 19, 22, ...")
print(f"  = arithmetic progression with common difference 3 = KEY2!")

# ============================================================
section("HOOK LENGTH FORMULA AND TOURNAMENTS", 5)
# ============================================================

print("The hook length formula: f^lambda = n! / prod(hook lengths)")
print("counts the number of standard Young tableaux of shape lambda.")
print()

# For rectangular partitions
print("SYT counts for rectangle shapes (a x b):")
for a in range(1, 7):
    for b in range(a, 7):
        n = a * b
        # Hook lengths for a x b rectangle
        hooks = []
        for i in range(a):
            for j in range(b):
                hook = (a - i) + (b - j) - 1
                hooks.append(hook)
        f_lambda = factorial(n)
        for h in hooks:
            f_lambda //= h
        if f_lambda in [1, 2, 3, 5, 6, 7, 8, 12, 14, 18, 21, 24, 30, 42, 48, 56, 120, 132]:
            star = " *"
        else:
            star = ""
        if a <= 4 and b <= 5 and f_lambda > 0:
            print(f"  ({a}x{b}), n={n:2d}: f^lambda = {f_lambda}{star}")

print()
print("SYT counts for staircase shapes (n, n-1, ..., 1):")
for n in range(1, 8):
    shape = list(range(n, 0, -1))
    total = sum(shape)
    # Compute hook lengths
    hooks = []
    for i in range(len(shape)):
        for j in range(shape[i]):
            # arm = cells to the right in row i
            arm = shape[i] - j - 1
            # leg = cells below in column j
            leg = sum(1 for i2 in range(i+1, len(shape)) if shape[i2] > j)
            hook = arm + leg + 1
            hooks.append(hook)
    f_lambda = factorial(total)
    for h in hooks:
        f_lambda //= h
    print(f"  staircase({n}) = {shape}, |lambda|={total}: f^lambda = {f_lambda}")

# ============================================================
section("CATALAN NUMBERS AND THEIR (2,3) CONTENT", 6)
# ============================================================

print("Catalan numbers C_n = C(2n,n)/(n+1):")
for n in range(15):
    cn = comb(2*n, n) // (n+1)
    star = ""
    if cn == 2: star = " = KEY1"
    elif cn == 5: star = " = KEY1+KEY2"
    elif cn == 14: star = " = dim(G2)"
    elif cn == 42: star = " = f(9) = 2*H_forb_2"
    elif cn == 132: star = " = dim(so(12))"
    elif cn == 429: star = " = 3*11*13"
    elif cn == 1430: star = " = 2*5*11*13"
    elif cn == 4862: star = " = 2*11*13*17"
    print(f"  C_{n:2d} = {cn:8d}{star}")

print()
print("Catalan numbers at tournament-meaningful indices:")
print(f"  C_2 = 2 = KEY1")
print(f"  C_3 = 5 = KEY1+KEY2")
print(f"  C_4 = 14 = dim(G2)")
print(f"  C_5 = 42 = f(9) = 2*H_forb_2")
print(f"  C_6 = 132 = dim(so(12))... ")
print()
print(f"  C_3/C_2 = 5/2 = (KEY1+KEY2)/KEY1 = |BI|/|BO|")
print(f"  C_4/C_3 = 14/5")
print(f"  C_5/C_4 = 42/14 = 3 = KEY2")
print(f"  C_6/C_5 = 132/42 = 22/7 ~ pi!")

# ============================================================
section("THE PARTITION FUNCTION MOD 2 AND MOD 3", 7)
# ============================================================

print("p(n) mod 2 — distribution:")
counts = {0: 0, 1: 0}
for n in range(100):
    counts[partition(n) % 2] += 1
print(f"  p(n) even: {counts[0]}/100, p(n) odd: {counts[1]}/100")
print()

# The pattern of p(n) mod 2
print("p(n) mod 2 for n=0..49:")
row = ""
for n in range(50):
    row += str(partition(n) % 2)
    if (n+1) % 10 == 0:
        print(f"  n={n-9:2d}-{n}: {row}")
        row = ""

print()
print("p(n) mod 3 for n=0..49:")
for n in range(50):
    if n % 10 == 0:
        vals = [partition(n+i) % 3 for i in range(min(10, 50-n))]
        print(f"  n={n:2d}-{n+9}: {vals}")

# ============================================================
section("SUMMARY", 8)
# ============================================================

print("="*70)
print("PARTITION — TOURNAMENT DICTIONARY")
print("="*70)
print()
print("Connection                           Value              Meaning")
print("-"*70)
print(f"p(5) = {partition(5)}                             5=KEY1+KEY2        H_forb_1 partitions")
print(f"p(6) = {partition(6)}                            6=h(G2)            H_forb_1 + 4")
print(f"p(7) = {partition(7)}                            7=H_forb_1         next forbidden")
print(f"p(8) = {partition(8)}                            8=rank(E8)         largest exceptional")
print(f"p(10) = {partition(10)}                           10=V(Petersen)     h(E8) partitions")
print(f"p(12) = {partition(12)}                           12=h(E6)           77 = 7*11")
print(f"p(30) = {partition(30)}                        30=h(E8)           5604 = 2^2*3*467")
print()
print("Ramanujan congruences mod (KEY1+KEY2, H_forb_1, 11)")
print("Pentagonal recurrence coefficient = KEY2")
print("Catalan: C_2=KEY1, C_3=KEY1+KEY2, C_4=dim(G2), C_5=f(9)")
print("Euler's pentagonal numbers grow at rate KEY2*k^2/2")
print()
print("p(7) = 15 = C(6,2) = C(h(G2),2)")
print("The number of partitions of H_forb_1 equals the number of edges of K_{h(G2)}!")
