"""
ten_eleven_deep.py -- kind-pasteur-2026-03-14-S63

Deep exploration of I(Omega, 10) and I(Omega, 11) in tournament theory,
and positional number theory implications.

KEY CONNECTIONS:
  10 = 2 * 5    = 1010 in binary (the "2 interleaved with 2")
  11 = prime    = 1011 in binary (the "2 interleaved with 3")

  In the alpha basis: I(Omega, x) = 1 + a1*x + a2*x^2
  I(10) = 1 + 10*a1 + 100*a2
  I(11) = 1 + 11*a1 + 121*a2

  In the beta basis (centered at x=-1): I(y-1) = b0 + b1*y + b2*y^2
  I(10) = b0 + 11*b1 + 121*b2   (y = x+1 = 11)
  I(11) = b0 + 12*b1 + 144*b2   (y = x+1 = 12)

POSITIONAL NUMBER THEORY:
  10 and 11 are "10" and "11" in decimal.
  In base b, "10" = b and "11" = b+1.

  For base 2: "10"=2, "11"=3 — the OCF pair!
  For base 3: "10"=3, "11"=4
  For base 7: "10"=7, "11"=8 — the Mersenne transition!
  For base 10: "10"=10, "11"=11

  The pair (b, b+1) for b=2 is the FUNDAMENTAL pair.

  I(Omega, b) and I(Omega, b+1) determine alpha via Vandermonde:
    I(b) = 1 + b*a1 + b^2*a2
    I(b+1) = 1 + (b+1)*a1 + (b+1)^2*a2

    Subtracting: I(b+1) - I(b) = a1 + (2b+1)*a2
    From I(b): a1 = (I(b)-1-b^2*a2)/b

    So: I(b+1)-I(b) = (I(b)-1-b^2*a2)/b + (2b+1)*a2
                     = (I(b)-1)/b + a2*(-(b^2/b) + 2b+1)
                     = (I(b)-1)/b + a2*(b+1)

    Therefore: a2 = (I(b+1) - I(b) - (I(b)-1)/b) / (b+1)
              = (b*I(b+1) - (b+1)*I(b) + 1) / (b*(b+1))

  This is the VANDERMONDE EXTRACTION for general b!

  At b=2: a2 = (2*I(3) - 3*I(2) + 1) / 6 = (2*I(3) - 3*H + 1) / 6
  (Confirmed by opus HYP-867)

  At b=10: a2 = (10*I(11) - 11*I(10) + 1) / 110

The key insight: ANY pair of consecutive evaluations (I(b), I(b+1))
suffices to extract a1 and a2 (for n <= 7 where higher alpha vanish).
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
from math import comb, factorial

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_ham_cycles_exact(A_sub, k):
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask_size in range(2, k+1):
        for mask in range(1 << k):
            if bin(mask).count('1') != mask_size:
                continue
            if not (mask & 1):
                continue
            for v in range(k):
                if not (mask & (1 << v)):
                    continue
                if v == 0:
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(k):
                    if (prev_mask & (1 << u)) and A_sub[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << k) - 1
    total_cycles = 0
    for v in range(1, k):
        if A_sub[v][0] and dp.get((full, v), 0):
            total_cycles += dp[(full, v)]
    return total_cycles

def get_alpha_1_2(A, n):
    cycles = []
    for size in range(3, n+1, 2):
        for combo in combinations(range(n), size):
            verts = list(combo)
            sub = np.zeros((size, size), dtype=int)
            for a in range(size):
                for b in range(size):
                    sub[a][b] = A[verts[a]][verts[b]]
            if size <= 5:
                c = int(np.trace(np.linalg.matrix_power(sub, size))) // size
            else:
                c = count_ham_cycles_exact(sub, size)
            for _ in range(c):
                cycles.append(frozenset(combo))
    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_1, alpha_2

def rand_bits(total_bits):
    if total_bits <= 30:
        return np.random.randint(0, 1 << total_bits)
    bits = 0
    remaining = total_bits
    shift = 0
    while remaining > 0:
        chunk = min(remaining, 30)
        bits |= int(np.random.randint(0, 1 << chunk)) << shift
        shift += chunk
        remaining -= chunk
    return bits

# ============================================================
# PART 1: I(Omega, x) at many x values — the "spectrum"
# ============================================================
print("=" * 70)
print("PART 1: I(Omega, x) spectrum at n=5 and n=7")
print("=" * 70)

for n in [5, 7]:
    total_bits = n*(n-1)//2
    np.random.seed(42)

    n_samp = min(2**total_bits, 200) if n <= 5 else 100

    # Collect evaluations
    evals = {x: [] for x in [-1, 0, 1, 2, 3, 5, 7, 10, 11]}

    for trial in range(n_samp):
        if n <= 5 and trial < 2**total_bits:
            bits = trial
        else:
            bits = rand_bits(total_bits)

        A = bits_to_adj(bits, n)
        a1, a2 = get_alpha_1_2(A, n)

        for x in evals:
            evals[x].append(1 + a1*x + a2*x**2)

    print(f"\n  n={n}: I(Omega, x) statistics:")
    print(f"  {'x':>4} {'min':>8} {'max':>8} {'mean':>10} {'mod_x':>12} {'always_1_mod_x':>15}")
    for x in sorted(evals.keys()):
        vals = evals[x]
        mod_dist = Counter(v % abs(x) for v in vals) if x != 0 else "N/A"
        always_1 = all(v % abs(x) == 1 for v in vals) if x not in [0, 1, -1] else "N/A"
        print(f"  {x:>4} {min(vals):>8} {max(vals):>8} {np.mean(vals):>10.1f} "
              f"{str(dict(mod_dist) if x not in [0] else 'N/A'):>12} "
              f"{str(always_1):>15}")

# ============================================================
# PART 2: Vandermonde extraction from consecutive pairs
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Vandermonde extraction from (I(b), I(b+1))")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

# For each base b, verify: a2 = (b*I(b+1) - (b+1)*I(b) + 1) / (b*(b+1))
print(f"\n  Vandermonde extraction verification at n={n}:")
for trial in range(5):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)

    print(f"\n  Trial {trial}: a1={a1}, a2={a2}")
    for b in [2, 3, 5, 7, 10]:
        Ib = 1 + a1*b + a2*b**2
        Ib1 = 1 + a1*(b+1) + a2*(b+1)**2

        # Vandermonde extraction
        a2_ext = (b*Ib1 - (b+1)*Ib + 1) / (b*(b+1))
        a1_ext = (Ib - 1 - b**2 * a2_ext) / b

        print(f"    b={b:>2}: I({b})={Ib:>8}, I({b+1})={Ib1:>8}, "
              f"a2_ext={a2_ext:.1f}, a1_ext={a1_ext:.1f}, "
              f"match={'YES' if abs(a2_ext-a2)<0.01 and abs(a1_ext-a1)<0.01 else 'NO'}")

# ============================================================
# PART 3: The positional representation of H
# ============================================================
print("\n" + "=" * 70)
print("PART 3: H in different bases")
print("=" * 70)

n = 7
total_bits = n*(n-1)//2
np.random.seed(42)

def to_base(num, base):
    if num == 0:
        return "0"
    digits = []
    while num > 0:
        digits.append(num % base)
        num //= base
    return ''.join(str(d) for d in reversed(digits))

print(f"\n  H in bases 2, 3, 6, 10:")
print(f"  {'H':>5} {'base2':>12} {'base3':>10} {'base6':>8} {'base10':>6} {'a1':>4} {'a2':>3}")
for trial in range(20):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    a1, a2 = get_alpha_1_2(A, n)

    print(f"  {H:>5} {to_base(H,2):>12} {to_base(H,3):>10} {to_base(H,6):>8} {H:>6} "
          f"{a1:>4} {a2:>3}")

# ============================================================
# PART 4: H in base 2: the bits encode OCF structure
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Base-2 representation of H and OCF structure")
print("=" * 70)

print("""
H = 1 + 2*a1 + 4*a2

In binary: H = ...b3 b2 b1 1
  bit 0 = 1 (always, Redei)
  bit 1 = a1 mod 2 (cycle parity)
  bits 2+ encode carries from a1 and a2

H - 1 = 2*a1 + 4*a2 = 2*(a1 + 2*a2)
(H-1)/2 = a1 + 2*a2

In binary: (H-1)/2 = ...c2 c1 c0
  c0 = a1 mod 2
  c1 = (a1//2 + a2) mod 2 (carries!)
  c2 = ((a1//2 + a2)//2 + ...) mod 2
""")

np.random.seed(42)
print(f"  Base-2 digits of (H-1)/2 vs a1, a2:")
print(f"  {'a1':>4} {'a2':>3} {'(H-1)/2':>8} {'base2':>10} {'a1%2':>4} {'(a1//2+a2)%2':>12}")
for trial in range(15):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2
    half = (H-1)//2

    print(f"  {a1:>4} {a2:>3} {half:>8} {to_base(half,2):>10} {a1%2:>4} {(a1//2+a2)%2:>12}")

# ============================================================
# PART 5: H in base 3: the 3-adic tower
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Base-3 representation and the beta basis")
print("=" * 70)

print("""
In the beta basis: H = b0 + 3*b1 + 9*b2

This is NOT the base-3 representation of H (because b0, b1, b2 can be
negative or > 2). But it's a "signed ternary" representation.

The base-3 representation of H IS meaningful though:
  H mod 3 = last base-3 digit = b0 mod 3
  (H - (H mod 3)) / 3 mod 3 = second-to-last digit

The base-3 digits of H encode the 3-adic structure.
""")

np.random.seed(42)
print(f"  Base-3 digits of H:")
print(f"  {'H':>5} {'base3':>10} {'b0':>5} {'b1':>5} {'b2':>3} {'H%3':>3} {'b0%3':>4}")
for trial in range(15):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2
    b0 = 1 - a1 + a2
    b1 = a1 - 2*a2
    b2 = a2

    print(f"  {H:>5} {to_base(H,3):>10} {b0:>5} {b1:>5} {b2:>3} {H%3:>3} {b0%3:>4}")

# ============================================================
# PART 6: The "10" and "11" of each base
# ============================================================
print("\n" + "=" * 70)
print("PART 6: I(Omega, b) and I(Omega, b+1) for various bases")
print("=" * 70)

print("""
For any base b:
  "10" in base b = b
  "11" in base b = b + 1

  I(Omega, b) = 1 + a1*b + a2*b^2    = "1 a2 a1 1" in mixed-radix
  I(Omega, b+1) = 1 + a1*(b+1) + a2*(b+1)^2

The DIFFERENCE: I(b+1) - I(b) = a1 + (2b+1)*a2
The RATIO: I(b+1)/I(b)

At b=2: I(3)/I(2) = (1+3*a1+9*a2)/(1+2*a1+4*a2)
  For a2=0: I(3)/I(2) = (1+3*a1)/(1+2*a1) -> 3/2 as a1 -> inf

This 3/2 ratio is the k-nacci convergence ratio!
""")

np.random.seed(42)
ratios = defaultdict(list)

for trial in range(200):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)

    for b in [2, 3, 5, 7, 10]:
        Ib = 1 + a1*b + a2*b**2
        Ib1 = 1 + a1*(b+1) + a2*(b+1)**2
        if Ib > 0:
            ratios[b].append(Ib1/Ib)

print(f"  I(b+1)/I(b) ratios at n=7:")
for b in [2, 3, 5, 7, 10]:
    r = ratios[b]
    # Theoretical: for large a1, ratio -> (b+1)/b
    print(f"    b={b:>2}: mean={np.mean(r):.4f}, min={min(r):.4f}, max={max(r):.4f}, "
          f"limit={(b+1)/b:.4f}")

# ============================================================
# PART 7: The (2,3) Vandermonde determinant = 1
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Vandermonde determinants and the specialness of (2,3)")
print("=" * 70)

print("""
The Vandermonde determinant of (b, b+1) for extracting (a1, a2) is:
  V = | 1  b    b^2   |
      | 1  b+1  (b+1)^2 |

  V_22 = det | b    b^2    | = b*(b+1)^2 - b^2*(b+1)= b*(b+1)*(b+1-b) = b*(b+1)
              | b+1  (b+1)^2|

So the Vandermonde factor is b*(b+1).

For b=2: V = 2*3 = 6. The extraction formula divides by 6.
For b=3: V = 3*4 = 12. Divides by 12.
For b=10: V = 10*11 = 110. Divides by 110.

The (2,3) pair has the SMALLEST Vandermonde determinant: 6 = 3!

This means:
  a2 = (2*I(3) - 3*H + 1) / 6

  For this to be an INTEGER, we need 2*I(3) - 3*H + 1 = 0 mod 6.

  Since H is odd (Redei): 3*H is odd. 2*I(3) is even. So 2*I(3)-3*H+1 is even.
  Also I(3) = 1 mod 3 (proved), so 2*I(3) = 2 mod 3, and 3*H = 0 mod 3.
  So 2*I(3) - 3*H + 1 = 2-0+1 = 3 = 0 mod 3.

  Therefore 2*I(3) - 3*H + 1 is divisible by BOTH 2 and 3, hence by 6. CHECK!
""")

# Verify
np.random.seed(42)
for trial in range(10):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2
    I3 = 1 + 3*a1 + 9*a2

    vand = 2*I3 - 3*H + 1
    print(f"  a1={a1:>3}, a2={a2:>2}: H={H:>4}, I(3)={I3:>5}, "
          f"2*I(3)-3*H+1={vand:>5}, /6={vand//6:>3}, "
          f"match_a2={vand//6 == a2}")

# ============================================================
# PART 8: Base-6 representation (the 2*3 base)
# ============================================================
print("\n" + "=" * 70)
print("PART 8: Base-6 representation — the natural base for H")
print("=" * 70)

print("""
Since H mod 6 is determined by (H mod 2, H mod 3) via CRT,
and H mod 2 = 1 always, H mod 6 encodes ONLY the topology (H mod 3).

In base 6:
  H mod 6 = last digit = 1, 3, or 5 (always odd)

  The base-6 representation of H is "the natural representation"
  because 6 = 2*3 captures both fundamental primes simultaneously.

  Base-6 digit 0 = Redei + topology
  Base-6 digit 1 = (H//6) mod 6 = higher structure
  etc.
""")

np.random.seed(42)
# Distribution of base-6 digits
digit_0 = Counter()
digit_1 = Counter()

for trial in range(500):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2

    digit_0[H % 6] += 1
    digit_1[(H // 6) % 6] += 1

print(f"  Base-6 digit 0 distribution: {dict(sorted(digit_0.items()))}")
print(f"  Base-6 digit 1 distribution: {dict(sorted(digit_1.items()))}")

# ============================================================
# PART 9: I(Omega, 10) and I(Omega, 11) — decimal structure
# ============================================================
print("\n" + "=" * 70)
print("PART 9: I(Omega, 10) and I(Omega, 11)")
print("=" * 70)

np.random.seed(42)

print(f"\n  I(10) and I(11) at n=7:")
print(f"  {'a1':>4} {'a2':>3} {'H':>5} {'I(10)':>8} {'I(11)':>8} "
      f"{'I(10)%10':>8} {'I(11)%11':>8} {'last_digit_I10':>14}")

for trial in range(20):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2
    I10 = 1 + 10*a1 + 100*a2
    I11 = 1 + 11*a1 + 121*a2

    print(f"  {a1:>4} {a2:>3} {H:>5} {I10:>8} {I11:>8} "
          f"{I10%10:>8} {I11%11:>8} {I10%10:>14}")

# Check: I(10) mod 10
# I(10) = 1 + 10*a1 + 100*a2 = 1 mod 10 (since 10*a1 = 0 mod 10)
# I(11) = 1 + 11*a1 + 121*a2 = 1 + a1 + a2 mod 11 (since 11=0, 121=0 mod 11)

print(f"\n  I(10) mod 10 = 1 always (trivially, since 10|10*a1 and 10|100*a2)")
print(f"  I(11) mod 11 = (1 + a1 + a2) mod 11 (since 11*a1 = 0 mod 11, 121*a2 = 0 mod 11)")

# I(11) mod 11 = I(Omega, 0) mod 11 = 1 mod 11? NO!
# I(11) mod 11 = 1 + a1 + a2 mod 11 = I(Omega, 1) mod 11
# Because 11 = 0 mod 11, so I(11) = I(0) = 1 mod 11? No:
# I(Omega, x) = 1 + a1*x + a2*x^2
# I(11) mod 11 = I(Omega, 11 mod 11) mod 11 = I(Omega, 0) mod 11 = 1 mod 11
# Wait: polynomial evaluation: I(11) = 1 + 11*a1 + 121*a2
# mod 11: = 1 + 0 + 0 = 1. Always!

# Similarly I(10) mod 10: I(10) = 1 + 10*a1 + 100*a2 = 1 mod 10. Always!

# So BOTH I(10) and I(11) have trivial mod structure.
# The interesting mod structure is at SMALL primes.

print(f"\n  I(11) mod 11 = 1 always (since 11 | 11*a1 and 11 | 121*a2)")
print(f"  In fact: I(b) mod b = 1 for ALL b (since b | b*a1 and b | b^2*a2)")
print(f"  This is just I(Omega, 0) = 1 mod b.")

# ============================================================
# PART 10: The LAST DECIMAL DIGIT of H
# ============================================================
print("\n" + "=" * 70)
print("PART 10: Last decimal digit of H")
print("=" * 70)

# H mod 10 = H mod 2 * H mod 5 via CRT
# H mod 2 = 1 (Redei)
# H mod 5 = (1 + 2*a1 + 4*a2) mod 5

# Since 2 mod 5 = 2, 4 mod 5 = 4 = -1:
# H mod 5 = 1 + 2*a1 - a2 mod 5

np.random.seed(42)
h_mod10 = Counter()
h_mod5 = Counter()

for trial in range(500):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)
    H = 1 + 2*a1 + 4*a2

    h_mod10[H % 10] += 1
    h_mod5[H % 5] += 1

print(f"  H mod 10: {dict(sorted(h_mod10.items()))}")
print(f"  H mod 5:  {dict(sorted(h_mod5.items()))}")
print(f"  Possible last digits: {sorted(h_mod10.keys())}")
print(f"  (H is odd, so last digit in {{1,3,5,7,9}})")

# H mod 5 = I(Omega, 2) mod 5 = I(Omega, 2 mod 5) mod 5 = I(Omega, 2) mod 5
# Since 2 has order 4 mod 5: 2^1=2, 2^2=4, 2^3=3, 2^4=1
# I(Omega, 2) mod 5 = 1 + 2*a1 + 4*a2 mod 5 = 1 + 2*a1 - a2 mod 5

# CRT: H mod 10 from (H mod 2, H mod 5)
# H = 1 mod 2, H = r mod 5
# r=0: H=5 mod 10, r=1: H=1 mod 10, r=2: H=7 mod 10,
# r=3: H=3 mod 10, r=4: H=9 mod 10

print(f"\n  CRT (mod 2=1, mod 5=r) -> mod 10:")
for r in range(5):
    for x in range(10):
        if x % 2 == 1 and x % 5 == r:
            print(f"    r={r}: H = {x} mod 10")

# ============================================================
# PART 11: Binary representations and self-similarity
# ============================================================
print("\n" + "=" * 70)
print("PART 11: Binary self-similarity of 10 and 11")
print("=" * 70)

print("""
Binary representations:
   2 = 10      (the original "10")
   3 = 11      (the original "11")
  10 = 1010    ("10" concatenated with "10")
  11 = 1011    ("10" concatenated with "11")

This is NOT a coincidence in the sense that:
  10 = 2 * 5 = 2 * (4 + 1) = 8 + 2 = 1010_2
  11 = 2 * 5 + 1 = 8 + 2 + 1 = 1011_2

But it IS structurally meaningful:
  10_decimal encodes "the pair (2,2)" in binary
  11_decimal encodes "the pair (2,3)" in binary

And (2,3) is the FUNDAMENTAL pair of tournament theory!

Moreover:
  10 + 11 = 21 = 3 * 7
  10 * 11 = 110 = 2 * 5 * 11
  11 - 10 = 1

  21 = 3 * 7: product of first two odd primes (3-cycles and 7-vertex tournaments)
  110 = 2 * 55 = 2 * 5 * 11: the Vandermonde factor for base 10
""")

# ============================================================
# PART 12: The H-spectrum at n=11 (Paley)
# ============================================================
print("\n" + "=" * 70)
print("PART 12: n=11 — the first prime with I(Omega, 11) structure")
print("=" * 70)

# At n=11, Paley tournament T_11 has H = 95095
# Let's compute I(Omega, x) for T_11 if we can find a1, a2

# H(T_11) = 95095 (known from memory)
# For Paley T_11: c3 = 55 (known)
# c5 = 594 (known)
# But we need the FULL a1 (all odd cycles) and a2 (disjoint pairs)

# From the known cycle counts:
# c3 = 55 (on C(11,3)=165 triples, each has P=1/3 for Paley)
# c5 = 594 (on C(11,5)=462 5-sets)
# c7 = 3960
# c9 = 11055
# c11 = 5505

# a1 = c3 + c5 + c7 + c9 + c11 = 55 + 594 + 3960 + 11055 + 5505 = 21169

a1_T11 = 55 + 594 + 3960 + 11055 + 5505
H_T11 = 95095

# a2 from H = 1 + 2*a1 + 4*a2 + 8*a3 + ...
# At n=11, we have alpha_3 and higher (three disjoint 3-cycles = 9 vertices < 11)
# So H = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + ...
# Can't determine a2 from H alone without knowing higher alpha

# But we can compute: lower bound for a2
# 2*a1 = 42338
# H - 1 = 95094
# H - 1 - 2*a1 = 95094 - 42338 = 52756
# This = 4*a2 + 8*a3 + 16*a4 + ...
# So 4*a2 <= 52756, a2 <= 13189

print(f"  Paley T_11:")
print(f"    H = {H_T11}")
print(f"    a1 = {a1_T11} (total directed odd cycles)")
print(f"    H - 1 - 2*a1 = {H_T11 - 1 - 2*a1_T11} (= 4*a2 + 8*a3 + ...)")
print(f"    a2 <= {(H_T11 - 1 - 2*a1_T11)//4} (upper bound)")

# I(Omega, x) = 1 + a1*x + a2*x^2 + a3*x^3 + ...
# At x=-1: I(-1) = 1 - a1 + a2 - a3 + ...
# At x=2: I(2) = H = 95095
# At x=3: I(3) = 1 + 3*a1 + 9*a2 + 27*a3 + ...

# From Vandermonde (b=2):
# a2_eff = (2*I(3) - 3*H + 1) / 6
# But we don't know I(3) for T_11 without computing a3, a4, etc.

# Let's focus on what we CAN say:
# H mod 3 = I(Omega, -1) mod 3
# H_T11 mod 3 = 95095 mod 3 = ?
print(f"    H mod 3 = {H_T11 % 3}")
print(f"    H mod 6 = {H_T11 % 6}")
print(f"    H mod 10 = {H_T11 % 10}")
print(f"    H mod 11 = {H_T11 % 11}")
print(f"    H mod 12 = {H_T11 % 12}")

# H/|Aut| = 95095/55 = 1729 (the taxicab number!!)
print(f"    H/|Aut| = {H_T11 // 55} (= 1729 = 12^3 + 1^3 = 10^3 + 9^3)")
print(f"    1729 = Hardy-Ramanujan taxicab number!")
print(f"    1729 = 7 * 13 * 19")
print(f"    1729 mod 6 = {1729 % 6}")
print(f"    1729 mod 10 = {1729 % 10}")

# ============================================================
# PART 13: The taxicab connection: 1729 = 10^3 + 9^3 = 12^3 + 1^3
# ============================================================
print("\n" + "=" * 70)
print("PART 13: 1729 and the 10-11 connection")
print("=" * 70)

print("""
H(T_11) / |Aut(T_11)| = 95095 / 55 = 1729

1729 is the SMALLEST number expressible as sum of two cubes in two ways:
  1729 = 12^3 + 1^3 = 10^3 + 9^3

The 10 and 9 appear here! And:
  12 = 2^2 * 3 (contains 2 and 3)
  10 = 2 * 5
   9 = 3^2
   1 = 1

  12^3 = (4*3)^3 = 64 * 27 = 2^6 * 3^3
   1^3 = 1
  10^3 = (2*5)^3 = 8 * 125 = 2^3 * 5^3
   9^3 = 3^6

  1729 = 2^6 * 3^3 + 1 = 2^3 * 5^3 + 3^6

The factorization: 1729 = 7 * 13 * 19
  7 = tournament order for first nontrivial Paley
  13 = next prime after 11
  19 = next Paley prime after 11

  ALL THREE factors are primes p = 1 mod 6 (quadratic residue primes)!
  7 = 1 mod 6, 13 = 1 mod 6, 19 = 1 mod 6

  These are exactly the primes where Paley tournaments T_p are defined
  with QR structure.
""")

# Check: is 7*13*19 all == 1 mod 6?
for p in [7, 13, 19]:
    print(f"  {p} mod 6 = {p % 6}")

# ============================================================
# PART 14: I(Omega, x) mod x and mod (x+1)
# ============================================================
print("\n" + "=" * 70)
print("PART 14: Universal residues of I(Omega, x)")
print("=" * 70)

print("""
For ANY graph G and ANY positive integer x:
  I(G, x) mod x = I(G, 0) mod x = 1 mod x

This is because I(G, x) = 1 + a1*x + a2*x^2 + ..., and all terms
after 1 are divisible by x.

Similarly:
  I(G, x) mod (x-1) = I(G, 1) mod (x-1) = (1 + a1 + a2 + ...) mod (x-1)

  I(G, x) mod (x+1) = I(G, -1) mod (x+1) = (1 - a1 + a2 - ...) mod (x+1)

At x = 2 (H):
  H mod 2 = I(G, 0) mod 2 = 1 (Redei)
  H mod 1 = 0 (trivial)
  H mod 3 = I(G, -1) mod 3 (topology!)

At x = 10:
  I(10) mod 10 = 1 (trivial)
  I(10) mod 9 = I(G, 1) mod 9 (total independent sets mod 9)
  I(10) mod 11 = I(G, -1) mod 11 (Euler char mod 11)

At x = 11:
  I(11) mod 11 = 1 (trivial)
  I(11) mod 10 = I(G, 1) mod 10 (total indep sets mod 10)
  I(11) mod 12 = I(G, -1) mod 12 (Euler char mod 12)

So evaluating at x=10 and taking mod 11 gives I(G,-1) mod 11,
and evaluating at x=11 and taking mod 12 gives I(G,-1) mod 12.

The PAIR (10, 11) gives:
  I(G, -1) mod 11 (from I(10) mod 11)
  I(G, -1) mod 12 (from I(11) mod 12)

  CRT: I(G, -1) mod 132 (since gcd(11,12)=1)!
""")

# Verify at n=7
np.random.seed(42)
print(f"  Verification at n=7:")
ok_10_11 = 0
ok_11_12 = 0
for trial in range(200):
    bits = rand_bits(total_bits)
    A = bits_to_adj(bits, n)
    a1, a2 = get_alpha_1_2(A, n)

    I10 = 1 + 10*a1 + 100*a2
    I11 = 1 + 11*a1 + 121*a2
    I_neg1 = 1 - a1 + a2

    if I10 % 11 == I_neg1 % 11:
        ok_10_11 += 1
    if I11 % 12 == I_neg1 % 12:
        ok_11_12 += 1

print(f"    I(10) mod 11 = I(-1) mod 11: {ok_10_11}/200")
print(f"    I(11) mod 12 = I(-1) mod 12: {ok_11_12}/200")

print("\n\nDone.")
