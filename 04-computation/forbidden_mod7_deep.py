#!/usr/bin/env python3
"""
Deep analysis of the mod-7 forbidden pattern at n=6.
opus-2026-03-14-S84

At n=6: missing H = {7, 21, 35, 39}. The values 7, 21, 35 are ALL odd multiples
of 7 up to max_H=45. Is H ≡ 0 (mod 7) always forbidden at n=6?
What about H mod 7 distribution? What about n=5, n=7?

Also: 39 = 3*13. What makes 13 special? 13 is the 7th prime.
And 7 = 2^3 - 1 (Mersenne), 13 = 2^4 - 3 (near-Mersenne).

Plan:
1. H mod p distribution at n=5, 6 (exhaustive)
2. Which primes p have "H never ≡ 0 mod p" property?
3. Forbidden H structure across n=3..7
4. Arithmetic progressions in achievable H
5. Connection to quadratic residues mod p
"""

from itertools import permutations
from collections import Counter, defaultdict
from fractions import Fraction
import sys

def compute_all_H(n):
    """Exhaustively compute H for all 2^C(n,2) tournaments on n vertices."""
    m = n * (n - 1) // 2
    N = 1 << m
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    all_perms = list(permutations(range(n)))

    H_values = []
    for bits in range(N):
        if bits % 10000 == 0 and N > 10000:
            print(f"  n={n}: {bits}/{N} ({100*bits/N:.1f}%)", file=sys.stderr)

        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(arcs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1

        H = 0
        for p in all_perms:
            valid = True
            for i in range(n-1):
                if adj[p[i]][p[i+1]] != 1:
                    valid = False
                    break
            if valid:
                H += 1
        H_values.append(H)

    return H_values

# ============================================================
# Part 1: H mod p distribution
# ============================================================
print("=" * 70)
print("PART 1: H mod p DISTRIBUTION")
print("=" * 70)

# Cache H values
H_cache = {}
for n in [3, 4, 5, 6]:
    print(f"\nComputing n={n}...")
    H_cache[n] = compute_all_H(n)
    dist = Counter(H_cache[n])
    achievable = sorted(dist.keys())
    print(f"  Achievable H: {achievable}")
    print(f"  Max H = {max(achievable)}, #distinct = {len(achievable)}")

# Check H mod p for small primes
print("\n" + "=" * 70)
print("PART 2: WHICH PRIMES FORBID H ≡ 0?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    H_vals = H_cache[n]
    achievable = sorted(set(H_vals))
    print(f"\nn={n}: achievable = {achievable}")

    for p in [2, 3, 5, 7, 11, 13]:
        residues = sorted(set(h % p for h in achievable))
        missing_residues = sorted(set(range(p)) - set(residues))

        # Count tournaments per residue class
        res_counts = Counter()
        for h in H_vals:
            res_counts[h % p] += 1

        zero_present = 0 in residues
        info = ""
        if not zero_present:
            info = " *** H ≢ 0 (mod p) ***"
        if missing_residues:
            info += f" missing residues: {missing_residues}"

        print(f"  mod {p:2d}: residues = {residues}{info}")

# ============================================================
# Part 3: Arithmetic progression analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 3: ARITHMETIC PROGRESSIONS IN ACHIEVABLE H")
print("=" * 70)

for n in [3, 4, 5, 6]:
    achievable = sorted(set(H_cache[n]))
    print(f"\nn={n}: {achievable}")

    # Check which APs are fully contained
    max_H = max(achievable)
    ach_set = set(achievable)

    # Find longest AP
    best_ap = (0, 0, 0)  # (length, start, step)
    for d in range(1, max_H + 1):
        for a in achievable:
            length = 0
            while a + length * d <= max_H and (a + length * d) in ach_set:
                length += 1
            if length > best_ap[0]:
                best_ap = (length, a, d)

    print(f"  Longest AP: start={best_ap[1]}, step={best_ap[2]}, length={best_ap[0]}")
    print(f"    AP: {[best_ap[1] + i*best_ap[2] for i in range(best_ap[0])]}")

    # Check gaps structure
    gaps = []
    all_odd = list(range(1, max_H + 1, 2))
    for h in all_odd:
        if h not in ach_set:
            gaps.append(h)

    if gaps:
        print(f"  Missing odd values: {gaps}")
        # Factor each
        for g in gaps:
            factors = []
            temp = g
            for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
                while temp % p == 0:
                    factors.append(p)
                    temp //= p
            if temp > 1:
                factors.append(temp)
            print(f"    {g} = {'×'.join(map(str, factors))}, mod 6 = {g%6}, mod 8 = {g%8}, mod 12 = {g%12}")

# ============================================================
# Part 4: The 39 = 3*13 mystery
# ============================================================
print("\n" + "=" * 70)
print("PART 4: WHY IS 39 FORBIDDEN AT n=6?")
print("=" * 70)

# 39 = 3*13. The achievable set includes:
# 37 and 41, but not 39. Gap of 4 in odd integers.
# Is 39 special because 13 = H_forb_1 + h(G2) = 7 + 6?
# Or because 39 = 3*(2*7-1) = 3*13?

achievable_6 = sorted(set(H_cache[6]))
ach_set_6 = set(achievable_6)

print(f"Achievable at n=6: {achievable_6}")
print(f"Missing: {sorted(set(range(1, 46, 2)) - ach_set_6)}")

# How many tournaments achieve each neighboring value?
dist_6 = Counter(H_cache[6])
for h in [33, 35, 37, 39, 41, 43, 45]:
    print(f"  H={h}: {dist_6.get(h, 0)} tournaments")

# Check: is 39 near any forbidden value at other n?
print(f"\n39/3 = 13 (prime)")
print(f"39/7 = {39/7:.4f} (not integer)")
print(f"39 = 40 - 1 = 5*8 - 1")
print(f"39 = 42 - 3 = 6*7 - 3")
print(f"39 ≡ {39 % 7} mod 7, {39 % 6} mod 6, {39 % 5} mod 5, {39 % 8} mod 8")

# ============================================================
# Part 5: Forbidden values at n=5 and n=7 (sample)
# ============================================================
print("\n" + "=" * 70)
print("PART 5: n=5 DETAILED + n=7 SAMPLED")
print("=" * 70)

achievable_5 = sorted(set(H_cache[5]))
print(f"\nn=5 achievable: {achievable_5}")
missing_5 = sorted(set(range(1, max(achievable_5)+1, 2)) - set(achievable_5))
print(f"n=5 missing odd: {missing_5}")

# n=7: sample 100k tournaments
import random
random.seed(42)
n = 7
m = n * (n - 1) // 2  # 21
arcs7 = [(i, j) for i in range(n) for j in range(i+1, n)]
perms7 = list(permutations(range(n)))
print(f"\nn=7: sampling 100000 tournaments (m={m})...")

H_sample_7 = Counter()
for trial in range(100000):
    if trial % 20000 == 0:
        print(f"  trial {trial}/100000", file=sys.stderr)

    bits = random.randint(0, (1 << m) - 1)
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs7):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    H = 0
    for p in perms7:
        valid = True
        for i in range(n-1):
            if adj[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            H += 1
    H_sample_7[H] += 1

achievable_7_sample = sorted(H_sample_7.keys())
max_H_7 = max(achievable_7_sample)
print(f"\nn=7 sampled achievable: {achievable_7_sample}")
print(f"n=7 max H seen = {max_H_7}")
print(f"n=7 #distinct seen = {len(achievable_7_sample)}")

# Check mod 7
missing_odd_7 = sorted(set(range(1, max_H_7+1, 2)) - set(achievable_7_sample))
print(f"n=7 missing odd (in sample): {missing_odd_7}")
print(f"n=7 multiples of 7 found: {[h for h in achievable_7_sample if h % 7 == 0]}")

# Check if H ≡ 0 mod 7 ever at n=7
mod7_zero = [h for h in achievable_7_sample if h % 7 == 0]
print(f"n=7 H ≡ 0 mod 7: {mod7_zero}")
for h in mod7_zero[:10]:
    print(f"  H={h}: {H_sample_7[h]} times in sample")

# ============================================================
# Part 6: Quadratic residue connection
# ============================================================
print("\n" + "=" * 70)
print("PART 6: QUADRATIC RESIDUES AND ACHIEVABLE H")
print("=" * 70)

# QR mod p for small p
for p in [5, 7, 11, 13]:
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    qnr = set(range(1, p)) - qr
    print(f"\nmod {p}: QR = {sorted(qr)}, QNR = {sorted(qnr)}")

    for n in [5, 6]:
        ach = set(H_cache[n])
        ach_mod = set(h % p for h in ach)
        qr_only = ach_mod & qr
        qnr_only = ach_mod & qnr
        missing = set(range(p)) - ach_mod
        print(f"  n={n}: achievable residues mod {p} = {sorted(ach_mod)}")
        if missing:
            in_qr = missing & qr
            in_qnr = missing & qnr
            print(f"    missing in QR: {sorted(in_qr)}, missing in QNR: {sorted(in_qnr)}")

# ============================================================
# Part 7: Score sequence → H relationship at n=6
# ============================================================
print("\n" + "=" * 70)
print("PART 7: SCORE SEQUENCES AND FORBIDDEN H")
print("=" * 70)

# Group n=6 tournaments by score sequence
n = 6
m_6 = 15
arcs6 = [(i, j) for i in range(n) for j in range(i+1, n)]

# For each score sequence, what H values are achievable?
score_to_H = defaultdict(set)
score_to_count = Counter()

for idx, H in enumerate(H_cache[6]):
    # Reconstruct score sequence
    bits = idx
    scores = [0] * n
    for k, (i, j) in enumerate(arcs6):
        if (bits >> k) & 1:
            scores[i] += 1
        else:
            scores[j] += 1
    score_seq = tuple(sorted(scores))
    score_to_H[score_seq].add(H)
    score_to_count[score_seq] += 1

print(f"\nn=6: {len(score_to_H)} distinct score sequences")
for seq in sorted(score_to_H.keys()):
    h_vals = sorted(score_to_H[seq])
    h_range = f"{min(h_vals)}-{max(h_vals)}"
    # Check if 39 is achievable for this score seq
    has_39 = 39 in score_to_H[seq]
    near_39 = {37, 41} & score_to_H[seq]

    info = ""
    if max(h_vals) >= 35:
        if not has_39 and near_39:
            info = " [39 GAP!]"
        elif has_39:
            info = " [has 39]"

    print(f"  {seq}: {score_to_count[seq]:5d} tours, H ∈ {{{', '.join(map(str, h_vals))}}}{info}")

# ============================================================
# Part 8: Number-theoretic structure of gaps
# ============================================================
print("\n" + "=" * 70)
print("PART 8: GAP STRUCTURE — NUMBER THEORY")
print("=" * 70)

# At n=6, gaps = {7, 21, 35, 39}
# 7 = 7, 21 = 3*7, 35 = 5*7, 39 = 3*13
# Sum of gaps = 102 = 2*3*17
# Product of gaps = 7*21*35*39 = 7^3 * 3^2 * 5 * 13

gaps_6 = [7, 21, 35, 39]
print(f"Gaps at n=6: {gaps_6}")
print(f"Sum = {sum(gaps_6)}")
print(f"Product = {7*21*35*39}")
print(f"  = 7^3 * 3^2 * 5 * 13 = {7**3 * 9 * 5 * 13}")

# Mean of achievable H
ach_6 = sorted(set(H_cache[6]))
print(f"\nMean of achievable set = {sum(ach_6)/len(ach_6):.2f}")
print(f"Median of achievable set = {ach_6[len(ach_6)//2]}")

# Check: gap positions within odd integers 1..45
all_odd_6 = list(range(1, 46, 2))  # 23 odd values
gap_positions = [all_odd_6.index(g) for g in gaps_6]
print(f"\nGap positions (0-indexed in odd integers 1..45): {gap_positions}")
print(f"  7 is the 3rd odd integer")
print(f"  21 is the 10th odd integer")
print(f"  35 is the 17th odd integer")
print(f"  39 is the 19th odd integer")
print(f"Position differences: {[gap_positions[i+1]-gap_positions[i] for i in range(len(gap_positions)-1)]}")

# Are positions 3, 10, 17 an AP? (step 7!)
print(f"\n3, 10, 17 form AP with step 7! ({10-3}, {17-10})")
print(f"Position 19 breaks the pattern (would need 24 for AP)")

# So the 7-multiples gaps are at positions 3, 10, 17 = AP(3, 7)
# The non-7 gap (39) is at position 19

# ============================================================
# Part 9: Generating function and divisibility
# ============================================================
print("\n" + "=" * 70)
print("PART 9: GENERATING FUNCTION ANALYSIS")
print("=" * 70)

dist_6 = Counter(H_cache[6])
# Q_6(q) = sum_{H achievable} count(H) * q^H
# Check: is Q_6 divisible by (1 - q^7) or related cyclotomic?

# Build polynomial coefficients
max_h = max(dist_6.keys())
coeffs = [0] * (max_h + 1)
for h, c in dist_6.items():
    coeffs[h] = c

print(f"Q_6(q) = sum c_h * q^h, max h = {max_h}")
print(f"Q_6(1) = {sum(coeffs)} = 2^15 = {2**15}")

# Evaluate at roots of unity
import cmath
for k in [2, 3, 5, 7]:
    omega = cmath.exp(2j * cmath.pi / k)
    val = sum(c * omega**h for h, c in enumerate(coeffs))
    print(f"  Q_6(ω_{k}) = {val.real:.6f} + {val.imag:.6f}i, |Q_6| = {abs(val):.6f}")

# Check if Q_6 vanishes at 7th root of unity
omega7 = cmath.exp(2j * cmath.pi / 7)
for j in range(1, 7):
    val = sum(c * (omega7**j)**h for h, c in enumerate(coeffs))
    print(f"  Q_6(ω_7^{j}) = {val.real:.4f} + {val.imag:.4f}i, |Q_6| = {abs(val):.4f}")

# ============================================================
# Part 10: Comparison with theoretical max and Gamma function
# ============================================================
print("\n" + "=" * 70)
print("PART 10: GAMMA FUNCTION AND FACTORIAL CONNECTIONS")
print("=" * 70)

import math

# max_H at n = n!/2^(n-1) * something? No, max H values:
for n in [3, 4, 5, 6]:
    ach = sorted(set(H_cache[n]))
    max_h = max(ach)
    mean_h = Fraction(sum(H_cache[n]), len(H_cache[n]))
    n_fact = math.factorial(n)
    print(f"n={n}: max_H={max_h}, n!={n_fact}, max_H/n! = {max_h/n_fact:.6f}, mean_H = {float(mean_h):.6f}")
    print(f"  max_H * 2^(n-1) = {max_h * 2**(n-1)}, n! = {n_fact}, ratio = {max_h * 2**(n-1) / n_fact:.6f}")

    # Check: is mean_H = n!/2^(n-1)?
    expected_mean = Fraction(n_fact, 2**(n-1))
    print(f"  mean_H = {mean_h} = {float(mean_h):.6f}, n!/2^(n-1) = {expected_mean} = {float(expected_mean):.6f}, match = {mean_h == expected_mean}")

# ============================================================
# Part 11: Is H ≡ 0 mod 7 forbidden for ALL n?
# ============================================================
print("\n" + "=" * 70)
print("PART 11: IS H ≡ 0 MOD 7 FORBIDDEN FOR ALL n?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    ach = set(H_cache[n])
    mod7_zero = [h for h in sorted(ach) if h % 7 == 0]
    print(f"n={n}: H ≡ 0 mod 7 values: {mod7_zero}")

# n=3: achievable = {1, 3}, no multiples of 7
# n=4: achievable = {1, 3, 5}, no multiples of 7
# n=5: achievable = {1, 3, 5, 9, 11, 13, 15},
# 15 is not 0 mod 7. But is any achievable value 0 mod 7?

for n in [3, 4, 5, 6]:
    ach = set(H_cache[n])
    has_7k = any(h % 7 == 0 for h in ach)
    print(f"n={n}: any H ≡ 0 mod 7? {has_7k}")

print(f"\nn=7 (sampled): any H ≡ 0 mod 7? {any(h % 7 == 0 for h in H_sample_7.keys())}")
if any(h % 7 == 0 for h in H_sample_7.keys()):
    mod7_zero_7 = sorted([h for h in H_sample_7.keys() if h % 7 == 0])
    print(f"  Values: {mod7_zero_7}")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)
print("""
KEY FINDINGS:
1. H mod 7 behavior across n — is H ≡ 0 mod 7 truly forbidden?
2. The 7-multiple gaps at n=6 sit at positions 3, 10, 17 in odd integers — AP with step 7!
3. Score sequences constrain which H values are achievable
4. Q_6(ω_7^j) values reveal cyclotomic divisibility structure
5. The 39 = 3*13 outlier breaks the pure mod-7 pattern
""")
