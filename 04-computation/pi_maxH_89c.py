#!/usr/bin/env python3
"""
pi_maxH_89c.py — Does the Paley tournament achieve maximum H?

opus-2026-03-14-S89c

OEIS A038375: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095
These are the maximum H values for tournaments on n=1..11 vertices.

H(P_3) = 3 = A038375(3) ✓
H(P_7) = 189 = A038375(7) ✓
H(P_11) = 95095 = A038375(11) ✓

Hypothesis: The Paley tournament P_p achieves max H(T) among all
tournaments on p vertices, for p ≡ 3 mod 4.
"""

import math
import itertools

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def paley_tournament(p):
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and is_qr(j - i, p):
                adj[i][j] = 1
    return adj

def cyclic_tournament(n):
    """Regular cyclic tournament: i→j if (j-i) mod n ∈ {1,...,⌊n/2⌋}"""
    adj = [[0]*n for _ in range(n)]
    forward = set(range(1, n//2 + 1))
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in forward:
                adj[i][j] = 1
    return adj

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

# A038375 values from OEIS
a038375 = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661,
           9: 3357, 10: 15745, 11: 95095}

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: Paley vs Max H")
print("=" * 70)

print("\n  n | max H (A038375) | H(Paley) | H(Cyclic) | Paley=Max? | Cyclic=Max?")
print("  " + "-" * 72)

for p in [3, 7, 11, 19]:
    if p % 4 != 3:
        continue
    adj_paley = paley_tournament(p)
    H_paley = count_H(adj_paley, p)
    adj_cyclic = cyclic_tournament(p)
    H_cyclic = count_H(adj_cyclic, p)
    max_h = a038375.get(p, "?")
    paley_is_max = H_paley == max_h if isinstance(max_h, int) else "?"
    cyclic_is_max = H_cyclic == max_h if isinstance(max_h, int) else "?"
    print(f"  {p:2d} | {max_h:>15} | {H_paley:>8} | {H_cyclic:>9} | {str(paley_is_max):>10} | {str(cyclic_is_max):>11}")

# For non-Paley n:
for n in [4, 5, 6, 8, 9, 10]:
    adj_cyclic = cyclic_tournament(n) if n % 2 == 1 else None
    H_cyclic = count_H(adj_cyclic, n) if adj_cyclic else "N/A"
    max_h = a038375.get(n, "?")
    cyclic_is_max = H_cyclic == max_h if isinstance(max_h, int) and isinstance(H_cyclic, int) else "?"
    print(f"  {n:2d} | {max_h:>15} | {'N/A':>8} | {str(H_cyclic):>9} | {'N/A':>10} | {str(cyclic_is_max):>11}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: Verify Exhaustively at n=3..7")
print("=" * 70)

for n in range(3, 8):
    if n <= 6:
        max_found = 0
        max_count = 0
        for adj in all_tournaments(n):
            h = count_H(adj, n)
            if h > max_found:
                max_found = h
                max_count = 1
            elif h == max_found:
                max_count += 1

        print(f"\n  n={n}: max H = {max_found} (achieved by {max_count} tournaments)")
        print(f"    A038375({n}) = {a038375[n]}")
        print(f"    Match: {max_found == a038375[n]}")

    # Check Paley and cyclic
    if n % 4 == 3:
        H_paley = count_H(paley_tournament(n), n)
        print(f"    H(Paley P_{n}) = {H_paley}, is max: {H_paley == a038375[n]}")
    if n % 2 == 1:
        H_cyclic = count_H(cyclic_tournament(n), n)
        print(f"    H(Cyclic C_{n}) = {H_cyclic}, is max: {H_cyclic == a038375[n]}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: Score Sequences of Max-H Tournaments")
print("=" * 70)

for n in range(3, 7):
    max_h = a038375[n]
    print(f"\n  n={n}: max H = {max_h}")
    for adj in all_tournaments(n):
        h = count_H(adj, n)
        if h == max_h:
            scores = sorted(sum(adj[i][j] for j in range(n)) for i in range(n))
            print(f"    scores = {scores}")
            break  # just show one

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: Ratios max H(n+1)/max H(n)")
print("=" * 70)

print("\n  n | max H | ratio(n)/ratio(n-1)")
prev = 1
for n in range(2, 12):
    curr = a038375[n]
    ratio = curr / prev
    print(f"  {n:2d} | {curr:>10} | {ratio:.6f}")
    prev = curr

# Known: for regular tournaments, H ≈ n!/(2^{n-1}) × constant
# The ratio should grow roughly as n/2 if max H ~ n!/2^{n-1} × c
# But max H / Mean H = max H × 2^{n-1} / n!

print("\n  Max H / Mean H:")
for n in range(3, 12):
    max_h = a038375[n]
    mean_h = math.factorial(n) / 2**(n-1)
    print(f"  n={n:2d}: max/mean = {max_h/mean_h:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: H(P_p) / max H — Is Paley Always Optimal?")
print("=" * 70)

# We know H(P_3)=3=max, H(P_7)=189=max, H(P_11)=95095=max
# What about H(P_19)? We computed H(P_19)=1172695746915
# A038375 doesn't give us n=19, but we can compare to cyclic

for p in [3, 7, 11, 19]:
    if p % 4 != 3:
        continue
    H_paley = count_H(paley_tournament(p), p)
    H_cyclic = count_H(cyclic_tournament(p), p)
    mean_h = math.factorial(p) / 2**(p-1)
    print(f"  P_{p:2d}: H(Paley) = {H_paley}, H(Cyclic) = {H_cyclic}")
    print(f"         Paley/Cyclic = {H_paley/H_cyclic:.6f}")
    print(f"         Paley/Mean = {H_paley/mean_h:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: The Max H / Mean H Ratio — Approaching e?")
print("=" * 70)

# Let r_n = max_H / mean_H
# r_3 = 3/1.5 = 2
# r_4 = 5/3 = 1.667
# r_5 = 15/7.5 = 2
# r_6 = 45/22.5 = 2
# r_7 = 189/78.75 = 2.4
# r_8 = 661/315 = 2.098
# r_9 = 3357/1417.5 = 2.368
# r_10 = 15745/7087.5 = 2.221
# r_11 = 95095/38981.25 = 2.440

print("\nMax H / Mean H sequence:")
ratios = []
for n in range(3, 12):
    max_h = a038375[n]
    mean_h = math.factorial(n) / 2**(n-1)
    r = max_h / mean_h
    ratios.append(r)
    print(f"  n={n:2d}: r = {r:.6f}")

# Is this converging? To what?
# If max H ~ c·n!/2^{n-1}, then r → c
# The alternation suggests parity-dependent behavior
print(f"\n  Avg of odd n: {sum(ratios[i] for i in range(0,len(ratios),2))/((len(ratios)+1)//2):.4f}")
print(f"  Avg of even n: {sum(ratios[i] for i in range(1,len(ratios),2))/(len(ratios)//2):.4f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 7: The Alon Bound — Maximum H ≤ c·n^{3/2}·n!/2^n")
print("=" * 70)

# Noga Alon (1990) proved: max H(T) ≤ c·n^{3/2}·n!/2^n for some constant c
# The lower bound from Paley: H(P_p) ≥ c'·n!/2^{n-1} (constant × mean)
# So max H / mean H = O(n^{3/2}) (Alon upper bound) and Ω(1) (Paley lower bound)

print("\nAlon bound analysis: max H ≤ c·n^{3/2}·n!/2^n")
print("log(max H / (n!/2^{n-1})) / log(n):")
for n in range(3, 12):
    max_h = a038375[n]
    mean_h = math.factorial(n) / 2**(n-1)
    r = max_h / mean_h
    if r > 0:
        exponent = math.log(r) / math.log(n) if n > 1 else 0
        print(f"  n={n:2d}: max/mean = {r:.6f}, log ratio / log n = {exponent:.4f}")

# If max/mean ~ n^α, what is α?
# From data: ratios hover around 2-2.5, so α ≈ 0
# This means max H ~ Θ(mean H) — the Paley bound is TIGHT!

print("\n  CONCLUSION: max H / mean H appears BOUNDED")
print("  The Paley tournament gives max/mean ≈ 2.4 for large p ≡ 3 mod 4")
print("  This suggests max H ~ c · n!/2^{n-1} with c ≈ e (Euler's number?)")

# Check: is c approaching e = 2.718...?
print(f"\n  e = {math.e:.6f}")
print(f"  Max odd-n ratio (n=11): {a038375[11] / (math.factorial(11)/2**10):.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 8: Paley H via Formula")
print("=" * 70)

# For Paley P_p: H = IP(G(P_p), 2) = 1 + 2·Σ_cycles + 4·Σ_pairs + ...
# The number of 3-cycles in P_p = p(p²-1)/24
# For regular tournaments generally: t₃ = C(p,3) - p·C((p-1)/2, 2) = p(p²-1)/24

# Is there a closed form for H(P_p)?
# Values: H(P_3)=3, H(P_7)=189, H(P_11)=95095, H(P_19)=1172695746915

# Check if H(P_p) = (2p-1)!! × something
for p, H in [(3, 3), (7, 189), (11, 95095), (19, 1172695746915)]:
    # (2p-1)!!
    double_fact = 1
    for k in range(1, 2*p, 2):
        double_fact *= k
    print(f"  H(P_{p}) = {H}, (2p-1)!! = {double_fact}, ratio = {H/double_fact:.6f}")

# Hmm: H(P_3)/5!! = 3/15 = 1/5
# H(P_7)/13!! = 189/135135 ≈ 0.001399
# Doesn't seem clean.

# Check: H(P_p) / (p!/2^{(p-1)/2})
for p, H in [(3, 3), (7, 189), (11, 95095), (19, 1172695746915)]:
    val = math.factorial(p) / 2**((p-1)//2)
    print(f"  H(P_{p}) / (p!/2^{{(p-1)/2}}) = {H/val:.6f}")

# Check: C(p-1, (p-1)/2)
for p, H in [(3, 3), (7, 189), (11, 95095), (19, 1172695746915)]:
    c = math.comb(p-1, (p-1)//2)
    print(f"  H(P_{p}) / C(p-1,(p-1)/2) = {H/c:.6f}")
    # Also H/C(p-1,(p-1)/2) * something?

# Check: H = C(2k, k) where p = 2k+1?
for p, H in [(3, 3), (7, 189), (11, 95095)]:
    k = (p-1)//2
    ck = math.comb(2*k, k)
    print(f"  p={p}, k={k}: C(2k,k) = {ck}, H/C(2k,k) = {H/ck:.6f}")

# The central binomial coefficient C(2k,k) ~ 4^k/√(πk) by Stirling
# So C(2k,k) involves π!

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 9: H(P_p) and the Central Binomial Coefficient")
print("=" * 70)

# k = (p-1)/2
# C(2k,k) for k=1,3,5,9: 2, 20, 252, 48620
# H(P_p):                 3, 189, 95095, 1172695746915
# H/C(2k,k):             1.5, 9.45, 377.36, ...
# Hmm, 1.5 = 3/2, 9.45 = 189/20
# 189/20 = 189/20. Is 189 = 20 × 9 + 9 = 189? 20*9=180, so 189/20 = 9.45. Not clean.

# Let me try: H(P_p) / p
for p, H in [(3, 3), (7, 189), (11, 95095), (19, 1172695746915)]:
    print(f"  H(P_{p})/p = {H//p}, H mod p = {H%p}")

# H(P_3)/3 = 1, H(P_7)/7 = 27, H(P_11)/11 = 8645, H(P_19)/19 = ...
# 1, 27, 8645, ...
# 27 = 3³. Is 8645 = something^something? 8645 = 5 × 1729 = 5 × 7 × 247 = 5 × 7 × 13 × 19

# Try H(P_p) as a product
for p, H in [(3, 3), (7, 189), (11, 95095)]:
    x = H
    factors = []
    for d in range(2, x+1):
        while x % d == 0:
            factors.append(d)
            x //= d
        if x == 1:
            break
    print(f"  H(P_{p}) = {H} = {' × '.join(str(f) for f in factors)}")

# H(P_3) = 3
# H(P_7) = 189 = 3³ × 7
# H(P_11) = 95095 = ?

print("\n  FACTORIZATIONS:")
print(f"  H(P_3) = 3")
print(f"  H(P_7) = 189 = 3³ × 7 = 27 × 7")
print(f"  H(P_11) = 95095 = 5 × 7 × 11 × 13 × 19")

# Wait: 5 × 19 = 95, 7 × 13 = 91, 11... let me check
print(f"  Check: 5 × 7 × 11 × 13 × 19 = {5*7*11*13*19}")
# = 95095. YES!
# So H(P_11) = 5 × 7 × 11 × 13 × 19 = product of 5 consecutive primes!?
# Actually: 5, 7, 11, 13, 19... these are NOT consecutive (missing 17)
# But: 5 = C(5,2)+1? No, 5 is (p-1)/2 + ... hmm

# 5 × 7 × 11 × 13 × 19 — are these related to p=11?
# 5 = (11-1)/2, 7 = prime, 11 = p, 13 = p+2, 19 = 2p-3
# Not obvious.

# Actually let me double-check with a different factorization
# 95095 = 95095
# /5 = 19019
# /7 = 2717
# /11 = 247  wait let me redo
# 95095/5 = 19019
# 19019/7 = 2717
# 2717/11 = 247
# 247/13 = 19
# So 95095 = 5 × 7 × 11 × 13 × 19 ✓

# These are: 5, 7, 11, 13, 19 — the 5 primes > 3 and ≤ 19
# But 17 is missing!
# Actually: 5 = C(5,1), 7, 11, 13, 19 — no clear pattern from p=11

# Let me check H(P_7) = 189 = 3³ × 7
# = 27 × 7
# 27 = 3³, and 7 = p

# H(P_7) = 7 × 27 = 7 × 3³
# H(P_11) = 5 × 7 × 11 × 13 × 19
# H(P_3) = 3

# Try: H(P_p) = C(p, (p-1)/2) × something?
for p, H in [(3, 3), (7, 189), (11, 95095)]:
    k = (p-1)//2
    c = math.comb(p, k)
    print(f"\n  p={p}: C(p,k) = C({p},{k}) = {c}")
    print(f"  H/C(p,k) = {Fraction(H, c)}")

from fractions import Fraction

# p=3: C(3,1)=3, H/3 = 1
# p=7: C(7,3)=35, H/35 = 189/35 = 27/5... hmm
# p=11: C(11,5)=462, H/462 = 95095/462 = ...
print(f"  95095/462 = {Fraction(95095, 462)}")

# 95095/462 = 95095/462. GCD(95095, 462)?
import math as m
g = m.gcd(95095, 462)
print(f"  GCD = {g}, simplified = {95095//g}/{462//g}")
# = 95095/462 = ... let me compute
# 462 = 2 × 3 × 7 × 11
# 95095 = 5 × 7 × 11 × 13 × 19
# GCD = 7 × 11 = 77
# 95095/77 = 1235, 462/77 = 6
# So H(P_11)/C(11,5) = 1235/6
print(f"  H(P_11)/C(11,5) = {Fraction(95095, 462)}")

print("\n" + "=" * 70)
print("END — opus-2026-03-14-S89c max H")
print("=" * 70)
