#!/usr/bin/env python3
"""
es_elliptic_s90ct.py — Erdos-Straus conjecture via elliptic/rational curve perspective.

4/n = 1/a + 1/b + 1/c  with 1 <= a <= b <= c, integers.

Clearing denominators: 4abc = n(bc + ac + ab).

For each prime n up to 1000:
  1. Find ALL ordered decompositions (a,b,c) with a<=b<=c
  2. Count solutions per prime
  3. Classify by n mod 42 and correlate
  4. Fix a, get rational curve in (b,c): 4bc = n(bc + ac + ab)
     => bc(4 - n) = n*a*(b + c)  => parametric family
  5. Find hardest primes (fewest solutions)
  6. Compare systematic vs greedy for p = 1 mod 4
  7. Print summary statistics

opus-2026-03-16-S90ct
"""

import sys
from math import gcd, ceil, isqrt
from collections import defaultdict

# ---------------------------------------------------------------------------
# 1. Find all ES decompositions for a given n
# ---------------------------------------------------------------------------
def es_decompositions(n):
    """
    Find all (a, b, c) with a <= b <= c and 4/n = 1/a + 1/b + 1/c.

    Since 1/a >= (4/n)/3, we have a >= ceil(3n/4).
    Also 1/a <= 4/n, so a <= n (since 1/a >= ... but we need 1/a > 0 leftover).
    Actually: 1/a = 4/n - 1/b - 1/c < 4/n, so a > n/4.
    And 1/a >= 1/b >= 1/c, so 4/n = 1/a+1/b+1/c <= 3/a => a <= 3n/4.

    For b: 1/b = (4/n - 1/a) - 1/c <= (4/n - 1/a), and 1/b >= 1/c,
    so 2/b >= 4/n - 1/a => b <= 2/(4/n - 1/a) = 2an/(4a - n).
    Also b >= a.
    """
    solutions = []
    # a ranges: a > n/4 (strict) and a <= 3n/4
    a_min = n // 4 + 1  # strictly > n/4
    a_max = 3 * n // 4 + 1  # a <= 3n/4, +1 for ceiling issues

    for a in range(a_min, a_max + 1):
        r = 4 * a - n  # 4/n - 1/a = (4a - n)/(na), need r > 0
        if r <= 0:
            continue
        # 1/b + 1/c = r / (n*a)
        # b >= a, and 1/b >= 1/c => 2/b >= r/(na) => b <= 2na/r
        b_max = 2 * n * a // r + 1
        b_min = a
        for b in range(b_min, b_max + 1):
            # 1/c = r/(na) - 1/b = (rb - na) / (nab)
            num = r * b - n * a
            if num <= 0:
                continue
            den = n * a * b
            if den % num != 0:
                continue
            c = den // num
            if c >= b:
                solutions.append((a, b, c))
    return solutions


# ---------------------------------------------------------------------------
# 2. Greedy algorithm (Egyptian fraction style)
# ---------------------------------------------------------------------------
def greedy_es(n):
    """Greedy: 4/n = 1/a + remainder, repeat. Returns (a, b, c) or None."""
    # 1/a where a = ceil(n/4)
    a = ceil(n / 4)
    # remainder = 4/n - 1/a = (4a - n) / (na)
    num = 4 * a - n
    den = n * a
    if num <= 0:
        return None
    g = gcd(num, den)
    num, den = num // g, den // g
    # 1/b where b = ceil(den/num)
    b = ceil(den / num)
    # remainder = num/den - 1/b = (num*b - den) / (den*b)
    num2 = num * b - den
    den2 = den * b
    if num2 <= 0:
        if num2 == 0:
            # exact: 4/n = 1/a + 1/b, set c = infinity — not valid 3-term
            return None
        return None
    if den2 % num2 != 0:
        return None  # greedy doesn't terminate in 3 steps
    c = den2 // num2
    result = tuple(sorted([a, b, c]))
    return result


# ---------------------------------------------------------------------------
# 3. Systematic formula for p = 1 mod 4
# ---------------------------------------------------------------------------
def systematic_1mod4(p):
    """
    For p = 1 mod 4: 4/p = 1/p + 1/((p+1)/4) + ...
    Actually: 4/p = 1/ceil(p/4) + remainder.

    Known identity: 4/p = 1/p + 1/((p+1)/4 * p) ... let's use:
    If p = 4k+1: 4/p = 1/(k+1) + 1/(p(k+1))
    Wait that's only 2 terms for 4/(4k+1).
    Check: 1/(k+1) + 1/(p(k+1)) = (p+1)/(p(k+1)) = (4k+2)/(p(k+1)) = 2/(p)... no.

    Identity: 4/(4k+1) = 1/(k+1) + (4(k+1) - (4k+1)) / ((4k+1)(k+1))
                        = 1/(k+1) + 3/((4k+1)(k+1))
    So we need 3/((4k+1)(k+1)) = 1/b + 1/c.

    Simplest: b = (4k+1)(k+1), c = (4k+1)(k+1)/2... only if divisible.
    Use: 3/m = 1/m + 1/m + 1/m (but need distinct? no, repeats ok if a<=b<=c).
    Actually 3/m = 1/ceil(m/3) + ... greedy on remainder.

    Let's just return the canonical: 4/(4k+1) = 1/(k+1) + 1/((4k+1)(k+1)) + 1/((4k+1)(k+1))
    Check: 1/(k+1) + 2/((4k+1)(k+1)) = (4k+1+2)/((4k+1)(k+1)) = (4k+3)/((4k+1)(k+1))
    That's not 4/(4k+1).

    Better: just find the decomposition with smallest max(c).
    """
    # Find decomposition minimizing c (the largest denominator)
    sols = es_decompositions(p)
    if not sols:
        return None
    return min(sols, key=lambda t: t[2])


# ---------------------------------------------------------------------------
# 4. Rational curve analysis: fix a, curve in (b,c)
# ---------------------------------------------------------------------------
def curve_genus_info(n, a):
    """
    Fixing a: 4/n - 1/a = 1/b + 1/c => (4a-n)/(na) = 1/b + 1/c.
    Let R = 4a - n, D = na. Then 1/b + 1/c = R/D.
    => Rbc = D(b+c). This is a RATIONAL curve (genus 0, conic).
    Rbc - Db - Dc = 0 => (Rb - D)(Rc - D) = D^2.
    So solutions come from divisors of D^2 = (na)^2.
    Number of solutions = number of divisor pairs (d1,d2) of D^2 with
    d1 = Rb - D, d2 = Rc - D, d1*d2 = D^2, d1 > 0, d2 > 0, b <= c.
    """
    R = 4 * a - n
    if R <= 0:
        return 0, []
    D = n * a
    D2 = D * D
    count = 0
    points = []
    # d1 * d2 = D^2, d1 <= d2 (since b <= c)
    for d1 in range(1, isqrt(D2) + 1):
        if D2 % d1 != 0:
            continue
        d2 = D2 // d1
        # b = (d1 + D) / R, c = (d2 + D) / R
        bnum = d1 + D
        cnum = d2 + D
        if bnum % R != 0 or cnum % R != 0:
            continue
        b = bnum // R
        c = cnum // R
        if b >= a and c >= b:
            count += 1
            points.append((a, b, c))
    return count, points


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------
def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return [p for p in range(2, limit + 1) if is_prime[p]]


def main():
    LIMIT = 1000
    primes = sieve_primes(LIMIT)
    print(f"Erdos-Straus Conjecture: 4/n = 1/a + 1/b + 1/c")
    print(f"Analyzing all {len(primes)} primes up to {LIMIT}")
    print("=" * 72)

    # Storage
    sol_count = {}      # prime -> number of solutions
    all_sols = {}       # prime -> list of solutions
    mod42_counts = defaultdict(list)  # residue -> list of (prime, count)
    hardest = []        # (count, prime)

    for p in primes:
        sols = es_decompositions(p)
        sol_count[p] = len(sols)
        all_sols[p] = sols
        mod42_counts[p % 42].append((p, len(sols)))
        hardest.append((len(sols), p))

    # Verify conjecture holds
    zero_sol = [p for p in primes if sol_count[p] == 0]
    print(f"\nPrimes with NO solution: {zero_sol if zero_sol else 'NONE (conjecture holds for all primes <= 1000)'}")

    # ---------------------------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------------------------
    counts = sorted(sol_count.values())
    print(f"\nSolution count statistics:")
    print(f"  Min: {counts[0]}, Max: {counts[-1]}, Median: {counts[len(counts)//2]}")
    print(f"  Mean: {sum(counts)/len(counts):.1f}")

    # ---------------------------------------------------------------------------
    # Hardest primes (fewest solutions)
    # ---------------------------------------------------------------------------
    hardest.sort()
    print(f"\n{'='*72}")
    print(f"TOP 20 HARDEST PRIMES (fewest ES decompositions):")
    print(f"{'Prime':>8} {'#Sol':>6} {'mod4':>6} {'mod42':>6} {'mod12':>6}  Sample decomposition")
    print(f"{'-'*72}")
    for cnt, p in hardest[:20]:
        sample = all_sols[p][0] if all_sols[p] else None
        sample_str = f"4/{p} = 1/{sample[0]} + 1/{sample[1]} + 1/{sample[2]}" if sample else "NONE"
        print(f"{p:>8} {cnt:>6} {p%4:>6} {p%42:>6} {p%12:>6}  {sample_str}")

    # Easiest primes
    print(f"\nTOP 10 EASIEST PRIMES (most decompositions):")
    print(f"{'Prime':>8} {'#Sol':>6} {'mod4':>6} {'mod42':>6}")
    print(f"{'-'*50}")
    for cnt, p in sorted(hardest, reverse=True)[:10]:
        print(f"{p:>8} {cnt:>6} {p%4:>6} {p%42:>6}")

    # ---------------------------------------------------------------------------
    # Mod-42 classification
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"SOLUTION COUNT BY n mod 42:")
    print(f"{'Residue':>8} {'#Primes':>8} {'AvgSol':>10} {'MinSol':>8} {'MaxSol':>8}")
    print(f"{'-'*50}")
    for r in sorted(mod42_counts.keys()):
        entries = mod42_counts[r]
        if not entries:
            continue
        cnts = [c for _, c in entries]
        avg = sum(cnts) / len(cnts)
        print(f"{r:>8} {len(entries):>8} {avg:>10.1f} {min(cnts):>8} {max(cnts):>8}")

    # ---------------------------------------------------------------------------
    # Mod-4 classification
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"SOLUTION COUNT BY n mod 4:")
    mod4_data = defaultdict(list)
    for p in primes:
        if p <= 3:
            continue
        mod4_data[p % 4].append(sol_count[p])
    for r in sorted(mod4_data.keys()):
        cnts = mod4_data[r]
        if cnts:
            print(f"  p = {r} mod 4: {len(cnts)} primes, avg={sum(cnts)/len(cnts):.1f}, "
                  f"min={min(cnts)}, max={max(cnts)}")

    # ---------------------------------------------------------------------------
    # Rational curve analysis: for hardest primes, show divisor structure
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"RATIONAL CURVE ANALYSIS (genus-0 conic for each fixed a):")
    print(f"For the 10 hardest primes, show how solutions distribute across 'a' values.\n")
    for cnt, p in hardest[:10]:
        if cnt == 0:
            continue
        a_min = p // 4 + 1
        a_max = 3 * p // 4 + 1
        print(f"  n = {p} ({cnt} solutions total, a ranges [{a_min}, {a_max}]):")
        for a in range(a_min, min(a_max + 1, a_min + 30)):
            R = 4 * a - p
            if R <= 0:
                continue
            D = p * a
            n_div, pts = curve_genus_info(p, a)
            if n_div > 0:
                pts_str = "; ".join(f"({t[0]},{t[1]},{t[2]})" for t in pts[:3])
                extra = f" +{n_div-3} more" if n_div > 3 else ""
                print(f"    a={a}: D^2/(na)^2 has {n_div} valid divisor pairs. "
                      f"R=4a-n={R}. {pts_str}{extra}")

    # ---------------------------------------------------------------------------
    # Greedy vs systematic for p = 1 mod 4
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"GREEDY vs OPTIMAL for p = 1 mod 4:")
    print(f"{'Prime':>8} {'Greedy_c':>12} {'Optimal_c':>12} {'Ratio':>8} {'Greedy wins?':>14}")
    print(f"{'-'*60}")
    greedy_wins = 0
    optimal_wins = 0
    ties = 0
    primes_1mod4 = [p for p in primes if p > 3 and p % 4 == 1]
    for p in primes_1mod4[:30]:  # first 30
        g = greedy_es(p)
        o = systematic_1mod4(p)
        if g and o:
            gc = g[2]
            oc = o[2]
            ratio = gc / oc if oc > 0 else float('inf')
            winner = "GREEDY" if gc <= oc else ("TIE" if gc == oc else "OPTIMAL")
            if gc < oc: greedy_wins += 1
            elif gc == oc: ties += 1
            else: optimal_wins += 1
            print(f"{p:>8} {gc:>12} {oc:>12} {ratio:>8.2f} {winner:>14}")
        elif g:
            print(f"{p:>8} {g[2]:>12} {'N/A':>12} {'---':>8} {'GREEDY':>14}")
        elif o:
            print(f"{p:>8} {'N/A':>12} {o[2]:>12} {'---':>8} {'OPTIMAL':>14}")

    print(f"\n  Summary (first 30 primes = 1 mod 4):")
    print(f"    Greedy better: {greedy_wins}, Optimal better: {optimal_wins}, Ties: {ties}")

    # Full comparison
    g_better = o_better = tie = 0
    for p in primes_1mod4:
        g = greedy_es(p)
        o = systematic_1mod4(p)
        if g and o:
            if g[2] < o[2]: g_better += 1
            elif g[2] == o[2]: tie += 1
            else: o_better += 1
    print(f"  Full (all {len(primes_1mod4)} primes = 1 mod 4 up to {LIMIT}):")
    print(f"    Greedy better: {g_better}, Optimal better: {o_better}, Ties: {tie}")

    # ---------------------------------------------------------------------------
    # Correlation: solution count vs prime size
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"SOLUTION COUNT GROWTH:")
    brackets = [(5, 100), (100, 200), (200, 400), (400, 600), (600, 800), (800, 1001)]
    for lo, hi in brackets:
        bucket = [sol_count[p] for p in primes if lo <= p < hi]
        if bucket:
            print(f"  Primes in [{lo},{hi}): n={len(bucket)}, "
                  f"avg={sum(bucket)/len(bucket):.1f}, "
                  f"min={min(bucket)}, max={max(bucket)}")

    # ---------------------------------------------------------------------------
    # Mod-42 heatmap: which residues are HARDEST?
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"MOD-42 DIFFICULTY RANKING (by average solution count, ascending):")
    ranked = []
    for r in mod42_counts:
        entries = mod42_counts[r]
        cnts = [c for _, c in entries]
        if len(entries) >= 2:  # need enough data
            ranked.append((sum(cnts)/len(cnts), r, len(entries)))
    ranked.sort()
    print(f"{'Residue':>8} {'AvgSol':>10} {'#Primes':>8}  Note")
    print(f"{'-'*50}")
    for avg, r, n in ranked[:15]:
        note = ""
        if r % 2 == 0: note = "(even residue)"
        if r % 3 == 0: note = "(div by 3)"
        if r % 7 == 0: note = "(div by 7)"
        if gcd(r, 42) == 1: note = "(coprime to 42)"
        print(f"{r:>8} {avg:>10.1f} {n:>8}  {note}")

    # ---------------------------------------------------------------------------
    # Discovery: primes needing large denominators
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"PRIMES REQUIRING LARGEST DENOMINATORS (max c across all solutions):")
    large_c = []
    for p in primes:
        if all_sols[p]:
            min_c = min(s[2] for s in all_sols[p])
            large_c.append((min_c, p))
    large_c.sort(reverse=True)
    print(f"{'Prime':>8} {'MinMax_c':>12} {'mod4':>6} {'mod24':>6} {'#Sol':>6}")
    print(f"{'-'*50}")
    for mc, p in large_c[:15]:
        print(f"{p:>8} {mc:>12} {p%4:>6} {p%24:>6} {sol_count[p]:>6}")

    # ---------------------------------------------------------------------------
    # Key discovery check: Bott slot (mod 8) correlation
    # ---------------------------------------------------------------------------
    print(f"\n{'='*72}")
    print(f"BOTT SLOT (mod 8) ANALYSIS:")
    mod8_data = defaultdict(list)
    for p in primes:
        if p > 3:
            mod8_data[p % 8].append(sol_count[p])
    for r in sorted(mod8_data.keys()):
        cnts = mod8_data[r]
        print(f"  p = {r} mod 8: n={len(cnts)}, avg={sum(cnts)/len(cnts):.1f}, "
              f"min={min(cnts)}, max={max(cnts)}")

    print(f"\n{'='*72}")
    print(f"DONE. Erdos-Straus verified for all primes <= {LIMIT}.")
    print(f"Total decompositions found: {sum(sol_count.values())}")
    print(f"Key insight: the conic (Rb-D)(Rc-D)=D^2 reduces ES to divisor counting.")


if __name__ == "__main__":
    main()
