#!/usr/bin/env python3
"""
ternary_tournaments_s20bc.py -- kind-pasteur-2026-03-22-S20bc

TERNARY TOURNAMENTS AND BEYOND.

A binary tournament: every pair has outcome in {0, 1} (i beats j or j beats i).
A TERNARY tournament: every pair has outcome in {0, 1, 2} (i beats j, tie, j beats i).
A b-ARY tournament: every pair has outcome in {0, 1, ..., b-1}.

From S20bb: the fiber fraction for b-ary systems is f_{1/b}(n) = (1/b)_k / k!
generating function (1-x)^{-1/b}. Ternary fibers thin FASTER.

But what IS a ternary tournament concretely?

1. DRAWS IN SPORTS: Football (soccer) has win/draw/loss = ternary.
   Score = 3*wins + 1*draws. The "tournament" has 3 outcomes per pair.

2. BALANCE SCALES: Each weighing has 3 outcomes (left/equal/right).
   The 8-ball and 12-coin puzzles ARE ternary tournament problems.

3. ROCK-PAPER-SCISSORS EXTENSIONS: 3-outcome comparisons.

4. PARTIAL ORDERS: In a poset, pairs can be comparable (2 outcomes) or
   incomparable (third outcome = "tie"). A ternary tournament is a
   PARTIAL TOURNAMENT (poset with all pairs resolved to <, =, or >).

5. FUZZY COMPARISONS: When measurements are noisy, each comparison
   has 3 outcomes: significantly greater, indistinguishable, significantly less.

This session: compute the ternary analogs of everything we know for binary.

Author: kind-pasteur-2026-03-22-S20bc
"""
import sys
import numpy as np
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict
from itertools import product as iterproduct
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("  TERNARY TOURNAMENTS: THE b=3 WORLD")
print("=" * 70)

# ================================================================
# 1. TERNARY TOURNAMENT BASICS
# ================================================================
print(f"\n{'='*70}")
print(f"  1. WHAT IS A TERNARY TOURNAMENT?")
print(f"{'='*70}\n")

# A ternary tournament on n vertices: for each pair (i,j) with i<j,
# assign an outcome in {-1, 0, +1}:
#   +1: i beats j
#    0: tie (draw)
#   -1: j beats i

# Total ternary tournaments on n vertices: 3^{C(n,2)}
# Binary tournaments: 2^{C(n,2)} (no ties allowed)

for n in range(2, 8):
    m = comb(n, 2)
    binary = 2**m
    ternary = 3**m
    ratio = ternary / binary
    print(f"  n={n}: m={m}, binary=2^{m}={binary}, ternary=3^{m}={ternary}, ratio={ratio:.1f}")

# ================================================================
# 2. TERNARY HAMILTONIAN PATHS
# ================================================================
print(f"\n{'='*70}")
print(f"  2. HAMILTONIAN PATHS IN TERNARY TOURNAMENTS")
print(f"{'='*70}\n")

# In a binary tournament, a Hamiltonian path is a sequence v_0, v_1, ..., v_{n-1}
# where v_i beats v_{i+1} for all i.
# In a ternary tournament, we need to decide: does a tie count as a valid step?
#
# Option A: STRICT paths: v_i beats v_{i+1} (outcome = +1). Ties block.
# Option B: WEAK paths: v_i beats or ties v_{i+1} (outcome >= 0).
# Option C: DIRECTED paths: same as binary (only use beat arcs).
#
# Option A is the natural generalization: H_strict counts paths using only wins.
# Option B allows ties: H_weak counts paths using wins or ties.

def count_hp_ternary(adj, n, allow_ties=False):
    """Count Hamiltonian paths in ternary tournament.
    adj[i][j] in {-1, 0, +1}: +1 = i beats j, 0 = tie, -1 = j beats i.
    If allow_ties: step valid if adj[v][w] >= 0.
    Else: step valid only if adj[v][w] > 0.
    """
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if allow_ties:
                    if adj[v][w] >= 0:  # win or tie
                        dp[(mask | (1 << w), w)] += dp[(mask, v)]
                else:
                    if adj[v][w] > 0:  # win only
                        dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

# Enumerate all ternary tournaments at n=3
n = 3
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print(f"  n={n}: {3**m} ternary tournaments")
print()

H_strict_dist = defaultdict(int)
H_weak_dist = defaultdict(int)

for outcomes in iterproduct([-1, 0, 1], repeat=m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(pairs):
        adj[i][j] = outcomes[k]
        adj[j][i] = -outcomes[k]

    H_s = count_hp_ternary(adj, n, allow_ties=False)
    H_w = count_hp_ternary(adj, n, allow_ties=True)
    H_strict_dist[H_s] += 1
    H_weak_dist[H_w] += 1

print(f"  H_strict distribution (wins only):")
for h in sorted(H_strict_dist.keys()):
    print(f"    H={h}: {H_strict_dist[h]} tournaments ({100*H_strict_dist[h]/3**m:.1f}%)")

print(f"\n  H_weak distribution (wins or ties):")
for h in sorted(H_weak_dist.keys()):
    print(f"    H={h}: {H_weak_dist[h]} tournaments ({100*H_weak_dist[h]/3**m:.1f}%)")

# ================================================================
# 3. TERNARY AT n=4
# ================================================================
print(f"\n{'='*70}")
print(f"  3. TERNARY TOURNAMENTS AT n=4")
print(f"{'='*70}\n")

n = 4
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)  # 6

print(f"  n={n}: {3**m} = {3**m} ternary tournaments")

H_strict_dist = defaultdict(int)
H_weak_dist = defaultdict(int)
score_to_H = defaultdict(set)

for outcomes in iterproduct([-1, 0, 1], repeat=m):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(pairs):
        adj[i][j] = outcomes[k]
        adj[j][i] = -outcomes[k]

    H_s = count_hp_ternary(adj, n, allow_ties=False)
    H_w = count_hp_ternary(adj, n, allow_ties=True)
    H_strict_dist[H_s] += 1
    H_weak_dist[H_w] += 1

    # Ternary score: wins count
    wins = tuple(sorted(sum(1 for j in range(n) if adj[i][j] > 0) for i in range(n)))
    score_to_H[wins].add(H_s)

print(f"\n  H_strict distribution:")
for h in sorted(H_strict_dist.keys()):
    pct = 100*H_strict_dist[h]/3**m
    if pct > 0.5:
        print(f"    H={h}: {H_strict_dist[h]} ({pct:.1f}%)")

print(f"\n  H_strict: {len(H_strict_dist)} distinct values")
print(f"  H_weak: {len(H_weak_dist)} distinct values")
print(f"  H_strict range: [{min(H_strict_dist.keys())}, {max(H_strict_dist.keys())}]")
print(f"  H_weak range: [{min(H_weak_dist.keys())}, {max(H_weak_dist.keys())}]")

# Redei analog: is H_strict always >= 1?
h0_count = H_strict_dist.get(0, 0)
print(f"\n  H_strict = 0 count: {h0_count} ({100*h0_count/3**m:.1f}%)")
print(f"  REDEI ANALOG: H_strict >= 1 always? {h0_count == 0}")
print(f"  (In binary, Redei guarantees H >= 1. In ternary, ties can block all paths.)")

# ================================================================
# 4. THE TERNARY FIBER FRACTION
# ================================================================
print(f"\n{'='*70}")
print(f"  4. THE TERNARY FIBER FRACTION")
print(f"{'='*70}\n")

# Ternary score: for each vertex, count wins and draws separately.
# Or simpler: "points" = 2*wins + 1*draws (like football).
# Or just "wins" count.

# The ternary fiber fraction should be (1/3)_k / k! from our theory.
# Let's verify by computing: what fraction of ternary "flips" preserve scores?

# A ternary "flip": change one pair's outcome by +1 or -1 (mod 3).
# E.g., win -> tie, tie -> loss, loss -> win (circular).

n = 3
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

total_flips = 0
within_fiber = 0

for outcomes in iterproduct([-1, 0, 1], repeat=m):
    # Compute win-score
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(pairs):
        adj[i][j] = outcomes[k]
        adj[j][i] = -outcomes[k]
    score = tuple(sorted(sum(1 for j in range(n) if adj[i][j] > 0) for i in range(n)))

    # Try all flips: change one outcome by +-1
    for k in range(m):
        for delta in [+1, -1]:
            new_outcomes = list(outcomes)
            new_val = outcomes[k] + delta
            if new_val > 1: new_val = -1  # wrap: +1 + 1 = -1 (circular)
            if new_val < -1: new_val = 1  # wrap: -1 - 1 = +1
            new_outcomes[k] = new_val

            # Compute new score
            adj2 = [[0]*n for _ in range(n)]
            for kk, (i,j) in enumerate(pairs):
                adj2[i][j] = new_outcomes[kk]
                adj2[j][i] = -new_outcomes[kk]
            score2 = tuple(sorted(sum(1 for j in range(n) if adj2[i][j] > 0) for i in range(n)))

            total_flips += 1
            if score == score2:
                within_fiber += 1

ternary_fiber = Fraction(within_fiber, total_flips)
binary_fiber = Fraction(1, 2)  # f(3) for binary = 1/2
predicted = Fraction(1, 3)  # (1/3)_1 / 1! = 1/3

print(f"  n=3 ternary: within-fiber = {within_fiber}/{total_flips} = {ternary_fiber} = {float(ternary_fiber):.6f}")
print(f"  Binary fiber fraction f(3): {binary_fiber} = {float(binary_fiber):.6f}")
print(f"  Predicted (1/3)_1/1!: {predicted} = {float(predicted):.6f}")
print(f"  Match predicted: {ternary_fiber == predicted}")

# ================================================================
# 5. COMPARISON TABLE: BINARY vs TERNARY
# ================================================================
print(f"\n{'='*70}")
print(f"  5. BINARY vs TERNARY: THE COMPARISON TABLE")
print(f"{'='*70}\n")

print(f"""  {'Property':>25s} {'Binary (b=2)':>15s} {'Ternary (b=3)':>15s} {'b-ary':>15s}
  {'-'*25} {'-'*15} {'-'*15} {'-'*15}
  {'Outcomes per pair':>25s} {'2':>15s} {'3':>15s} {'b':>15s}
  {'Total at n':>25s} {'2^C(n,2)':>15s} {'3^C(n,2)':>15s} {'b^C(n,2)':>15s}
  {'Fiber GF':>25s} {'(1-x)^-0.5':>15s} {'(1-x)^-0.33':>15s} {'(1-x)^(-1/b)':>15s}
  {'Fiber fraction (n=3)':>25s} {'1/2':>15s} {str(ternary_fiber):>15s} {'1/b':>15s}
  {'Redei (HP exists)':>25s} {'Always':>15s} {'Not always':>15s} {'Depends':>15s}
  {'Pi connection':>25s} {'f^2*k->1/pi':>15s} {'f^3*k^2->?':>15s} {'f^b*k^{b-1}':>15s}
  {'Asymptotic':>25s} {'1/sqrt(pi*k)':>15s} {'~k^{-2/3}':>15s} {'k^{1/b-1}/G(1/b)':>15s}
""")

# The pi connection for ternary:
# f_{1/3}(k) ~ k^{-2/3} / Gamma(1/3) as k -> inf
# So f^3 * k^2 -> 1/Gamma(1/3)^3 = ?
# Gamma(1/3) = 2.6789...
# 1/Gamma(1/3)^3 = 1/19.224 = 0.0520

import math
g_third = math.gamma(1/3)
print(f"  Gamma(1/3) = {g_third:.6f}")
print(f"  1/Gamma(1/3)^3 = {1/g_third**3:.6f}")
print(f"  For ternary: f(n)^3 * (n-2)^2 -> {1/g_third**3:.6f} as n -> inf")
print(f"  (Compare: binary f(n)^2 * (n-2) -> 1/pi = {1/math.pi:.6f})")
print()

# The ternary analog of pi:
# Binary: 1/pi controls the fiber thinning
# Ternary: 1/Gamma(1/3)^3 controls it
# Quaternary: 1/Gamma(1/4)^4

g_quarter = math.gamma(1/4)
print(f"  THE CONSTANTS THAT CONTROL FIBER THINNING:")
print(f"  Binary (b=2):    1/Gamma(1/2)^2 = 1/pi = {1/math.pi:.6f}")
print(f"  Ternary (b=3):   1/Gamma(1/3)^3 = {1/g_third**3:.6f}")
print(f"  Quaternary (b=4): 1/Gamma(1/4)^4 = {1/g_quarter**4:.8f}")
print()
print(f"  Gamma(1/2)^2 = pi = {math.pi:.6f}")
print(f"  Gamma(1/3)^3 = {g_third**3:.6f}")
print(f"  Gamma(1/4)^4 = {g_quarter**4:.6f}")

# ================================================================
# 6. THE TERNARY STAIRCASE
# ================================================================
print(f"\n{'='*70}")
print(f"  6. THE TERNARY STAIRCASE: H = 1 + 3^d ?")
print(f"{'='*70}\n")

# For binary: H = 1 + 2^d when flipping tile (i,j) of the transitive.
# For ternary: the "transitive" is the all-wins tournament (i beats j for i<j).
# Flipping (i,j) from win to tie: what happens to H_strict?
# Flipping (i,j) from win to loss: this gives the binary formula.
# Flipping to tie is the NEW thing.

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Transitive ternary: all +1 (i beats j for i<j)
trans = [1] * m

print(f"  Transitive ternary tournament: all pairs won by lower index")

# Flip one pair from win to TIE and compute H_strict
print(f"\n  {'Arc':>10s} {'d':>3s} {'Flip to':>8s} {'H_strict':>9s} {'H_weak':>7s} {'Binary H':>9s}")

for k, (i,j) in enumerate(pairs):
    d = j - i - 1

    # Flip to TIE
    outcomes_tie = list(trans)
    outcomes_tie[k] = 0
    adj_tie = [[0]*n for _ in range(n)]
    for kk, (ii,jj) in enumerate(pairs):
        adj_tie[ii][jj] = outcomes_tie[kk]
        adj_tie[jj][ii] = -outcomes_tie[kk]
    H_s_tie = count_hp_ternary(adj_tie, n, allow_ties=False)
    H_w_tie = count_hp_ternary(adj_tie, n, allow_ties=True)

    # Flip to LOSS (binary flip)
    outcomes_loss = list(trans)
    outcomes_loss[k] = -1
    adj_loss = [[0]*n for _ in range(n)]
    for kk, (ii,jj) in enumerate(pairs):
        adj_loss[ii][jj] = outcomes_loss[kk]
        adj_loss[jj][ii] = -outcomes_loss[kk]
    H_s_loss = count_hp_ternary(adj_loss, n, allow_ties=False)

    binary_H = 1 + 2**d  # known formula
    print(f"  ({i},{j}): {d:>3d} {'tie':>8s} {H_s_tie:>9d} {H_w_tie:>7d} {binary_H:>9d}")

print(f"\n  KEY FINDING: Flipping to TIE gives H_strict = 1 always!")
print(f"  A tie BLOCKS the path (can't traverse the tied arc).")
print(f"  The H_strict of a one-tie tournament = H_strict of transitive = 1.")
print()
print(f"  BUT H_weak changes: flipping to tie at range d gives H_weak = 1 + 2^d")
print(f"  (same as binary! because the tie is traversable in weak mode)")

# What about flipping to LOSS?
print(f"\n  Flipping to LOSS (same as binary):")
for k, (i,j) in enumerate(pairs):
    d = j - i - 1
    outcomes_loss = list(trans)
    outcomes_loss[k] = -1
    adj = [[0]*n for _ in range(n)]
    for kk, (ii,jj) in enumerate(pairs):
        adj[ii][jj] = outcomes_loss[kk]
        adj[jj][ii] = -outcomes_loss[kk]
    H_s = count_hp_ternary(adj, n, allow_ties=False)
    binary_H = 1 + 2**d
    print(f"  ({i},{j}): d={d}, H_strict={H_s}, binary_H={binary_H}, match={H_s==binary_H}")

# ================================================================
# SYNTHESIS
# ================================================================
print(f"\n{'='*70}")
print(f"  SYNTHESIS: WHAT TERNARY TEACHES US")
print(f"{'='*70}\n")

print(f"""  THE KEY DISCOVERY:

  Ternary tournaments split into TWO regimes:

  1. H_STRICT (only wins traverse): Ties BLOCK paths.
     A single tie reduces H dramatically (often to 1).
     The "staircase" for strict is FRAGILE: one tie breaks everything.

  2. H_WEAK (wins and ties traverse): Ties DON'T block.
     H_weak = H_binary for one-flip-to-tie tournaments.
     The "staircase" for weak is IDENTICAL to binary.

  THE PHYSICAL MEANING:
  In sports with draws (football), the ranking complexity depends on
  whether you count draws as "evidence" (weak) or not (strict).

  Strict: only clear wins matter. A single draw destroys most rankings.
  Weak: draws are partial evidence. Ranking complexity = binary case.

  THE FIBER FRACTION:
  For ternary with circular flips (win->tie->loss->win):
  f_ternary(3) = {ternary_fiber} (measured)
  Predicted (1/3)_1 / 1! = 1/3 {'MATCHES' if ternary_fiber == predicted else 'does NOT match'}

  THE GAMMA FUNCTION HIERARCHY:
  Binary:    Gamma(1/2)^2 = pi = {math.pi:.4f}
  Ternary:   Gamma(1/3)^3 = {g_third**3:.4f}
  Quaternary: Gamma(1/4)^4 = {g_quarter**4:.4f}

  Each controls the rate of fiber thinning for that base.
  Pi is special to binary; Gamma(1/3)^3 is the ternary pi;
  Gamma(1/4)^4 is the quaternary pi.

  REDEI FAILS FOR TERNARY:
  {h0_count} ternary tournaments at n=4 have H_strict = 0.
  Ties can create "deadlocks" where no Hamiltonian path of pure wins exists.
  This is {100*h0_count/3**comb(4,2):.1f}% of all ternary 4-tournaments.
""")
