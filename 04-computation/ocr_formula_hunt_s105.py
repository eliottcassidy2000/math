#!/usr/bin/env python3
"""
ocr_formula_hunt_s105.py — Hunting for an exact OCR formula
opus-2026-03-21-S105

KNOWN EXACT VALUES (from S103/S104):
  OCR(3) = 1
  OCR(4) = 1
  OCR(5) = 129/133
  OCR(6) = 460807/480480

  1 - OCR(5) = 4/133
  1 - OCR(6) = 19673/480480

The goal: find a pattern, a formula, or at least understand WHY
these specific fractions appear.

APPROACH:
1. Factor the numerators and denominators
2. Connect to known tournament quantities
3. Derive OCR from the OCF (H = 1 + 2α₁ + 4α₂ + ...)
4. Use the c₃ formula (c₃ = C(n+1,3)/4 - S₂/2) to understand the residual
5. Compute the EXACT Var(H) and Var(H|scores) as rational functions of n
"""

from math import comb, factorial, gcd, log, sqrt
from fractions import Fraction
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

print("=" * 72)
print("  HUNTING FOR AN EXACT OCR FORMULA")
print("  opus-2026-03-21-S105")
print("=" * 72)
print()

# =============================================================================
# PART 1: FACTOR THE KNOWN OCR FRACTIONS
# =============================================================================

print("=" * 72)
print("  PART 1: FACTORING THE OCR FRACTIONS")
print("=" * 72)
print()

def factorize(n):
    if n == 0: return {0: 1}
    factors = {}
    d = 2
    m = abs(n)
    while d * d <= m:
        while m % d == 0:
            factors[d] = factors.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        factors[m] = factors.get(m, 0) + 1
    return factors

# The exact values
ocr_data = {
    3: {'var_H': Fraction(3, 4), 'var_cond': Fraction(0), 'ocr': Fraction(1)},
    5: {'var_H': Fraction(285, 16), 'var_cond': Fraction(15, 28), 'ocr': Fraction(129, 133)},
    6: {'var_H': Fraction(585, 4), 'var_cond': Fraction(59019, 9856), 'ocr': Fraction(460807, 480480)},
}

for n, d in sorted(ocr_data.items()):
    res = 1 - d['ocr']
    print(f"  n = {n}:")
    print(f"    Var(H) = {d['var_H']} = {float(d['var_H']):.6f}")
    print(f"    Var(H|s) = {d['var_cond']} = {float(d['var_cond']):.6f}")
    print(f"    OCR = {d['ocr']} = {float(d['ocr']):.10f}")
    print(f"    1-OCR = {res}")
    if res != 0:
        print(f"    Factorization of 1-OCR:")
        print(f"      numerator {res.numerator}: {factorize(res.numerator)}")
        print(f"      denominator {res.denominator}: {factorize(res.denominator)}")
    print()

print("""
  OBSERVATIONS:
  1-OCR(5) = 4/133 = 2²/(7·19)
  1-OCR(6) = 19673/480480 = 19673/(2⁵·3·5·7·11·13)

  The denominator of OCR(5) is 133 = 7 × 19.
  The denominator of OCR(6) is 480480 = 2⁵ × 3 × 5 × 7 × 11 × 13.

  Is 480480 related to tournament quantities?
  C(6,2) = 15. 15! = 1307674368000. 480480 = ... let's check.
""")

print(f"  480480 = {factorize(480480)}")
print(f"  133 = {factorize(133)}")
print(f"  C(5,2)·(C(5,2)-1)·...·1 products:")
for k in range(1, 15):
    p = 1
    for i in range(k):
        p *= (comb(5,2) - i) if comb(5,2) - i > 0 else 1
    if p == 133 or p % 133 == 0:
        print(f"    Product of {k} terms from C(5,2)=10 down: {p}")

# Check if 133 relates to n=5 tournament counting
print(f"\n  n=5: C(5,2) = 10, 2^10 = 1024")
print(f"  1024/133 = {1024/133:.4f}")
print(f"  133 = 7 × 19")
print(f"  7 appears in tournament theory as the SPLIT prime")
print(f"  19 = C(5,2) + C(5,1) = 10 + 9? No. 19 is the 8th prime.")
print()

# Connection to the ambiguous class
print(f"  At n=5, the residual comes ENTIRELY from score class (1,2,2,2,3).")
print(f"  This class has 280 tournaments = 5!·C(5,2)/... hmm")
print(f"  280 = 2³ × 5 × 7 = C(5,2) × 28 = 10 × 28")
print(f"  280/1024 = 35/128")
print()
print(f"  Var(H|sc) = 96/49")
print(f"  Var(H|scores) = (35/128) × (96/49) = 35×96/(128×49) = 3360/6272 = 15/28")
print(f"  15 = C(6,2) = C(5+1, 2). Interesting!")
print(f"  28 = C(8,2) = C(7+1, 2). Hmm.")
print(f"  Or: 28 = 4 × 7.")
print(f"  So 15/28 = C(6,2) / (4 × 7).")
print()

# =============================================================================
# PART 2: EXACT Var(H) FORMULA
# =============================================================================

print("=" * 72)
print("  PART 2: EXACT Var(H) — CAN WE DERIVE IT?")
print("=" * 72)
print()

print("""
  Var(H) = E[H²] - E[H]².

  We know E[H] = n!/2^{n-1} (each permutation is a HP with probability 1/2^{n-1}).

  For E[H²], we need the second moment. H² = (Σ_σ Π_arcs)² = Σ_{σ,τ} Π_{arcs of σ ∩ τ}.

  Two Hamiltonian paths σ and τ share some arcs. If they share k arcs:
    E[Π_{arcs of σ} × Π_{arcs of τ}] = 1/2^{2(n-1)-k}
  because 2(n-1)-k DISTINCT arc choices are made.

  So: E[H²] = Σ_{σ,τ} 1/2^{2(n-1)-overlap(σ,τ)}
            = Σ_{k=0}^{n-1} N(k) / 2^{2(n-1)-k}

  where N(k) = #{pairs (σ,τ) with exactly k common arcs}.

  This is a COMBINATORIAL quantity depending on the overlap structure
  of Hamiltonian paths in K_n.
""")

# Compute N(k) for small n by direct enumeration
for n in [3, 4, 5]:
    perms = list(permutations(range(n)))
    paths = []
    for perm in perms:
        arcs = frozenset((perm[i], perm[i+1]) for i in range(n-1))
        paths.append(arcs)

    overlap_counter = Counter()
    for i, p1 in enumerate(paths):
        for j, p2 in enumerate(paths):
            k = len(p1 & p2)
            overlap_counter[k] += 1

    # N(k) from our data
    N_k = {}
    for k, count in sorted(overlap_counter.items()):
        N_k[k] = count

    # Verify E[H²]
    E_H2_from_overlap = sum(N_k.get(k, 0) / (2 ** (2*(n-1) - k))
                            for k in range(n))
    E_H = factorial(n) / (2 ** (n-1))
    Var_H = E_H2_from_overlap - E_H**2

    print(f"  n = {n}:")
    print(f"    N(k) overlap distribution:")
    for k in range(n):
        print(f"      k={k}: N(k) = {N_k.get(k, 0)} "
              f"(fraction of all (σ,τ) pairs with k shared arcs)")
    print(f"    E[H] = {E_H}")
    print(f"    E[H²] from overlap = {E_H2_from_overlap}")
    print(f"    Var(H) = {Var_H}")
    var_exact = ocr_data.get(n, {}).get('var_H', None)
    if var_exact is not None:
        print(f"    Var(H) exact (from S104) = {float(var_exact)}")
    else:
        print(f"    Var(H) exact (from S104) = not available")
    print()

# =============================================================================
# PART 3: WHAT DETERMINES THE RESIDUAL?
# =============================================================================

print("=" * 72)
print("  PART 3: ANATOMY OF THE RESIDUAL")
print("=" * 72)
print()

print("""
  The residual Var(H|scores) comes from score classes where H is not unique.

  KEY INSIGHT FROM OCF: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
  where α_k = #{independent sets of size k in Ω(T)}.

  Since c₃ is determined by scores (c₃ = C(n+1,3)/4 - S₂/2),
  the residual comes from the HIGHER odd cycle counts (c₅, c₇, ...)
  and the ARRANGEMENT of cycles (e_cyc, etc.).

  More precisely: α₁ = c₃ + c₅ + c₇ + ... (total odd cycles in T).
  c₃ is determined by scores. c₅ is NOT.
  So: α₁ = c₃(scores) + c₅ + higher = known + unknown.
  The unknown part c₅ + higher contributes the residual.

  AT n=5: α₁ = c₃ + c₅ (no 7-cycles at n=5).
  c₃ is determined by scores: c₃ = 5 - S₂/2.
  c₅ is NOT determined by scores for the class (1,2,2,2,3).
  The residual IS the variance of c₅ within this class.
""")

# Compute: how much of α₁ is explained by c₃?
print("  Variance decomposition of α₁ at n=5:")
print()

n = 5
m = comb(n, 2)
alpha1_by_score = defaultdict(list)
c3_by_score = defaultdict(list)
c5_by_score = defaultdict(list)
H_by_score = defaultdict(list)

for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    scores_list = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1:
                adj[i][j] = 1; scores_list[i] += 1
            else:
                adj[j][i] = 1; scores_list[j] += 1
            idx += 1
    score_key = tuple(sorted(scores_list))

    c3 = sum(1 for a,b,c in combinations(range(n), 3)
             if (adj[a][b] and adj[b][c] and adj[c][a]) or
                (adj[a][c] and adj[c][b] and adj[b][a]))

    # Count 5-cycles (all directed 5-cycles through all 5 vertices)
    c5 = 0
    for perm in permutations(range(n)):
        if all(adj[perm[k]][perm[(k+1)%n]] for k in range(n)):
            c5 += 1
    c5 //= n  # each cycle counted n times

    alpha1 = c3 + c5
    H = sum(1 for perm in permutations(range(n))
            if all(adj[perm[k]][perm[k+1]] for k in range(n-1)))

    alpha1_by_score[score_key].append(alpha1)
    c3_by_score[score_key].append(c3)
    c5_by_score[score_key].append(c5)
    H_by_score[score_key].append(H)

# Check: is c₃ constant within each score class?
print(f"  Is c₃ constant within each score class at n=5?")
for sc in sorted(c3_by_score.keys()):
    vals = set(c3_by_score[sc])
    print(f"    {sc}: c₃ values = {sorted(vals)} {'✓ CONSTANT' if len(vals) == 1 else '✗ VARIES'}")

print()

# Check: is c₅ the ONLY source of variance?
print(f"  Is c₅ the source of H variance within score classes?")
for sc in sorted(c5_by_score.keys()):
    c5_vals = set(c5_by_score[sc])
    H_vals = set(H_by_score[sc])
    if len(H_vals) > 1:
        print(f"    {sc}: c₅ = {sorted(c5_vals)}, H = {sorted(H_vals)}")
        # Verify H = 1 + 2*alpha1 at n=5 (where alpha2=0 always)
        for c5_v in sorted(c5_vals):
            c3_v = c3_by_score[sc][0]  # constant
            alpha1_v = c3_v + c5_v
            H_pred = 1 + 2 * alpha1_v
            matching_H = [H_by_score[sc][i] for i in range(len(c5_by_score[sc]))
                         if c5_by_score[sc][i] == c5_v]
            actual_H = matching_H[0] if matching_H else '?'
            print(f"      c₅={c5_v}: α₁=c₃+c₅={c3_v}+{c5_v}={alpha1_v}, "
                  f"H_pred=1+2·{alpha1_v}={H_pred}, H_actual={actual_H} "
                  f"{'✓' if H_pred == actual_H else '✗'}")

print()

# =============================================================================
# PART 4: THE c₅ DISTRIBUTION WITHIN SCORE CLASSES
# =============================================================================

print("=" * 72)
print("  PART 4: c₅ DISTRIBUTION — THE KEY TO THE RESIDUAL")
print("=" * 72)
print()

print("""
  Since the residual comes entirely from the VARIANCE OF c₅ within
  score classes, we need to understand what determines c₅ given scores.

  At n=5: c₅ = number of directed 5-cycles (Hamiltonian cycles).
  This is the number of CYCLIC Hamiltonian paths.

  CLAIM: c₅ is determined by the tournament's REGULAR decomposition.
  A tournament can be decomposed as: near-transitive part + cycle part.
  The more "cyclic" the tournament, the more 5-cycles it has.
""")

# For the ambiguous class (1,2,2,2,3):
sc_ambig = (1, 2, 2, 2, 3)
c5_list = c5_by_score[sc_ambig]
H_list = H_by_score[sc_ambig]

c5_counter = Counter(c5_list)
print(f"  Score class {sc_ambig}:")
print(f"  c₅ distribution:")
for c5_val in sorted(c5_counter.keys()):
    count = c5_counter[c5_val]
    H_for_c5 = [H_list[i] for i in range(len(c5_list)) if c5_list[i] == c5_val]
    H_val = H_for_c5[0]
    print(f"    c₅ = {c5_val}: {count} tournaments ({count/len(c5_list)*100:.1f}%), H = {H_val}")

# The exact Var(c₅) within this class
c5_arr = [float(x) for x in c5_list]
var_c5 = sum((x - sum(c5_arr)/len(c5_arr))**2 for x in c5_arr) / len(c5_arr)
print(f"\n  E[c₅ | sc] = {sum(c5_arr)/len(c5_arr):.4f}")
print(f"  Var(c₅ | sc) = {var_c5:.6f}")
print(f"\n  Since H = 1 + 2·(c₃ + c₅) and c₃ is constant within the class:")
print(f"  Var(H | sc) = 4 × Var(c₅ | sc) = {4 * var_c5:.6f}")
print(f"  Exact Var(H | sc) = 96/49 = {96/49:.6f}")
print(f"  Match: {abs(4 * var_c5 - 96/49) < 0.001}")
print()

# =============================================================================
# PART 5: Var(c₅ | scores) AS THE KEY QUANTITY
# =============================================================================

print("=" * 72)
print("  PART 5: THE OCR DECOMPOSITION THEOREM")
print("=" * 72)
print()

print("""
  THEOREM (OCR decomposition at n=5):

  OCR(5) = 1 - Var(H|scores)/Var(H)
         = 1 - 4·Var(c₅|scores) / Var(H)

  Because H = 1 + 2·α₁ = 1 + 2·(c₃ + c₅) and c₃ is determined by scores,
  Var(H|scores) = 4·Var(c₅|scores).

  This generalizes: for any n,
    Var(H|scores) = 4·Var(c₅|scores) + 16·Var(α₂|scores, c₅) + ...
                  ≈ 4·Var(c₅|scores)  [dominant term]

  So the OCR is controlled by HOW MUCH c₅ VARIES within score classes.

  OCR ≈ 1 - 4·Var(c₅|scores) / Var(H)

  To get an exact formula, we need:
  (a) Var(H) as a function of n  [E[H²] via overlap counting]
  (b) Var(c₅|scores) as a function of n  [harder — depends on cycle structure]
""")

# Verify at n=5
total_var_c5_weighted = 0
for sc, c5_vals in c5_by_score.items():
    c5_arr_sc = [float(x) for x in c5_vals]
    var_c5_sc = sum((x - sum(c5_arr_sc)/len(c5_arr_sc))**2 for x in c5_arr_sc) / len(c5_arr_sc)
    weight = len(c5_vals) / (2**m)
    total_var_c5_weighted += weight * var_c5_sc

print(f"  Verification at n=5:")
print(f"    Var(c₅|scores) = {total_var_c5_weighted:.6f}")
print(f"    4 × Var(c₅|scores) = {4*total_var_c5_weighted:.6f}")
print(f"    Var(H|scores) = {float(Fraction(15,28)):.6f}")
print(f"    Match: {abs(4*total_var_c5_weighted - float(Fraction(15,28))) < 0.001}")
print()

# =============================================================================
# PART 6: THE Var(H) EXACT FORMULA
# =============================================================================

print("=" * 72)
print("  PART 6: EXACT Var(H) FORMULA")
print("=" * 72)
print()

print("""
  From Part 2, Var(H) = E[H²] - E[H]² where:
    E[H] = n! / 2^{n-1}
    E[H²] = Σ_k N(k) / 2^{2(n-1)-k}

  The N(k) (overlap structure of Hamiltonian paths) is key.

  Let's compute E[H] and Var(H) as exact fractions.
""")

for n in [3, 4, 5, 6]:
    E_H = Fraction(factorial(n), 2**(n-1))
    print(f"  n={n}: E[H] = {E_H} = {float(E_H):.4f}")

    # E[H²] from direct computation
    m = comb(n, 2)
    if n <= 6:
        t0 = time.time()
        H_vals = []
        for mask in range(2**m):
            adj_flat = [0] * (n*n)
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if (mask >> idx) & 1: adj_flat[i*n+j] = 1
                    else: adj_flat[j*n+i] = 1
                    idx += 1

            # DP for H
            dp = [0] * ((1 << n) * n)
            for v in range(n): dp[(1 << v) * n + v] = 1
            for S in range(1, 1 << n):
                for v in range(n):
                    if not (S & (1 << v)): continue
                    val = dp[S * n + v]
                    if val == 0: continue
                    for u in range(n):
                        if S & (1 << u): continue
                        if adj_flat[v * n + u]:
                            dp[(S | (1 << u)) * n + u] += val
            full = (1 << n) - 1
            H = sum(dp[full * n + v] for v in range(n))
            H_vals.append(H)

        elapsed = time.time() - t0
        E_H2 = Fraction(sum(h*h for h in H_vals), 2**m)
        Var_H = E_H2 - E_H**2
        print(f"    E[H²] = {E_H2} = {float(E_H2):.4f}")
        print(f"    Var(H) = {Var_H} = {float(Var_H):.6f}  ({elapsed:.1f}s)")

        # Try to express Var(H) in terms of n
        # Known: Var(H) at n=3: 3/4, n=4: 3, n=5: 285/16, n=6: 585/4
        print()

print()

# Sequence of Var(H): 3/4, 3, 285/16, 585/4
print("  Var(H) sequence:")
vars_H = [Fraction(3,4), Fraction(3), Fraction(285,16), Fraction(585,4)]
for i, (v, n_val) in enumerate(zip(vars_H, [3,4,5,6])):
    ratio_prev = v / vars_H[i-1] if i > 0 else None
    per_arc = v / comb(n_val, 2)
    print(f"    n={n_val}: Var(H) = {v} = {float(v):.4f}, "
          f"per arc = {float(per_arc):.4f}"
          + (f", ratio to prev = {float(ratio_prev):.4f}" if ratio_prev else ""))

print()

# =============================================================================
# PART 7: THE RESIDUAL FRACTION — PATTERN HUNTING
# =============================================================================

print("=" * 72)
print("  PART 7: PATTERN IN THE RESIDUAL FRACTION")
print("=" * 72)
print()

# 1-OCR values
residuals = {
    3: Fraction(0),
    4: Fraction(0),
    5: Fraction(4, 133),
    6: Fraction(19673, 480480),
}

print("  Residual 1-OCR as fraction:")
for n_val in sorted(residuals.keys()):
    r = residuals[n_val]
    if r > 0:
        print(f"    n={n_val}: 1-OCR = {r} = {float(r):.10f}")
        print(f"      num = {r.numerator} = {factorize(r.numerator)}")
        print(f"      den = {r.denominator} = {factorize(r.denominator)}")

print()

# Check: denominator of 1-OCR(5) = Var(H) numerator when cleared?
# Var(H) = 285/16, Var(H|s) = 15/28
# 1-OCR = (15/28)/(285/16) = 15×16/(28×285) = 240/7980 = 4/133
# So denominator 133 = 7980/60 = (28×285)/(15×16) * 4 ... hmm
# Actually: 133 = 7 × 19. And 285 = 3 × 5 × 19. So 19 divides both!
print("  Connection to Var(H):")
print(f"    n=5: Var(H) = 285/16 = (3·5·19)/2⁴")
print(f"    n=5: Var(H|s) = 15/28 = (3·5)/(2²·7)")
print(f"    1-OCR(5) = (15/28)/(285/16) = (15·16)/(28·285)")
print(f"             = 240/7980 = 4/133 = 2²/(7·19)")
print()
print(f"    The 19 in the denominator comes from Var(H) = 285/16 = 3·5·19/16.")
print(f"    The 7 comes from Var(H|s) = 15/28 (the class has 280 = 2³·5·7 tournaments).")
print()

# Can we predict 1-OCR(7)?
# From S103: Var(H) at n=7 = 1611.91..., Var(H|s) = 67.53...
# 1-OCR(7) ≈ 67.53/1611.91 ≈ 0.04189
print(f"  n=7 (from S103 float computation):")
print(f"    1-OCR(7) ≈ 67.5281/1611.9141 ≈ {67.5281/1611.9141:.10f}")
print(f"    If this is a fraction p/q, what could it be?")
# Try to find a rational approximation
from fractions import Fraction
target = Fraction(67528, 1611914)  # approximate
approx = target.limit_denominator(10000000)
print(f"    Best rational approximation: {approx} = {float(approx):.10f}")
print(f"    num = {approx.numerator}, den = {approx.denominator}")
print(f"    num factors: {factorize(approx.numerator)}")
print(f"    den factors: {factorize(approx.denominator)}")
print()

# =============================================================================
# PART 8: THE KEY INSIGHT — OCR AND TOTAL ODD CYCLES
# =============================================================================

print("=" * 72)
print("  PART 8: THE INSIGHT — OCR = f(RESIDUAL CYCLE VARIANCE)")
print("=" * 72)
print()

print("""
  SUMMARY OF WHAT WE KNOW:

  1. H = 1 + 2α₁ + 4α₂ + ... (OCF)
  2. α₁ = c₃ + c₅ + c₇ + ... (total odd cycles)
  3. c₃ = C(n+1,3)/4 - S₂/2 (PROVED, determined by scores)
  4. c₅, c₇, ... are NOT determined by scores
  5. Var(H|scores) = 4·Var(α₁|scores) + higher order terms
                   = 4·Var(c₅ + c₇ + ...|scores) + h.o.t.
  6. At n=5: Var(H|scores) = 4·Var(c₅|scores) EXACTLY

  THE FORMULA FOR OCR:

    OCR(n) = 1 - 4·Var(c₅ + c₇ + ... | scores) / Var(H)
           ≈ 1 - 4·Var(c₅ | scores) / Var(H)  [dominant term]

  This is EXACT if α₂ = 0 always (true at n=5, false at n≥6).

  FOR A RIGOROUS FORMULA, WE NEED:
  (a) A formula for Var(c₅|scores) in terms of n
  (b) A formula for Var(H) in terms of n (the overlap counting)
  (c) The higher-order correction from Var(α₂|scores, α₁)

  WHAT THIS TELLS US ABOUT SCALING:
  - Var(H) grows as (n!)²/2^{2n} ∼ exponentially
  - Var(c₅|scores) grows polynomially (bounded by C(n,5)²)
  - So Var(c₅)/Var(H) → 0 as n → ∞ IF Var(H) grows fast enough

  BUT: the number of higher odd cycle counts (c₇, c₉, ...) also grows.
  The total Var(c₅ + c₇ + ...|scores) might grow as fast as Var(H).
  This is the OPEN QUESTION: does the sum of cycle residuals stay bounded
  relative to Var(H)?
""")

# Compute the contributions at n=5 and n=6
print("  Contribution breakdown:")
print(f"    n=5: Var(c₅|scores)/Var(H) = {float(total_var_c5_weighted/Fraction(285,16)):.6f}")
print(f"    n=5: 4×Var(c₅|scores)/Var(H) = {4*float(total_var_c5_weighted/Fraction(285,16)):.6f}")
print(f"    n=5: 1-OCR = {float(Fraction(4,133)):.6f}")
print(f"    Match: {'YES' if abs(4*total_var_c5_weighted/float(Fraction(285,16)) - float(Fraction(4,133))) < 0.001 else 'NO'}")
print()

# =============================================================================
# PART 9: WHAT WOULD A FORMULA LOOK LIKE?
# =============================================================================

print("=" * 72)
print("  PART 9: CANDIDATE FORMULA STRUCTURES")
print("=" * 72)
print()

print("""
  Based on the analysis, the OCR formula must have the structure:

    1 - OCR(n) = Var(residual cycle counts | scores) / Var(H)

  where the numerator counts the variance in cycle arrangements
  NOT explained by the score sequence.

  CANDIDATE 1 (simplest):
    1 - OCR(n) ≈ C(n,5) · c₅_var / Var(H)
  where c₅_var is the average within-class variance of c₅.
  This would give polynomial/exponential → 0 as n → ∞.

  CANDIDATE 2 (exact at n=5):
    1 - OCR(5) = 4/133 = 4/(7·19)
  The factor 4 = 2² comes from the H = 1 + 2α₁ coefficient squared.
  The factor 7 = forbidden prime = the SPLIT channel.
  The factor 19 = appears in Var(H) = 285/16 = 3·5·19/16.

  CANDIDATE 3 (structural):
    1 - OCR(n) = Σ_{ambiguous score classes} P(sc) · Var(c₅|sc) · 4 / Var(H)

  This is the LAW OF TOTAL VARIANCE applied to the OCF decomposition.
  It is not a "formula" in closed form, but it identifies the
  EXACT MECHANISM of the residual: the variance of 5-cycles within
  score classes, weighted by class probability.

  THE HONEST CONCLUSION:
  A closed-form OCR(n) formula does not yet exist.
  But the EXACT MECHANISM is identified:
    OCR residual = variance of higher odd cycle counts within score classes.
  At n=5, this is purely Var(c₅|scores).
  At n≥6, it includes c₇ terms and α₂ corrections.
""")

# =============================================================================
# PART 10: A PRACTICAL SPEEDUP — COMPUTING OCR WITHOUT FULL ENUMERATION
# =============================================================================

print("=" * 72)
print("  PART 10: SPEEDUP — OCR WITHOUT FULL ENUMERATION")
print("=" * 72)
print()

print("""
  Current: OCR(n) requires enumerating all 2^{C(n,2)} tournaments.
  At n=7: 2^21 ≈ 2M tournaments, 520 seconds.
  At n=8: 2^28 ≈ 268M tournaments, ~14 hours (infeasible in one session).

  SPEEDUP IDEA 1: Compute within each score class independently.
  Since OCR = 1 - Σ P(sc)·Var(H|sc) / Var(H), we can:
  (a) Enumerate tournaments BY SCORE CLASS
  (b) Within each class, compute only Var(H|sc)
  (c) Skip classes that are known to be unambiguous

  Many score classes have Var(H|sc) = 0 (H is unique).
  At n=5: 8/9 classes are unambiguous → only need to enumerate 280/1024 tournaments.
  At n=6: 13/22 classes unambiguous → only need ~78% of tournaments.

  SPEEDUP IDEA 2: Use the c₅ formula directly.
  Since Var(H|scores) ≈ 4·Var(c₅|scores), we can:
  (a) For each score class, count the DISTRIBUTION of c₅ values
  (b) Compute Var(c₅|sc) from the distribution
  (c) This avoids computing H entirely — c₅ is much cheaper

  c₅ computation: O(n!) per tournament (check all cyclic permutations)
  vs H computation: O(2^n · n²) per tournament (bitmask DP)
  At n=7: c₅ via O(5040) vs H via O(128·49) = O(6272). Similar cost.
  At n=8: c₅ via O(40320) vs H via O(256·64) = O(16384). c₅ slower!

  So for n≥8, the DP approach for H is actually FASTER per tournament.
  But the score-class decomposition still helps by skipping unambiguous classes.
""")

print()
print("=" * 72)
print("  SUMMARY OF FINDINGS")
print("=" * 72)
print()

print("""
  1. EXACT OCR FRACTIONS:
     OCR(5) = 129/133 (4/133 residual, factors: 2²/(7·19))
     OCR(6) = 460807/480480

  2. THE MECHANISM:
     The residual comes from Var(c₅ + c₇ + ...|scores).
     At n=5: purely Var(c₅|scores) = 15/(28·4) = 15/112.
     The factor 4 amplification: Var(H|scores) = 4·Var(c₅|scores).

  3. THE DECOMPOSITION THEOREM (exact at n=5):
     OCR(n) = 1 - 4·Var(α₁ - c₃|scores) / Var(H)
     where c₃ = C(n+1,3)/4 - S₂/2 is the score-determined part.

  4. THE OPEN QUESTION:
     Does Σ_{k≥5,odd} Var(c_k|scores) grow slower than Var(H)?
     If yes: OCR → 1. If no: OCR plateaus.
     The trend (1-OCR)×n = 0.15, 0.25, 0.29 is ambiguous.

  5. SPEEDUP PATH:
     Score-class decomposition allows skipping unambiguous classes.
     At n=7: 96.6% of tournaments are in ambiguous classes (not much savings).
     But the DOMINANT class contributes most of the variance — focus there.

  NO CLOSED-FORM FORMULA YET.
  But the exact mechanism (c₅ variance within score classes) is identified.
  This is the right target for a formal proof.
""")

print("opus-2026-03-21-S105: HUNTING FOR AN EXACT OCR FORMULA")
