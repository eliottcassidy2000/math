"""
recursion_correction.py — opus-2026-04-04-S14

THE EXACT FIBER BUNDLE RECURSION:

a(n) = a(n-1) * 2^{n-1} / n * (1 + epsilon(n))

where epsilon(n) → 0 exponentially fast.

What IS epsilon(n)? It comes from the non-identity permutations in Burnside.
The dominant correction is from the cycle type (1^{n-3}, 3^1):
- Number of such permutations: C(n,3) * 2 = n(n-1)(n-2)/3
- Contributing fix count: 2^{e(n,3)} where e(n,3) is the orbit count

Let's compute epsilon(n) exactly and see if there's a formula.
"""
from math import factorial, comb, log2
from fractions import Fraction

# Known values of a(n)
a = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
    15: 31021002160355166848
}

print("EXACT RECURSION CORRECTION epsilon(n)")
print("=" * 60)

print("\n  a(n) = a(n-1) * 2^{n-1} / n * (1 + epsilon(n))")
print(f"  {'n':>3s} {'a(n-1)*2^(n-1)/n':>25s} {'a(n)':>25s} {'epsilon':>15s} {'epsilon*2^n':>15s}")

for n in range(4, 16):
    pred = Fraction(a[n-1]) * Fraction(2**(n-1)) / Fraction(n)
    actual = Fraction(a[n])
    if pred != 0:
        epsilon = actual / pred - 1
        eps_scaled = float(epsilon) * 2**n
        print(f"  {n:3d} {float(pred):25.1f} {a[n]:25d} {float(epsilon):15.10f} {eps_scaled:15.6f}")

# The epsilon(n) sequence
print(f"\n  epsilon(n) * 2^n:")
eps_scaled_list = []
for n in range(4, 16):
    pred = Fraction(a[n-1]) * Fraction(2**(n-1)) / Fraction(n)
    actual = Fraction(a[n])
    epsilon = actual / pred - 1
    eps_2n = float(epsilon) * 2**n
    eps_scaled_list.append((n, eps_2n))
    print(f"    n={n}: epsilon*2^n = {eps_2n:.6f}")

# Can we predict epsilon from the structure?
# The Burnside formula: a(n) = (1/n!) Σ_{λ: all odd} (n!/z_λ) 2^{f(λ)}
# The dominant term: λ = (1^n): gives 2^{C(n,2)}
# The ratio: a(n) / (2^{C(n,2)}/n!) = 1 + correction

# The NEXT term: λ = (1^{n-3}, 3): z_λ = 3 * (n-3)!
# Its contribution: n! / (3 * (n-3)!) * 2^{f(1^{n-3},3)} / n!
# = 1/(3 * (n-3)!) * n! * 2^{f(1^{n-3},3)} / n! ... simplifying
# = C(n,3) * 2! / 3 * 2^{f} ...

# Actually for λ = (1^{n-3}, 3):
# z_λ = 3 * 1^{n-3} * (n-3)! = 3 * (n-3)!
# n! / z_λ = n! / (3 * (n-3)!) = n * (n-1) * (n-2) / 3
# f(λ) = ? (number of arc orbits)

# For a permutation with one 3-cycle (a,b,c) and n-3 fixed points:
# Arc orbits:
# - Arcs between fixed points: C(n-3, 2) arcs, each a self-orbit → C(n-3,2) orbits
# - Arcs between a fixed point and a cycle vertex: 3*(n-3) arcs.
#   The 3-cycle rotates cycle vertices, so arcs from fixed point f to cycle vertex
#   are in orbits of size 1 (since f is fixed, but the cycle vertex rotates).
#   Wait: arc (f, a) maps to (f, b) under the permutation. So arcs from f to the
#   3-cycle vertices form ONE orbit of size 3.
#   Similarly arcs TO f from cycle vertices form another orbit of size 3.
#   Total: (n-3) * 2 orbit groups, each of size 3, giving (n-3)*2/... hmm.
#
#   Actually: for each fixed point f, the arcs f→a, f→b, f→c form one orbit (size 3),
#   and arcs a→f, b→f, c→f form another orbit (size 3). But in a tournament,
#   for each ordered pair (f, cycle_vertex) there's one arc direction.
#   The orbit acts on ORDERED pairs, so (f,a), (f,b), (f,c) → one orbit, and
#   (a,f), (b,f), (c,f) → another orbit.
#   For tournaments: exactly one of f→a or a→f exists. The 3-cycle preserves the
#   tournament, so all three arcs f→a, f→b, f→c have the SAME direction (all toward
#   the cycle or all away). So for each f, 1 orbit of size 3 (either all f→cycle or all cycle→f).
#   But wait: the 3-cycle acts on tournament arcs by σ(u→v) = (σu→σv). So if f→a exists,
#   then f→b must exist too (since σ(f→a) = f→b). So all 3 arcs go the same way.
#   So: (n-3) arcs from f to cycle → 1 bit each (all same) → (n-3) orbits of size 1
#   No, each fixed point contributes 1 orbit between itself and the cycle (all 3 arcs
#   in the same direction). So (n-3) orbits from fixed-to-cycle connections.

# - Arcs within the 3-cycle: 3 directed arcs (a→b, b→c, c→a or reverse).
#   The 3-cycle permutes these cyclically: a→b maps to b→c maps to c→a.
#   So they form 1 orbit of size 3. In a tournament, either all go clockwise
#   or all counter-clockwise. 1 orbit.

# Total orbits: C(n-3,2) + (n-3) + 1
# f(1^{n-3}, 3) = C(n-3,2) + (n-3) + 1 = (n-3)(n-4)/2 + (n-3) + 1
#                = (n-3)(n-2)/2 + 1 = C(n-2,2) + 1

# So the contribution of λ = (1^{n-3}, 3):
# = [n(n-1)(n-2)/3] * 2^{C(n-2,2) + 1} / n!
# = [n(n-1)(n-2)/3] * 2^{C(n-2,2) + 1} / n!

# Relative to the main term 2^{C(n,2)} / n!:
# ratio = [n(n-1)(n-2)/3] * 2^{C(n-2,2)+1} / 2^{C(n,2)}
# = [n(n-1)(n-2)/3] * 2^{C(n-2,2)+1 - C(n,2)}
# = [n(n-1)(n-2)/3] * 2^{(n-2)(n-3)/2 + 1 - n(n-1)/2}
# = [n(n-1)(n-2)/3] * 2^{-(n-1) + 1}  [since C(n,2)-C(n-2,2) = n-1... wait]

# C(n,2) = n(n-1)/2, C(n-2,2) = (n-2)(n-3)/2
# C(n,2) - C(n-2,2) = [n(n-1) - (n-2)(n-3)]/2 = [n²-n-n²+5n-6]/2 = [4n-6]/2 = 2n-3
# So C(n-2,2) - C(n,2) = -(2n-3)

# ratio = [n(n-1)(n-2)/3] * 2^{-(2n-3)+1} = [n(n-1)(n-2)/3] * 2^{-2n+4}
# = n(n-1)(n-2) / (3 * 2^{2n-4})
# = n(n-1)(n-2) / (3 * 4^{n-2})

print(f"\n{'='*60}")
print(f"LEADING CORRECTION FROM 3-CYCLE PERMUTATIONS")
print(f"{'='*60}")

print(f"\n  The correction from cycle type (1^{{n-3}}, 3):")
print(f"  ratio = n(n-1)(n-2) / (3 * 4^{{n-2}})")

for n in range(4, 16):
    correction = n * (n-1) * (n-2) / (3 * 4**(n-2))
    pred_eps = correction  # leading correction to epsilon
    actual_pred = Fraction(a[n-1]) * Fraction(2**(n-1)) / Fraction(n)
    actual_eps = float(Fraction(a[n]) / actual_pred - 1)
    print(f"  n={n:3d}: predicted correction = {correction:.10f}, "
          f"actual epsilon = {actual_eps:.10f}, "
          f"ratio = {actual_eps/correction:.4f}" if correction > 1e-20 else "")

# Actually the epsilon in a(n)/[a(n-1)*2^{n-1}/n] is DIFFERENT from the
# correction to 2^{C(n,2)}/n!. Let me recompute.

# a(n) = (1/n!) [2^{C(n,2)} + C(n,3)*2*2^{C(n-2,2)+1} + ...]
# a(n-1) = (1/(n-1)!) [2^{C(n-1,2)} + ...]
# a(n)/a(n-1) = (n-1)!/n! * [2^{C(n,2)} + corrections] / [2^{C(n-1,2)} + corrections]
# = (1/n) * 2^{C(n,2)-C(n-1,2)} * [1 + corrections/2^{C(n,2)}] / [1 + corrections/2^{C(n-1,2)}]
# C(n,2) - C(n-1,2) = n-1
# So a(n)/a(n-1) ≈ 2^{n-1}/n * (1 + delta_n - delta_{n-1}) where delta_n = correction/main_term

# delta_n = C(n,3)*2*2^{C(n-2,2)+1} / 2^{C(n,2)} = n(n-1)(n-2)/3 * 2^{C(n-2,2)+1-C(n,2)}
# = n(n-1)(n-2)/3 * 2^{-(2n-4)} = n(n-1)(n-2) / (3 * 2^{2n-4})

# So epsilon(n) ≈ delta_n - delta_{n-1}
# = n(n-1)(n-2)/(3*4^{n-2}) - (n-1)(n-2)(n-3)/(3*4^{n-3})
# = (n-1)(n-2)/(3*4^{n-3}) * [n/4 - (n-3)]
# = (n-1)(n-2)/(3*4^{n-3}) * [n/4 - n + 3]
# = (n-1)(n-2)/(3*4^{n-3}) * [(-3n+12)/4]
# = (n-1)(n-2)(12-3n) / (12 * 4^{n-3})
# = -(n-1)(n-2)(n-4) / (4^{n-2})   for n > 4

print(f"\n  Refined: epsilon ≈ delta_n - delta_{{n-1}} = -(n-1)(n-2)(n-4) / 4^(n-2)")
for n in range(5, 16):
    if n in a and n-1 in a:
        predicted_eps = -(n-1)*(n-2)*(n-4) / 4**(n-2)
        actual_pred = Fraction(a[n-1]) * Fraction(2**(n-1)) / Fraction(n)
        actual_eps = float(Fraction(a[n]) / actual_pred - 1)
        ratio = actual_eps / predicted_eps if abs(predicted_eps) > 1e-30 else 0
        print(f"  n={n}: predicted = {predicted_eps:.12f}, actual = {actual_eps:.12f}, "
              f"ratio = {ratio:.4f}")
