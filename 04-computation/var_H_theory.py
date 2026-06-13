#!/usr/bin/env python3
"""
var_H_theory.py — kind-pasteur-2026-03-21-S16

THEORETICAL COMPUTATION OF Var(H) FOR RANDOM TOURNAMENTS.

E[H^2] = sum_{sigma, tau in S_n} P(both sigma and tau are HPs of random T)
       = sum_{sigma, tau} 2^{-|E(sigma) union E(tau)|}

where E(sigma) = set of directed arcs used by sigma as a Hamiltonian path.

|E(sigma)| = n-1 always. |E(sigma) union E(tau)| depends on the overlap.

For two permutations sigma, tau:
  arcs(sigma) = {(sigma(0)->sigma(1)), (sigma(1)->sigma(2)), ...}
  arcs(tau) = {(tau(0)->tau(1)), (tau(1)->tau(2)), ...}

  Each arc is an ORDERED pair (u, v).
  Two arcs can agree: (u,v) in both (same direction).
  Two arcs can conflict: (u,v) in one and (v,u) in other.
  Or be unrelated.

  |union| = |arcs(sigma)| + |arcs(tau)| - |arcs(sigma) ∩ arcs(tau)|
  = 2(n-1) - |agree|

  where |agree| = #{arcs shared (same direction)}.

  But WAIT: in a tournament, arc (u,v) means u->v. If sigma uses (u,v) and
  tau uses (v,u), these are DIFFERENT arcs but involve the SAME PAIR {u,v}.
  P(both present) = 0 in a tournament! So if sigma uses (u,v) and tau uses (v,u),
  P(both HPs) = 0.

  Actually: P(sigma is HP) = prod_{k} P(sigma(k) -> sigma(k+1)) = (1/2)^{n-1}
  assuming all arcs are independent (which they are for a random tournament).

  P(sigma AND tau both HPs) = prod_{pairs {u,v}} P(arc direction consistent with both)

  For each unordered pair {u,v}:
  - If neither sigma nor tau uses {u,v}: P = 1 (any direction ok)
  - If only sigma uses (u,v): P = 1/2
  - If only tau uses (u,v): P = 1/2
  - If both use same direction (u,v): P = 1/2 (same constraint)
  - If sigma uses (u,v) and tau uses (v,u): P = 0! (impossible)

  So P(both) = (1/2)^{# of distinct unordered pairs used by sigma or tau}
             * I(no conflicting arcs)

  where "conflicting" means sigma uses (u,v) and tau uses (v,u) for same {u,v}.

  E[H^2] = sum_{sigma, tau: no conflict} 2^{-|pairs(sigma) union pairs(tau)|}

Let me compute this for small n.
"""
from itertools import permutations
from fractions import Fraction
from math import factorial, gcd

for n in range(3, 8):
    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")

    perms = list(permutations(range(n)))
    n_perms = len(perms)

    E_H2 = Fraction(0)
    n_compatible = 0
    n_total_pairs = 0

    # Pre-compute arcs for each permutation
    perm_arcs = []
    perm_pairs = []
    for sigma in perms:
        arcs = set()
        pairs = set()
        for k in range(n-1):
            arcs.add((sigma[k], sigma[k+1]))
            pairs.add(frozenset({sigma[k], sigma[k+1]}))
        perm_arcs.append(arcs)
        perm_pairs.append(pairs)

    for i in range(n_perms):
        for j in range(n_perms):
            n_total_pairs += 1
            # Check for conflicting arcs
            conflict = False
            for u, v in perm_arcs[i]:
                if (v, u) in perm_arcs[j]:
                    conflict = True
                    break

            if not conflict:
                n_compatible += 1
                # Count distinct unordered pairs used
                union_pairs = perm_pairs[i] | perm_pairs[j]
                n_pairs = len(union_pairs)
                E_H2 += Fraction(1, 2**n_pairs)

    E_H = Fraction(factorial(n), 2**(n-1))
    Var_H = E_H2 - E_H**2

    print(f"  E[H] = {E_H} = {float(E_H):.6f}")
    print(f"  E[H^2] = {E_H2} = {float(E_H2):.6f}")
    print(f"  Var(H) = {Var_H} = {float(Var_H):.6f}")
    print(f"  Compatible pairs: {n_compatible} / {n_total_pairs} = {n_compatible/n_total_pairs:.4f}")

    # Compare with exact computation
    # At n=5: Var(H) should be 285/16
    # At n=6: Var(H) should be 585/4
    known = {3: Fraction(3, 4), 4: Fraction(3), 5: Fraction(285, 16), 6: Fraction(585, 4)}

    if n in known:
        print(f"  KNOWN Var(H) = {known[n]} = {float(known[n]):.6f}")
        if Var_H == known[n]:
            print(f"  MATCH ✓")
        else:
            print(f"  MISMATCH ✗ (difference = {Var_H - known[n]})")

    # Factor Var(H)
    p, q = Var_H.numerator, Var_H.denominator
    print(f"  Var(H) = {p}/{q}")

    # Check OCR
    # Var(c3) from Rao: c3 = C(n,3) - sum C(si,2)
    # For random tournament, E[si] = (n-1)/2, Var(si) = (n-1)/4
    # c3 = C(n,3) - sum si*(si-1)/2
    # Var(c3) from exact computation (previously validated)
    # For the OCR formula: OCR = 4*Var(c3) / Var(H) assuming Cov(eps,c3)=0
    # But we showed Cov(eps,c3) != 0 at n>=5.
    # Still, the exact R^2 = Cov(H,c3)^2 / (Var(c3)*Var(H))

# Factor analysis
print("\n\nVAR(H) FACTORIZATION:")
known_var = {
    3: Fraction(3, 4),        # Wait, earlier said 5/4? Let me recheck.
    4: Fraction(3, 1),
    5: Fraction(285, 16),
    6: Fraction(585, 4),
    7: Fraction(206325, 128),
}

for n, v in known_var.items():
    p, q = v.numerator, v.denominator
    # Factor p
    factors = []
    pp = p
    for f in range(2, 1000):
        while pp % f == 0:
            factors.append(f)
            pp //= f
    if pp > 1:
        factors.append(pp)
    print(f"  n={n}: Var(H) = {p}/{q}, numerator factors: {factors}")

# Can we find the formula for Var(H)?
# Var(H) values: 3/4, 3, 285/16, 585/4, 206325/128
# Numerators: 3, 3, 285, 585, 206325
# 285 = 3*5*19
# 585 = 5*9*13
# 206325 = 3*5^2*2751 = 3*5^2*3*917 = 9*25*917 = 225*917
# Actually: 206325 / 3 = 68775, /5 = 13755, /5 = 2751, /3 = 917
# 917 = 7*131
# So 206325 = 3^2 * 5^2 * 7 * 131
#
# 285 = 3 * 5 * 19
# 585 = 3^2 * 5 * 13
# 206325 = 3^2 * 5^2 * 7 * 131
#
# Pattern: the PRIME factor appears as:
# n=5: 19 (the OCR denominator!)
# n=6: 13 (the OCR denominator!)
# n=7: 131 (the OCR denominator!)
#
# So the OCR denominator is the largest prime factor of the Var(H) numerator!

print("\n\n*** KEY PATTERN ***")
print("The OCR denominator appears as the largest prime factor of Var(H) numerator!")
print("  n=5: Var(H) = 285/16 = (3*5*19)/16,     OCR denom = 19 ✓")
print("  n=6: Var(H) = 585/4 = (3^2*5*13)/4,       OCR denom = 13 ✓")
print("  n=7: Var(H) = 206325/128 = (3^2*5^2*7*131)/128, OCR denom = 131 ✓")
print()
print("If this pattern holds, OCR(n) = 1 - k/(largest prime factor of Var(H)·2^m)")
print("The denominator IS the largest prime factor of the numerator of Var(H).")
