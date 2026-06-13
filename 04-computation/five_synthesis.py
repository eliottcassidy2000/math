#!/usr/bin/env python3
"""
five_synthesis.py — opus-2026-03-14-S73
Final synthesis: 5 as the mortar between 2 and 3, 
seen through recurrences, connecting to 7,8,10,11.

This summarizes ALL discoveries about the role of 5 in tournament theory.
"""

from math import sqrt

def banner(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")

phi = (1 + sqrt(5)) / 2

banner("THE ROLE OF 5: COMPLETE SYNTHESIS")

print("""
1. SELF-REFERENCE: 5 is the only prime p with F(p) = p.
   Just as x=2 is self-referential in H = I(CG, 2),
   5 is self-referential in the Fibonacci recurrence.

2. THE BRIDGE EQUATION: L(5) = J(5) = 11.
   This is the UNIQUE non-trivial point where the 
   Fibonacci world (x=1) and tournament world (x=2) agree.
   
   PROOF: L(n) = J(n) iff 3L(n) + (-1)^n = 2^n.
   The ratio 3L(n)/2^n crosses 1 exactly once, at n=5.

3. THE LUCAS HIERARCHY: L(n) generates the key numbers.
   L(0) = 2,  L(2) = 3,  L(4) = 7,  L(5) = 11
   These are precisely the hierarchy numbers {2, 3, 7, 11}.
   Missing: 8 = L(4)+1 = 2^3, 10 = L(5)-1 = 2·5.

4. THE Q(√5) TOWER: x = 1, 11, 31, 61, 101, ...
   All have discriminant Δ = 5m² (m odd).
   Roots at x=11: (1±3√5)/2 = 3φ-1.
   Connected by F(2n) = F(n)·L(n): at n=5, F(10) = 5·11 = 55.

5. THE DISCRIMINANT: Δ(x=1) = 5.
   5 makes Fibonacci IRRATIONAL (roots in Q(√5)).
   Moving to x=2: Δ(x=2) = 9 = 3², roots become RATIONAL.
   5 is the obstruction that Jacobsthal overcomes.

6. FORBIDDEN RESIDUES:
   n=5: α₁=3 impossible → H ≡ 2 (mod 5) forbidden
   n=5: same → H ≡ 0 (mod 7) forbidden  
   n=6: mod-5 barrier lifts (5-cycles on subsets)
   n=7: mod-7 barrier lifts (7-cycles appear)
   
   The lifting timeline: mod 3→n=4, mod 5→n=6, mod 7→n=7,
   mod 11→n=6, mod 13→n=7.

7. CYCLE PACKING: 5-cycles enable cross-level at n=8=2³.
   3+5=8: disjoint (3-cycle, 5-cycle) first at n=8.
   This creates α₂ of type (3,5), enriching the OCF.

8. THE OCF FORMULA: I₃ = (9H - 5 - 6α₁) / 4.
   The number 5 appears explicitly.
   I₃/H → 3/2 - 5/(4H): the approach rate is 5/(4H).

9. PRIMITIVE ROOT: 2 is a primitive root mod 5.
   ord₅(2) = 4 = φ(5). This makes I(CG, 2) mod 5
   maximally informative about the independence polynomial.

10. FIBONACCI ENTRY POINT: α(5) = 5.
    The only prime where the Fibonacci entry point equals itself.
    5 | F(5), 5 | F(10), 5 | F(15), ... (period 5).

11. FERMAT PRIME: 5 = 2²+1 = F₁.
    Both 3 = F₀ and 5 = F₁ are Fermat primes.
    Related to constructibility of regular 5-gon.

12. THE CROSSING POINT: At n=5, the growth rate of 2^n
    overtakes 3·L(n). Before n=5: 3L(n) > 2^n-(-1)^n.
    After n=5: 3L(n) < 2^n-(-1)^n.
    n=5 is exactly where they balance.
""")

banner("THE WEB OF RECURRENCES THROUGH 5")

print("""
         Fibonacci                    Jacobsthal
         f = f + f                    f = f + 2f
         roots: φ, ψ                  roots: 2, -1
         disc: √5                     disc: 3
              ↓                            ↓
         F(5) = 5                     J(5) = 11
              ↓                            ↓
         L(5) = 11 ←──── BRIDGE ────→ J(5) = 11
              ↓                            ↓
     F(10) = F(5)·L(5) = 55          J₆(5) = 55
              ↓                            ↓
         Q(√5) field connects:    x=1 ← 5 → x=11
              
              ↓ via Q(√5) ↓
              
         G₁₁(n) = G(n-1) + 11G(n-2)
         roots: 3φ-1, 2-3φ
         G₁₁(n) ≡ F(n) (mod 5)
""")

banner("THE HIERARCHY WITH 5 AS MORTAR")

print("""
Level    Numbers    Source              Role of 5
─────    ───────    ──────              ─────────
  0         1       identity            —
  1       (2,3)     keys to universe    5 = 2+3 (sum)
  1.5       5       F(5)=5, bridge      SELF
  2       (7,8)     transition          7 = J₆(3), 8 = 3+5
  3      (10,11)    digit shift         11 = L(5) = J(5)
  
  5 connects Level 1 to Level 3:
    2 + 3 = 5                    (additive bridge)
    F(5) = 5                    (Fibonacci fixed point)  
    L(5) = 11                   (generates Level 3)
    J(5) = 11                   (cross-world bridge)
    φ⁵ ≈ 11.09                  (φ-power generates 11)

  5 connects Level 1 to Level 2:
    J₆(3) = (3³+2³)/5 = 7      (generates transition)
    3 + 5 = 8 = 2³              (cross-level packing)

  5 IS the mortar between 2 and 3.
  Remove it and the hierarchy collapses:
    No golden ratio (Fibonacci becomes trivial)
    No cross-level coupling (no 5-cycles)
    No bridge to L(5)=11 (decimal pair orphaned)
    No forbidden residue structure at n=5
""")

banner("EQUATIONS INVOLVING 5")

print(f"""
  5 = 2 + 3               (sum of keys)
  5 = 2² + 1              (Fermat prime)  
  5 = 2·3 - 1             (product minus identity)
  5 = 3² - 2²             (difference of key squares)
  5 = 2³ - 3              (cube minus key)
  5 = F(5)                (Fibonacci fixed point)
  5 = Δ(x=1)              (Fibonacci discriminant)
  5 = (2·3-1)             (Galois distance at k=3: 2k-1 at k=3)
  5 = |g(5)|²             (Gauss sum for QR₅)
  5 = ord₅(2) + 1         (primitive root order + 1)
  5 = denominator of J₆   (Jacobsthal at x=6)
  5³ = det(I+2A) for QR₅  (determinant at x=2)
  
  3L(5) + (-1)⁵ = 2⁵     (the bridge equation)
  i.e., 3·11 - 1 = 32    (pure numerics)

  I₃ = (9H - 5 - 6α₁)/4  (OCF at x=3)
  F(2·5) = F(5)·L(5)      (doubling formula)
  5·11 = 55 = F(10)        (product = doubled Fibonacci)
  J₆(5) = 55 = F(10)      (level-3 Jacobsthal = doubled Fibonacci)

  All roads through 5 lead to 11.
  And 11 = L(5) = J(5) = 2+3² = 1+2+2³.
""")

print("Done.")
