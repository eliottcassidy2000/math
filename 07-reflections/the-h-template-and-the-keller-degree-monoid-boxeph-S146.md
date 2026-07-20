# The H-template: tournament spectral gaps and the Keller-degree monoid
**boxeph-2026-07-20-S146** · Owner brief: sharpen the JC₂ map; find tournament
analogies; think H ∈ {7, 21} impossibility.

## 1. The tournament laws (verified this session, exhaustive n ≤ 6 + 60k at n = 7)
h(T) = #Hamiltonian paths. Verified: (R) Rédei — h odd on every tournament tested;
(M) multiplicativity h(T) = Π h(strong components), 300/300; (C) the CONDENSATION-
MONOID LAW: the attainable set is the multiplicative monoid on strong-attainable
values — matched exactly at n ≤ 6. Data: attainable = {1,3,5,9,11,13,15,17,19,23,
25,...}; strong minima grow 3, 5, 9, 15 at n = 3..6. **The {7, 21} impossibility
mechanism**: 7 is prime and below the strong-growth curve at every n (min-strong
exceeds it permanently), so 7 is never strong-attainable and never a product;
21 = 3·7 needs a 7 or a strong 21 — also below the curve. Nuance: 35 appears at
n = 7 as a NEW STRONG value (not 5·7) — gaps close only from the strong side;
{7, 21} never do. Permanent spectral gaps = min-strong growth + multiplicativity.

## 2. The Keller-degree monoid (the analogy, made exact)
Keller ∘ Keller = Keller (dets multiply) and cover-degrees multiply: **the
attainable cover-degrees of Keller self-maps of ℂⁿ form a multiplicative monoid**;
"strong tournament" ↔ composition-irreducible Keller map; the condensation
(transitive order of strong components) ↔ the composition series. JC₂ says the
ℂ²-monoid is {1}. The H-template strategy: prove degree gaps for IRREDUCIBLES
(the strong analogue) and let the monoid propagate them.

## 3. The Euler ledger (dim 2, clean case) — the gap-proving instrument
THEOREM (clean case: A(F) smooth components, pairwise disjoint, no fiber jumps
over the bad finite sets). For a Keller counterexample F: ℂ²→ℂ² of cover degree d:
   Σᵢ (d − kᵢ)·χ(Cᵢ) = d − 1,
kᵢ = #persisting sheets over component Cᵢ, with kᵢ ≤ #Fix(local monodromy gᵢ)
(étale-section lemma, S145) and the ghost meridians generating the monodromy.
Proof: χ-multiplicativity of the covering V → U plus additivity of algebraic χ:
1 − χ(D) = d(1 − χ(C)), χ(D) = Σ kᵢχ(Cᵢ) in the clean case. ∎
CONSEQUENCES (new sharpenings of the JC₂ map):
- **d = 2 with all-smooth asymptotic components is IMPOSSIBLE**: ghosts have
  Fix = 0 (term 2χ = 2), silents have kᵢ ≤ 1 (term ≥ 1), ≥1 ghost required to
  generate ℤ/2 ⟹ LHS ≥ 2 > 1 = d − 1. Together with the S144 equivariant no-go
  and Wang, degree 2 is cornered from three sides — the first solid step toward
  the ODD-DEGREE CONJECTURE (Rédei analogue: Keller counterexample cover-degrees
  are odd; compositions preserve oddness, so it suffices for irreducibles).
- **ℤ/3-monodromy d = 3 needs a χ = 0 ghost**: 3χ_g ≤ 2 forces χ_g ≤ 0 with a
  persisting silent making up the ledger — a singular parametric ghost, highly
  constrained.
- **S₃-monodromy d = 3 minimal profile**: one transposition-ghost is not
  transitive, so ≥ 2 are needed; the ledger 2(χ₁ + χ₂) = 2 forces χ-profile
  (1, 0): ONE SMOOTH LINE plus ONE NODAL PARAMETRIC CURVE (or silent-laden
  variants). Any degree-3 JC₂ counterexample now has a nearly explicit shape.

## 4. The dictionary (established ↔ established, conjecture flagged)
h multiplicative over condensation ↔ deg multiplicative over composition ·
strong tournament ↔ irreducible Keller map · min-strong growth ⟹ {7,21} gaps ↔
ledger/Ghost exclusions ⟹ degree gaps · Rédei odd ↔ ODD-DEGREE CONJECTURE ·
tilings×|Aut| = H ↔ fibers×deck = d (Galois case — which the Ghost theorem
forbids in all-ghost form) · S₃ cycle types ↔ the two 3-vertex tournaments.

## 5. Honest scope + named next
Ledger proved in the clean case; the general-defect bookkeeping (fiber jumps,
sheet intersections) is the named residual alongside S145's silent-component
question. Next lever: a "min-irreducible-degree growth" theorem — any lower
bound on irreducible Keller degrees beyond 3 would create permanent monoid gaps
exactly as min-strong growth creates {7, 21}.
Files: h_spectrum_keller_monoid_boxeph_S146.py + .out (all tournament laws frozen).
