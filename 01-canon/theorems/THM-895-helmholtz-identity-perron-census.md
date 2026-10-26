---
id: THM-895
title: Corrected unit-flow Helmholtz identity and Perron census — potential d/n, gradient energy x/n, residual energy C(n,2) − x/n
status: PROVED — corrected unit-flow Hodge identity and λ = 0 iff transitive; FINITE-EXACT — independent integral checks through n5; VERIFIED NUMERICALLY — floating-point Perron correlations and spreads in labelled n5/n6 censuses
source: mac-mini-2026-07-16-S115 (owner: Helmholtz + arborescences + Kendall–Wei self-compositions)
related:
  - THM-866 (the walk = gradient-energy ascent, now literal), THM-855 F3/F5 (per-flip and drift laws = gradient-side mechanics)
  - klein cont.9 HYP-7023/7024 (arborescences = determinant shadow of H — the curl side's tree census; complementary)
  - opus THM-894 (Kendall–Wei on the LRC resonance side)
script: 04-computation/helmholtz_perron_toothpick_macmini_S115.py -> 05-knowledge/results/helmholtz_perron_toothpick_macmini_S115.out
---

# THM-895 — the Helmholtz identity and the Perron census

**Correction, 2026-10-26.** The former potential d/(2n) was half the
least-squares solution for the stated unit flow. Its energy calculation
was internally consistent but never tested the normal equations. A single
arc K2 refutes that normalization. The correction and independent exact
controls are in the [information audit](../../05-knowledge/results/rank_geometry_metagraph_20261026/information/RESULT.md).
The numerical Perron census is not a machine-exact spectral proof.

**(i) Helmholtz/HodgeRank.** Put a_ij=1 if i beats j and -1 otherwise,
d_i=sum_j a_ij=2s_i-(n-1), and x=sum_i d_i². Assume n>=1 and sum_i φ_i=0.
The complete-graph Laplacian is nI-J, so the least-squares normal equations
are nφ=d. Thus **φ_i=d_i/n**. The residual a_ij-(φ_i-φ_j) has zero
divergence and is orthogonal to the gradient. Consequently

    gradient energy = sum_(i<j)(φ_i-φ_j)² = x/n,
    residual energy = C(n,2)-x/n.

The harmonic part is zero on the **filled clique complex**: all triangles
are included. The graph considered only as a one-dimensional complex has
a cycle space; completeness alone must not be used to remove that space.
For K2 the potential is (1/2,-1/2), gradient energy1 and residual0. ∎

So: **x is the gradient energy of the tournament flow** (up to n), the transitive
tournament is the max-gradient/min-curl point, the regular floor is max-curl, and the
tie-splitting +8 walk of THM-866 is gradient-energy ascent in steps of 8/n.
For transitive tournaments the residual energy is (n-1)(n-2)/6, generally
nonzero: an orientation induced by a potential is not the same as a unit
arc array equal to a potential difference. In general the residual energy
is (n-1)(n-2)/6 + 8c3(T)/n. This follows from the classical score formula
for the number c3 of cyclic triangles.

**(ii) Perron census (Kendall–Wei).** λ(T) is the spectral radius of the
adjacency matrix A, equivalently lim_(k->infinity) ||A^k||^(1/k) for any
matrix norm. It is a growth rate, not the generally divergent limit of
unnormalized iterated score vectors.
- **λ = 0 ⟺ T transitive** (λ = 0 ⟺ A nilpotent ⟺ acyclic ⟺ transitive); labeled count
  at λ = 0 equals n! exactly (n = 5, 6).
- corr(λ, −x) = **0.9313** (n = 5), **0.9218** (n = 6) — and identically equal to
  corr(λ, cycle energy), which is the mechanism: λ is a curl-side invariant.
- **λ splits x-levels**: at n = 6 all 7 interior levels carry positive λ-spread; at n = 5
  only levels x = 8, 16 split (the other levels are single-orbit-thin). The
  joint invariant (x,λ) therefore strictly refines x in these censuses.
  This does not say λ alone determines x: C3⊕C3 and C3⊕TT3 both have
  λ=1 but x=54 and62 respectively.

**Reading.** The Helmholtz split organizes the invariant zoo: gradient side = scores, x,
the walk, majorization (THM-866/868-869 territory — fully solvable/graded); curl side =
λ, odd cycles, H's higher digits, arborescence discord (klein cont.9) — where the
monodromy and the lawlessness live (THM-865/466). The 0.93 correlation with its exact
complement-form is the quantitative version of "the two sides see the same tournament
from opposite ends".
