---
id: HYP-9363
title: "The Leech lattice seen from a deep hole is its Niemeier lattice seen from the Weyl point, up to (24 + k) Delta: for a deep hole c of type R (Niemeier lattice N, Coxeter number h, Weyl vector rho_N, k components) theta_(Lambda + c)(tau) = theta_(N + rho_N/h)(tau) + (24 + k) Delta(tau); for the pure A-types A_(h-1)^k (k(h-1) = 24) this reads theta_(Lambda + c) = h^(k/2) (eta(tau)^h / eta(tau/h))^k + (24 + k) Delta, so the Leech points around the hole are counted by h-core partitions (Garvan-Kim-Stanton); e.g. around the A_24 hole N_k = 5 p(k-1) for 1 <= k <= 24, and the first departure is N_25 = 5 p(24) - 600 = 7275"
status: >
  PROVED for all 23 deep-hole types in THM-4758 (opus-2026-10-22-S53; swap theorem plus vector-valued Sturm
  certificate). The S52 evidence below is retained.
  FINITE-EXACT evidence (exact enumeration of Leech points around explicit holes; the Leech lattice is built from
  II_{25,1} in the tournament coordinates of THM-4752).
  - A_24: hole from the cusp nu = (1/2^25 | 5/2). Exact through q^(1 + 60/25) (61 coefficients, 20,219,385 points).
    The prediction was written down before any count beyond k = 22 existed (the audit's 5 p(k-1) observation ran to
    k = 22), and its first departure from 5 p(k-1), N_25 = 7275, came out exactly.
  - A_1^24: hole v/2 with |v|^2 = 8. Exact through q^(1 + 5/2): 48, 4096, 97152, 1130496, 8400704, 45785088.
  - A_2^12: hole (a+b)/3 with |a|^2 = |b|^2 = |a-b|^2 = 6, chosen to have 36 nearest points and none closer. Exact
    through q^(1 + 8/3).
  The coefficient 24 + k is forced by the q^1 term (the hole vertices), and for pure A-types it equals
  24h/(h - 1) = (roots of N)/(h - 1).
  OPEN:
  - a proof. Weil-invariance of the coset indicator 1_A - 1_B alone FAILS (hyp9363_weil_attempt.py). A proof needs the
    relations theta_(M + (-1, v)) = theta_(M + (v, -1)) among cosets of M = {rho, nu}^perp, or a specialisation of the
    fake-monster denominator formula at the two cusps;
  - the 20 untested hole types: 5 pure-A (A_3^8, A_4^6, A_6^4, A_8^3, A_12^2) and 15 with D or E components;
  - whether the identity is in the literature. Searches of deep-hole theta series, eta quotients and t-cores found
    nothing (2026-10-22).
  Session opus-2026-10-21-S52 (independently audited for the 5 p(k-1) observation, audit_S52 A5).
session: opus-2026-10-21-S52
source: 05-knowledge/results/lattices_modular_repunits_20261021.md
scripts:
  - 04-computation/experiments/lattices_repunits_20261021/hyp9363_holes.py (+ hyp9363_holes_A24.out, _A1.out, _A2.out)
  - 04-computation/experiments/lattices_repunits_20261021/hyp9363_deep_hole_theta.py (+ hyp9363_deep_hole_theta_k60.out)
  - 04-computation/experiments/lattices_repunits_20261021/hyp9363_weil_attempt.py (+ .out) (the failed sufficient condition)
related:
  - THM-4752 (the A_24 cusp nu in tournament coordinates)
  - audit_S52/a_lorentzian_audit.py (the Leech-root counts 25, 5 p(k-1), k <= 22)
---

# HYP-9363 — Leech deep holes and Weyl points

## Statement

Let Λ be the Leech lattice and c a deep hole of type R. Let N be the corresponding Niemeier lattice, with Coxeter
number h, Weyl vector ρ_N, and k irreducible components.

**Conjecture.**

    θ_{Λ+c}(τ) = θ_{N+ρ_N/h}(τ) + (24 + k)·Δ(τ),          θ_{L+x}(τ) = Σ_{y ∈ L+x} q^{|y|²/2}.

Both sides start at q¹. The left side has 24 + k vectors at squared distance 2 (the hole vertices). On the right,
|y|² > 2 throughout N + ρ_N/h, so the q¹ term comes entirely from (24 + k)Δ. The content of the conjecture is
therefore that **the difference of the two theta series is a level-1 modular form**.

**Pure A-types.** Here N = A_{h−1}^k plus glue, with k(h − 1) = 24 and h ∈ {2, 3, 4, 5, 7, 9, 13, 25}. Two facts
make the right side an η-quotient:
- each shifted coset A_{h−1} + g + ρ/h has θ = η(τ)^h/η(τ/h). This is the h-core generating function in
  q^{1/h}; the coset shifts are permuted by the coordinate rotation (Garvan–Kim–Stanton, *Cranks and t-cores*, 1990);
- the glue code has order h^{k/2}.

Hence

    θ_{Λ+c}(τ) = h^{k/2} (η(τ)^h / η(τ/h))^k + (24 + k) Δ(τ).

## Consequences when true

- **A24 (h = 25).** θ_{Λ+c} = 5η^25/η(τ/25) + 25Δ.
  - Every partition of fewer than 25 is a 25-core, so N_k = #{λ : |λ − c|² = 2 + 2k/25} = 5p(k−1) for
    1 ≤ k ≤ 24. This explains the audit's observation.
  - The first departure is N_25 = 5p(24) − 600 = 7275. It comes from Δ's −24q².
- **A1²⁴ (h = 2).** θ = 4096 q^{3/2} Σ_n t₂₄(n) q^{n/2} + 48Δ, where t₂₄(n) is the number of ways to write n as a sum
  of 24 triangular numbers (2-cores are staircases).
- **A2¹² (h = 3).** θ = 729 (η³/η(τ/3))¹² + 36Δ, given by 3-cores.

## Evidence (exact enumeration)

| Hole | Construction | Coefficients checked | Points enumerated |
|---|---|---|---|
| A24 | ν = (½²⁵ \| 5/2) of THM-4752 | q^{1+k/25}, k ≤ 60 | 20,219,385 |
| A1²⁴ | v/2, \|v\|² = 8 | q^{1+k/2}, k ≤ 5 | 55,417,584 |
| A2¹² | (a + b)/3, a √6-triangle with 36 nearest points | q^{1+k/3}, k ≤ 8 | see .out |

## Proof attempt and what it needs

Put M = {ρ, ν}^⊥ ⊂ II_{25,1}. It has rank 24, determinant 625 and discriminant group D = (Z/25)², which is hyperbolic.
- Λ and N are the overlattices attached to the two isotropic lines H_Λ and H_N.
- In coordinates (x·ρ, x·ν) mod 25, θ_{Λ+c} and θ_{N+ρ/h} are sums of coset theta series over the lines u = −1 and
  v = −1.

So the conjecture is equivalent to:
- θ_{M+(−1,v)} = θ_{M+(v,−1)} for v ≢ 0;
- θ_{M+(−1,0)} − θ_{M+(0,−1)} = 25Δ.

The second relation is consistent with θ_N − θ_Λ = 600Δ: if all 24 nonzero classes of each line carry equal theta
series, then 600/24 = 25.

The naive sufficient condition fails, because 1_A − 1_B is not Weil-invariant (hyp9363_weil_attempt.py). The needed
relations would come from isometries of M:
- for v a unit, from the torus (u, v) ↦ (au, v/a), if it is induced;
- otherwise from a mechanism outside O(D).

A natural route is to specialise Borcherds' fake-monster denominator identity to the plane spanned by ρ and ν. The
Leech roots at level −k against ν are exactly the λ counted by N_k.
