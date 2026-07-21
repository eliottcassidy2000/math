---
id: THM-1775
title: "THE MOMENT-NULLCONE TEMPLATE — three pillars of this project (tournament transitivity, TNC, GMC(2)) are ONE structure at three levels of a rational ⊂ algebraic ⊂ holonomic complexity hierarchy, and the tournament case is Cayley–Hamilton. The shared template: an invariant φ(Xᵐ) is the projection of the m-th power onto the trivial component of a symmetry (a 'moment'); its generating function F(t) = Σₘ φ(Xᵐ)tᵐ obeys a finite-order linear recurrence; the NULLCONE {φ(Xᵐ)=0 ∀m} is detected at finite DEPTH = recurrence order, and equals the locus where F collapses to its trivial value. INSTANCES, each verified: (T) TOURNAMENT — φ = trace, X = adjacency A, F(t) = Σ tr(Aᵐ)tᵐ = −t·d/dt·log det(I−tA) is RATIONAL (poles 1/λᵢ), nullcone = {tr(Aᵐ)=0 ∀m} = nilpotent = TRANSITIVE, detection depth = n, recurrence = CAYLEY–HAMILTON; verified transitive ⟺ tr(Aᵏ)=0 for k≤n over all tournaments n=3..7, and the char-poly recurrence tr_k = −Σcᵢtr_{k−i} for all n. (Λ) TNC — φ = constant term, X = Laurent Λ, F ALGEBRAIC of degree D, nullcone = one-sided (F≡1), depth D (THM-1710). (G) GMC(2) — φ = Gaussian E, X = P(Z,Z̄), F HOLONOMIC (polar/Laplace bridge THM-1645), nullcone = charge-one-sided, depth K (THM-1740). THE UNIFICATION IS EXACT, not analogical: Cayley–Hamilton is precisely the RATIONAL (finite-matrix) case of the holonomic recurrence, so the tournament trace-nullcone and the GMC moment-nullcone are the SAME theorem at the two ends of the {rational, algebraic, holonomic} spectrum; the detection depth is the arithmetic complexity of F. And the nullcone is always the MOST DEGENERATE object — transitive tournament, one-sided support, charge-one-sided P — the collapse of F to its trivial value"
status: >
  TEMPLATE: a structural unification.  Its content is the single fact that all three
  functionals are trivial-isotypic projections whose power-generating-functions are
  recurrence-governed, plus the exact observation that Cayley–Hamilton is the rational case
  of the holonomic recurrence.  (T) is PROVED here and machine-verified (transitive ⟺
  depth-n trace vanishing, all tournaments n=3..7; Cayley–Hamilton recurrence, n=5,6,7).
  (Λ) and (G) are my prior theorems THM-1710 / THM-1740, cited.  The hierarchy placement
  (rational/algebraic/holonomic) is exact for (T) and (Λ) and is the established class of the
  GMC generating function for (G).
  What is NOT claimed: a new proof of TNC or GMC.  This is a frame that makes precise WHY the
  detection-depth machinery of the last several sessions is not three coincidences but one
  phenomenon, and it lands the recent GMC/TNC work back onto the project's tournament core.
renumbered: "claimed THM-1750; renumbered to THM-1775 by first-pusher rule — the arborescence-ranking THM-1750 was pushed 8s earlier (18:06:24 vs 18:06:32)."
source: kind-pasteur-2026-07-20-S128c128 (owner: look for more creative unifying frames; think abstract; explore past concepts for connections)
depends_on:
  - THM-1710    # TNC detection depth D (the algebraic level)
  - THM-1740    # GMC(2) bounded finite Gröbner (the holonomic level)
related: [THM-895, THM-1555, THM-1670, THM-1720, THM-133]
reflection: 07-reflections/the-moment-nullcone-template-and-the-complexity-hierarchy.md
script: 04-computation/moment_nullcone_template_kps_S128c128.py (+ .out)
---

# THM-1775 — the moment-nullcone template

Three things this project keeps proving separately are one thing.

## The template

Let a symmetry `G` act, and let `φ(Y) = ⟨Y⟩_G` be the projection onto the trivial component.
For an object `X`, form the **moment sequence** `φ(X^m)` and its generating function
`F(t) = Σ_{m≥0} φ(X^m) t^m`. Then:

1. `F(t)` obeys a **finite-order linear recurrence** (its arithmetic class — rational,
   algebraic, or holonomic — bounds the order);
2. the **nullcone** `𝒩 = {X : φ(X^m) = 0 ∀ m ≥ 1}` is detected at **finite depth** = that
   order: `φ(X^m) = 0` for `m ≤ depth` already forces it for all `m`;
3. `𝒩` is exactly the locus where `F` **collapses to its trivial value**.

## The three instances

| | symmetry `G` | `φ(X^m)` | `X` | nullcone `𝒩` | `F(t)` class | depth | recurrence |
|---|---|---|---|---|---|---|---|
| **(T) tournament** | trace | `tr(A^m)` | adjacency `A` | transitive (nilpotent) | **rational** | `n` | Cayley–Hamilton |
| **(Λ) TNC** | circle `U(1)` | `CT(Λ^m)` | Laurent `Λ` | one-sided | **algebraic** | `D` | THM-1710 |
| **(G) GMC(2)** | `U(1)`×radial | `E[P^m]` | `P(Z,Z̄)` | charge-one-sided | **holonomic** | `K` | THM-1740 |

## (T) verified

`F(t) = Σ_m tr(A^m) t^m = −t·(d/dt) log det(I − tA)` is **rational** with poles `1/λ_i`, so
`tr(A^m)` obeys the **Cayley–Hamilton recurrence** `tr_k = −Σ_{i=1}^n c_i tr_{k−i}` (`c_i` the
char-poly coefficients). Hence:

> **`A` is transitive `⟺` `tr(A^m) = 0` for `m = 1,…,n`.**

*Proof.* `tr(A^m)=0` for `m=1..n` `⟺` (Newton) all power sums vanish `⟺` char poly `= x^n`
`⟺` `A` nilpotent `⟺` `A` acyclic `⟺` `A` transitive. ∎

Machine-verified: the depth-`n` implication and `transitive ⟺ nullcone` over **all**
tournaments `n = 3..7`; the Cayley–Hamilton recurrence for `n = 5,6,7`.

This is `THM-895`'s `λ = 0 ⟺ transitive` read through the template: the transitive tournament
is the **trace nullcone**, the tournament analogue of a one-sided Laurent polynomial and a
charge-one-sided Gaussian.

## Why the unification is exact, not analogical

**Cayley–Hamilton is the rational (finite-dimensional) case of the holonomic recurrence.** A
matrix is a rank-1-summable "period": `Σ tr(A^m) t^m = Σ_i λ_i t/(1−λ_i t)`, rational of degree
`n`. A diagonal of a rational function (TNC's `CT(Λ^m)`) is algebraic; a period integral
(GMC's `E[P^m]`, via the polar/Laplace bridge THM-1645) is holonomic. Same recurrence
phenomenon, three rungs of

> **rational ⊂ algebraic ⊂ holonomic (D-finite)**,

and **the detection depth is the arithmetic complexity of `F`**. The tournament trace-nullcone
(THM-895) and the GMC moment-nullcone (THM-1740) are literally the two ends of one theorem.

## The nullcone is always the most degenerate object

- tournament: transitive `⟺ det(I − tA) = 1 ⟺ F ≡ 0` (all `λ_i = 0`);
- TNC: one-sided `⟺ F ≡ 1` (small-branch sum trivial, THM-1670);
- GMC: charge-one-sided `⟺` `E`-series trivial past `m=0` (charge threshold).

In each, `𝒩` is where the generating function loses all its structure and becomes constant.
"Maximally ordered / one-sided / extremal" is one notion across the three.

## Resonances with the wider project (in the reflection)

The roots-of-unity thread runs through all three (regular-tournament spectra on the line
`Re = −½` are Gauss sums for the circulant/Paley case, THM-1555; TNC's singular indices are
roots-of-unity exponents, THM-1720/1725; GMC's `φ` is a circle average). Burnside's
`A000568 = (1/n!)Σ_σ Fix(σ)` is a *second* tournament projection (onto `S_n`-invariants), so
the tournament count and the tournament trace are two moments of the same object. And the
"extremal = nullcone" reading suggests LRC's tight AP is the combinatorial analogue of a
nullcone member. These are developed, and kept honest as suggestive, in
`07-reflections/the-moment-nullcone-template-and-the-complexity-hierarchy.md`.

## Named next

- **Place `H` itself (not just `tr(A^m)`) in the template.** `H` = Hamiltonian-path count
  (OCF `= I(Ω,2)`); is it a moment of a projection whose GF is rational/algebraic, with its
  own detection depth? THM-133's spectral-OCF chain is the lead.
- **The LRC rung.** LRC's "for all `t`" is an *extremal* (min-max), not a moment sum — but its
  finite-certificate structure (BV/comb bounds) is a detection depth of a different flavor.
  Ask whether the tight-AP locus is a nullcone in a moment functional built from the
  resonance matrix (THM-894's Kendall–Wei on speeds).
- **A single Lean interface** `MomentNullcone` with `zeros_propagate` (already proved,
  TNCDetectionDepth.lean) as the shared engine and the three instances as specializations.
