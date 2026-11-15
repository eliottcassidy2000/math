---
id: THM-5212
title: "The bottleneck-free growth bound in the floor-unit languages: the span-2 floor units compose by an exact 3-dimensional Hankel law whose bottleneck-free words are exactly the alternating ones (growth ((9 + sqrt 113)/2)^(1/4) = 1.7700005); in the span-3 language of the units 0.2, 1.2, 0.1.3, 0.3.1, 1.2.3 the bottleneck-free words are the paths of an explicit 9-state plane automaton, and an exact rational certificate gives H >= C * 1.7699^n for all of them"
status: >
  (A) PROVED (exact linear algebra). (B) PROVED, computer-assisted: exact rational verification of invariant sectors and
  edge inequalities. (C) PROVED from (A) and (B).
  OPEN: the full span-3 alphabet (13 floor units and non-unit letters) and arbitrary tournaments. For those the premium
  remains HYP-9742 with the periodic evidence of THM-5202 (D).
session: mac-mini-2026-11-15-growthbound (wave 29 of chanpak_waves_20261113)
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 29)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w29_unit_language_rep.py, w29_plane_automaton.py,
    w29_plane_cycles.py, w29_lp_reachable_sectors.py, w29_exact_growth_certificate.py, w29_receipts.txt
related:
  - THM-5202 (bottleneck factorization), HYP-9742 (NF2 fails at large primes), THM-4880 (ribbon transfer)
---

# THM-5212 - the bottleneck-free growth bound

Notation as in THM-5202.
- A *unit word* is the first vertex (letter 0) followed by a sequence of floor units.
- It is *bottleneck-free* if its transfer has rank ≥ 2 in the minimal representation of its unit language.
  Equivalently it contains no rank-one factor, since rank never increases under products.
- By THM-5202 (C), only bottleneck-free words can have prime counts outside bounded end effects.

## (A) Span 2: an exact composition law

Let A = 0.2 and B = 1.2, the two span-2 3-units. The language {A, B}* has a minimal representation of dimension 3. In
the Hankel coordinates

    z(W) = (H(W), H(WA), H(WB)),   with z(0) = (1, 3, 1),

the units act by the exact composition laws

    H(WAA) = 3 H(WA),
    H(WAB) = 4 H(WA) − 3 H(W),
    H(WBA) = (7 H(WB) + 4 H(WA) − 4 H(W))/3,
    H(WBB) = (5/3) (H(WA) + H(WB) − H(W)).

- *Rank structure.* M_A and M_B have characteristic polynomial x²(x − 3). The products AA and BB have rank one
  (characteristic polynomial x²(x − 9)), so they are bottlenecks. AB and BA have rank 2 and characteristic polynomial
  x(x² − 9x − 8).
- *Consequence.* The bottleneck-free unit words are exactly the alternating ones. Their counts satisfy the recurrence
  with characteristic polynomial x² − 9x − 8, so they grow exactly like ρ*^n with ρ* = ((9 + √113)/2)^(1/4).
- *Reading the laws.* A compatible repetition multiplies (H(WAA) = 3·H(WA)). An incompatible junction is affine in two
  earlier counts. That affine room is the only way a count can avoid factoring.

## (B) Span 3: a five-unit language and an exact certificate

The units 0.2 and 1.2 (3-units) and 0.1.3, 0.3.1, 1.2.3 (5-units) generate a language with minimal dimension 4.
- *Single units.* Each has characteristic polynomial x³(x − λ), λ ∈ {3, 5}.
- *Compatible pairs* have characteristic polynomial x³(x − λλ'): eigenvalues 9, 15, 25.
- *Incompatible pairs* have a second eigenvalue: x² − 9x − 8, x² − 17x − 16, x² − 19x − 2, x² − 25x + 10,
  x² − 47x + 22, x² − 33x − 32.

**The plane automaton.** Record the image plane of the current product. Every unit maps a plane either onto another
plane from a finite list or to rank one. The rank-≥2 words are exactly the paths of a 9-state automaton, with recurrent
part {2, …, 8}. Its 31 simple cycles have minimum growth 1.7700005 at (0.2, 1.2); the next are 1.7804 and 1.7805.

**Certificate (exact).** For each recurrent state s, with the plane written in a rational basis:
- *Invariant sectors.* There is a rational sector K_s. Its rays are snapped to the exact fixed directions, e.g.
  (−2, 1), where an edge map fixes a boundary. Every edge map T_u sends K_s into K_t, checked exactly on the extreme
  rays; there are 0 failures.
- *Expanding functionals.* There is a rational functional φ_s, positive on K_s, with

      φ_t(T_u x) ≥ (17699/10000)^|u| · φ_s(x)   for every edge s → t and every x in K_s

  (exact, on both extreme rays).
- *Reachability.* Every actual state vector of a word ending in s lies in K_s. The sectors were built from all
  vectors of paths of depth ≤ 8 and are invariant.
- *Positivity of the readout.* The count H is a positive functional on every K_s.

Hence every bottleneck-free word in this language satisfies H(W) ≥ C·1.7699^|W|. The words (0.2·1.2)^k reach ρ*, so
the bound is sharp up to the margin 1.7699 < ρ* = 1.7700005.

## (C) Consequence

In both languages a prime count satisfies H ≥ C·1.7699^n (≥ C·ρ*^n in span 2). That is exponentially above the
floor 5^(n/3) = 1.7099759^n. This is the bottleneck-free growth bound of HYP-9742, proved inside the floor-unit
languages that contain every floor word.

**Scope.** Not covered:
- the other eight span-3 floor units (a 13-unit computation did not finish);
- non-unit letters;
- non-ribbon tournaments.

The periodic evidence of THM-5202 (D) covers the full span-3 alphabet up to period 6.

**Composition view.** The floor units form a small "composition algebra" of tournament blocks:
- compatible pairs multiply exactly (nilpotent remainder, an exact coalescence);
- incompatible pairs carry one extra eigenvalue, which is the additive room primes need.

The concurrent packet chanpak_symmetry_interfaces_20261115 studies substitution interfaces. The two views are
complementary.
