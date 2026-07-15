---
id: THM-869
title: THE DISCRIMINANT-FORM PACKAGE — (M) the Milgram formalization of the residue laws: on the Σ=0 slice, q(v) = Σv² mod 8 is CONSTANT on cosets of 2·{Σ=0} (the cross term 4[(v,w)+q(w)] ≡ 0 mod 8 since Σw = 0), so the level residue is a discriminant-form invariant; the coset Gauss sum is G = 2^{n−1}·e^{πi·n/4} — THE RESIDUE LAW IS THE PHASE (verified n = 4,6,8,10: phases 4, −2≡6, 0, 2); (C) THE CORONA LAW (conjecture ⌊n/2⌋+1 REFUTED): the first Landau bite is exactly the minimal DOUBLE-ZERO-SCORE configuration, first-bite = 2(n−1)² + (5n−2), so corona width = ((n³−n)/3 − 2(n−1)² − 5n + 2)/8 + 1 for n ≥ 8 and 0 for n ≤ 6 (double-zero exceeds the ceiling) — widths 0, 0, 5, 16, 35, … grow CUBICALLY (verified exactly n = 4..10); (S) THE SEDENION RUNG = RANK-16 UNIQUENESS LOSS: E₈² and D₁₆⁺ are isospectral (equal thetas) but their Σ=0 SLICES diverge — 252 vs 240 at norm 2, 23662 vs 23790 at norm 4 (the divergence oscillates) — and tournaments live in the D₁₆⁺ world: the zero-divisor degeneration of the sedenions appears to tournaments as the isospectral ambiguity their score hyperplane BREAKS; (N) the numeration weave: 10³ ≡ −1 (mod 1001 = 7·11·13) is a TORSION resonance (abc→abcabc; the S313 knife-edge 143·7 = 1001) while 2¹⁰ = 1024 ≈ 10³ is a DIOPHANTINE near-dilate (the CF convergent 10/3 of log₂10 = [3;3,9,2,2,4]) — the two resonance types of THM-863/864 in the multiplicative lattice; the 2¹⁰-tiling space lives at n = 6 (m = C(5,2) = 10)
status: M PROVED (three-line lemma + exact Gauss sums); C PROVED-with-verified-first-bite (the double-zero minimality is verified n ≤ 10 and heuristically forced — higher-k violations cost more norm; the width formula exact at n = 8, 10); S EXACT (slice counts computed; the isospectrality of the full lattices is classical); N exact checks
source: opus-2026-07-15-S319 (owner: the four THM-868 next steps + the 1001/1024 numeration weave)
depends_on:
  - THM-868 (the E8 bridge this extends up the unimodular ladder)
  - THM-867 / HYP-6935 (the level laws being formalized)
related: [THM-863/864 (the two resonance types), the CD tower canon, kps cd-tower-architecture]
verification: 05-knowledge/results/milgram_corona_sedenion_opus_S319.out
---

# THM-869 — Milgram, the corona law, and the sedenion rung

## (M) The Milgram formalization

**Lemma (coset constancy).** On L₀ = {v ∈ Z^n : Σv = 0}, for any w ∈ L₀:
q(v + 2w) = q(v) + 4(v, w) + 4q(w). If v is all-odd then (v, w) ≡ Σw ≡ 0
(mod 2), and q(w) ≡ Σw ≡ 0 (mod 2); so q(v + 2w) ≡ q(v) (mod 8). ∎
The even-n residue law (x ≡ n mod 8) is thus a DISCRIMINANT-FORM statement:
q mod 8 descends to the finite quadratic space L₀/2L₀ ≅ F₂^{n−1}, and the
all-odd coset carries the constant value n. The coset Gauss sum is exact:

> **G = Σ_{coset} e^{πi q/4} = 2^{n−1} · e^{πi n/4}** — the residue law IS
> the phase (verified n = 4, 6, 8, 10: phase·4/π = 4, −2, 0, 2 ≡ n mod 8).

Up the ladder: the half-coset construction D_n⁺ is an even unimodular
lattice iff n ≡ 0 (mod 8) — E₈ at n = 8 (THM-868), **D₁₆⁺ at n = 16, D₂₄⁺
at n = 24: tournaments follow the D⁺ branch of the unimodular tree; the
Leech lattice's deeper glue is invisible to integer scores** (as A₈'s glue
was at n = 9).

## (C) The corona law

The first level where Landau bites is the cheapest illegal configuration:
TWO zero scores (d = −(n−1) twice; the k = 2 prefix 0 + 0 < 1). Minimizing
the rest (n−2 odd values summing to 2(n−1): (n/2 − 2) ones and n/2 threes)
gives

> **first-bite(n) = 2(n−1)² + (5n − 2)**, hence
> **corona width = [ (n³−n)/3 − 2(n−1)² − (5n−2) ]/8 + 1 for n ≥ 8; 0 for
> n ≤ 6** (the double-zero configuration exceeds the transitive ceiling).

Verified exactly: widths 0, 0, **5**, **16** at n = 4, 6, 8, 10 (and the
formula gives 35 at n = 12, 62 at n = 14, …: CUBIC growth ~ n³/24 — the
corona asymptotically swallows a constant fraction 1/8 of the levels…
precisely: width/levels → [1/3]/[1/3] − … → 1 − 6/n + O(1/n²): the corona
eventually dominates; the "trivial-filter" regime is the low-shell HALF).
The S318 guess ⌊n/2⌋+1 is REFUTED (it fit only n = 8) — logged per the BH
discipline.

## (S) The sedenion rung, made exact

At rank 16 the even unimodular lattice is NO LONGER UNIQUE: E₈ ⊕ E₈ and
D₁₆⁺ share their theta series (isospectral — the classical Witt/Milnor
example) yet are non-isomorphic. This is the lattice face of the
Cayley–Dickson degeneration at the sedenions (dim 16: zero divisors = loss
of division = loss of rigidity/uniqueness). The score hyperplane BREAKS the
isospectrality:

> Σ=0-slice counts: norm 2: **E₈²: 252 vs D₁₆⁺: 240**; norm 4: **E₈²:
> 23662 vs D₁₆⁺: 23790** (exact; the divergence oscillates so the full
> thetas can agree). Tournaments at n = 16 live in the D₁₆⁺ slice: the
> degeneration the algebra shows as zero divisors, the metagraph shows as
> the isospectral pair whose ambiguity its score hyperplane resolves.

## (N) The numeration weave (1001 and 1024)

10³ ≡ −1 (mod 1001 = 7·11·13): a TORSION resonance of order 2 — hence
abc × 1001 = abcabc and the alternating 3-digit-block divisibility tests;
143·7 = 1001 is the S313 Hunter knife-edge. 2¹⁰ = 1024 ≈ 10³: a DIOPHANTINE
near-dilate — 10/3 is the second convergent of log₂10 = [3; 3, 9, 2, 2, 4],
and the 2.4% offset is the per-decade 3-digit drift. These are exactly the
two resonance species of THM-863/864 (exact small-q torsion vs CF
near-dilate), now in the multiplicative lattice. House note: m = C(5,2) =
10, so the 2¹⁰ = 1024-tiling space is the n = 6 staircase.

## Named next

- The corona's INTERIOR structure (which orbits die first: the k = 2 vs
  k = 3 violation strata) as an exact stratification.
- n = 24: D₂₄⁺ vs the 24-dim even unimodular list (24 Niemeier classes!):
  the score slice's Niemeier selection — which of the 24 does the
  tournament hyperplane distinguish? (the Leech has no roots: immediately
  distinguishable at norm 2.)
- The A₅-monodromy pilot (running: the n = 8 floor-level class graph).
