---
id: HYP-9661
title: "Sources on lines. (A) Line density: for every odd N >= 129 the line D + 2U = N carries at least N/64 one-vertex sources (D = H(B), U = U_B(y)), and the share of its N/6 lattice points that are sources is bounded below and increases with N; near the diagonal U = D the share tends to 1. (B) The share at a fixed offset h is governed by the number of strong tournaments with H = h. (C) Remark on the head sheet. (D) Sources as a local-global problem: no source family with equidistributed states has a local obstruction, and the two missing values are not local."
status: >
  OPEN (opus-2026-11-10-S91, 2026-11-11). (A) is FINITE-EXACT as a lower bound for odd N in [129, 16384) (THM-5129
  (F); four pre-registered tests at rooted order 10, all confirmed) and a conjecture beyond; it implies Chan-Pak
  completeness with the finite base below 129 and is a quantitative form of HYP-9453 (1), (2) and of HYP-9454 (1)
  (klein: "the expected number of pairs on a line of height v grows linearly in v"), which are first. No proof
  route is known; THM-5129 (D) shows that single-base windows are not available at small h. (B) NUMERICAL
  (correlation 0.68 on 23 <= h <= 73). (C) is a remark: "head value" alone is vacuous (PROVED), the question at
  non-sink roots is OPEN. (D) is a reading of THM-5128 (A), (C) with one PROVED statement and one ANALOGY; the cited theorems on thin
  orbits are quoted from the literature and were not re-read in this session.
session: opus-2026-11-10-S91
date: 2026-11-11
source: 05-knowledge/results/s91_f1_local_densities_20261110/README.md
scripts:
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/line_density_order10.py (+ line_density_order9.out, line_density_order10.out), PREREGISTRATION.md
  - 04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_supply.py, head_parity_check.py, census10.py (+ .out)
related:
  - 01-canon/theorems/THM-5129-* (finite supply, permanent omissions, line counts), THM-5128-* (local densities of F1)
  - 05-knowledge/hypotheses/HYP-9453-*, HYP-9454-* (klein: near-join sweep, lines 3H + 2 sigma), HYP-9658-* (seeds of the 2-adic lines)
  - 01-canon/theorems/THM-4865-*, THM-5126-*, THM-5130-*; 05-knowledge/results/session261_creative_chanpak_20261111/, session262_interface_supply_20261111/
  - 05-knowledge/hypotheses/HYP-9659-* (tilted local limit law; the bivariate version is what (A) asks for)
---

# HYP-9661 - sources on lines

Notation of THM-5129. A one-vertex source is `(D, U) = (H(B), U_B(y))`; it gives the head value `D - E + U` and the
closed values `D + 2^t U`, `t >= 1`. Sinks appended to `B` double `U`, so the closed values of all sources are the
numbers `D + 2U`. `U >= D`, hence `D <= N/3` on the line `D + 2U = N`, which has about `N/6` lattice points with
`D` odd. `r(N)` is the number of sources on the line.

## (A) Line density (CONJECTURE; implies Chan-Pak)

1. `r(N) >= N/64` for every odd `N >= 129`.
2. The share `r(N) / (N/6)`, averaged over a dyadic block, is increasing in the block and bounded below by `0.27`.
3. For every `eps > 0` the share of the lattice points with `D <= U <= (1 + eps) D` that are sources tends to 1.

Data (THM-5129 (F)): shares `0.296, 0.334, 0.381, 0.464, 0.518, 0.588, 0.616` on `[2^k, 2^(k+1))`, `k = 7..13`,
from rooted tournaments of order `<= 10`; minimum of `r` over a block `5, 7, 24, 63, 146, 352, 794`; top decile of
`3D/N` on `[4096, 8192)`: `0.94` (odd `U`), `0.84` (even `U`). The two-step family misses `7, 11, 13, 21, 47, 67`
and nothing else below 382,688 (THM-5129 (C), THM-4865 (e)).

Why a count and not a window. For a fixed offset `h` the slopes are finitely many odd parts times powers of two
(THM-5129 (A)): at most `|Delta(h)|` per octave. So the line is not swept by one base; it is hit because the number
of bases `h <= N/3` is linear in `N` and each has a fixed chance. The statement to prove is a bivariate local law
for `(H(B), U_B(y))`, the rooted analogue of HYP-9659.

Refutation: an odd `N >= 129` with fewer than `N/64` sources among all rooted tournaments, or a dyadic block whose
saturated share falls below the previous one.

## (B) Supply follows fibre size (NUMERICAL)

The covered share of `[h, 3h/2]` and the number of strong classes with `H = h` have correlation `0.68` on
`23 <= h <= 73`; the poorest offsets `35, 39, 47, 63` have two strong classes each. So a uniform supply at offset
`h` is a statement about the multiplicity of `h` in the strong spectrum (HYP-9453 (3): at least two strong
realisers for `h >= 15`; the atlas of THM-4810). Question: is `|Delta(h)| >= c h` for all large attained `h`?

## (C) The head sheet (REMARK; one PROVED triviality, one question)

Every attained number is a head value for a trivial reason: if `y` is a sink of `B` then `E = U = D` and
`H(B + z) = D`. So "head value" is not a constraint; the informative heads are at roots that are not sinks, where
`H(B + z) = D + 2 U_even` with `U_even >= 1` (THM-5129 (E)). The four numbers that no tail reaches are such heads:
`11 = 5 + 2 * 3`, `13 = 5 + 2 * 4`, `47 = 17 + 2 * 15`, `67 = 33 + 2 * 17`. The head is a different sheet from
the tails: it counts only the qualifying paths with an even number of vertices, and in the join family F1 it is
the sheet `s = 0`, whose favoured residue class is `-E = -1` instead of `0` (THM-5128 (C)). Question: is every
attained odd `N >= 9` of the form `D + 2 U_even` with `D <= N/3`? (An earlier helper run of this session found
every odd `n` in `[1025, 3069]` to be a head over a base with path count in `(n/5, n/3]`; not re-run.)

## (D) Sources as a local-global problem

**PROVED (THM-5128 (A)).** Let two independent states be glued by a nondegenerate bilinear form of rank `r` (for
F1: `N = x_1 x_2 + 2^s y_1 y_2 - [s = 0]`, `r = 2`). Modulo an odd prime the value is uniform except for one class,
whose excess is to first order the mass of the states that vanish modulo `p`. A family with equidistributed states
therefore has a positive local density in every class: no local obstruction, at any prime.

**Consequence.** In the language of sources: supply is uniform at every prime, the only local structure is the
mass of `p`-divisible sources (`p | D` and `p | U`), and it distinguishes the head (class `-E`) from the tails
(class `0`). The missing values `7 = N(1, 1, 1) = 3 + 2^2` and `21 = N(1, 2, 2) = 5 + 2^4` are invisible to every
such test. In the one-vertex interface they sit in finite heads: the offset 3 has no dyadic slope at all, and the
cone `5 + 2^m` is supplied only from `m = 6` (THM-5129 (B)). Those heads also contain attained numbers (`11, 13,
19, 37`), reached through other interfaces, so a finite head is where the uniform supply stops, not an obstruction.

**ANALOGY (three matching ingredients; the fourth is missing).** The one-apex states are the orbit of `(psi, h,
floor(h/2), 1)` under two integer matrices (one per binary digit), F1 is a bilinear form on two orbits, and the
twisted operators `T_k` have a spectral gap (`rho(T_k) < 2`). Orbit, bilinear form and congruence gap are the three
inputs of the density-one theorems of Bourgain and Kontorovich for integers represented by thin orbits (Zaremba;
Apollonian curvatures). What is missing here is the minor-arc bound. Counter-indication, and the reason no density
argument can be read as evidence for completeness: for Apollonian packings the local-global conjecture is false
(Haag, Kertzer, Rickards, Stange), by reciprocity obstructions that no congruence sees. The repository has the same
kind of obstruction on sub-families: a self-join `x^2 + 2^s y^2 - [s = 0]` is never `7 (mod 8)`, and F1 misses the
primes `4q - 1` with `q` prime at `7, 43, 67, 163, ...` (HYP-9657 (B)). Question: is `{7, 21}` the trace of a
reciprocity-type obstruction on the whole spectrum that happens to be finite?
