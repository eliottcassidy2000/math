---
id: THM-883
title: THE VERBLUNSKY WALL COORDINATES — kps's OPUC/Blaschke dictionary built into exact machinery: (V) the tight-locus measure (uniform on the primitive 14th roots; moments = Ramanujan sums c₁₄(k)/6) has EXACT Verblunsky vector α = (1/6, −1/5, 1/4, −1/3, 1/2, −1), i.e. α_n = (−1)ⁿ/(6−n) — the alternating harmonic countdown to the boundary — with Φ₆^{OPUC} = Φ₁₄^{cyclotomic} and termination |α₅| = 1 (TIGHT = ATOMIC = the orbit reaching the circle, kps HYP-3793 made exact); the heptagon (primitive 7th roots) has α_n = −1/(6−n) and THM-882's half-turn is EXACTLY the Verblunsky sign law α_n ↦ (−1)^{n+1}α_n (verified in exact rationals); (F) THE CRITICALITY FLOW: for the AP {1..13}, the safe-set measure's α₅(δ) blows up to the boundary AS δ → 1/14 (1−|α₅|: 0.79 at δ = 1/16 → 1.5e−4 at 5/72 → 3.9e−9 at threshold−1e−4; empirical rate ~ (1/14−δ)^3.5): THE MAXIMIN IS THE VERBLUNSKY TERMINATION TIME; (W) the wall coordinates of loose species at δ = 1/14: the THM-777 floor family {1..13}∖{6} has the near-period-2 boundary-hugging signature (0.49, −1.00, 0.51, −1.00, 0.47, −1.00, …) — the two-effective-atom deep-well geometry read instantly; GAP: an α₅ ≈ −0.92 spike; near-AP: spread small-α profile — the α-vector is a discriminating coordinate system for safe-set geometry
status: V PROVED-EXACT (rational Gram-Schmidt from Ramanujan moments; the cyclotomic identification and termination exact; the closed form α_n = (−1)ⁿ/(6−n) exact for q = 14, α_n = −1/(6−n) for q = 7; general-q conjecture: uniform-on-primitive-q-th-roots has α_n = μ-structured/(φ(q)−n) — the prime case is the known one-point-deleted Geronimus family); F/W computed (float OPUC on exact arc moments)
source: opus-2026-07-16-S322 (owner: build the OPUC/Blaschke coordinate and the Verblunsky wall coordinates); builds kps's a-dictionary-of-loop-maps (the Szegő recursion = Blaschke composition; s(2,7) = 1/14 the Dedekind quantum) and HYP-3793
depends_on:
  - THM-882 (the cyclotomic face this coordinatizes)
  - THM-873 (the tight locus)
related: [kps loop-map dictionary + HYP-3793 (tight = atomic, now exact), THM-877 (Ramanujan truncation — the moments ARE Ramanujan sums), Simon OPUC (finite-gap ↔ isospectral tori: the safe sets' α-asymptotics live on tori — the D₁₆⁺/E₈² isospectrality thread's OPUC cousin)]
verification: 05-knowledge/results/verblunsky_wall_coordinates_opus_S322.out
---

# THM-883 — the Verblunsky wall coordinates

## (V) The tight locus, coordinatized exactly

The uniform measure on LRC(14)'s tight locus (= the primitive 14th roots of
unity, THM-873/882) has moments c_k = c₁₄(k)/6 — **Ramanujan sums** — hence
rational OPUC data. Exact rational Gram–Schmidt gives:

> **α = (1/6, −1/5, 1/4, −1/3, 1/2, −1)**, i.e. **α_n = (−1)ⁿ/(6−n)** —
> the alternating harmonic countdown; Φ₆^{OPUC} = Φ₁₄ (the cyclotomic);
> |α₅| = 1 (termination: the measure is atomic with 6 atoms).

The heptagon (primitive 7th roots): α_n = −1/(6−n) (all negative), and the
half-turn z ↦ −z acts on Verblunsky coefficients as α_n ↦ (−1)^{n+1}α_n —
carrying the heptagon vector EXACTLY to the tight-locus vector: **THM-882's
cyclotomic transport is a Verblunsky sign flip.** (General-q conjecture: the
one-point-deleted/primitive-roots family has α_n = −1/(φ(q)−n)-structured;
the prime case is the classical Geronimus-type example; q = 14 adds the
Möbius sign alternation.)

## (F) The criticality flow — the maximin as termination time

For S = {1..13}, the normalized Lebesgue measure on the safe set G(δ) flows,
as δ ↑ 1/14, from the full circle (all α = 0) to the 6-atom tight-locus
measure. Computed: 1 − |α₅(δ)| = 0.93, 0.89, 0.79, 1.5e−4, 3.8e−7, 3.9e−9 at
δ = 1/28, 1/20, 1/16, 5/72, ~1/14 − 1e−4, ~1/14 − 1e−5. **Tightness is the
Verblunsky orbit hitting the boundary, and M(S) is its termination time** —
kps's verified metaphor (HYP-3793: "tight = atomic = the orbit reaching the
circle") is now a computed blow-up curve with empirical rate ~ (1/14 − δ)^3.5.

## (W) Wall coordinates of the loose species (δ = 1/14)

| family | arcs | μ(G) | α-signature |
|---|---|---|---|
| floor {1..13}∖{6} | 4 | 0.0082 | (0.49, −1.00, 0.51, −1.00, 0.47, −1.00, …) — near-period-2, boundary-hugging |
| GAP {1..7}∪{10..16} | 8 | 0.0516 | spike at α₅ ≈ −0.92 |
| near-AP {21..32} | 38 | 0.144 | spread, small |

The floor family's signature is the classical two-band/two-atom pattern:
its 4 arcs pair into two clusters and the α's say so at a glance — the
THM-777 deep-well geometry read directly from the coordinate. Loose species
are separated by their α-profiles: **the Verblunsky vector is the wall
coordinate system** — walls (arc endpoints) in, spectral shape out, with
Simon's finite-gap theory placing each species on an isospectral torus of
dimension #gaps − 1 (the OPUC cousin of the S319 slice-isospectrality
thread).

## Named next

- Prove the general-q closed form (the primitive-roots Verblunsky family).
- The α-flow's blow-up exponent (~3.5): exact value and mechanism (six arcs
  shrinking linearly — the rate should be derivable from the arc widths).
- The selector re-coordinatization: run THM-797/803/817's signed-complement
  families through the α-map — the period-2 signature may classify the
  two-sheet survivors (the "deep" families should all hug the boundary).
- Verblunsky coefficients as the LRC continuation state (the frontier's
  "retain owner/metric stalks" — the α-vector is a canonical compressed
  state; test whether it is Markov for comb insertion where exact endpoints
  were not, THM-840).
