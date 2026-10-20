---
id: HYP-9338
title: "Carmichael and Collatz transfer dictionary: conditional analytic corollaries and an open carry-count heuristic with clock, period, and ensemble guards"
status: >
  OPEN / REPAIRED 2026-10-20: the Collatz carry-count heuristic needs an explicit ensemble, primitive-period correction,
  clock-gap bounds and a positive-cycle guard; scalar affine characters are blind to carries and sparse global equidistribution is obstructed.
  Carmichael asymptotic consequences remain conditional here on the supplied paper's stated analytic inputs (S),(A), not independently proved in this session.
  The generic labelled Apollonius parity-cube group now has an independent squareclass proof; finite-field distribution estimates are separate.
  The original dictionary and conjecture below are retained for provenance; the dated scope notice controls their use.
session: mac-mini-2026-10-19-apollonius
source: 05-knowledge/results/apollonius_enumerative_20261019_galois_parity_cube_ford_descartes_windows.md
related:
  - THM-4718, THM-4719, THM-4717, THM-4714, THM-4520, THM-4591, THM-4641/THM-4651 (Carmichael rows), 05-knowledge/results/carmichael_successor_blocks_20261014b.md (codex's product-one compiler)
---

## Scope repair, 2026-10-20

The Carmichael exponent-one and rational-base-pseudoprime corollaries retain the supplied paper's precise analytic inputs `(S),(A)`; those inputs have not been independently established in this session. The old product-one compiler and the paper's Proposition 5.2 share a mechanism but have different conductor hypotheses. See the [reviewed carry and inheritance packet](../results/carmichael_collatz_carry_20261019/INHERITANCE_AND_BRIDGES.md).

The Collatz carry heuristic remains **OPEN** after repair. Every scalar character of its affine group erases translations, so the scalar-character product argument does not directly transfer. A changing-level radius near `1/2` is not a fixed-clock or fixed-composition gap. The new finite-clock permutation estimates and their 3-smooth almost-invariant vectors give the properly typed comparison.

The positive expression `binom(q,p)/(q*|2^q-3^p|)` cannot be a deterministic relative asymptotic for an integer count when it tends to zero. Any average must specify its ensemble and normalization. Sparse occupancy forbids global total-variation equidistribution when `binom(q,p)/|2^q-3^p|` tends to zero; near-critical entropy alone does not establish that regime without a clock-gap estimate. Dividing labelled words by `q` requires primitive period `q`; otherwise use the proved Burnside/Mobius correction. Positive integer cycles require `D=2^q-3^p>0`. None of these results proves convergence or eventual cycle exclusion.

The [generic Apollonius squareclass proof](../results/apollonius_squareclass_code_20261019/APOLLONIUS_SQUARECLASS_CODE.md) separately closes the full-group gap for labelled input circles with their radii in the base field. A finite-field census still does not prove specialization error terms.

### Original status, retained for provenance

> Dictionary; rows typed. (a) THM-4718 PROVED; the Carmichael transfer CITED (Howe-Huffman-Kohlmeyer-Kruer, Erdos problem
> 1057, Lean-verified per the paper's own account); the Collatz analogy is a reading. (b) Ford/Descartes PROVED elementary
> (473 Farey pairs checked); Kontorovich-Oh and HKRS CITED from memory, to be checked. (c) elementary. (d) the corollary for
> A073631 is PROVED from the paper (all constructed primes exceed y^delta > 3); the carry-count conjecture is OPEN and
> consistent with THM-4591 (no integer cycle of 3x+1 up to period 301,993).

## Historical dictionary and conjecture below: the scope repair above controls their use


# HYP-9338 — enumerative geometry, Apollonius, Ford–Descartes and the Carmichael transfer

| object | count / structure | mechanism | status |
|---|---|---|---|
| Apollonius over `C` | 8 = 2³ | Bézout for three quadrics | KNOWN |
| Apollonius over `F_q` | mean 4; counts 0/4/8 with weights 1/8, 3/4, 1/8 | Galois group = parity cube (THM-4718) | PROVED / FINITE-EXACT |
| Apollonius over `R` | `2·#{e : Δ_e > 0}`; "anything but seven" | signs of pairwise tangency defects | PROVED / CITED (Pedoe) |
| Ford circles | Descartes' rule = mediant; parent = anti-mediant | Soddy degeneration with the line | PROVED (elementary) |
| Apollonian packings | `N(X) ~ c X^1.3057`; local-global false | thin group orbit; quadratic reciprocity | CITED |
| Carmichael numbers | `C(X) = X^(1-o(1))` | subgroup pruning, Fejér energy, Bernoulli products, Korselt | CITED (the attached paper) |
| `(3/2)`-pseudoprimes A073631 | counting exponent 1 | every constructed Carmichael number is coprime to 6 | PROVED corollary |
| Collatz carries mod a clock | expected `C(q,p)/(q|D|)` solutions | character energy ↔ THM-4520's spectral radius 1/2 | OPEN (conjecture) |
| 3-smooth windows | complete iff `M < 3m`, `M < min(2^(q−1), 3^p)` | lattice points in a strip | PROVED (THM-4719) |

## The conjecture

For shape `(q, p)` with `D = 2^q − 3^p`, the number of necklaces `w` with `D | B_w` (the integer cycles of `3x+1` of that shape) is `C(q,p)/(q |D|) · (1 + o(1))` on average over shapes near the critical slope, hence `→ 0` exponentially since `C(q,p) ≤ 2^(q H(p/q))` with `H(log_3 2) = 0.951 < 1`; the Fourier proof would need the carry characters `Σ_w e(a B_w/D)` to cancel like the paper's `exp(−ρ(1−ρ)E_χ)`, which THM-4520 supplies on 3-adic levels but not yet modulo `D`.
