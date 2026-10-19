---
id: HYP-9338
title: "Enumerative-geometry dictionary: (a) Apollonius over F_q averages 4 rational circles (THM-4718), and the same Chebotarev-character style counts the Carmichael numbers of the Erdos-1057 paper (random subset products in a finite abelian group hitting the identity; Fejer energy) and the Collatz carries modulo a clock (THM-4520's spectral radius 1/2 is the analogue of the paper's character-energy lower bound); (b) Ford circles are the Soddy degeneration of Apollonius with the Stern-Brocot mediant as Descartes' rule and Minkowski's ? as the binary address; integral Apollonian packings have curvature growth exponent 1.3057 (Kontorovich-Oh) and the local-global conjecture fails by quadratic reciprocity (Haag-Kertzer-Rickards-Stange 2023), while the Ford sub-packing grows linearly (3/pi^2 fractions per unit curvature); (c) THM-4719's complete windows are the lattice-point count of a strip (an Ehrhart-type intersection), and the number of 3-smooth numbers below X is the lattice-point count under a line, ~ (log X)^2/(2 log 2 log 3); (d) leverage of the Carmichael paper: every Carmichael number it constructs is coprime to 6, hence a (3/2)-pseudoprime (A073631, THM-4713's note), so A073631 also has counting exponent 1; codex's product-one compiler is Proposition 5.2 of the paper; conjecture: the Collatz-carry analogue (necklaces of shape (q,p) with 2^q - 3^p | B_w) has expected count C(q,p)/(q |2^q - 3^p|) -> 0 exponentially at the critical slope"
status: >
  Dictionary; rows typed. (a) THM-4718 PROVED; the Carmichael transfer CITED (Howe-Huffman-Kohlmeyer-Kruer, Erdos problem
  1057, Lean-verified per the paper's own account); the Collatz analogy is a reading. (b) Ford/Descartes PROVED elementary
  (473 Farey pairs checked); Kontorovich-Oh and HKRS CITED from memory, to be checked. (c) elementary. (d) the corollary for
  A073631 is PROVED from the paper (all constructed primes exceed y^delta > 3); the carry-count conjecture is OPEN and
  consistent with THM-4591 (no integer cycle of 3x+1 up to period 301,993).
session: mac-mini-2026-10-19-apollonius
source: 05-knowledge/results/apollonius_enumerative_20261019_galois_parity_cube_ford_descartes_windows.md
related:
  - THM-4718, THM-4719, THM-4717, THM-4714, THM-4520, THM-4591, THM-4641/THM-4651 (Carmichael rows), 05-knowledge/results/carmichael_successor_blocks_20261014b.md (codex's product-one compiler)
---

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
