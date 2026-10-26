---
id: HYP-9418
title: "After THM-4807 and THM-4808: (1) untraced arcs become rare: the share of strong n-classes that carry one tends to 0, and the endpoint relation alone is not a reconstruction carrier; (2) the branch-character height sum has a Chebyshev-type cusp bias: S(X) = (1/2 + 1/sqrt 2) sqrt X + E(X), where E(X) = O(sqrt X) changes sign infinitely often and has logarithmic mean 0; (3) the naive density of the critical-level twin sieve exceeds the true twin count by a ratio increasing to the Mertens constant 4 e^(-2 gamma) (pre-registered values confirmed at 10^8, 10^9, 10^10); (4) the Berggren non-backtracking graphs G_q are nearly Ramanujan: exceptions at 3% of primes below 1000 and 11% in 1000..4000, far below the random-regular-graph 31%, with limiting statistics open; (5) the signed operator of the branch character behaves the same way."
status: >
  OPEN, with NUMERICAL and FINITE-EXACT evidence.
  (1) FINITE-EXACT: shares 1/6, 2/35, 9/353, 33/6008, 167/178133 for n = 5..9 and 1072/9355949 at n = 10 (audit).
  (2) S(X) to 10^10. Both pre-registered predictions (2fe36c8088) CONFIRMED: S/sqrt X = 1.2421 and 1.2194 at 2e9 and
  5e9.
  (3) All four pre-registered predictions (2fe36c8088) CONFIRMED: rho = 1.11158, 1.13123, 1.14480.
  (4) FINITE-EXACT for primes q < 1000 (exact Sturm census in the audit, audit_thm4808/spectra_exact.out). For
  1000 < q < 4000 all three pre-registered predictions are CONFIRMED (non-Ramanujan fraction 0.1099, worst ratio
  1.0105, max |mu| 3.5004), and the audit reproduced them with independent code.
  (5) FINITE-EXACT (audit, exact Sturm census): B_chi is Ramanujan for every prime q < 1000 except q = 251, 479, 863
  and 947. NUMERICAL: 42 exceptions among the 382 primes 1000 < q < 4000.
session: mac-mini-2026-10-25-pyth (2026-10-26)
source: 05-knowledge/results/berggren_nonbacktracking_branch_character_sieve_levels_20261026.md
scripts:
  - 04-computation/experiments/berggren_sieve_parity_20261026/PREREGISTRATION.md: predictions committed before the runs (2fe36c8088 and the (5)-entry commit)
  - 04-computation/experiments/berggren_sieve_parity_20261026/parity_height_big.c (+ .out), parity_height_remainder.c (+ .out), ray_pairing.py (+ .out)
  - 04-computation/experiments/berggren_sieve_parity_20261026/twin_sieve_ratio.c (+ .out)
  - 04-computation/experiments/berggren_sieve_parity_20261026/berggren_ramanujan_scan.py (+ .out, _1000_4000.out), berggren_ihara.gp (+ .out)
  - 04-computation/experiments/hamilton_path_recovery_20261026/audit_thm4807/audit_census_c_n10strong.out
  - 04-computation/experiments/berggren_sieve_parity_20261026/audit_thm4808/   # REPORT.md; exact spectral census, full-resolution E(X), independent sieve
related:
  - 01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md
  - 01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md
  - 05-knowledge/hypotheses/HYP-9417-legality-language-injectivity-goldbach-dependent-decay-and-perfect-hermits.md   # (1), (1b) closed by THM-4807
  - 05-knowledge/results/growing_sieve_path_20261026/sieve/README.md
---

# HYP-9418 — After THM-4807 and THM-4808

## (1) Untraced arcs and smaller legal-path observers

**Status of injectivity.** All-order injectivity of `T ↦ HP(T)`, HYP-9417 (1), is now proved four ways:
- THM-4807, elementarily;
- guard_reconstruction G1, Session70 G4/G5 and growing_sieve_path P1–P4, each relative to the Busch–Jacobson–Reid
  separator.

S61's strong form, HYP-9417 (1b), is THM-4807 (b).

**What remains is quantitative.**
- **Untraced arcs become rare.** Strong `n`-classes that carry an untraced arc number 1, 2, 9, 33, 167, 1072 for
  `n = 5..10`. Their share 0.17, 0.057, 0.025, 0.0055, 0.00094, 0.00011 tends to 0, and labelled tournaments behave
  the same way. Every untraced arc needs two cut vertices (THM-4807 (a)), and almost every tournament is 2-strong.
  - Untraced-arc totals in strong classes are 1, 2, 11, 39, 200, 1207.
  - The per-tournament maximum `(n² − 4n + 3)/8` is Busch's Theorem 4.26.
  - **Open:** the growth rate of the class counts.
- **Smallest reconstructing observers.**
  - The traced core reconstructs.
  - The endpoint relation alone does not. It first collides at `n = 4`, and its numbers of distinct values 1, 2, 8,
    44, 614, 15343 for `n = 1..6` are not in OEIS.
  - **Open:** characterise the image of the endpoint relation, and the smallest set of named guards whose answers
    determine every tournament. Session70 shows that one-step elementary relocations fail at order 4.

## (2) The branch character's cusp bias

THM-4808 (f) proves `S(X) = 1 + N_U(X) + N_D(X) + E(X)`. Here `N_U + N_D = (1/2 + 1/√2)√X + O(1)` counts the two root
rays, and `E(X) = Σ_w χ(w) β_w(X)` has bounded local terms.

**Conjecture.**
- `E(X) = O(√X)`;
- `E(X)` changes sign infinitely often;
- `(1/log Y) ∫_2^Y E(X) X^(−3/2) dX → 0`.

Equivalently, `S(X)/√X` has logarithmic mean `1/2 + 1/√2 = 1.2071…`. This is a Chebyshev-type bias: of the same order as
the fluctuation, and located on the two root cusps. By analogy with Chebyshev's bias for primes (Rubinstein–Sarnak),
the natural form is that `E(X)/√X` has a limiting logarithmic distribution with mean 0. The uniform bound `O(√X)` is
the strongest item, and the data cannot separate it from `O(√X log log X)`.

**Evidence.**

| X | `E(X)/√X` |
|---|---|
| 10⁸ | 0.0090 |
| 2·10⁸ | 0.0849 |
| 5·10⁸ | −0.0107 |
| 10⁹ | −0.0538 |
| 2·10⁹ | 0.0350 |
| 5·10⁹ | 0.0123 |
| 10¹⁰ | 0.0870 |

At every sampled height, `|E|/√X ≤ 0.30`; from `X = 10³` on, `|E|/√X ≤ 0.21`. The sign changes 10 times among the 22 sampled heights from 100 to `10⁹`.

At full resolution (every integer `X ≤ 10⁹`, audit) the bound is larger.
- `max |E|/√X = 0.368`, at `X = 4,222,837`.
- The decade maxima lie between 0.28 and 0.37, with no upward trend.
- `E` changes sign 13,538 times on `[100, 10⁹]`.
- The logarithmic mean of `E/√X` is −0.025 over `[5, 10⁹]` and −0.008 over `[10⁴, 10⁹]`.

**Pre-registered (committed `2fe36c8088`, then run):**
- `S/√X ∈ [1.00, 1.42]` at `2·10⁹` and `5·10⁹`: observed 1.24205 and 1.21936, CONFIRMED.
- `S > 0` at both heights: CONFIRMED.

**Mechanism.** Each non-root ray pair contributes one unit exactly when `X` falls between the interleaved heights.
`E` is therefore a signed lattice-point discrepancy. A proof of `O(√X)` would need cancellation among about `X/2π`
local terms with signs `χ(w)`. A transfer-operator or twisted-Eisenstein route is suggested; none is claimed.

## (3) The twin sieve at its exact level

Define the Mertens ratio

    ρ(X) = (X/2) ∏_{3≤q≤√X} (1 − 2/q) / π₂(X).

The exact level-`√X` sieve (growing_sieve_path W1–W3) decides twin-ness. Its naive independent-residue density
overcounts by `ρ(X) → 4e^(−2γ) = 1.2609…`. This follows from the Hardy–Littlewood and Mertens asymptotics, with
`π₂ ~ 2C₂ Li₂`. It is a heuristic for the limit, not a theorem.

**Pre-registered (committed `2fe36c8088`, then run):** `ρ(10⁸) ∈ [1.100, 1.130]`, `ρ(10⁹) ∈ [1.120, 1.145]`,
`ρ(10¹⁰) ∈ [1.130, 1.160]`, increasing, and below 1.2609. All four were CONFIRMED.

| X | ρ(X) | HL–Mertens model | `π₂(X)` |
|---|---|---|---|
| 10⁸ | 1.11158 | 1.11417 | 440,312 |
| 10⁹ | 1.13123 | 1.13178 | 3,424,506 |
| 10¹⁰ | 1.14480 | 1.14555 | 27,412,679 |

The `π₂` values agree with OEIS A007508, and the audit reproduced them with an independent segmented sieve. These
predictions concern quantities already fixed by published data, so they check the arithmetic, not the model.

## (4) Near-Ramanujan Berggren graphs

`G_q` is the non-backtracking graph of THM-4808 (G), with about `q/4` vertices and degree 4.

**Conjecture.**
- `G_q` is Ramanujan (all nontrivial `|μ| ≤ 2√3`) for most but not all primes `q`.
- The exceptions have positive proportion.
- The excess `max|μ| − 2√3` is `O(q^(−2/3+ε))`. That is the Tracy–Widom scale of random 4-regular graphs; compare
  Huang–McKenzie–Yau on the Ramanujan property of random regular graphs.
- Uniform expansion is a theorem (Selberg plus the Brooks–Burger transfer, Lubotzky 1994 Ch. 4; see THM-4808 (G));
  the Ramanujan property is not expected for all `q`.

**Evidence.** Below 1000 there are 5 exceptions (331, 569, 719, 809, 887) among the 166 primes `5 ≤ q < 1000`, with
worst ratio 1.0133. The census is exact (audit: Sturm counts on `charpoly(A)` for every prime `q < 1000`). The
violating eigenvalue is negative at 331, 569 and 809 and positive at 719 and 887, so both edges matter.

**Pre-registered for 1000 < q < 4000** (entry (5) of PREREGISTRATION.md, committed before the run):
- the non-Ramanujan fraction lies in `[0.03, 0.35]`;
- the worst ratio is below 1.03;
- all `max|μ| ≤ 3.8`.

**Outcome:** 42 of the 382 primes are exceptional, a fraction of 0.1099; the worst ratio is 1.01047 (at
`q = 3539`); the largest nontrivial `|μ|` is 3.5004. All three predictions are CONFIRMED.

The exceptional fraction rose from 0.03 below 1000 to 0.11.
- For random 4-regular graphs, Huang–McKenzie–Yau (arXiv:2412.20263) give a limiting Ramanujan probability of
  about 0.69, so a non-Ramanujan fraction near 0.31. Both observed fractions are well below that.
- A slow approach to 0.31 is possible, but the Tracy–Widom picture does not predict the rise. Whether `G_q` follows
  random-graph statistics is open.
- In `1000 < q < 4000` the worst case (`q = 3539`) is on the negative edge.

## (5) The signed operator

`B_χ = P_U − P_A + P_D` is the branch-signed operator on `P¹(F_q)`. It satisfies a twisted Ihara–Bass identity
(audit, exact for every prime `q ≤ 401`):

    charpoly(B_χ)(x) = (x − 1)^α (x + 1)^β · x^n F_{A_χ}(x + 3/x).

Here `A_χ` is the adjacency twisted by `χ(ι) = χ(Sι) = −1`, `χ(S) = 1`, acting on the `χ`-equivariant functions;
the orbits `{0, ∞}` and `{±1}` drop out. Exactly, for all primes `q ≤ 211` (audit: all primes `q < 251`):
- its eigenvalues are `±1`, together with conjugate pairs of modulus exactly `√3`. The multiplicities of `±1` are
  `(q+9)/4` and `(q+1)/4` for `q ≡ 3 (mod 4)`, and `(q+7)/4` and `(q−1)/4` for `q ≡ 1 (mod 4)`;
- at `q = 251`, a degree-62 irreducible factor has two real roots off the circle, `−1.99990` and `−1.50007`. Their
  product is 3, and the twisted-adjacency eigenvalue is `μ = −3.49998`. The roots are irrational; neither is 2;
- below 1000 the exceptions are exactly `q = 251, 479, 863, 947` (audit, exact). Among the 382 primes
  `1000 < q < 4000` there are 42 (numerical, audit), with worst ratio 1.0098 at `q = 1567`.

**Conjecture.** As in (4), the signed version is nearly Ramanujan. Consequently the branch parity of depth-`k` nodes
decorrelates from every residue class mod `q` at rate about `3^(−k/2)`, up to non-normality constants. This is the
quantitative form of THM-4808 (c).
