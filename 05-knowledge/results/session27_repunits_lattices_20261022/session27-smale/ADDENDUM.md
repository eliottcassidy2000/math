# Smale, binary repunits, and modular forms: bounded source check

Status: CITED classical topology + FINITE-EXACT arithmetic; no new theorem or H-rarity implication. This checks the Kervaire–Milnor lead attributed to concurrent S52 / MSG1913 at a08f25ec; it does not claim that lead or its prior correction as new work.

## What the primary-author source establishes

Milnor's 2011 survey, Theorem 5 (p. 806), gives the cyclic order of the subgroup of homotopy spheres bounding parallelizable manifolds. With standard Bernoulli indexing, for integer m ≥ 2:

|bP_(4m)| = 2^(2m−2) (2^(2m−1)−1) numerator(|4 B_(2m)/m|).

Milnor uses positive old-index numbers B_m = |B_(2m)|; the displayed list on the same page fixes the convention. No additional parity multiplier is inserted. Here bP_(4m) lies inside Θ_(4m−1), not necessarily equal to it. Table 1 gives |Θ_15| = 2·8128 = 16256. Footnote 2 (p. 805) identifies Smale's h-cobordism theorem as the link between the sphere-structure and h-cobordism classifications in the dimensions used here. It is not the order formula itself.

Source: John Milnor, [Differential Topology Forty-six Years Later](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/milnor2011.pdf), Notices AMS 58 (2011), 804–809. The downloaded PDF and p. 806 rendering were directly inspected; the AMS URL returned HTTP 403.

## Exact arithmetic and scope

| m | B_(2m) | 2^(2m−1)−1 | numerator of absolute 4 B_(2m)/m | bP_(4m) order |
|---|---|---|---|---|
| 2 | −1/30 | 7 | 1 | 28 |
| 3 | 1/42 | 31 | 2 | 992 |
| 4 | −1/30 | 127 | 1 | 8128 |

Thus 7, 31, 127 genuinely occur as classical topological order factors. In contrast, 63 = 2^6−1 is not in this odd-exponent factor subfamily. This does not say that 63 cannot divide other topological orders. It gives no equality between tournament Hamiltonian-path counts H and bP orders, and no proof of missing or rare H-values.

## Precisely limited modular-forms connection

Using standard Bernoulli numbers, the constant-one Eisenstein series of weight 2m ≥ 4 is

E_(2m)(q) = 1 − (4m/B_(2m)) Σ_(n≥1) σ_(2m−1)(n) q^n.

This follows by rescaling William Stein's q-coefficient-one convention, E_k = −B_k/(2k) + Σ σ_(k−1)(n)q^n, in [Modular Forms of Level 1, Definition 2.9 and equation (3)](https://wstein.org/books/modform/modform/level_one.html). Bernoulli arithmetic is genuinely shared, but the two normalizations are different: 4B_(2m)/m versus 4m/B_(2m). The modular formula alone supplies neither the Mersenne factor nor a map from H-counts to topology.

## Existing-repository reconciliation

This agrees with the prior correction in `05-knowledge/results/collatz_last_dip_cst_alphabet_exotic_20261005.md`, §3.2–3.3, and `01-canon/MISTAKES.md`, MISTAKE-563. The historical `surgery_exotic_23.py` output double-counted an odd-index factor; do not reuse its bP table. This addendum independently confirms the convention, not a new discovery.

Reproduce the three rows with `python check_arithmetic.py`; output is in `check_arithmetic.out`. The script computes Bernoulli numbers by exact rational recurrence and asserts the three orders. No repository files were changed.
