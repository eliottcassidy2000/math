---
id: THM-5131
title: "The mod-4 race of the Hamiltonian-path count: triangle term and square term. (T) For every tournament R on m vertices the sum of (-1)^(c_3) over its 2^m one-vertex extensions is 2^ceil(m/2); hence E (-1)^(c_3) = 2^(-floor((n-1)/2)) over labelled n-tournaments: an even number of cyclic triangles leads at every order by an exact margin. (S) Self-converse tournaments of even order have H = 1 mod 4 (KNOWN), more precisely H = 1 + 4s mod 8 with s the number of orbits of odd cycles that meet their image; there are at most (n-1)!! 2^(n^2/4) labelled ones, so the lead D_n = P(H = 1 mod 4) - P(H = 3 mod 4) splits as SC_n/2^C(n,2) + D_n^ns. FINITE-EXACT: SC_n = 48, 7040, 6737920, 31596429312 at n = 4, 6, 8, 10, carrying 150, 86, 60, 37 per cent of D_n; three-strata split and class race for n <= 10; strong-only race for n <= 11 through the free monoid of ordered sums; exact residues mod 8 at n = 10, 11 (HYP-9355's pull persists at 5 * 10^-4). NUMERICAL (pre-registered): D_12 4^12 = 989 +- 99, positive, half the plateau value 2 * 10^3 of orders 8-11 (the registered plateau FAILED); rank law at order 11 and deletion-pair correlation as predicted. Chan-Pak OPEN"
status: >
  PROVED: (T) (1)-(3); (S) (2)-(4); Proposition 5.1 (strong race); the deck remarks.
  KNOWN, credited:
  - H = 1 + 2 c_odd mod 4 and the quadratic form of c_3 mod 2 (THM-466, THM-478);
  - (S) (1): a fixed-point-free involutive anti-automorphism forces c_odd even (THM-791), c_3 even (THM-552); the
    odd-order balance of self-converse classes (THM-4922 (4)); the self-converse values (THM-4930);
  - exact D_n for n <= 11, the residues mod 8 at n = 10, 11, |D_n| <= 1/2, the vertex-cube quadratic form and the
    gate law (THM-5127 and its two audits);
  - H multiplies over ordered sums (classical; THM-1370 uses it), the sign coordinates of the same-day codex note
    spectrum_multiplicity_signed_tables_20261111;
  - the question answered in section 5 is HYP-9355 (opus S50);
  - the Chebyshev / Mertens taxonomy and the owner's prompt on prime powers (opus S32,
    prime_power_bias_collatz_collector_20261014.md); the atlas of exceptional sets (kind-pasteur S15,
    root_change_local_factors_20261107/ATLAS.md); the crevasses (THM-1370, THM-4889, THM-4986).
  NEW here: (T); the mod-8 refinement (S) (2); the count (S) (3) and the reading of self-converse tournaments as the
  square term of the race; the exact splits and races listed under FINITE-EXACT.
  FINITE-EXACT: sections 4-5 of the source note (programs and outputs in the packet).
  NUMERICAL: the order-12 estimates of section 7 (pre-registered R1 FAILED, R2-R4 CONFIRMED).
  NOT CLAIMED: the Chan-Pak conjecture; any part of HYP-9663; novelty of (T) outside the repository.
  AUDIT: AUDITSTATUS
session: mac-mini-2026-11-11-race
source: 05-knowledge/results/mod4_race_square_term_20261111/README.md
related:
  - HYP-9663
  - THM-5132
  - THM-5127
  - HYP-9355
  - THM-466 (the-2-adic-digits-of-H-are-odd-cycle-collection-counts)
  - THM-478 (digit-one-degree-law)
  - THM-552 (c3-parity-dichotomy-self-converse)
  - THM-791 (H-companion-laws-to-the-transitivity-flow)
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - THM-1980
  - THM-4860
  - THM-4889
  - THM-4922
  - THM-4930
  - THM-4986
verification:
  - 05-knowledge/results/mod4_race_square_term_20261111/scrace2.c
  - 05-knowledge/results/mod4_race_square_term_20261111/race_tables.py
  - 05-knowledge/results/mod4_race_square_term_20261111/c3_bias.py
  - 05-knowledge/results/mod4_race_square_term_20261111/strong_signed.py
  - 05-knowledge/results/mod4_race_square_term_20261111/out/
---

# THM-5131 - the mod-4 race of the Hamiltonian-path count: triangle term and square term

Full statements, proofs and tables: [source note](../../05-knowledge/results/mod4_race_square_term_20261111/README.md).
The Chan-Pak conjecture (every odd number except 7 and 21 is a Hamiltonian path count) is **OPEN**.

## Notation

`H(T)`: number of Hamiltonian paths; `c_k`: number of directed `k`-cycles; `c_odd = c_3 + c_5 + ...`;
`chi(T) = (-1)^(c_odd(T)) = (-1)^((H(T)-1)/2)`; `D_n = E chi` over labelled tournaments on `n` vertices;
`R + v_X`: add a vertex `v` with out-neighbourhood `X`; `SC_n`: number of labelled self-converse tournaments.

## (T) The triangle race (PROVED)

For a tournament `R` on `m >= 0` vertices:

1. `c_3(R + v_X) = c_3(R) + e_R(X, V \ X)`, and `e_R(X, V \ X) = sum_{a in X} s_R(a) - C(|X|, 2)`.
2. `sum_X (-1)^(c_3(R + v_X)) = 2^ceil(m/2)`, for every `R`.
3. `E (-1)^(c_3) = 2^(-floor((n-1)/2))` over labelled tournaments on `n >= 1` vertices.

Mechanism: the triangle part of the vertex-cube form is `sigma . x + e_2(x)` with `sigma` the score parities; its
Gauss sum is `2^((m+1)/2) cos((p - q - 1) pi/4)` (`p`, `q` the numbers of even and odd scores), and
`(-1)^(c_3(R))` times it is invariant under arc reversal. For `n <= 4`, `chi = (-1)^(c_3)`, so `D_3 = D_4 = 1/2`.
In score form: `c_3 = C(q, 2) + floor(m/4) (mod 2)`, `q` the number of vertices of odd score (THM-5132 (1)). The
next layer - pentagons, and with them `D_5 = 13/32` and `D_6 = 1/4` - is THM-5132.

## (S) The square term

1. (KNOWN: THM-791, THM-552.) A self-converse tournament has an involutive anti-automorphism, with no fixed vertex
   when the order is even; then `c_odd` is even and `H = 1 (mod 4)`.
2. (PROVED.) With `psi` such an anti-automorphism, `H = 1 + 4s (mod 8)`, `s` the number of orbits `{C, psi(C)}` of
   odd cycles with `V(C) cap V(psi(C))` non-empty.
3. (PROVED.) `2^(n^2/4) <= SC_n <= (n-1)!! 2^(n^2/4)` for even `n`; `SC_n <= n!! 2^((n^2-1)/4)` for odd `n`.
4. (PROVED.) For even `n`: `D_n = SC_n/2^C(n,2) + D_n^ns` with `0 < SC_n/2^C(n,2) <= (n-1)!! 2^(-n(n-2)/4)`.

## Proposition 5.1 (strong race; PROVED)

With exponential generating functions over labelled tournaments, `T_chi(x) = 1/(1 - S_chi(x))`, where the
coefficients are the sums of `chi` over all and over strong tournaments.

## Deck remarks (PROVED)

For even `n`, `chi(T) = prod_v chi(T - v)`; for odd `n`, `prod_v chi(T - v) = 1`.

## FINITE-EXACT (all classes of orders 3-10, labelled weights; `out/race_tables.out`)

- `SC_n = 8, 48, 704, 7040, 380416, 6737920, 958578688, 31596429312` at `n = 3, ..., 10`. Self-converse tournaments
  are the majority of labelled tournaments exactly for `n <= 5`.
- Share of `D_n` carried by self-converse tournaments at even `n`: `150%, 85.9%, 60.4%, 37.1%` at `n = 4, 6, 8, 10`
  (`3/4, 55/256, 1645/65536, 1928493/2147483648` against `1/2, 1/4, 1361/32768, 650471/268435456`).
- Three strata in units of `4^(-n)` (self-converse / non-trivial automorphism group and not self-converse / rigid
  and not self-converse): `192, -64, 0`; `256, -80, 240`; `880, -216, 360`; `295.5, -7.0, 1417.5`;
  `1645, -104.2, 1181.2`; `79.3, -73.0, 1965.7`; `941.6, -54.1, 1653.3` at `n = 4, ..., 10`.
- Class race (classes with `H = 1` minus classes with `H = 3 mod 4`): `0, 0, 0, 8, 36, 264, 1268, 22600` at
  `n = 3, ..., 10`; at `n = 4, 6` the classes that are not self-converse favour `H = 3` (`-2, -4`).
- Signed counts of strong labelled tournaments: `-2, 24, 176, 4080, 123952, 8135680, 332488064, 75937800960,
  15259291222528` at `n = 3, ..., 11` (two routes for `n <= 10`); the strong-only lead is positive for `4 <= n <= 11`.
- Relative masses of `H = 1, 3, 5, 7 (mod 8)`: `1.00330, 0.99603, 1.00155, 0.99912` at `n = 10`;
  `1.00052, 0.99948, 1.00044, 0.99956` at `n = 11`.

## NUMERICAL (Monte Carlo with exact averaging over the last vertex; pre-registered)

`N = 14,400,000` random labelled bases of order 11 (`biasmc.c`; no failure of THM-5127 (C) (1)-(3) or of the parity
law in any sample): `D_12 = (5.89 +- 0.59) * 10^(-5)`, i.e. `D_12 4^12 = 989 +- 99`, positive at ten standard errors
(R1, registered as `[1500, 2500]`: FAILED; the values `1706, 2722, 1972, 2541, 1988` at `n = 7, ..., 11` were a
plateau). Rank law of the pair form at order 11: `0.83869 / 0.15992 / 0.00140` for ranks `10 / 8 / <= 6`, against
`0.83898 / 0.15965 / 0.00137` for a random alternating matrix with an odd number of edges (R2 CONFIRMED);
`rho_11 = 1.0172` (R3 CONFIRMED); corank-1 bias frequency `0.50399` (R4 CONFIRMED); `E_R abs(b_R) = 0.014543`.
