---
id: THM-5131
title: "The mod-4 race of the Hamiltonian-path count: triangle term and square term. (T) For every tournament R on m vertices the sum of (-1)^(c_3) over its 2^m one-vertex extensions is 2^ceil(m/2) (new); hence E (-1)^(c_3) = 2^(-floor((n-1)/2)) over labelled n-tournaments, a closed form recorded and verified for n <= 10 in THM-555 and proved here for every n. (S) Self-converse tournaments of even order have H = 1 mod 4, more precisely H = 1 + 4s mod 8 (both KNOWN), and there are at most (n-1)!! 2^(n^2/4) labelled ones (corollary of a known count); new is their weighing against the lead D_n = P(H = 1 mod 4) - P(H = 3 mod 4) = SC_n/2^C(n,2) + D_n^ns: they carry 150, 86, 60, 37 per cent of D_n at n = 4, 6, 8, 10 (SC_8 = 6737920, SC_10 = 31596429312). FINITE-EXACT: three-strata split and class race for n <= 10; strong-only race for n <= 11; exact residues mod 8 at n = 10, 11 (HYP-9355's pull persists at 5 * 10^-4). NUMERICAL (pre-registered): D_12 4^12 = 989 +- 99 (audit: 1039 +- 125; pooled 1008 +- 77), positive, half the plateau value 2 * 10^3 of orders 8-11 (the registered plateau FAILED; rate OPEN); rank law at order 11 and deletion-pair correlation as predicted. Chan-Pak OPEN"
status: >
  PROVED: (T) (1)-(3); (S) (4); Proposition 5.1 (strong race); the deck remarks.
  KNOWN, credited:
  - the closed form (T) (3): THM-555 (kind-pasteur S20, tiling ensemble, verified n <= 10) with the Gauss-sum mechanism
    of HYP-2707; the classical score formula c_3 = C(n,3) - sum C(s_v,2), of which (T) is an elementary corollary;
  - H = 1 + 2 c_odd mod 4 and the quadratic form of c_3 mod 2 (THM-466, THM-478, THM-559);
  - (S) (1): a self-converse tournament has an anti-involution (Session29 converse-congruences section 1, THM-4930);
    with no fixed vertex it forces c_odd even (THM-791), c_3 even (THM-552);
  - (S) (2), H = 1 + 4s mod 8: Session29 converse-congruences section 3 and Session130 hamiltonian/PROOFS Theorem 1
    (same statement, same proof);
  - (S) (3): corollary of the count of labelled tournaments with a fixed anti-involution (Session29 half-path
    section 5, THM-4743 (6), THM-4930); SC_5, SC_6, SC_7 and the bounds at these orders are in the S312 ledgers;
  - the odd-order balance of self-converse classes (THM-4922 (4)); the self-converse values (THM-4930);
  - exact D_n for n <= 11, |D_n| <= 1/2, the vertex-cube quadratic form and the gate law (THM-5127 and its audits);
    D_4..D_7 also in THM-466; the class split 32 : 24 at n = 6 in THM-1980;
  - T = 1/(1 - S) for labelled tournaments (THM-4816, THM-4817) and the multiplicativity of H over ordered sums
    (classical; THM-1370), of which Proposition 5.1 is a one-line consequence; the sign coordinates of the same-day
    codex note spectrum_multiplicity_signed_tables_20261111;
  - the question answered in section 5 is HYP-9355 (opus S50);
  - the Chebyshev / Mertens taxonomy and the owner's prompt on prime powers (opus S32,
    prime_power_bias_collatz_collector_20261014.md); the atlas of exceptional sets (kind-pasteur S15,
    root_change_local_factors_20261107/ATLAS.md); the crevasses (THM-1370, THM-4889, THM-4986).
  NEW here: the pointwise form (T) (2) and the proof of (T) (3) for every n in the uniform ensemble (which also
  proves THM-555's verified line for n >= 4); the reading of self-converse tournaments as the square term of the
  race and the split (S) (4) with its shares; SC_8, SC_9, SC_10; the exact splits and races listed under FINITE-EXACT.
  FINITE-EXACT: sections 4-5 of the source note (programs and outputs in the packet; every entry recomputed by the
  audit with its own programs; the residues mod 8 at n = 11 are the audit's census).
  NUMERICAL: the order-12 estimates of section 7 (pre-registered R1 FAILED, R2-R4 CONFIRMED; an independent
  estimator of the audit agrees).
  NOT CLAIMED: the Chan-Pak conjecture; any part of HYP-9663; novelty outside the repository.
  AUDIT: PASS WITH CORRECTIONS (05-knowledge/results/mod4_race_square_term_20261111/audit/AUDIT.md): no mathematical
  error; MAJOR M1 ((S) (2) and (S) (3) typed new, both known) and M2 ((T) (3) typed "not found", in THM-555 since
  June); 22 minor; all applied. MISTAKE-760.
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

## (T) The triangle race (PROVED; the value in (3) KNOWN as a verified formula, THM-555)

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

1. (KNOWN: Session29 converse-congruences section 1, THM-4930, THM-791, THM-552.) A self-converse tournament has
   an involutive anti-automorphism, with no fixed vertex when the order is even; then `c_odd` is even and
   `H = 1 (mod 4)`.
2. (KNOWN: Session29 converse-congruences section 3; Session130 `hamiltonian/PROOFS.md` Theorem 1.) With `psi` such
   an anti-automorphism, `H = 1 + 4s (mod 8)`, `s` the number of orbits `{C, psi(C)}` of odd cycles with
   `V(C) cap V(psi(C))` non-empty.
3. (Corollary of the KNOWN count per anti-involution: Session29 half-path section 5, THM-4743 (6), THM-4930.)
   `2^(n^2/4) <= SC_n <= (n-1)!! 2^(n^2/4)` for even `n`; `SC_n <= n!! 2^((n^2-1)/4)` for odd `n`.
4. (PROVED; the weighing is new.) For even `n`: `D_n = SC_n/2^C(n,2) + D_n^ns` with
   `0 < SC_n/2^C(n,2) <= (n-1)!! 2^(-n(n-2)/4)`.

## Proposition 5.1 (strong race; PROVED, a one-line consequence of THM-4816 / THM-4817)

With exponential generating functions over labelled tournaments, `T_chi(x) = 1/(1 - S_chi(x))`, where the
coefficients are the sums of `chi` over all and over strong tournaments.

## Deck remarks (PROVED)

For even `n`, `chi(T) = prod_v chi(T - v)`; for odd `n`, `prod_v chi(T - v) = 1`.

## FINITE-EXACT (all classes of orders 3-10, labelled weights; `out/race_tables.out`)

- `SC_n = 8, 48, 704, 7040, 380416, 6737920, 958578688, 31596429312` at `n = 3, ..., 10` (orders 5-7 were on file).
  Labelled self-converse tournaments are the majority of labelled tournaments exactly for `n <= 5`.
- Share of `D_n` carried by self-converse tournaments at even `n`: `150%, 85.9%, 60.4%, 37.1%` at `n = 4, 6, 8, 10`
  (`3/4, 55/256, 1645/65536, 1928493/2147483648` against `1/2, 1/4, 1361/32768, 650471/268435456`).
- Three strata in units of `4^(-n)` (self-converse / non-trivial automorphism group and not self-converse / rigid
  and not self-converse): `192, -64, 0`; `256, -80, 240`; `880, -216, 360`; `295.5, -7.0, 1417.5`;
  `1645, -104.2, 1181.2`; `79.3, -73.0, 1965.7`; `941.6, -54.1, 1653.3` at `n = 4, ..., 10`.
- Class race (classes with `H = 1` minus classes with `H = 3 mod 4`): `0, 0, 0, 8, 36, 264, 1268, 22600` at
  `n = 3, ..., 10`; at `n = 4, 6` the classes that are not self-converse favour `H = 3` (`-2, -4`).
- Signed counts of strong labelled tournaments: `-2, 24, 176, 4080, 123952, 8135680, 332488064, 75937800960,
  15259291222528` at `n = 3, ..., 11` (two routes for `n <= 10`); the strong-only lead is positive for `4 <= n <= 11`.
- Relative masses of `H = 1, 3, 5, 7 (mod 8)`: `1.00330, 0.99603, 1.00155, 0.99912` at `n = 10` (this packet and the
  first audit of THM-5127); `1.00051, 0.99949, 1.00044, 0.99956` at `n = 11` (census of the audit of this packet;
  labelled counts `9011754468980736, 9002583394936832, 9011181642270720, 9003277512775680` of `2^55`).

## NUMERICAL (Monte Carlo with exact averaging over the last vertex; pre-registered)

`N = 14,400,000` random labelled bases of order 11 (`biasmc.c`; no failure of THM-5127 (C) (1)-(3) or of the parity
law in any sample): `D_12 = (5.89 +- 0.59) * 10^(-5)`, i.e. `D_12 4^12 = 989 +- 99`, positive at ten standard errors
(R1, registered as `[1500, 2500]`: FAILED; the values `1706, 2722, 1972, 2541, 1988` at `n = 7, ..., 11` were a
plateau; the rate of decay is OPEN). Independent estimator of the audit, 9,000,000 bases: `1039 +- 125`; pooled
`1008 +- 77`. Rank law of the pair form at order 11: `0.83869 / 0.15992 / 0.00140` for ranks `10 / 8 / <= 6`, against
`0.83898 / 0.15965 / 0.00137` for a random alternating matrix with an odd number of edges (R2 CONFIRMED);
`rho_11 = 1.0172` (R3 CONFIRMED); corank-1 bias frequency `0.50399` (R4 CONFIRMED); `E_R abs(b_R) = 0.014543`. (R2
passes its registered tolerance; the observed fractions are about three standard errors from the reference law.)
