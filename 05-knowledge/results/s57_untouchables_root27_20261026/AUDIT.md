# Independent audit of S57 (THM-4790, THM-4792, HYP-9400 and the results note)

Two independent referee agents, briefed only with file paths and claims, were run on the pre-audit commit
`c13c51552a` (2026-10-26). They edited no repository files; their scratch code lived in separate temporary directories.
Their reports are condensed here, each finding followed by the response. All corrections are applied in the same
commit as this file; the ledger entry is MISTAKE-672.

## Referee A: THM-4790 (tc fibres) and THM-4792 (sheet tables)

### THM-4790: no MAJOR finding

**OK** (each independently re-derived or recomputed):
- **(a) Reduction.** H and Aut factor over the strong components. Busch's Theorem 1 and Corollary 1 were checked in
  the EJC original. `L_n = 3, 3, 5, 7, 11` holds even under the weaker bound `5^((n-1)/3)`.
- **(b) The group lemma, case by case.**
  - The whole subgroup lattice of `Z_3 wr Z_3` was brute-forced: orders `1, 3, 9, 27, 81` occur `1, 22, 22, 4, 1`
    times. Every order-9 subgroup contains `(3)(3)` or `(3)(3)(3)`. Every order-27 subgroup, and the group itself,
    contains `(3)(3)(3)`.
  - The `n_5 in {1, 81}` step was checked, together with the forced 9-point orbit and the primitive/imprimitive Jordan
    argument.
  - So was the `Z_5^2` argument.
- **The candidate filter is necessary** (Lagrange). With the filter removed entirely, no hit.
- **Invariant enumeration covers all classes.** An odd-order `sigma` never reverses a pair.
- **Independent re-implementation.** It used different conjugates, a pull Held–Karp, and an automorphism counter on
  refined invariants. It reproduced every number in the table, including the `|Aut|`-computation counts 12, 6,240,
  1,660, 69,616 and 57,394.
  - Validation against brute force: 150 random and 304 symmetric tournaments.
  - Agreement with the session's `ham`/`aut_count` on 25,674 strong invariant tournaments with
    `|Aut| in {3, ..., 81}`.
- **Lemma-bypassing runs**, all with no hit:
  - `n = 9`, type `(3)`, `2^22`, minimum `tc` 25;
  - `n = 11`, type `(5)`, `2^23`;
  - `n = 11`, type `(3)(3)`, `2^25`.
  So `n <= 11` holds with only Cauchy and the order-9 argument.
- **The graph-chart identity.** `tc(T) = #{G : T_G iso T}` was checked numerically on 158 symmetric tournaments.
- **Prior art.** Nothing on `H/|Aut|` fibres. The H version is MathOverflow 232751, Grinberg's notes and THM-1370.

**MINOR**, with responses:
1. Only strong invariant tournaments were enumerated. The title and (c) now say so.
2. The self-tests exercise `|Aut| > 1` only on C3, T5 and P7. The audit's cross-validation is now recorded in the
   status.
3. "133,000" in the README should be 134,922. Fixed.
4. Name Sylow's theorem in the status. Done.
5. Cite MathOverflow 232751 directly. Done.

### THM-4792: one MAJOR finding (prior art)

**MAJOR.** The classification in (c), that `-1` lies in `<2>` mod `m` iff every `p | m` has `ord_p(2)` even with a
common 2-adic valuation, is P. Moree, "On the divisors of `a^k + b^k`", Acta Arith. 80 (1997) 197–212, Theorem 1
(`a = 2`, `b = 1`).
- Moree's Lemma 1 is the same compatibility argument.
- His Propositions 3–4 give the prime-power reduction.
- His sec. 5 gives the union-over-classes description and the prime densities 7/24, 7/24, 8/24, 1/24, ...
- The moduli are OEIS A014657.

*Response:* (c)2 is retyped KNOWN with the proof marked as reproduced. The README's "classification left open by
kind-pasteur" now credits Moree.

**MINOR**, with responses:
1. The mod-8 valuations are the second supplementary law. Now labelled classical.
2. (c)3, the `+-1/4` coset, is not in kind-pasteur's note; it is a trivial restatement. Status corrected.
3. The source of `U_-(x) = -U_+(-x)` is the 2026-10-05 note; THM-4474 (E) is negation as an isomorphism of parity
   graphs. Corrected.
4. The "Syracuse not conjugate" control was vacuous, since the unit shift maps odd to even. Removed from the script;
   the THM now states the parity reason.
5. "15 is the least such modulus" is trivial: it is the least odd number with two distinct prime factors.
   Rephrased.

**OK:**
- every proof line by line;
- the classification on all 50,000 odd `m <= 100,001`, four criteria with zero mismatches (12,743 alternating);
- the mod-8 classes for all primes below `2·10^6`;
- the literal table shift for `m <= 1201`;
- the four identities for `|x| <= 10^6`;
- the root-27 orbit `8 -> 12 -> 18 -> 27 -> 40` against `7 -> 11 -> 17 -> 26 -> 39`;
- faithful representation of the 2026-10-05 Prop. 1 and kind-pasteur's A.3–A.4.

## Referee B: untouchables (README secs. 2–4, 8; HYP-9400)

(The first run was cut off by a network error and resumed with its context intact.)

**MAJOR 1: the failure criterion for PP's assumption was wrong.**
- The note read PP sec. 3.1 correctly. But "membership in `T_a` needs every small prime power to divide `sigma(m)`" is
  sufficient, not necessary.
- `s(n)` leaves `T_a` at `l` iff `v_l(s(n)) != v_l(n)`, capped at `v_l(A_y)`. For `l` not dividing `n` this is
  `sigma(n) = n (mod l)`.
- On `(2^26, 2^27]` the literal per-prime failure is `0.177, 0.241, 0.176, 0.150, 0.124` for `l = 3..13`. The note's
  `l` not dividing `sigma(n)` shares are `0.145, 0.539, 0.582, 0.827, 0.806`, a different quantity.
- So "fails at 5 for half of all inputs, at 13 for four fifths" was wrong, and so was "below 10% only after `log X`
  exceeds `10^12`": the real figure is about 240.
- The joint failure is larger than the note claimed: `gcd(s(n), A_y) != gcd(n, A_y)` for 62% of `n` at `y = 13` and
  82% at `y = 53`.
- The referee proposed relabelling `l` not dividing `sigma(n)` as the failure of `s` to act as `n -> -n mod l`.

*Response:* reproduced independently by the session (`pp_class_preservation_20261026.py`): `0.623` / `0.820` joint,
per prime `0.177, 0.241, 0.175, 0.150, 0.124`. Sec. 2.4 now reports two departures: (a) the literal class
assumption, per prime and joint, and (b) residue neutrality (`f_l`). The HYP-9400 title and (3) are rewritten, and the
`10^12` sentence is removed.

**MAJOR 2: the rate in HYP-9400 (4)(a) contradicts the data.**
- `(spread - 1)` fitted against `log k` over `W_24..W_33` gives exponents `-0.76, -0.71, -1.09, -0.75` (mod 5, 7, 11,
  13), not `-1/(l-1)`.
- The clause "mod-13 spread above 1.10 for `k < 10^30`" is unsupported: the fit extrapolates to `1.088`.
- Mod 3 is non-monotone: it falls to `1.0095` at `W_18`, where the argmax flips, then rises to `1.0923`. That
  contradicts "decrease monotonically" and "fall slowly from 1.08".
- The `W_40` intervals are consistent: the extrapolations `1.392` and `1.174` fall inside.

*Response:*
- the `10^30` clause is removed;
- the measured exponents are stated, with the thin-core exponent kept only as an unobserved heuristic for the eventual
  regime;
- "monotonically" is restricted to `l = 5, 7, 11, 13`, and the mod-3 behaviour is described;
- the `W_40` ranges are kept.

**MINOR**, with responses:
- **Framing.** PP sec. 3.2 already says `Delta` is effectively computable and that their printed numbers are
  "intermediate values" of (3.4). Credited. What is new is the evaluation itself.
- **Kill/serve rule.** It is the classical residue exclusion: PP 2016 sec. 2, proof of Thm. 1.3; Pomerance–Yang 2014
  Prop. 2.3. Retyped as a classical mechanism. What is new is the exact core decomposition of the data.
- **Mod-7 ordering.** It has three discordant pairs, not one near-tie. Now "rank correlation 0.77"; mod 5 still has
  one near-tie.
- **Overclaims.** "The same failure explains the per-number law" and "accounts for the measured residue and abundancy
  structure" are softened. The deficient-end excess is now marked unexplained, with the class failure only as a
  candidate.
- **Monte Carlo tail.** The correction above the primes is `2.07e-5`, not `1.4e-5`. Fixed. The MC checks `Delta` only
  to about `1e-5`; the `+-1e-8` rests on the window averages plus the referee's deterministic convolution
  (`0.17865689139`, agreeing to `1e-11` across cutoffs).
- **Partial-sum fit.** The converged constant is `0.1621075`, with residual `6–7e-8` at `y = 10^6`. Fixed.
- **`W_34`** is the partial window `(2^33, 10^10]`. Stated.
- **Citation.** For the density of `n` with `l` not dividing `sigma(n)`: Ramanujan's manuscript sec. 19, Watson
  (Math. Z. 39, 1935), Rankin (Proc. Glasgow Math. Assoc. 5, 1961), Scourfield (Acta Arith. 10, 1964) and
  Ciolan–Languasco–Moree (JMAA 519, 2023). Adopted.

**OK** (independently recomputed):
- `Delta = 0.17865689139`, by deterministic convolution and by window averages on `[1e7, 2e7)` and `[5e7, 6e7)`.
- `S(y)/log y` for `y = 10^3..10^8`, to 10 digits.
- The formula reproduces all ten of PP's intermediate values (`10^6` to `2·10^10`) and CGJM's Table 2, to 6 digits.
- The kill/serve rule holds on every preimage in `[1e5, 1e6]` with `P > 13`.
- The `10^6` core JSON was reproduced exactly; eight `10^8` core rows agree with exact prime counts.
- All `W_33` numbers: local density `0.168696`, mod-5 shares and spread `1.4533`, mod-13 `1.2014`, mod-3 `1.0923`,
  x-bins `0.22406` and `0.50386`. Also the mod-15 extremes and the U-P5 moments.
- The arrays agree on their overlap, and with brute force for every even `k <= 2·10^6`; they equal A057978 at every
  power of ten.
- The decade fits: `0.1773–0.1778` (`1/log^2`) and `0.1845–0.1866` (`1/log`).
- Odd Goldbach dependence: `443,732` / `56,266` at `10^6`, and `8,250` at `10^5`, by an independent brute force.
  `10^7` and `10^8` were not recomputed.
- The F1 statements, and `785` of `3,661` prime powers with `q + 1` untouchable (first `q = 289`).
- Prior art: no source states the limit of PP's constant, and no prior report of nonzero-residue share variation was
  found. Searched:
  - PP 2016;
  - CGJM;
  - OEIS A005114 and A057978;
  - erdosproblems.com 418 and its 9 forum comments;
  - Booker (arXiv:1610.07471);
  - Pomerance–Yang;
  - targeted searches for 0.1786, 0.1787 and 0.17866.

**Found by the session before the audit report:** README sec. 3 listed only 9, 11, 39 as single-source odd numbers.
3 (`= s(4)`) and 7 (`= s(8)`) also have exactly one source, so the complete list below `10^8` is `3, 7, 9, 11, 39`.
Brute force to 2000 confirms it; fixed.
