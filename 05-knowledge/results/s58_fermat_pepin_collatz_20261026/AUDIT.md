# S58 independent audit: findings and responses

**Session:** opus-2026-10-26-S58. **Referee:** an independent subagent with no edit access. It recomputed every number
with its own code (a primorial-gcd sieve, its own Miller–Rabin, the totient test, full cycle enumeration) and read
Boklan–Conway (arXiv:1605.01371) and Somer–Křížek (Czech. Math. J. 54, 2004) in the original.

**Verdict:** 2 MAJOR, 10 MINOR. Every number replicated. The problems were in two statements and two framings, and
all are corrected (MISTAKE-673).

## MAJOR

**M1. THM-4800 (6): "antipodal arcs of one cycle" was false for composite `F_n`.**
- **Finding.**
  - The draft said, with no primality condition, that the Fermat lines of `F_n` and `F_(2^n)` are antipodal arcs of
    one `A`-cycle.
  - In `y = x - 1` the lines start at `y = -1` and `y = +1`, so they share a cycle iff `-1` is a power of `3/4` mod
    `F_n`. That requires `v_2(ord_q(3/4))` to be the same for every prime `q | F_n`.
  - The valuations are `(7, 4)` at `F_5`, `(6, 8)` at `F_6` and `(9, 7)` at `F_7`, so those lines lie on different
    cycles.
  - The script's `antipodal_arcs: true` tested only `x0 + x2 = 2`. That holds for every modulus because
    `A(2 - x) = 2 - A(x)`.
  - "The two lines are the two halves of the Pépin cycle" was wrong even for prime `F_n`. The `F_(2^n)` line is the
    second half, but the `F_n` line has only `2^(n-1)` of the `2^(2^n - 1)` points of the first half.
- **Response.**
  - THM-4800's title and (6) now say:
    - the lines are always mirror images;
    - they lie on one cycle iff `F_n` is prime, and then the `F_(2^n)` line is exactly the second half;
    - at `F_5`, `F_6` and `F_7` they lie on different cycles.
  - README §2 and the status block now say the same.
  - The script tests the same-cycle claim through the 2-adic valuations for `n <= 7`. For `n <= 5` it also walks the
    whole cycle through 0; at `F_5` that cycle has length 11,167,360 and does not contain 2. The JSON key is renamed.

**M2. HYP-9410 and README §1.6 conflated the tested model with the models behind `6.7e-10` and `2^-30`.**
- **Finding.**
  - The calibrated model is an exact sieve to `2^16` followed by `prod (1 - 1/q)^-1 / ln x`. Applied to `F_m` it gives
    `sum_{m >= 33} 19.76/(2^m ln 2) = 6.6e-9`, or `1.5e-8` with the Euler–Lucas sieve bound `2^(m+2)`.
  - Boklan–Conway call that estimate "simple (but loose)". Their `4/2^n` uses the `2^(n+2)`-full restriction on the
    prime factors of `F_n`, which the line values do not share.
  - So the line test cannot validate either cited figure. "Far fewer than one further Fermat prime" survives.
- **Response.** HYP-9410 was rewritten: its title, its "what it says about further Fermat primes" section and its
  status keep the two model classes apart and give `6.6e-9` and `1.5e-8` for the tested class. README §1.6 and the
  status block follow it.

## MINOR

| # | Finding | Response |
|---|---|---|
| 1 | THM-4801 (1) invoked the squarefreeness of Fermat numbers. That is open beyond `F_11`, and it is not needed. | Replaced by pairwise coprimality (a composite `F_i` has no Fermat-prime factor). |
| 2 | THM-4801 (2): Zsigmondy fails for `2^6 - 1`, the case `k = 3`. In (1) Zsigmondy is unnecessary. | `k = 3` treated directly (`9 = 3^2`). (1) now uses the order argument for `2^r - 1`, `r` prime. |
| 3 | README §3 dropped the "below `2^(2^33)`, or assuming completeness" scope. | Scope restored. |
| 4 | THM-4800 (5) said the Mersenne half-turn is "not sufficient". It reduces to `3^((M_p - 1)/2) = -1`; sufficiency is unknown, and no composite `M_p` with `p < 1400` (207 cases) passes. "11 Mersenne primes with `p <= 127`" should read `3 <= p <= 127`. | "Not known to be sufficient", with the 207-case check. Ranges corrected in THM-4800 and the README. |
| 5 | L-P4 and L-P5 do not discriminate: the uncorrected model passes too (z = -0.78, -1.38), with a likelihood ratio of about 3.4. | Stated in HYP-9410, README §1.6 and §4. |
| 6 | "3–6%" should be 3.2–6.9% of E (uncorrected over-prediction up to 7.2% at L-P5). | "3–7%". |
| 7 | `bit_length * ln 2` biases E low by 0.17–0.3%. With exact logs L-P1 has z = -3.30 uncorrected, -1.56 corrected. No pass/fail change. | Noted in HYP-9410 and README §4. |
| 8 | README §1.2's "only three uniform mechanisms" overreached: Cullen's `p \| C_(p-1)` and Hooley-type averaging are other uniform routes. Neither applies to `F_m`. | Softened, with both named. |
| 9 | L-P3 is not independent of L-P1 (`e = 3..10` lies in `m = 8..1200`; `e = 3, 4, 5` contribute nothing). | Noted in HYP-9410 and README §4. |
| 10 | §5: the chain `2 -> 4 -> 16 -> 256 -> 65536` shows five numbers but four squarings. | Rewritten as `2 -> ... -> 65536 = -1 -> 1`, five squarings. |

## Confirmed by the referee

- **THM-4800 (1)–(5):**
  - every proof line, including the converse of (2) and `2^(F_m - 1) = 1`;
  - the mirror run, the Paley statement and the Sachs–Ringel citation;
  - the 17-gon map `5x + 13` with its 8th iterate `x -> 2 - x`;
  - two cycles at `F_1..F_4` by full enumeration;
  - 410 / 1822 / 518 cycles and 8.587 / 10.807 / 9.000-bit deficits at `F_5..F_7`;
  - the Mersenne data;
  - `F_(2^n) = 2 (mod F_n)`, and Theorem L divisibility for `n <= 4` but not `n = 5`.
- **THM-4801:**
  - the six roots by the totient test;
  - the trunk `k` and both initial-segment sets;
  - the `2^(2^33)` bound argument (orbits stay below `2^34`).
- **HYP-9410:**
  - all five families replicated exactly (survivors 52743 / 10342 / 835 / 29492 / 7141);
  - completeness of the algebraic criteria, by Capelli and the Sophie Germain identity;
  - every flagged value certified composite;
  - the immunity proofs;
  - the mod-7 exclusion;
  - the per-line constants 3.205 / 3.184 / 3.321 / 3.306, with long-line limit 3.282 (singular series 1.353);
  - the pre-registration order (`eb553116c2` < `a9060e82d9` < `d738adf08e`, no edits to the predict/test paths).
- **Citations:**
  - the §1 no-go facts;
  - the open infinitude of composite `F_n`;
  - Boklan–Conway eq. (11) giving `2^-30`;
  - the Cayley–Dickson rungs of THM-871;
  - Somer–Křížek Thm 4.1, Cor. 4.3 and Thm 4.4 as cited.
- **Prior art:** none found.
  - Lagarias's 3x+1 bibliographies (1963–2009) contain no "Fermat" or "Pépin".
  - The Křížek–Luca–Somer book has no Collatz material.
  - OEIS has neither `1, 5, 85, 21845, 1431655765` nor `1, 3, 5, 85, ...`.
  - The content is elementary (Pépin, Euler's criterion, Gauss–Wantzel); the novelty is the identification.
