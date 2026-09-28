# The Lehmer five and Collatz: two 2-adic engines with opposite memory (aliquot drivers persist for free, Collatz shadows are priced), the sequence of 276 read through its drivers, perfect numbers at the Mersenne indices, and D21/D22 closed

**Session:** opus, `collatz-posets-zeta5-20260927` (ninth note), 2026-09-27.
**Owner's directive:** "pursue D21 and D22; my reference to the Lehmer 5 was
the concept of aliquot sums, and the family around 276 that appears
infinite, and its relation as with Collatz, looking for as creative useful
connections between those concepts as possible."
**Correction.** The eighth note read "the Lehmer 5" as Emma Lehmer's simplest
quintic; the owner meant the **Lehmer five** `276, 552, 564, 660, 966`, the
five smallest numbers whose aliquot sequences are not known to terminate or
enter a cycle (named after D. H. Lehmer; Wikipedia's aliquot-sequence page
read today). The eighth note's Brocard–conductor coincidence stands as an
exact numerical fact under its own name; its label is corrected there and
the misreading is logged in MISTAKES.
**Inherits (cited):** the fifth note (price sheet: a Collatz cycle shadow of
depth `K` costs `K` bits of 2-adic precision and yields `c_w K` bits of
growth; Terras's memoryless valuations), the fourth note (family dichotomy),
the seventh and eighth notes (Catalan parity at the Mersenne indices; the
Seidel doubling law and the graph tower), HYP-2220 (the repo's one aliquot
touchpoint: even perfect numbers as triangular rows `C(2^p, 2)`, the LRC
`14/21` shadow), the Kuratowski note's typing discipline. Classical facts
from memory unless stated: Catalan 1888 / Dickson 1913 (every aliquot
sequence terminates or cycles), Guy–Selfridge 1975 (drivers; the
counter-conjecture that many sequences diverge), Erdős 1973 (untouchable
numbers have positive density), Euclid–Euler.

**Status: PROVED elementary (Propositions 1–4, Theorem 5 for D21) +
FINITE-EXACT (the 276 sequence to 401 terms; the 2-adic transition matrices
for even `n <= 2·10^6` and for Collatz orbits below `2·10^5`; the 138
excursion; untouchables and perfect numbers by sieve; D21 to `n = 5`
exhaustive and random `n = 6, 7`; D22 to `|n| <= 60`) + CITED + ANALOGY and
NUMEROLOGY typed + DIRECTION. Collatz OPEN; the Catalan–Dickson question
OPEN. Independent audit OWED (as for the session's other notes).** Script
`04-computation/experiments/collatz_aliquot_lehmer_five_20260927.py`, output
beside it.

## 0. The answer in one paragraph

Collatz and the aliquot map are the same kind of machine, a multiplicative
step steered by the 2-adic part of the current number, with opposite
memory, and that one difference is why the two communities conjecture
opposite fates. Along a Collatz orbit the valuation `v_2(3m+1)` is
memoryless (Terras: geometric of ratio `1/2` whatever the previous
valuation; measured rows `0.52, 0.26, 0.10, ...` at every previous value),
and the only class with positive growth is `v = 1` (`+0.585` bits; every
other class descends). Along an aliquot sequence the valuation `v_2(s(n))`
is sticky (measured on even `n <= 2·10^6`: it repeats with probability
`0.90, 0.75, 0.57, 0.38` from `a = 1, 2, 3, 4`; along the sequence of 276 it
repeats in `91%` of the steps), and the growth classes are exactly the sticky
ones: the mean of `log_2(s(n)/n)` is `-0.35` bits in the class `a = 1` and
`+0.13, +0.32, +0.40, +0.45, +0.47` in the classes `a = 2, ..., 6`. The
mechanism is exact (Proposition 1): for `n = 2^a m` with `m` odd, `s(n) =
(2^(a+1) - 1) sigma(m) - 2^a m`, so the valuation `a` persists whenever
`v_2(sigma(m)) > a`, and `v_2(sigma(m))` is at least the number of prime
powers `p^e || m` with `e` odd, which grows with the size of `m`; the
driver locks in as the numbers grow. In Collatz the residue that produces a
growth pattern is consumed by the pattern (the fifth note's price sheet:
depth `K` buys `c_w K` bits and costs `K` bits), so growth is priced and
divergence is heuristically impossible; in the aliquot map the driver is
preserved by the algebra of `sigma`, growth is free, and divergence is
heuristically expected (Guy–Selfridge), with the Lehmer five as the smallest
open witnesses. The sequence of 276 (401 terms computed, 33 digits at the
end) shows the machine at work: driver `2^2·7` for most of the first 170
steps with mean growth `+0.25` bits per step, a fall from 28 to 15 digits
under the driver `2` (steps 176–224, growth `-1` bit per step), then a
recovery under drivers `2^3, 2^7, 2^4` and climb to 33 digits, i.e. a
completed excursion below the running maximum followed by regenerated
growth, which is the reset lane's configuration realized for free. The
cycles of the two maps are both "`2^k` minus something" problems (Euclid–
Euler's `2^p - 1` prime for perfect numbers; `2^A - 3^p | S_w` for Collatz
cycles), and the even perfect numbers sit at the same Mersenne indices
`2^p - 1` where the Catalan numbers are odd. The drivers `2^3·3 = 4!` and
`2^3·3·5 = 5!` are the Brocard factorials and `7!` is not a driver; that is
recorded as numerology. D21 and D22 are closed: `omega(D^2(G))` is the
largest four-part induced structure of Theorem 5, and Lehmer's conductor
polynomial meets the factorials only at the three Brocard cases.

## 1. The two engines (PROVED + FINITE-EXACT)

**Proposition 1 (the aliquot driver lock).** Let `n = 2^a m` with `m` odd,
`a >= 1`, and `s(n) = sigma(n) - n`. Then `s(n) = (2^(a+1) - 1) sigma(m) -
2^a m`, and
(a) if `v_2(sigma(m)) > a` then `v_2(s(n)) = a`; if `v_2(sigma(m)) < a` then
`v_2(s(n)) = v_2(sigma(m))`; if `v_2(sigma(m)) = a` then `v_2(s(n)) >= a + 1`;
(b) `v_2(sigma(m)) >= #{p^e || m : e odd}`, since `sigma(p^e) = 1 + p + ... +
p^e` is even for odd `p` and odd `e`;
(c) `s(n)/n = (2 - 2^(-a)) sigma(m)/m - 1`, whose mean over odd `m` is `(2 -
2^(-a)) (3/4) zeta(2) - 1 = 0.85, 1.16, 1.31, 1.39, ...` for `a = 1, 2, 3, 4`.

*Proof.* Multiplicativity of `sigma` and `sigma(2^a) = 2^(a+1) - 1`, which is
odd; (a) compares valuations of the two terms; (b) is the parity of
`sigma(p^e)`; (c) uses `E[sigma(m)/m] = zeta(2) prod_(p in S)(1 - p^(-2))` for
`m` coprime to `S = {2}`. ∎

**Proposition 2 (the Collatz engine has no lock).** Along a Syracuse orbit
the valuation `v_(l+1) = v_2(3 m_l + 1)` is determined by `m_l mod 2^(v+1)`,
and the map `n mod 2^k -> (v_1, ..., v_j)` (words with `d_j <= k`) is
Terras's bijection: the valuations of a uniformly random residue are
independent geometric variables of ratio `1/2`. The per-step growth is
`log_2 3 - v`: positive only for `v = 1`. (Fifth note, Theorem 1: a growth
pattern of depth `K` is available only to residues that encode it, and is
consumed by it.)

**FINITE-EXACT (script part 4).** Transition matrices `P(next valuation = b |
current = a)`:

| aliquot, even `n <= 2·10^6` | `b = 0` | 1 | 2 | 3 | 4 | 5 | `>= 6` | mean `log_2(s/n)` |
|---|---|---|---|---|---|---|---|---|
| `a = 1` | 0.001 | **0.895** | 0.052 | 0.026 | 0.013 | 0.007 | 0.007 | `-0.348` |
| `a = 2` | 0.001 | 0.111 | **0.747** | 0.070 | 0.035 | 0.018 | 0.017 | `+0.129` |
| `a = 3` | 0.002 | 0.118 | 0.145 | **0.566** | 0.085 | 0.042 | 0.042 | `+0.318` |
| `a = 4` | 0.003 | 0.125 | 0.150 | 0.172 | **0.382** | 0.084 | 0.084 | `+0.404` |
| `a = 5` | 0.004 | 0.134 | 0.154 | 0.176 | 0.169 | 0.235 | 0.128 | `+0.445` |

| Collatz, orbits of odd `n < 2·10^5` | `b = 1` | 2 | 3 | 4 | 5 | 6 | step growth |
|---|---|---|---|---|---|---|---|
| `a = 1` | 0.520 | 0.258 | 0.101 | 0.061 | 0.039 | 0.021 | `+0.585` |
| `a = 2` | 0.470 | 0.230 | 0.160 | 0.096 | 0.019 | 0.025 | `-0.415` |
| `a = 3` | 0.574 | 0.165 | 0.090 | 0.136 | 0.019 | 0.016 | `-1.415` |
| Terras | 0.5 | 0.25 | 0.125 | 0.0625 | 0.031 | 0.016 | — |

(the rare Collatz rows `a >= 4` are dominated by the small values near the
end of orbits and are not quoted). Reading: **Collatz grows only in the
class that cannot be held; the aliquot map grows only in the classes that
hold themselves.** This is the exact content of "drivers persist for free"
against "shadows are priced", and it is why Guy–Selfridge and the Collatz
conjecture point in opposite directions from the same 2-adic engine.

## 2. The sequence of 276 (FINITE-EXACT)

The Lehmer five and their factorizations: `276 = 2^2·3·23`, `552 = 2^3·3·23`,
`564 = 2^2·3·47`, `660 = 2^2·3·5·11`, `966 = 2·3·7·23`; `276` is itself
untouchable (no `m` has `s(m) = 276`), a leaf of its own graph. The
sequence of 276 was followed for 401 terms (8 seconds with sympy; the
literature has it beyond 200 digits): the 2-adic valuation is `2` on 219 of
the 400 steps, `1` on 60, `3` on 55, `4` on 39, higher on 27; it repeats
from one step to the next in `90.7%` of the steps; the mean growth is `0.247`
bits per step. The shape (script part 3, selected indices): driver `2^2·7`
(often with `3` or `5` alongside) through step ~170 with the number rising
to 28 digits; then the driver `2` from step ~176 to ~224 with growth `-1.0,
-1.0, -0.9, -0.87` bits per step (the number falls to 15 digits); then
drivers `2^3, 2^7·3^2, 2^4·3·11, 2^6·3^2, 2^4·5·7^2` with growth back to
`+1` bit per step, reaching 33 digits at step 400. The Guy–Selfridge drivers
in the standard list (`2, 2^2·7, 2^3·3, 2^3·3·5, 2^5·3·7, 2^9·3·11·31` and
the even perfect numbers) all pass the driver test `v | 2^(a+1) - 1`,
`2^(a-1) | sigma(v)` (script), and `2^3·3 = 24 = 4!`, `2^3·3·5 = 120 = 5!`
are among them while `7! = 2^4·315` is not (`315` does not divide `31`).
For scale: the sequence of 138 rises to `179,931,895,322` and returns to `1`
after 177 steps, an excursion exponent `log(peak)/log(138) = 5.26`, against
Collatz record excursions below `1.9`.

Reading against the reset lane. The fall of 276 from 28 to 15 digits under
the driver `2` and its recovery is a completed excursion below the running
maximum followed by regenerated growth; in Collatz such a configuration is
bought with a source of size `2^(precision)` (fifth note, Proposition 4),
here it is bought with nothing: the driver `2` was lost when the odd
cofactor acquired a prime power with `v_2(sigma) >= 2`, and the growth
drivers re-formed by Proposition 1(b) as the number's prime factors
accumulated. The same event that the Collatz thread has to price, the
aliquot map performs for free.

## 3. Cycles, leaves, indices (PROVED + CITED + ANALOGY typed)

**Proposition 3 (cycles as `2^k`-minus-something problems).** The even
perfect numbers are `2^(p-1)(2^p - 1)` with `2^p - 1` prime (Euclid–Euler),
i.e. the fixed points of `s` among even numbers are the cases where
`sigma(2^(p-1)) = 2^p - 1` is prime. The Collatz cycle points are `x_w =
S_w/(2^A - 3^p)`, integers exactly when `2^A - 3^p` divides the carry. Both
cycle families are governed by a power of two minus a small quantity being
of the right arithmetic kind, and both "no other cycles" questions (no odd
perfect number; no nontrivial Collatz cycle) are open.

**Proposition 4 (the Mersenne indices).** The even perfect numbers are the
triangular numbers `T_(2^p - 1)` (HYP-2220's observation, re-verified below
`4·10^5`: `6, 28, 496, 8128`), and the Catalan numbers `C_n` are odd exactly
at `n = 2^k - 1` (seventh note). The two families share the index set of
Mersenne numbers by two different mechanisms (Euclid–Euler needs `2^p - 1`
prime; Kummer needs nothing).

Leaves: the Syracuse tree's leaves are the multiples of 3 (density exactly
`1/3`); the aliquot graph's leaves are the untouchable numbers (`2, 5, 52,
88, 96, 120, 124, 146, ...`; 52 below 632 in the exact scan; Erdős: positive
density; `276` is one). Cycle types: the plus sheet of Collatz has the fixed
point `1`, the minus sheet has a fixed point `1`, a 2-cycle `{5, 7}` and a
7-cycle `{17, ..., 91}`; the aliquot map has fixed points (perfect
numbers), 2-cycles (amicable pairs, `220–284`) and longer cycles (sociable
numbers, `12496` of period 5, `14316` of period 28). The correspondence
fixed point / 2-cycle / longer cycle is a typing, ANALOGY. The Catalan of
the Catalan–Dickson conjecture is Eugène Catalan, the Catalan of the
numbers: onomastic only.

## 4. D21 and D22 closed (PROVED + FINITE-EXACT)

**Theorem 5 (D21: the doubled complete-split number).** For a graph `G`,
`omega(D^2(G)) = cs(D(G))` is the largest `|K| + |I| + |J| + |L|` over
pairwise disjoint sets with `K, L` cliques, `I, J` independent, `K–I`,
`K–J`, `K–L`, `I–J` complete and `I–L`, `J–L` anticomplete, or `|K| + |J| +
1` with `K` a clique completely joined to a nonempty independent set `J`
(the twin case `L = {u} ⊆ J`, `I = ∅`).

*Proof.* A complete split subgraph of `D(G)` has clique part `K` (plain) `∪
I'` (prime) with `K` a clique, `I` independent, `K–I` complete, `K ∩ I = ∅`,
and independent part `J` (plain) `∪ L'` (prime) with `J` independent, `L` a
clique and `J–L` anticomplete except at twins; the join conditions read
`K–J` complete, `K–L` complete, `I–J` complete, `I–L` anticomplete; a twin
`u ∈ J ∩ L` forces `L = {u}` (else `u` would be adjacent to the rest of `L`
against `J–L` anticomplete) and `I = ∅` (else `I–J` complete against `I–L`
anticomplete at `u`). Verified against brute-force `omega(D^2(G))` on all
labeled graphs with `2 <= n <= 5` (1098 graphs) and on random graphs with
`n = 6` (60) and `n = 7` (15). ∎

The seventh and eighth notes' towers are therefore governed by nested
compatible split structures; `cs(D(G)) < cs(G) + sep(G)` because the four
cross conditions are never all satisfiable at full size on the graphs
computed.

**D22 (closed).** For the Lehmer conductor polynomial `f(n) = n^4 + 5n^3 +
15n^2 + 25n + 25` and `|n| <= 60`, the only factorial values among `f(n) -
1`, `f(n)^2 - 1`, `f(n) + 1`, `f(n)^2 + 1` are `f(0) - 1 = 4!`, `f(-1)^2 - 1 =
f(-2)^2 - 1 = 5!`, `f(1)^2 - 1 = 7!`: the three Brocard cases (with the
duplicate `f(-2) = f(-1) = 11`) and nothing else.

## 5. What is useful for Collatz (DIRECTION; typed)

* **D24 (exact contrast as a lemma).** The statement "Collatz valuations
  are memoryless and its only growth class is the unstable one; aliquot
  valuations are sticky and its growth classes are the stable ones" is a
  precise reason why the two heuristics diverge, and it sharpens the fifth
  note's verdict: any Collatz variant in which the valuation persisted
  (a sticky `v = 1`) would diverge at once (the `-1` shadow for ever), so
  Terras's memorylessness is the whole of what keeps the `-1` shadow finite.
  A proof of the divergence half would have to use memorylessness *for
  every integer*, which is what the residue classes give on average and
  never pointwise.
* **D25 (the aliquot side).** Conversely, a proof of Guy–Selfridge for a
  single sequence would have to show a driver persists for ever, which by
  Proposition 1 is the statement that `v_2(sigma(m))` stays above `a` along
  the odd cofactors: a pointwise statement about prime factorizations of a
  specific sequence, as hard as the Collatz pointwise statement in the
  opposite direction. The two problems are dual obstructions: prove a
  pattern persists (aliquot) versus prove every pattern breaks (Collatz).
* **D26 (untouchables against multiples of 3).** The exact leaf density
  `1/3` of the Syracuse tree is the Collatz analogue of the untouchable
  density (Pollack–Pomerance's heuristic `~0.17`); the tree of 1 must
  contain every leaf's ancestors in both graphs. Whether the aliquot
  in-degree formula (finite except at `1`) has a Collatz analogue for the
  descent tree of the parallel note (in-degree `c_D ≈ 1.67`) is a cheap
  comparison not made here.

## 6. Reproduction

    cd 04-computation/experiments
    python3 collatz_aliquot_lehmer_five_20260927.py > collatz_aliquot_lehmer_five_20260927.out

Needs sympy for the factorizations of the 276 sequence; about ten minutes
(the divisor-sum sieves to `2·10^6` and `4·10^5`).
