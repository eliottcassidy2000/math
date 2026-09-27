# Pricing the approaches (D12): the adelic price sheet of a cycle shadow, the three-place multiplier identity, explicit double-excursion sources that stay above the source, and where the growth of real orbits actually comes from

**Session:** opus, `collatz-posets-zeta5-20260927` (fifth note), 2026-09-27.
**Owner's directive:** "keep going, pursue D12 pricing the approaches toward a
Collatz proof" (D12 of the fourth note: a potential `V(m) = log m + c(m)`
with `c` a function of the 2-adic distances to the negative cycle points,
against the reset lane's regenerated precision).
**Inherits (cited):** THM-4512 (exact-word classes mod `2^(A+1)`, the affine
map `Q_w`), the S12 note (growth families are 2-adic shadows of negative
cycles), the S13 note (Proposition 4: no prefix rank; the reset lane's
obstruction "two complete excursions, both above the source, arbitrarily
large combined growth, regenerated precision around the same negative
cycle"), the parallel spine note (2-adic source classes to 3-adic landing
classes; section 6.1), THM-4476 (reciprocal sums of divergent orbits
converge: the real-place price of approaches to small values), the fourth
note (Propositions 4–5: the trivial cycle's exact linearization; cycles
attract or repel their shadows by sign), the HARD-class note (the
supercritical positive-entropy words untouched by every mechanism).

**Status: PROVED elementary (Theorem 1, Propositions 2–4) + FINITE-EXACT
(regeneration depths and growth attribution for all odd `n < 2·10^5` and
the record orbits; five explicit double-excursion sources with exact
ledgers) + DIRECTION. Collatz OPEN. D12 is closed as a proof route: the
price of one approach is exact and polynomial in the current value, the
lane's obstruction is reproduced with explicit integers and shown to be
bought with a source of size `2^(precision)`, and a Lyapunov function would
have to price infinitely many future regenerations, which is the height
walk again; the useful by-product is that deep cycle shadows carry about one
percent of the growth of real orbits.** Script
`04-computation/experiments/collatz_pricing_approaches_20260927.py`, output
beside it. Independent audit OWED (as for the session's other notes).

Notation: Syracuse map `U(m) = (3m+1)/2^v`; a cycle word `w = (v_1, ...,
v_p)` with `A = sum v_i`, carry `S_w = sum_(t<p) 3^(p-1-t) 2^(v_1+...+v_t)`,
rational cycle point `x_w = S_w/(2^A - 3^p)` (`U^p(x_w) = x_w` along `w`),
multiplier `mu_w = 3^p/2^A`. The integer cycle points of the plus sheet are
`x = -1` (`w = (1)`), `-5` (`(1,2)`), `-17` (`(1,1,1,2,1,1,4)`) and `1`
(`(2)`).

## 0. The answer in one paragraph

Every rational cycle is exactly linear on its 2-adic neighbourhood: for `m
- x_w = 2^K u` with `u` odd the orbit follows `w` for `r = floor((K-1)/A)`
periods and `U^(pr)(m) - x_w = mu_w^r (m - x_w) = 3^(pr) 2^(K - Ar) u`
(Theorem 1). So an approach of depth `K` to a cycle point buys growth at
the fixed rate `c_w = (p log_2 3 - A)/A` bits per bit of precision (`0.585`
for `-1`, the universal maximum; `0.0566` for `-5`; `0.0086` for `-17`;
`-0.2075` for the trivial cycle), converts each spent bit of 2-adic precision
into ternary digits of 3-adic proximity, and, because a positive integer
congruent to a negative `x` modulo `2^K` is at least `2^K - |x|`, costs size:
`K <= log_2(m + |x|)`, so one shadow multiplies the value by at most `(m +
|x|)^(c_w)`. The per-step multiplier `3/2^v` has absolute values `3/2^v`, `2^v`
and `1/3` at the real, 2-adic and 3-adic places, product one (Proposition
2): real growth is exactly the excess of 3-adic contraction over 2-adic
expansion, and the 3-adic metric is the contracting metric the hedgehog
analogy asked for, uniformly by `1/3` per odd step within a valuation class
and useless across a change of word. The reset lane's obstruction is
reproduced with explicit integers (Proposition 4 and the ledgers): the word
`w^(r_1) + (2,2) + w^(r_2)` has a least source `n` of `75`–`174` bits for
`-5` and `35`–`55` bits for `-1`, its second approach is deeper than the
first (`K_2 = 2K_1 + 1`), the dip between the two excursions is smaller than
the second growth, the whole orbit stays above `n`, and the combined growth
is `0.043`–`0.047 log_2 n` for `-5` and `0.485`–`0.49 log_2 n` for `-1`:
"arbitrarily large combined growth" means arbitrarily large source. Along
real orbits (all odd `n < 2·10^5` and the five record orbits) about half of
all positive growth happens inside runs of two or more ones (shallow
shadows of `-1`), an eighth inside shallow shadows of `-5`, and about one
percent inside deep shadows of any cycle: the growth that a cycle label
could price is the shallow, unavoidable part, and the rest is generic
no-descent words. A potential that prices approaches must therefore price
every future regeneration, each bounded by `c_w log_2` of the value at that
time, and the sum of those prices along a divergent orbit is the height
walk it was meant to control.

## 1. The price sheet (PROVED)

**Theorem 1.** Let `w` be a cycle word with cycle point `x_w` and
multiplier `mu_w`, and let `m` be an integer with `m - x_w = 2^K u`, `u`
odd, `K >= A + 1`. Then:
(a) the valuation word of `m` begins with `w^r`, `r = floor((K-1)/A)`, and
`U^(pr)(m) - x_w = mu_w^r (m - x_w) = 3^(pr) 2^(K - Ar) u`;
(b) over the `r` periods the deviation from the cycle point loses exactly
`Ar` bits of 2-adic valuation and gains exactly `pr` ternary digits of
3-adic valuation;
(c) `log_2 U^(pr)(m) - log_2 m <= c_w K + O_w(1)` with `c_w = (p log_2 3 -
A)/A`, and for every word of `L` odd steps with `d` halvings the growth `L
log_2 3 - d` is at most `(log_2 3 - 1) d`, with equality only for the
all-ones word; so `c_(-1) = log_2 3 - 1 = 0.585` is the universal maximum
growth per bit of precision;
(d) if `x_w` is a negative integer and `m > 0`, then `K <= log_2(m + |x_w|)`,
hence `U^(pr)(m) <= C_w (m + |x_w|)^(1 + c_w)`.

*Proof.* (a) `U^p` is the affine map `y -> (3^p y + S_w)/2^A` on the exact
class of `w` (THM-4512), whose fixed point is `x_w`; the class condition is
`m = x_w mod 2^(A+1)` and persists for `r` periods while `K - A(r-1) >= A +
1`. (b) is (a) read at the two places. (c) `mu_w^r = 2^(r(p log_2 3 - A))` and
`r <= K/A`; the universal bound is `L log_2 3 - d <= (log_2 3 - 1) d` since `d
>= L`. (d) `m = x_w mod 2^K` with `m > 0 > x_w` forces `m >= 2^K - |x_w|`.
Checked exactly for `-1, -5, -17, 1`, `K <= 60`, `u in {1, 3, 7}` (script
part 1). ∎

| cycle | word | `x_w` | `mu_w` | `c_w` (bits of growth per bit of precision) |
|---|---|---|---|---|
| `-1` | `(1)` | `-1` | `3/2` | `+0.585` (the universal maximum) |
| `-5` | `(1,2)` | `-5` | `9/8` | `+0.0566` |
| `-17` | `(1,1,1,2,1,1,4)` | `-17` | `2187/2048` | `+0.0086` |
| trivial | `(2)` | `1` | `3/4` | `-0.2075` (descent) |

The place exchange in (a) is the parallel note's Proposition 4 for the
periodic words: the 2-adic source class `x_w + 2^K Z` lands, after `r`
periods, in `x_w + 3^(pr) Z`. For the trivial cycle this was the fourth
note's `T^(2j)(1 + 4^j t) = 1 + 3^j t`.

**Proposition 2 (the three-place multiplier).** For an odd step with
valuation `v`, the multiplier `3/2^v` has `|3/2^v|_oo = 3/2^v`, `|3/2^v|_2 =
2^v`, `|3/2^v|_3 = 1/3`, product one. Along a shadow of `r` periods the
real, 2-adic and 3-adic distances of two same-word points are multiplied by
`3^(pr)/2^(Ar)`, `2^(Ar)` and `3^(-pr)` respectively (checked on same-word
pairs, script part 2). In particular `U` is a 3-adic contraction by `1/3`
per odd step within a valuation class, and real growth over any segment is
exactly the excess of the 3-adic contraction over the 2-adic expansion:
`3^L/2^(d_L) = (3^(-L))^(-1) (2^(d_L))^(-1)`.

*Proof.* Product formula for the rational number `3/2^v`; within a
valuation class `U(m) - U(m') = (3/2^v)(m - m')`. ∎

This is the exact content of the hedgehog analogy's "hyperbolic metric":
the 3-adic metric is the one in which the map contracts uniformly, and the
contraction is worthless exactly where the words differ, which is where
the dynamics is.

## 2. Regeneration, measured (FINITE-EXACT)

Along every orbit of odd `n < 2·10^5` (script part 3) the depth `v_2(m_j -
x)` of the approach to `x in {-1, -5, -17}` never exceeds `log_2(m_j + |x|)`
(the size bound of Theorem 1(d), a theorem), while depths exceeding
`log_2 n` do occur (maximal excess `1.25`, `3.25`, `2.51` bits for the three
points; maximal depths `17, 18, 18` at `n = 77671, 174759, 174751`): precision
is regenerated, but never beyond the size of the value carrying it. The
tail frequencies `P(depth >= k)` for `k = 2, 4, ..., 12` are `0.500, 0.124,
0.027, 0.0040, 0.00095, 0.00021` (`-1`), `0.500, 0.129, 0.027, 0.0118,
0.00082, 0.00012` (`-5`), `0.500, 0.124, 0.030, 0.0049, 0.00091, 0.00020`
(`-17`), against the uniform law `2^(1-k)` for odd integers `0.5, 0.125,
0.031, 0.0078, 0.0020, 0.00049`: deep approaches along orbits are somewhat
rarer than uniform for `-1` and `-17` and show an excess at depth 8 for
`-5`; the note records the numbers and draws no conclusion from them.

## 3. Where the growth of real orbits comes from (FINITE-EXACT)

Attribute each positive height increment `log_2 3 - v` of an orbit to the
cycle point `x` whose shadow the current value is in, at a depth threshold
(script part 4). For all odd `n < 2·10^5` (total `1.18·10^6` bits of
positive growth):

| threshold | inside `-1` shadows | inside `-5` shadows | inside `-17` shadows |
|---|---|---|---|
| shallow (`depth >= A + 2`: `3, 5, 13`) | `52.0%` | `13.4%` | `0.00%` |
| deep (`8, 10, 16`) | `0.8%` | `0.2%` | `0.00%` |

The record orbits `27, 703, 6171, 77031, 837799` have shallow `-1` shares
`58, 60, 56, 63, 62%` and deep shares of at most `3.9%`. Reading: a shallow
`-1` shadow of depth 3 is a run of two ones, i.e. an ordinary pair of
consecutive `v = 1` steps; half of all growth lives in such runs, which is
the generic frequency of ones, not an approach to a cycle in any useful
sense. The growth that a cycle label could price, the deep approaches,
is about one percent. The rest, and the bulk of every excursion, is the
generic no-descent word: the positive-entropy part of `E_inf`, the
HARD-class words no mechanism reaches.

## 4. The reset lane's obstruction, priced (FINITE-EXACT + PROVED)

**Proposition 3 (explicit double excursions above the source).** For the
words `w^(r_1) + (2,2) + w^(r_2)` and `w^(r_1) + (2,2,2) + w^(r_2)` with
`w` the word of `-5` or of `-1`, the least positive member `n` of the exact
class (THM-4512, class mod `2^(A_total + 1)`) has the following ledgers
(script part 5; all verified by running the orbit):

| cycle | `(r_1, bridge, r_2)` | bits of `n` | `K_1` | growth 1 | dip | `K_2` | growth 2 | orbit above `n` | total growth / `log_2 n` |
|---|---|---|---|---|---|---|---|---|---|
| `-5` | `(8, (2,2), 16)` | `75.1` | `25` | `1.36` | `-0.83` | `51` | `2.72` | yes | `0.043` |
| `-5` | `(12, (2,2,2), 24)` | `114.7` | `37` | `2.04` | `-1.25` | `74` | `4.08` | yes | `0.042` |
| `-5` | `(16, (2,2,2), 40)` | `174.3` | `49` | `2.72` | `-1.25` | `123` | `6.80` | yes | `0.047` |
| `-1` | `(10, (2,2), 20)` | `34.5` | `11` | `5.85` | `-0.83` | `23` | `11.70` | yes | `0.485` |
| `-1` | `(16, (2,2,2), 32)` | `54.8` | `17` | `9.36` | `-1.25` | `34` | `18.72` | yes | `0.490` |

In each case the second approach is deeper than the first, the dip between
the excursions is smaller than the second growth, and every orbit value
stays above the source: exactly the lane's configuration. So no potential
of the form `log m + (prepaid growth of the current shadow)`, nor any
potential depending on the current 2-adic distances to the cycle points
(bounded lookahead), is non-increasing along these orbits: at the re-entry
it jumps by `c_w K_2`.

**Proposition 4 (source price).** In any such configuration with the second
entry value not above the first exit value, `K_1 <= log_2(n + |x|)`, `K_2 <=
(1 + c_w) log_2(n + |x|) + O(1)`, and the combined growth of the two
excursions is at most `c_w (2 + c_w) log_2(n + |x|) + O(1)`. More generally,
along any orbit every shadow event of depth `K` at value `m` satisfies `K <=
log_2(m + |x|)` and contributes growth at most `c_w log_2(m + |x|)`.

*Proof.* Theorem 1(c),(d) applied at each entry; the first exit value is at
most `(n + |x|)^(1 + c_w) C_w`. ∎

So "arbitrarily large combined growth" is a statement about arbitrarily
large sources: the growth is a fixed small fraction (`0.043` for `-5`,
`0.49` for `-1`) of `log_2 n`, and it is paid for by the bits of `n` that
encode the two approaches. This does not rescue a bounded-lookahead
potential (Proposition 3 refutes those), but it locates the lane's
obstruction precisely: the obstruction is to potentials that do not read
the source; a potential that reads `log_2 n` prices any finite number of
excursions polynomially.

## 5. Why pricing stops here

A Lyapunov function `V` for the divergence half must decrease across
every completed excursion of every orbit. Pricing gives, for each future
regeneration event `i` at value `m_i`, a prepaid growth of at most `c_(w_i)
log_2 m_i`; a potential that prepays all of them is `V(m) = log_2 m + sum_i
c_(w_i) log_2 m_i` over the orbit's future shadow events, which is finite
iff the orbit's future is, and it is the total stopping time in disguise.
Between the events, the growth is the generic no-descent words (section 3),
which no cycle label prices at all. And the one exact contraction available
(Proposition 2, the 3-adic metric) contracts only within a valuation class,
so it prices the shadows and nothing else. D12 therefore reduces to the
same transversality statement as the rank routes: the sum of the dips must
exceed the sum of the prepaid growths, and that is the height walk's
return.

What survives as tools: the price sheet (Theorem 1) with the universal
rate `log_2 3 - 1` and the per-cycle rates; the polynomial bound per shadow
and per pair of excursions (Proposition 4); the explicit sources of
Proposition 3, which are the smallest integers realizing the lane's
configuration for those words and can serve as hostiles for any proposed
rank; and the attribution table, which says that any argument built on
cycle labels addresses at most a few percent of the growth budget of real
orbits.

## 6. Directions (DIRECTION; none pursued)

* **D14. Price the generic words.** The bulk of growth is in no-descent
  words that shadow no cycle. The same accounting applies to the fixed
  points `x_w` of *every* no-descent word (`E_inf` is their closure, S13):
  an orbit following a no-descent word `w` of length `L` with `d` halvings
  is in the `2^(-(d+1))`-shadow of the 2-adic point `x_w`, and grows by `L
  log_2 3 - d <= 0.585 d` bits. So every excursion is a shadow of some point
  of `E_inf`, and its price is the same: at most `0.585` bits per bit of
  precision, at most `0.585 log_2(m + |x_w|)` in total only when `x_w` is a
  negative integer. For the rational and irrational points of `E_inf` the
  size bound fails (a positive integer can be 2-adically close to a
  non-integer 2-adic point without being large), which is exactly why the
  generic words are not priced by size. **Corrected in the sixth note:** for
  every fixed 2-adic point `x` the depth-`K` approach class has least positive
  member `rho_K(x) = x mod 2^K`, so the price is `K - z_K(x)` with `z_K` the
  zero run of `x`'s bits just below position `K`, a geometric fluctuation for
  generic points of `E_inf`; a small integer approaches deeply only the
  points that look like it, so the obstruction is the tautology that an
  integer's cheap growth is its own word, not a lack of size price.
* **D15. The `-1` rate as the worst case.** Since every word's growth is at
  most `0.585 d`, and the halvings `d` are the precision consumed, the only
  way to make growth cheap is to consume precision without ones, i.e. to
  descend. A potential of the form `log_2 m - 0.585 v_2(m + 1)` prepays the
  worst case exactly and is refuted by the `-1` ledgers above at the
  re-entry; whether a *scaled* prepayment `log_2 m - lambda v_2(m+1)` with
  `lambda < 0.585` has a smaller violation set on real orbits is a cheap
  computation not done here.

## 7. Reproduction

    cd 04-computation/experiments
    python3 collatz_pricing_approaches_20260927.py > collatz_pricing_approaches_20260927.out

Standard library only; about ten minutes (the attribution pass over `2·10^5`
orbits twice).
