# The shadow error of a Collatz orbit is its 3x−1 copy: an exact calculus, the Flatto–Lagarias–Pollington carry structure with adaptive halvings, the bounded-shadow interval map, and what a spreading theorem would have to say

**Session:** opus, `collatz-shadow-flp-20260927`, 2026-09-27.
**Owner's directive:** "pursue the Flatto–Lagarias–Pollington digit spreading
angle and other related ones as they appear."
**Inherits:** [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)
(thin divergence; the Bernstein series `R(d) = sum_l 2^(d_l)/3^(l+1)` converges
for divergent orbits; the two-places corollary `R(d) < n` on the minus sheet),
the S13 note [`collatz_directions_20260926.md`](collatz_directions_20260926.md)
(the no-descent fractal, the Mahler sibling), the Mahler frontier
([`MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md`](../../00-navigation/MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md)),
Mahler 1968 (Z-numbers, the map `g -> ceil(3g/2)`), Flatto–Lagarias–Pollington
1995 (the fractional parts `{xi (p/q)^n}` have spread at least `1/p`), Dubickas
and Bugeaud's book (Chapter 3) for the surrounding literature, all cited from
memory and not re-derived here.

**Status: PROVED elementary propositions (the shadow-error calculus, sections
1–3) + FINITE-EXACT (checked on the truncated tails of `27, 703, 871, 6171`,
on the rational cycle points, on the `−5` shadow families, on both sheets) +
HYP-9164 (the bounded-shadow conjecture, an intermediate pointwise statement)
+ DIRECTION marked (what a spreading theorem would be). Collatz OPEN; nothing
here is a proof of anything about all `n`.** Scripts:
`04-computation/experiments/collatz_shadow_error_20260927.py` and
`collatz_shadow_error_20260927_itineraries.py`, with `.out` files.

## 0. The answer in one paragraph

A divergent orbit `m_l` of `U(m) = (3m+1)/2^v` has a real shadow: with
`d_l = v_1 + ... + v_l`, `C_l = prod_(j<l) (1 + 1/(3 m_j))` and `xi = n C_inf`
(finite by THM-4476), `xi 3^l / 2^(d_l) = m_l + eta_l` where
`eta_l = sum_(k>=0) 2^(d_(l+k) - d_l)/3^(k+1)` is the real Bernstein value of
the tail word. The error obeys `2^(v_(l+1)) eta_(l+1) = 3 eta_l − 1`: **it is
the 3x−1 copy of the orbit, driven by the orbit's own halvings.** It is always
`>= 1`, a large next valuation needs a large error (`2^(v_(l+1)) <= 3 eta_l − 1`),
the orbit can never dip below `m_l/(3 eta_l)` after time `l`, on a periodic tail
the error equals `|x_w|` (the two copies cancel exactly on the negative cycles:
`m + eta = 0`), and on the minus sheet the same recursion holds with
`m − eta`. The integer parts `M_l = floor(xi 3^l/2^(d_l))` satisfy
`2^v M_(l+1) = 3 M_l + c_l` with carries `c_l in {−2^v + 1, ..., 2}` determined by
the fractional parts: this is exactly the carry structure of
Flatto–Lagarias–Pollington, with the rigid shift `q` replaced by the adaptive
`2^v`. FLP's theorem refutes confinement of `{xi (3/2)^n}` to an interval of
length `< 1/3`; the Collatz shadow has no confinement hypothesis to refute,
because its fractional parts `{eta_l}` are unconstrained and its integer part
`floor(eta_l)` is a free 3x−1 orbit. So the angle does not transfer as a
theorem; what it yields is the right object and the right intermediate
question. Bounded shadow error forces bounded valuations and turns the word
into the itinerary of an explicit interval map (`eta < 2`: valuations `<= 2`,
`f(eta) = (3 eta − 1)/2` on `[1, 5/3)`, `(3 eta − 1)/4` on `[5/3, 2)`); the
bounded-shadow conjecture (HYP-9164: no integer orbit has `sup eta_l < inf`)
is a pointwise statement strictly between Mahler's rigid-cut problem and the
full conjecture, and the trivial spreading bound `sup eta_l >= 5/3` is sharp
over words.

## 1. The shadow-error calculus (PROVED)

Notation: odd `m_0 = n`, `m_(l+1) = (3 m_l + 1)/2^(v_(l+1))`, `d_l = sum v`,
`S_(l-1) = sum_(t<l) 3^(l-1-t) 2^(d_t)`, so `2^(d_l) m_l = 3^l n + S_(l-1)`. For
a divergent orbit `C_l -> C_inf < inf` (THM-4476) and the real series
`R(d) = sum_(t>=0) 2^(d_t)/3^(t+1) = n (C_inf − 1)` converges.

**Proposition 1 (definition and recursion).** Put
`eta_l = sum_(k>=0) 2^(d_(l+k) − d_l)/3^(k+1)`, the real Bernstein value of the
tail word `sigma^l d`. Then `eta_l = m_l (C_inf/C_l − 1)`, `xi 3^l/2^(d_l) = m_l +
eta_l` with `xi = n C_inf = n + eta_0`, and
```text
2^(v_(l+1)) eta_(l+1) = 3 eta_l − 1 .
```
*Proof.* Splitting off `k = 0`: `eta_l = 1/3 + (2^(v_(l+1))/3) eta_(l+1)`, which is
the recursion; `m_l(C_inf/C_l − 1) = m_l (prod_(j>=l)(1 + 1/(3 m_j)) − 1)` obeys
the same recursion (expand one factor) and both tend to `0` relative to `m_l`,
so they agree (uniqueness below). Adding `2^v m_(l+1) = 3 m_l + 1` gives
`2^v x_(l+1) = 3 x_l` for `x_l = m_l + eta_l`, hence `x_l = x_0 3^l/2^(d_l)` and
`x_0 = n + eta_0 = n C_inf = xi`. ∎ (Uniqueness: two sequences obeying
`2^v y' = 3y − 1` differ by `y − y'` multiplied by `3^k/2^(d_(l+k) − d_l)`, which
tends to infinity along a divergent orbit, so a bounded-ratio solution is
unique.)

**Proposition 2 (the error is at least one).** `eta_l >= sum_k 2^k/3^(k+1) = 1`,
with equality iff every later valuation is `1`; for an integer orbit
`eta_l > 1`. The largest value over no-descent words of length `12` is `2.95`
(words hugging the critical line), and it grows like a third of the length of
a critical stretch.

**Proposition 3 (budget and dip bounds).** `2^(v_(l+1)) <= 3 eta_l − 1`; more
generally `2^(d_(l+k) − d_l) <= 3^(k+1) eta_l` for every `k >= 0`, hence
```text
m_(l+k) >= m_l / (3 eta_l)      for every k >= 0 .
```
*Proof.* Each term of the series is at most the sum; `2^(d_(l+k)) m_(l+k) >= 3^k
2^(d_l) m_l`. ∎ So a valuation `v` at time `l+1` needs `eta_l >= (2^v + 1)/3`, and
the future of the orbit never falls below the current value divided by
`3 eta_l`.

**Proposition 4 (periodic tails and the two copies).** If the tail from `l` is
`w^inf` with `3^p > 2^A`, then `eta_l = |x_w| = S_w/(3^p − 2^A)` exactly (the
geometric series), the real and 2-adic values of the tail agree, and on the
negative cycles `m + eta = 0`: the orbit and its 3x−1 copy cancel and the real
shadow is `0`. Checked for `(1), (1,2), (2,1)`, the two rotations of the `−17`
word, `(1,1,2)` (`19/11`) and `(1,2,1,2,2,1)` (`1213/217`), and on all points of
the three negative cycles. For a positive divergent orbit the shadow is
`m_l + eta_l > m_l + 1 > 0`: the copies do not cancel.

**Proposition 5 (minus sheet).** For `U^−(m) = (3m − 1)/2^v` the same
construction gives `xi^− 3^l/2^(d_l) = m_l − eta^−_l`, `eta^−_l` the real value of
the tail, with the same recursion `2^v eta^−_(l+1) = 3 eta^−_l − 1`; on the
3x−1 cycles `1`, `5 -> 7`, `17 -> 25 -> 37 -> 55 -> 41 -> 61 -> 91` one has
`m − eta^− = 0` (checked), and for orbits reaching the fixed point `1` the real
value equals the 2-adic value `n` (checked for `3, 11, 29`), consistent with
THM-4476's strict inequality `R(d) < n` being reserved for non-eventually-
periodic words. In THM-4476's language: on the minus sheet `R(d) <= R_2(d)`
with equality iff the word is eventually periodic; on the plus sheet
`R(d) > 0 > R_2(d)` trivially. The conjecture "no divergent orbit on either
sheet" reads: **a halving word whose real Bernstein series converges and whose
2-adic value is a nonzero integer is eventually periodic.**

**Checks on finite orbits.** For `27, 703, 871, 6171` (`41, 62, 65, 96` odd
steps) the truncated tails `eta_l^(K)` (partial sums to the last odd step
before the `1`-cycle) satisfy the recursion exactly, the sharp budget for
tails of at least `30` terms, the dip bound inside the truncation, and the
carry structure of section 2; the maximal truncated errors are `332`
(`27`, at `m = 3077`), `17293`, `9909`, `61981`: the error is large exactly
where the future falls far below the current point. Along the `−5` shadow
families `n_k = 4 8^k − 5` the partial error over the shadow is `5(1 − (8/9)^k)`
exactly, tending to `|−5| = 5`.

## 2. The Flatto–Lagarias–Pollington carry structure (PROVED), and why the theorem does not transfer (verdict)

Let `M_l = floor(x_l) = m_l + floor(eta_l)`. Then
```text
2^(v_(l+1)) M_(l+1) = 3 M_l + c_l,      c_l = 3 {eta_l} − 2^(v_(l+1)) {eta_(l+1)}  in  {−2^v + 1, ..., 2},
```
an integer carry determined by the fractional parts of the shadow error
(checked: carries in `{−24, ..., 2}` on the four orbits, always within the
bounds). In FLP the sequence is `g_n = floor(xi (p/q)^n)`, `q g_(n+1) = p g_n +
c_n`, the shift `q` is rigid, the carry ranges over `p` values, and the theorem
is that the fractional parts `{xi (p/q)^n}` cannot all lie in an interval of
length `< 1/p` (the carries would be too constrained). Mahler's Z-numbers are
the case `{.} < 1/2`, which FLP's bound `1/3` does not reach. In Collatz:

* the orbit itself has the rigid carry `+1` and the adaptive shift `2^v`;
* the integer parts of the shadow have the FLP carry structure with adaptive
  shift and carry range growing with `v`;
* nothing confines `{eta_l}`: the only known constraints are `eta_l >= 1` and
  the budget of Proposition 3, both one-sided;
* FLP's own theorem applied to `xi (3/2)^l = 2^(e_l)(m_l + eta_l)`, `e_l = d_l − l`,
  says the fractional parts of `2^(e_l) eta_l` spread over `1/3`, which the
  adaptive exponents satisfy for free.

Verdict: the FLP hypothesis has no Collatz counterpart and its conclusion is
vacuous here; the transferable content is the object (a real shadow whose
integer parts follow a `3x + carry` recursion with the orbit's halvings) and
the recognition that Mahler's problem is the rigid-cut case. In the
coordinates `N = (m + 1)/2` the Syracuse step reads
`N -> (3N − 1 + 2^(v−1))/2^v`: for `v = 1` it is Mahler's even branch `3N/2`,
for `v = 2` it is `(3N + 1)/4`, one halving more than Mahler's odd branch
`(3N + 1)/2`; the two problems share the first branch and differ by the
adaptive halvings.

## 3. The bounded-shadow regime (PROVED small statements) and HYP-9164

**Proposition 6.** If `eta_l < B` for all `l` then every valuation is at most
`log_2(3B − 1)` (`B = 2`: `v <= 2`; `B = 5`: `v <= 3`; `B = 17`: `v <= 5`) and
`m_(l+k) >= m_l/(3B)` for all `l, k` (no deep dips). For `B = 2` the word is
the itinerary of the map
```text
f(eta) = (3 eta − 1)/2  on [1, 5/3),      f(eta) = (3 eta − 1)/4  on [5/3, 2),
```
which maps `[1, 2)` into itself (`v = 2` is forced exactly when `eta >= 5/3`,
since `v = 1` would give `eta' >= 2` and `v = 2` below `5/3` would give
`eta' < 1`). *Proof.* From `eta' = (3 eta − 1)/2^v >= 1` and `eta' < 2`. ∎ The
primitive periodic words with `p <= 10` whose error stays in `[1, 2)` are eight:
`(1)`, `(1^a 2)` for `a = 4..9` (errors `211/179, 665/601, 2059/1931, ...`,
maxima `1.905, 1.809, ..., 1.691`) and `(1^5 2 1^3 2)` (max `1.9994`); in each
case the itinerary of `f` is the word. Their 2-adic points are rational and
not integers. The maxima of the family `(1^a 2)` decrease to `5/3`, so

**Proposition 7 (runs in the `B = 2` regime).** Along an `f`-itinerary every
`v = 2` step is followed by at least three `v = 1` steps: after `v = 2` the error
lies in `[1, 5/4)`, and `g(eta) = (3 eta − 1)/2` maps `[1, 5/4)` to `[1, 11/8)`,
then to `[1, 25/16)`, both below `5/3`, so the third image `[1, 1.84)` is the
first that can reach `5/3`. Hence the frequency of `v = 2` is at most `1/4`,
every block `(1, 1, 1, 2)` multiplies the value by at least `81/32`, and a
bounded-shadow (`B = 2`) divergent orbit would grow at least like
`2^(0.335 l)`: bounded shadow error is a fast-divergence regime, whose
starting classes are thin (`(1 − h(rho))` with `rho >= 0.8`). Exploration
(`collatz_shadow_error_20260927_itineraries.py`): for all `490` rational
starting errors in `[1, 2)` with denominators up to `40`, the exact itineraries
of length `80` stay in `[1, 2)`, are growth words (frequency of `v = 2` between
`0` and `0.188`, mean `0.152`), reproduce the starting error as the limit of
their Bernstein partial sums (deviation `< 4e-10`, the self-consistency of
Proposition 1), and none has a 2-adic point that stabilises to an integer
over the second half of the window.

**Trivial spreading bound.** Every integer orbit with convergent series has
`sup_l eta_l >= 5/3` (all errors below `5/3` force every valuation to be `1`, the
word `1^inf`, the point `−1`), and `5/3` is the infimum over words.

**HYP-9164 (bounded shadow).** No positive integer orbit, on either sheet, has
`sup_l eta_l < inf`; equivalently, every divergent orbit has unbounded shadow
error, i.e. unbounded valuations or unboundedly deep dips relative to the
current point. The sub-case "no divergent orbit has bounded valuations" is the
Mahler-shaped case: Mahler's Z-numbers are the rigid-cut problem (`v = 1` for
ever in the `N`-coordinates of section 2), bounded valuations are the
bounded-cut problem, and the conjecture is the adaptive-cut problem. The
hypothesis is implied by the conjecture, is pointwise, and is the statement
an FLP-type argument would have to prove first.

## 4. Two completions (CITED + restated)

Along an orbit the pair `(m_l, eta_l)` is `(∓ R_2, R)` of the tail word: the
2-adic value (a positive integer on the plus sheet up to sign, the orbit point)
and the real value (the error). The real shadow `xi 3^l/2^(d_l) = m_l ± eta_l`
is the discrepancy between the two completions of the same rational partial
sums, `0` exactly on cycles and positive on divergent orbits. The pointwise
question is whether a word can have a positive integer as one completion and
a finite positive real as the other without being periodic. THM-4476's
`R(d) < n` is the only inequality between the two completions proved so far.

## 5. What a spreading theorem would have to say (DIRECTION)

1. **Target statement.** For every odd `n >= 3` whose real series converges,
   `sup_l eta_l = inf` (HYP-9164). A quantitative version: for every `B` there
   is `L(B, n)` with `eta_l >= B` for some `l <= L`. The trivial case is
   `B = 5/3`.
2. **The mechanism to look for.** An orbit with `eta_l < B` for all `l` has its
   word coded by an interval map with `floor(log_2(3B − 1))` branches (the
   generalisation of `f`), and its orbit points are `m_l = M_l − floor(eta_l)`
   with `M_l` following the carry recursion. Both are deterministic given
   `eta_0`; the integrality of `n` is the extra condition. An FLP-type
   argument would count, for each `L`, the pairs (itinerary of length `L`,
   residue class of `n`) that are compatible and show the count is below `1`:
   the itineraries have entropy at most `log 2` per step, the residue
   condition costs `2^(−A)` per word, and the 2-adic point of the itinerary
   must be a positive integer, i.e. must have an eventually-zero expansion,
   which the itinerary does not see. This is the same wall as S13's
   Proposition 4 unless the interval map's arithmetic (the rational slopes
   `3/2^v` and the break points `(2^v + 1)/3`) forces a 2-adic property on
   its itineraries. That is the open question this session leaves: **do the
   itineraries of `f` (and of its `B`-generalisations) have 2-adic values
   with unbounded binary expansions?** For the eight periodic itineraries
   the 2-adic values are non-integer rationals, which is consistent.
3. **Related angle that appeared.** The error `eta_l` reads the shadowed cycle:
   `eta_l ≈ |x_w|` while the tail shadows `x_w`. So the sequence `eta_l` is a
   real-valued observable of the 2-adic address that says which negative
   rational point the orbit is currently imitating, and the budget
   `2^v <= 3 eta − 1` says the orbit can only take a deep halving while it
   imitates a point of large absolute value. This turns the reset lane's
   "regenerated precision around the same negative cycle" into a quantitative
   ledger: precision regenerated around `x_w` costs nothing on the 2-adic
   side and raises the error to `|x_w|` on the real side.

## 6. Next

Prove or refute HYP-9164 in the `B = 2` regime: decide whether some positive
integer has all `eta_l < 2`, i.e. a word that is an `f`-itinerary. The
periodic itineraries have non-integer rational points; the aperiodic ones
form a set of dimension at most `log_2` of the growth rate of `f`'s symbolic
system, and the question is again transversality with `Z^+`, now for a much
smaller fractal than `E_inf`.
