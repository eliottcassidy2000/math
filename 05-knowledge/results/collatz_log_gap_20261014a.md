# Logarithmic cycle gaps with an authenticated source prefix

**Status.** PROVED: the marked occupancy bound, monotone clock refinement, exact ordered-carry membership test, and endpoint-corrected logarithmic identity below. FINITE-EXACT: the explicitly bounded controls in §7. CITED: the cycle-context linear-form bound in Simons–de Weger, Lemma 12, and the classical Legendre criterion. CONDITIONAL: the consequences of any proposed irrationality-measure bound for `log(3)/log(2)`. OPEN: nontrivial integer cycles, divergent positive orbits, and universal Collatz convergence. This is a source-specific certificate compiler, not a new literature record for cycle lengths.

## 1. Inheritance and the coordinate that must survive

The closest mechanism is the ordered cycle gate in
[Pillai/convergent cycle gates](collatz_mod6_20260921_pillai_convergents_cycle_gates.md),
§§1–3: a word's multiplier and its carry are separate data. That note already
contains the hostile clock `8/5`: all 35 positive valuation words have gap 13,
but none has carry divisible by 13. Its signed/rational boundary includes
`5/7 -> 11/7 -> 5/7`, which irrationality of the logarithmic slope does not
exclude. The recent
[rational-approximant note](collatz_rational_approximants_riemann_sphere_20261013.md)
explicitly leaves its remembered Rhin constant/hypothesis unverified; §6
below supplies the directly checked primary source, without silently importing
the remembered coefficient 1.

The live concept board is: a native valuation word; its exact carry; the two
clocks (odd steps and shortcut steps); marked source states; a logarithmic gap;
and the endpoint needed when the path is not closed. The useful transfer is
from a hypothetical primitive cycle to a finite set of necessary clock
candidates. It preserves the source and observed states. It forgets the order
of unobserved states only in its explicitly labelled upper bound. The ordered
carry decoder restores the missing test for a complete proposed clock.

We use the user's shortcut map

\[
T(n)=\begin{cases}n/2&n\text{ even},\\(3n+1)/2&n\text{ odd},\end{cases}
\]

whose positive ROOT cycle is `1 <-> 2`. The accelerated odd map is
`U(n)=(3n+1)/2^a`, with `a=v_2(3n+1)`. One odd-map edge takes exactly `a`
shortcut steps: one odd step and `a-1` even steps. Thus a word of length `r`
and valuation cost `A` has `r` odd steps and `A` shortcut steps. The odd-map
ROOT is the fixed point 1, with word `(2)`.

## 2. The exact gap, and why arbitrary prefixes are different

For a native path `x_0,...,x_r` with valuations `a_i` and `A=sum a_i`,

\[
\frac{2^A}{3^r}
 =\frac{x_0}{x_r}\prod_{i=0}^{r-1}\left(1+\frac1{3x_i}\right),
\qquad
A-r\alpha=\log_2\frac{x_0}{x_r}
 +\sum_{i=0}^{r-1}\log_2\left(1+\frac1{3x_i}\right),
\tag{1}
\]

where `alpha=log_2(3)`. This follows by multiplying the actual edge equations.
For a positive cycle the endpoint quotient is one, so

\[
0<\Lambda=A\ln2-r\ln3
 =\sum_i\ln\left(1+\frac1{3x_i}\right)
 <\frac{r}{3m},\qquad m=\min_i x_i.
\tag{2}
\]

The endpoint quotient cannot be discarded on a nonclosed path. For every
integer `K>=2`, the source `2^K-1` has the exact prefix `1^(K-1)` and endpoint
`2*3^(K-1)-1`. These arbitrarily long actual prefixes have `A/r=1`, uniformly
far from `alpha`; the endpoint term in (1) balances the discrepancy. The
source changes with `K`. This neither constructs an infinite positive orbit
of ones nor excludes divergent orbits.

## 3. Sharp occupancy bound with retained observations

Every member of a positive odd cycle is a 3-unit: an odd output of `U` cannot
be divisible by 3. In a primitive cycle its states are distinct. Let `m>1`
be its asserted minimum, and authenticate an actual prefix of `ell` edges
with word `w` and distinct states

\[
M=\{x_0=m,x_1,\ldots,x_\ell\},\quad x_i\ge m.
\]

The API rejects a wrong valuation, a repeated state, a state below `m`, ROOT,
or a source divisible by 3. These checks authenticate the observations and
the necessary minimum conditions; they do not assert that a future cycle exists.

For a proposed primitive length `r>=ell+1`, complete `M` to `r` elements by
inserting the smallest unused odd 3-units at least `m`. Write the resulting
set as `S_r(M)` and define the exact rational number

\[
C_r(M)=\prod_{x\in S_r(M)}\frac{3x+1}{3x}.
\tag{3}
\]

**Marked occupancy theorem.** Every cycle compatible with these observations
satisfies

\[
\frac{2^A}{3^r}\le C_r(M)
 \le C_r(\varnothing)
 \le\left(1+\frac1{3m}\right)^r.
\tag{4}
\]

The first product is the maximum over all `r`-element sets of odd 3-units at
least `m` containing `M`. To prove it, sort the unmarked elements and replace
each by the corresponding smallest available element. Every factor
`(3x+1)/(3x)` strictly decreases with `x`. This also proves sharpness for the
relaxed occupancy problem, not attainability by an actual cycle.

Let `c=sum(w)`. The complete finite clock gate is

\[
\max\{\operatorname{bitlength}(3^r),\ c+r-\ell\}
 \le A\le\left\lfloor\log_2(3^r C_r(M))\right\rfloor.
\tag{5}
\]

The first lower bound is exactly `2^A>3^r`; the second gives every unobserved
edge at least one valuation bit. The upper bound is evaluated by integer
cross multiplication, never by a floating logarithm. An empty interval proves
that this source cannot be the minimum of a primitive cycle of this length.

**Refinement theorem.** Retaining additional actual states at the same source
can only decrease (3), and replacing one previously unknown valuation by its
actual value `a>=1` can only increase the lower bound in (5), by `a-1`.
Therefore the candidate-clock sets are nested under authenticated prefix
extension. This is the recursive invariant supplied by this packet.

For the 11-edge prefix of 91 the marked states are

`91,137,103,155,233,175,263,395,593,445,167,251`.

At `r=17`,

\[
C_{17}(M)=\frac{76631701190606848}{73779331906564965}
 <\frac{2^{27}}{3^{17}}.
\]

The unmarked distinct-state bound allows this clock, but the marked bound
excludes it. At `r=29`, `A=46` survives (5), illustrating the necessary-only
boundary. The ordered test in the next section rejects that survivor.

## 4. Carry decoding: the remaining finite test is not a logarithm

For a positive valuation word `w=(a_1,...,a_r)`, put `s_i=a_1+...+a_i`,
`s_0=0`. Its affine carrier is

\[
U_w(n)=\frac{3^r n+B_w}{2^A},\qquad
B_w=\sum_{i=0}^{r-1}3^{r-1-i}2^{s_i}.
\tag{6}
\]

For `r>1`, the first letter and remaining carry obey

\[
B_w-3^{r-1}=2^{a_1}B_{\mathrm{tail}},\qquad
a_1=v_2(B_w-3^{r-1}),
\tag{7}
\]

because every nonempty tail carry is odd. The right side is positive.
`decode_carry(r,A,B)` recursively applies (7), checks `a_1>=1` and that the
remaining cost can contain `r-1` positive letters, and terminates at `r=1`
with the required carry 1 and final letter equal to the remaining cost.
This proves a total iff membership test for the prescribed length, cost and
carry. It is the fixed-length affine version of the inherited
[translation decoder](translation_phase_decoder_20261005.md); no priority
claim is made for the decoding idea.

A proposed source `m` and clock require

\[
B=m(2^A-3^r).
\tag{8}
\]

`primitive_cycle_word(m,r,A)` decodes (8), then replays every exact valuation
and checks first return with no earlier repeated state. The result is a
literal primitive-cycle receipt or `None`. Its finite bound is the requested
length; it does not search until convergence. Source 1 returns `(2)` for
`(r,A)=(1,2)` and rejects the padded repetition `(2,2)` as nonprimitive.

The clock `8/5` has gap 13 but its 35 carries all fail divisibility by 13.
Conversely `(1,3)` has the positive rational fixed point `5/7` and visits
`11/7`: rational cycles survive perfectly valid logarithmic bounds. The
integer-source API rejects rational sources rather than changing the problem.

## 5. What an irrationality exponent would actually imply

Suppose a bound for the **ratio** is independently established:

\[
\left|\alpha-\frac pq\right|\ge c q^{-\mu}
\quad\text{for the relevant integers }p,q.
\]

Equations (1)–(2), applied to a cycle only, give

\[
m<\frac{r^\mu}{3c\ln2}.
\tag{9}
\]

If `mu(alpha)=2`, then for every positive `epsilon` one obtains such a bound
with exponent `2+epsilon` and a suitable positive constant, after absorbing
the finitely many exceptional denominators. Thus the consequence is a
polynomial relation between cycle minimum and cycle length. It is neither
nonexistence of all cycles nor control of nonperiodic trajectories.

Separate statements about `mu(log 2)` and `mu(log 3)` do not supply a
statement about their quotient. This is a genuine logical boundary for
irrationality measures: for any fixed nonzero Liouville number `L`, there
exist real `b` for which both `b` and `Lb` have irrationality exponent 2,
although `(Lb)/b=L` has infinite exponent. Indeed exponent 2 holds on a
full-measure set: the intervals for approximations of exponent `2+epsilon`
have summable total lengths on a bounded interval, while Dirichlet gives
the lower exponent 2; intersect over positive rational `epsilon`. Scaling
preserves null sets, so that set and its `L`-preimage intersect. This is not
a claim about either particular logarithm; it shows why a quotient transfer
needs its own theorem.

For a basic continued-fraction placement, (2) and `ln2>2/3` give
`0<A/r-alpha<1/(2m)`. Consequently `m>=r^2` puts the reduced `A/r` within
the classical Legendre threshold and makes it a convergent. More generally,
the positive-gap condition `C_r(M)-1<=1/(3r)` suffices. The factor is valid
here because the gap is positive; the inherited two-sided absolute-gap
criterion has different constants. Continued-fraction placement still leaves
(8) and native replay unpaid.

## 6. A primary-source linear-form bound, with its constant retained

[Simons–de Weger, Acta Arithmetica 117 (2005), pp.51–70](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/117/1/82263/theoretical-and-computational-bounds-for-m-cycles-of-the-3n-1-problem),
[original PDF](https://www.impan.pl/shop/en/publication/transaction/download/product/82263),
p.61, Lemma 12, gives in its cycle setting

\[
\Lambda>\exp[-13.3(0.46057+\ln r)].
\tag{10}
\]

Their `K` is our number `r` of odd steps, and their `K+L` is our shortcut
cost `A`. The proof explicitly applies Rhin's proposition at height
`H=K+L`, not at height `K`; it then uses their Lemma 8. This verifies the
cited cycle-context result, not a best current irrationality measure.

An exact convenient weakening follows without decimal numerical error:
`ln2>56/81`, from two terms of its atanh series, and
`(133/10)(46057/100000)<9(56/81)`. Thus (10) is greater than
`1/(512*r^(133/10))`. Combining with (2) yields the integer-only necessary
test

\[
(3m)^{10}<512^{10}r^{143}.
\tag{11}
\]

`rhin_minimum_possible` implements (11) and labels the cited dependency.
The elementary marked gate (5) does not depend on Rhin or on any claimed
irrationality exponent. No other claim in the incoming rational-approximant
note is audited or promoted here.

## 7. Exact controls and reproducibility

The script uses rational atanh-series bounds, including an explicit positive
tail bound, to certify the first 15 continued-fraction digits of `alpha`:

`1,1,1,2,2,3,1,5,2,23,2,2,1,1,55`.

It checks all 1,364 words over valuations `1..4` of lengths `1..5`; complete
carry membership on every clock of length at most 5 and cost at most 10
against its full composition list and a separate odd-carry test range;
1,060 independent occupancy subsets; and all marked-prefix refinements for
the four supplied 11-edge observations below and lengths `12..256`.

| Source | Repeated-minimum lengths | Distinct 3-unit lengths | Marked lengths | Least marked length |
|---:|---:|---:|---:|---:|
|31|227|99|98|17|
|47|206|85|84|17|
|71|182|70|70|17|
|91|161|64|63|29|

These count necessary candidate lengths in the declared finite interval,
not cycles or unproved convergence cases. They deliberately use truncated
observations of familiar sources to test the information refinement.
The actual carry test rejects the surviving `91,29,46` candidate. Additional
controls retain the rational-cycle hostile, reject ROOT padding, reject
boolean/float/source-state substitutions, and authenticate the Mersenne
endpoint formula for `2<=K<=64`.

Reproduce from the repository root:

```text
python -B 04-computation/experiments/collatz_log_gap_20261014a.py
python -B -O 04-computation/experiments/collatz_log_gap_20261014a.py
```

The deterministic output is saved as
[collatz_log_gap_20261014a.json](collatz_log_gap_20261014a.json).
Normal and optimized runs agree on **36,858 checks**. Its LF-normalized
SHA256 is `ed5fd538643087fc25ec3abd6999af94d4989a3ee9e4aec0deb104cd39ed2ea8`.
The next missing global step is not more digits of the continued fraction:
one would need a source-valid rule excluding all remaining clocks and also
nonperiodic escape, or independently grounding every source. Neither is
asserted by this packet.
