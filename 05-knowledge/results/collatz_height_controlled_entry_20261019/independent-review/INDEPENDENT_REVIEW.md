# Independent audit: height-controlled prescribed-word entry

Date: 19 October 2026.

## Verdict

**PASS.** No mathematical gap was found in the concrete theorem, the general
prescribed-word theorem, or the optional excursion estimate. The proof is an
unconditional finite-parity counting argument. It does not assume a random
infinite orbit, condition a density-one theorem on an entry event, establish
pointwise entry, or establish infinite recurrence.

The reviewed source is `HEIGHT_CONTROLLED_ENTRY.md` in the adjacent
`collatz_height_controlled_entry_20261019` result package. Its audited SHA-256
is recorded in `check-receipt.json`; the author subsequently added only the
recommended Korec clock attribution, and the final receipt is regenerated
against that frozen text. No author proof/checker file, checkout, or Git state
was edited by this review.

One nonblocking attribution improvement was identified: Korec already notes
that his generic power descent can occur within `log_2 y` steps. The author was
notified and agreed to credit that explicitly. The distinctive content of the
present argument must therefore be described as the prescribed word together
with an explicit joint exceptional-count estimate, not logarithmic time by
itself. This review makes no priority claim for that combination.

## 1. Exact word and lawful cylinder

The 32-bit string has 20 odd bits. Independent reconstruction by successive
binary lifting and direct integer iteration gives the residue

    r = 2962948091 modulo 4294967296.

Direct iteration from this residue gives endpoint `2405410912`. The resulting
affine identity is exactly

    T^32(y) = (3486784401 y + 15810005461)/4294967296.

The slope is less than one, and

    15810005461 / (4294967296 - 3486784401) < 20 < r.

Every positive integer in the cylinder is `r + t*2^32`, with integer `t>=0`.
Thus all positive lawful lifts descend. There is no exceptional small lift
hidden in the asserted all-lift descent. The proof uses occurrence of the full
word, so branch lawfulness at the selected time is guaranteed.

## 2. Sampling space and selection

For `X=2^m`, the interval `[X,2X)` contains exactly `2^m` integers. When `N<=m`,
each residue modulo `2^N` occurs exactly `2^(m-N)` times. The parity encoding is
a bijection on these residues: the two lifts of an `N`-bit source differ at the
next output parity because their endpoints differ by the odd integer `3^h`.
Consequently the first `N` bits, jointly and without conditioning, have the
uniform product distribution.

All candidate blocks are disjoint deterministic subsets of those bits. Their
avoidance probability is exactly `(1-2^-32)^k`. The prefix odd-count and carry
events are also estimated on the same unconditional space. The union bound
does not require independence between those three bad events.

The subsequent choice of an occurrence time creates no selection bias: outside
the bad events, every relevant prefix already satisfies both bounds. Neither
the binomial law nor a carry moment is reapplied at a selected stopping time.

## 3. Carry and large deviations

With `B_0=0`, averaging the updates `B/2` and `(3B+1)/2` under a fresh uniform
bit gives conditional mean `B+1/4`. Hence `E[B_j]=j/4` exactly. Markov's inequality
at threshold `X^(15/16)/2` therefore gives `j/(2X^(15/16))`, and summing over
the relevant prefixes is bounded by

    N(N+1)/(4X^(15/16)).

The use of a supremum over many prefixes is legitimate because this is a union
bound, not an assertion that their carries are independent.

For `z=11/9`, exponential Markov gives the odd-count tail bounded by `q^j`, with

    q^20 = 10^20/(9^9 11^11).

The strict tail `H_j>11j/20` causes no rounding defect: on that event,
`z^H_j>z^(11j/20)`. Summing the geometric tail gives `q^J/(1-q)` as stated.
Exact rational arithmetic verifies `q<999/1000<1`. High-precision displays are

    q = 0.995004154110973160587971738874504688...
    -(3/4)log_2(q) = 0.0054191595091437797607714718...

These decimals are not used to establish the necessary inequalities.

## 4. Height, full-word clock, and rounding

Since `3^5<2^8`, the prefix bound gives

    A_j <= 2^(-3j/25) <= 2^(-9m/100),

using `j>=ceil(3m/4)`. Because `x<2X`,

    T^j(x) < 2X^(91/100) + X^(15/16)/2 <= X^(15/16).

The final comparison is exactly equivalent to `11m/400>=2`, and hence holds
for integer `m>=73`. The declared `m>=512` is more than sufficient.

The block count `k=floor((N-J)/32)` ensures that the last inspected word ends
at `J+32k<=N`. This includes the entire descent certificate, rather than just
the start of its parity word. Therefore

    j+32 <= N <= (7/8)m <= (7/8)log_2 x.

The important nonempty-window issue is correctly handled by `m>=512`.
For example `m=513` has only one candidate block, so a casual replacement of
floors and ceilings by their real-valued asymptotics would be unsafe.

The floor loss used in the proof is valid:

    k >= m/256 - 2.

The independent checker tests a complete period `512<=m<768` and the exact
translation law

    (N(m+256),J(m+256),k(m+256)) = (N(m)+224,J(m)+192,k(m)+1).

That periodic identity extends the rounding check to every `m>=512`; it is
not an extrapolation from a large finite cutoff.

## 5. Exceptional exponent and arbitrary source cutoff

Writing `p=2^-32`, block avoidance is at most

    (1-p)^(-2) X^(-delta),
    delta = -log_2(1-p)/256
          = 1.312123496115341563262705484710350336...e-12.

There is no lost factor of 32, 8, or `log 2`: asymptotically there are `m/256`
independent candidates, and `X=2^m`. The Chernoff exponent exceeds `3/4000`,
while `delta<2^-38<3/4000`. The latter bound follows, for example, from
`-ln(1-p)<=p/(1-p)`, `ln 2>1/2`, and `p<1/2`. The carry term has a much larger
power saving; its polynomial factor does not alter the dominant exponent.

For an arbitrary cutoff `Y`, bound the partially included last dyadic interval
by its complete-interval estimate and sum. Since `alpha=1-delta>0`, the sum of
`2^(alpha m)` is geometric and is `O(Y^alpha)`. There is no extra factor of
`log Y`, and the finitely many intervals below `2^512` are absorbed effectively.
The theorem is relative to each original source `x`, not just to the dyadic
lower endpoint, because `X<=x` and `m<=log_2 x` were used in the correct
directions.

### Optional explicit constant

The argument even gives the coarse all-cutoff bound

    #{positive x<=Y failing the concrete conclusion} <= 5000 Y^(1-delta)

for every real `Y>=1`. This is an audit-derived convenience bound, not a
practical numerical success guarantee.

Indeed, for `m>=512`, `m^2<=2^(m/2)` makes the carry contribution at most
`(1/2)2^(-delta m)`. The Chernoff term is at most `1000*2^(-delta m)` and the
avoidance coefficient is less than 2. Thus the dyadic count is at most
`1003*2^(alpha m)`. Since `alpha>1/2`, the geometric-sum factor is less than 4.
The initial range costs less than `2Y^alpha`, since `512delta<1`. The resulting
constant is at most 4014, safely rounded to 5000. When this bound exceeds `Y`,
the trivial bound is of course stronger.

## 6. General word and optional excursion

For fixed word length `L` and `log_4 3<beta<1`, the parameters exist as stated.
The quantity `eta=1-a(1-p0 log_2 3)` tends to `log_4 3` as `a` tends up to one
and `p0` tends down to one half. Strict inequality leaves room to choose
`a<c<1`. Both the word and its endpoint fit inside `floor(cm)` steps for all
sufficiently large `m`. The avoidance and large-deviation terms have positive
exponents, and the carry term is `O(m^2 2^(-beta m))`. One can choose a positive
effective exceptional exponent below all three rates. An arbitrary word need
not descend, so the conditional endpoint statement correctly requires the
extra uniform all-lift descent hypothesis.

For the excursion addendum, `A_j` is a nonnegative mean-one martingale on the
finite parity filtration. The finite-horizon maximal inequality gives the
claimed multiplier bound. The separate carry union bound is
`N(N+1)/(8X^(1+epsilon))`. Together with `x<2X`, these imply the ceiling
`3X^(1+epsilon)<=3x^(1+epsilon)` through time `N`, including the certified word
endpoint. With `epsilon=1/100`, the new error exponents exceed `delta`, so the
same final exception order remains valid. This does not give a ceiling of `x`.

## 7. Primary-source scope

- [Korec, Theorem 1 and page 86](https://dml.cz/bitstream/handle/10338.dmlcz/133225/MathSlov_44-1994-1_8.pdf)
  proves natural-density-one power descent above the threshold `log_4 3` and
  explicitly observes that a `log_2 y` time bound is available. The present
  argument neither improves that threshold nor establishes convergence for
  every integer. The source was inspected directly.
- [Rozier–Terracol, Lemma 2.5](https://arxiv.org/html/2502.00948v5)
  states the exact mean carry `j/4`; the elementary derivation in the reviewed
  note is correct. The cited 2026 journal reference agrees with the arXiv
  bibliographic record.
- [Tao's abstract](https://arxiv.org/abs/1909.03562)
  states almost-bounded minimal values for logarithmic-density-almost-all
  sources. The stronger height conclusion there is not a natural-density
  conclusion, and the reviewed note maintains that distinction.

## 8. Reproducible independent controls

Run `python independent_checks.py` in this review directory. The checker uses
explicit exceptions rather than Python assertions, so its checks remain active
under `python -O`. It imports no production checker and makes no network calls.

The independent controls comprise direct binary lifting for the 32-bit word;
direct integer-orbit parity, mean-multiplier, mean-carry and embedded block
avoidance checks on all 32,766 sources in dyadic intervals of exponents 1 through
14; exact rational parameter inequalities; and a complete affine-periodic
rounding check. The receipt records 460,126 checks. These bounded controls
authenticate identities and transcription. The proof of the all-size claim is
the symbolic argument audited above, not a source census.
