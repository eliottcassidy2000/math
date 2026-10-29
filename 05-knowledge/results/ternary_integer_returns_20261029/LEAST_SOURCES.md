# A pointwise replacement target: the least surviving integer

**PROVED:** the exact source-lift recurrence and pruning rule below.
**FINITE-EXACT:** least positive sources for 1 through 40 consecutive expanding
returns. **OPEN:** their unbounded growth. This tests a pointwise arithmetic
quantity instead of replacing it by a Haar probability.

Let `mu_k` be the smallest positive odd integer whose first `k` applications
of the cheap return controller all strictly increase. The five failure words
are indexed `0,1,2,3,4` in the table of INTEGER_RETURNS. The source
`2^(k+1)-1` realizes `k` copies of letter zero, so `mu_k` always exists.

## Exact source/endpoint recurrence

A prefix `u` has least positive source `r`, source modulus `2^(D+1)`, odd clock
`J`, and actual endpoint `h=(3^J r+C)/2^D`. On the rest of its cylinder,

\[
 n=r+2^{D+1}t,\qquad U_u(n)=h+2\,3^Jt.
\]

To append a failure word `v` with `(j,d,c)` and target guard `r_v`, choose

\[
 t_0=\left(\frac{r_v-h}{2}\,3^{-J}\right)\bmod 2^d,
 \quad 0\le t_0<2^d,                                  \tag{15}
\]

then update

\[
 r'=r+2^{D+1}t_0,\quad
 h'=\frac{3^j(h+2\,3^Jt_0)+c}{2^d},\quad
 (J',D')=(J+j,D+d).                                    \tag{16}
\]

This retains the original source and its actual endpoint. In particular
`r'>=r`; every descendant cylinder has least positive representative at least
that of its ancestor. Given an authenticated upper bound for `mu_k`, any
prefix with `r` at least that bound may therefore be pruned. Start with the
Mersenne upper bound and traverse all five legal continuations. This searches
**every positive source**, not just a selected finite input interval.

The independent implementation instead updates the total affine numerator
and recomputes the original-source congruence. Both methods give identical
minima and witness words at every horizon through forty. Every witness is
also run through the native integer Collatz operations, with exact block ends.

Selected results:

| `k` | `mu_k` | Mechanism |
|---:|---:|---|
| 14 | 32767 | all valuation-one blocks |
| 15 | 56575 | seven zero letters, one `3111` block, seven zero letters |
| 24 | 10190847 | an inserted nontrivial failure beats the Mersenne candidate |
| 31 | 3556769791 | two nontrivial failure blocks |
| 32–35 | 3823844863 | one source survives four successive tested horizons |
| 40 | 406360948735 | twenty-three zero letters, one `3111` block, sixteen zero letters |

Thus every positive odd source below `406360948735` has a non-growing return
among its first forty cheap returns. For a source greater than one this first
non-growing return is a strict drop from its **current** input. It need not
fall below the original source, and this finite statement is not a convergence
certificate for that interval.

The tempting candidate `mu_k=2^(k+1)-1` first fails at `k=15`. Its failure
mechanism is a nontrivial block that changes later dyadic phases; forgetting
that block's carry incorrectly predicts the all-zero-letter answer. The
survivor is the Mersenne **upper bound**, not its equality. Through forty, all
minimizing witness words use only the first two letters; this observation is
not used as a filter and is not asserted beyond the finite computation.

## Compare with the incoming all-precision recurrence

The concurrent [exceptional Cantor recurrence](../session97_overnight_arithmetic_carriers_20261029/controllers/EXCEPTIONAL_CANTOR_RECURRENCE.md)
proves an all-level counting theorem for the same exact inverse charts:

\[
 N_L=N_{L-1}+N_{L-6}+3N_{L-11}\quad(L\ge11),
 \qquad 2^{-s}+2^{-6s}+3\,2^{-11s}=1.
\]

Its initial values are `(1,2,3,4,6,9,13,18,24,31,40)`, and
`0.441<s<0.442`. It bounds possible forever-expanding positive sources below
`X` by `(17/8)X^s`. This is stronger than a qualitative zero-measure claim,
but does not prove that any such positive source exists or that none exists.
The `N_L` count uses division-clock precision, whereas `mu_k` uses completed
return blocks and actual least positive representatives. Their different
clocks and quantifiers must not be identified.

In ordinary 2-adic binary coordinates, a positive integer has an eventually
zero digit tail. The remaining arithmetic problem is therefore to exclude
such tails from this twisted inverse-chart attractor. Its address coding
includes odd-unit multiplication and full affine carries; it is not ordinary
concatenation of five fixed binary strings. The source-lift recurrence above
retains precisely those operations when testing the least-integer boundary.

## A precise pointwise helper question

The following are equivalent:

1. Every positive odd source eventually leaves the five-letter all-expanding
   return language.
2. `mu_k -> infinity`.

Indeed the sets of surviving integers are nested, so their least elements
are nondecreasing. A bounded sequence of integer minima eventually stabilizes
at a source surviving every horizon. Conversely any infinite surviving
positive source bounds every minimum.

A lower bound `mu_k>=c alpha^k` for some `alpha>1` would give a logarithmic
source-height bound for this escape. It would still owe the accumulated growth
before that escape. The native controller's endpoint inequality
`R(n)+1 <= (3/2)(n+1)` bounds that accumulated growth once an escape-time bound
is supplied, but does not by itself pay it. This separates two research
obligations that an averaged renewal score cannot merge.

Reproduce `check_least_sources.py` normally and with `-O`.
[least_sources.json](least_sources.json) records the complete finite horizons,
exact endpoints, labelled words and visited-prefix counts. All five letters
remain enabled throughout both searches; there is no inferred exclusion of
the less-used three branches.
