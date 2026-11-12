# Exact correction needs the right missing information

**Session277, 12 November 2026. Completed research packet. Chan–Pak remains OPEN.**

The request was not to accept a good approximation as good enough: measure
its signed error and build corrections that really recover the integer. This
session does that in two different settings, then identifies precisely what
the resulting recursive families can and cannot supply to Chan–Pak.

## 1. A sharper answer to why A000568 approximations fail

A000568(n) counts ordinary tournament isomorphism classes. Its odd-cycle
Burnside formula is a sum of positive rational symmetry sectors. We prove:

> For every n>=64, **every individual sector exceeds181**. Therefore no proper
> subset of those sectors can recover the integer just by rounding its sum.

This is stronger than saying that a particular asymptotic approximation has a
large error. It covers adaptive sector orderings and arbitrary positive-sector
truncations. It does not forbid signed resummation or an exact aggregate
completion. The proof bounds every sector below by
2^((n-(1+log_2 n)^2)/2).

We recorded **1,004 signed-error rows**, all sizes0..250 with support cutoffs
0,3,6,9. All errors are undercounts (or zero when the sum is complete). At250,
retaining every sector moving at most nine vertices has relative error about
10^-561, yet misses about10^8316 in absolute terms. Fixed support corrections
remove only linearly many error bits while the integer has quadratic size.

[Proof, exact tail interval and complete error table](absolute/README.md).

## 2. A complete, certified correction architecture for A000568

Instead of guessing the omitted amount, combine a **rigorous rational interval**
with an **exact residue modulo2^b**, where2^b exceeds the interval's integer
width. Exactly one integer in the interval has that residue. The code checks
this uniqueness and does not read the known exact answer to choose its result.

The residue backend uses the inherited exact response recurrence, now compiled
safely at finite precision. One inverse of the odd part of n! replaces inverses
for individual sectors. A new precision-debt lemma explains why the exact
signed-gauge recurrence cannot simply be reduced modulo2^b: right shifts need
extra high bits. A backward recurrence specifies sufficient precision at each
state of any exact-divisible signed-shift DAG.

| n | full answer bits | residue bits, support cutoff3 | bits saved |
|---:|---:|---:|---:|
|64|1721|1544|177|
|100|4426|4132|294|
|150|10303|9846|457|
|200|18655|18034|621|
|250|29489|28703|786|

The compiled250 reconstruction succeeds under800MiB/90s resource limits in a
recorded14.42s run and agrees with the inherited exact cache. The Python250
attempt was killed by the resource limit and is explicitly retained as a failed
attempt, not a successful result. No speed advantage over Session270's full
GMP counter is claimed. The achievement is **guaranteed exact completion**, not
constant-time computation: most bits and a growing state space remain needed.
Explicit exact output alone costs Omega(n^2) bit operations.

[Algorithm, proof, code, timings and reproduction](engine/README.md).

## 3. Recursive families: a coordinate change can destroy finite recurrences

For transitive modules substituted into a fixed quotient tournament Q, the
inherited run/path-cover compiler becomes a rational multivariate walk series
under u_i=exp(x_i)-1. This produces a sharp classification of the
**multinomial-weighted Hamiltonian aggregate**, not ordinary class counts:

- Its sequence is P-recursive **if and only if Q is transitive**.
- A directed cycle forces a genuine positive walk-series pole and infinitely
  many translated poles after the exponential substitution. Thus a finite
  run-state recursion need not become a finite polynomial-coefficient recursion
  in vertex-size coordinates.
- Every regular quotient of the same order shares the aggregate. That exact
  collapse loses the individual graph and target information.

The regular case has an exact geometric-moment/pole expansion. Directed
intervals and an analytic tail bound turn it into an all-order terminating
integer-recovery algorithm with adaptive precision. For the C3 aggregate:

| n | leading-pole approximation minus exact, approximate display | conjugate pairs used to certify all digits |
|---:|---:|---:|
|20|-82.4145|2|
|50|-8.26745×10^23|3|
|100|+3.65940×10^77|6|
|250|+6.18379×10^291|15|

Unlike the positive Burnside truncation, this correction uses signed complex
conjugate contributions. Its genuine success at integer recovery shows why
one must classify the error mechanism rather than extrapolate a relative
accuracy slogan. The discrete geometric measure is the useful bridge from
older moment-transform work; no Gaussian-moment or Jacobian conjecture is
assumed or claimed proved here.

[Family theorem, compiler, certified intervals and controls](families/README.md).

## 4. What this contributes to the constructive Chan–Pak task

For every fixed nontrivial strong quotient Q and positive transitive modules,
we prove H(Q[TT_r1,...,TT_rq])>=2^max(r_i)+1. A prescribed target c therefore
has a **complete finite inverse search**, with each r_i<=floor(log_2(c-1)).
This is not an arbitrary numerical cutoff. The singleton quotient is excluded.

The productive C3 face has

    H(C3[TT_a,TT_b,1]) = 2 + (2^a-1)(2^b-1).

It provides an exact factorization inverse and literal witnesses, including23
and47. The full C3 transitive-inflation family realizes45 but cannot realize
25,27,29,31 at any size; these are family exclusions, not tournament-spectrum
exclusions. Ordinal sum composes actual witnesses multiplicatively.

The same bound proves that a finite bank of fixed strong quotients with
transitive modules reaches only polylogarithmically many targets below X.
To finish Chan–Pak, this route must add growing quotients, nontransitive
children or nonuniform reusable interfaces. The full path-cover compiler can
accept those richer children; a scalar aggregate cannot replace their state.

A separate bridge makes the census/support distinction quantitative: each
fixed-H class fiber grows at most polynomially with tournament size, using
inherited finite strong-core inventories. Eventually a single omitted
Burnside sector outweighs the entire fiber. Thus even a remarkably accurate
scalar census can hide all information about a prescribed target.
[Precise scope and comparison](bridges/CENSUS_VERSUS_SUPPORT.md).

Concurrent S93's independently published bound A(Y)>=cY^(3/4)/(log Y)^4 at
commit f90b22ab5daa626c21d92da4fac37d9d5052cbd8 is an inherited improvement in
spectrum abundance, not this session's result and not A000568. It does not
supply the missing prescribed-target selector.

## 5. Completion and remaining obligations

New mathematics passed independent reviews, with minor boundary/precision
repairs applied. [Recovery review](audit/RECOVERY_REVIEW.md) includes907 Python
and128 compiled independent checks. [Family review](audit/FAMILY_REVIEW.md)
includes direct permutation checks, aggregate controls, pole boundaries,
Mersenne-face checks and complete fixed-family target boxes. Author controls
include5,357 independent subset-DP inflation comparisons and18 directed
integer certificates. Exact finite observations, all-order proofs and
benchmarks are separately labelled throughout.

**Not proved:** Chan–Pak, an eventual all-odd threshold, a finite universal
alphabet of tournament structures, or a polynomial-time arbitrary-order exact
A000568 algorithm. The result is stronger error correction and a precise
classification of a compositional boundary, not a claim that every remaining
constructive obligation has disappeared.
