# Tournament tails as a cycle-count polynomial module

2026-10-21. **PROVED:** the normalized transfer, the exact change of state,
the terminal linear-functional fold, and their integer implementations.
**FINITE-EXACT:** the coefficient and hostile controls in `receipt.json`.
**MEASURED:** the separately recorded shared-host benchmarks. This is an
enumeration algorithm, not a tournament classification or a polynomial-time
all-order bound. In particular, the incoming degree-local scalar algorithm
remains an important ordinary-count baseline; fewer profile keys alone do not
prove faster execution or lower memory.

## 1. Inheritance and target

The closest mechanisms are [THM-514, unit-tail and prime-tail kernels](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md),
[THM-516, the divisor-profile quadratic form](../../../01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md),
[future-response enumeration](../tournament_future_boundary_20261020/README.md),
and [the packed common-odd-denominator implementation](../tournament_enumeration_20261021a/README.md).
The identity `F'(x)=F(2x)` and the ordinary unit-tail formula are inherited,
not new here. Reconstructed Session21 evidence is idea provenance, not an
executable performance or proof premise.

During this session, commit `622e445138` supplied the independently developed
[degree-local transfer](../degree_local_tournament_enumeration_20261021/session22-enumeration/README.md).
Its exact `local_scale.py` is frozen as [incoming_local.py](incoming_local.py)
for reproducible comparison. That method removes state-dependent divisions
but retains the response-state graph. This packet changes that graph into a
polynomial module and folds its final functional before signed excess is
distributed. The live concepts are future gcd responses, cycle count,
quadratic energy, integer denominator lattices, signed core excess, and
terminal linear functionals.

The hostile boundary is deleting cycle count because its universal quadratic
term can be removed algebraically. It can be moved into a polynomial payload;
it cannot be discarded. The least-used sidecar is the distinction between a
profile key and the several coefficient cells it owns.

## 2. Exact normalized state and transfer

For an odd-cycle history `a_1,...,a_t`, let actual mass be `M=sum a_i`,
`z=product ell^k k!`, and

\[
 e=\sum_i(a_i-1)/2+\sum_{i<j}\gcd(a_i,a_j),\qquad
 g=e-\binom t2.
\]

The Burnside weight is `2^e/z`. Define

\[
 S_j=\sum_i(\gcd(a_i,j)-1)=R_j-t.                 \tag{1}
\]

All these integers are nonnegative; `S_1=0` and `S_j<=M`. At committed
degree `d` and a fixed descending stage, replace the old scalar state
`(d,R_j)` by `(d,S_j:j>=3)` with polynomial payload

\[
 P(y)=\sum_{\text{histories in this key}}
            \frac{2^g}{z}\,y^t.                \tag{2}
\]

In core mode the inherited signed excess coefficient is included in each
summand. Committed degree `d` can then exceed actual mass `M`; neither `t`
nor the responses count excess as an actual cycle. Only still-affordable
future odd `j` are retained, exactly as in the inherited descending-order
mask. Once unaffordable, a length cannot become affordable later.

**Transfer theorem.** Adding `k` actual cycles of odd length `ell` acts by

\[
 P(y)\longmapsto
 \frac{2^{(\ell-1)k^2/2+kS_\ell}}{\ell^k k!}
       y^kP(y),\qquad
 S_j\longmapsto S_j+k(\gcd(\ell,j)-1).          \tag{3}
\]

For ordinary counting the degree increases by `ell*k`. For core counting,
each excess `s>=0` contributes additionally
`(-1)^s binom(k+s-1,s)` at degree `d+ell*(k+s)`, with `k=0` allowing only
`s=0`. Formula (3) is independent of the payload's exponent `t`.

**Proof.** The inherited energy increment is
`k(ell-1)/2+ell*binom(k,2)+kR_ell`. Subtract
`binom(t+k,2)-binom(t,2)=kt+binom(k,2)` and use `R_ell=t+S_ell`.
The result is exactly the exponent in (3). All merged histories at a key
therefore undergo the same polynomial transformation. Linearity proves
validity after aggregation, including signed cancellation. This proves the
algorithm for every finite bound, rather than inferring it from the checks.

The final readout is the linear functional

\[
 \mathcal G(y^t)=2^{\binom t2}.                 \tag{4}
\]

It restores the omitted energy and hence the exact original coefficients.
It is not evaluation at a fixed number or a ring homomorphism:
`G(y^2)=2` but `G(y)^2=1`.

## 3. The final fold and integer arithmetic

At the final permitted length (1 normally, 3 in the no-fixed observer), fix
a key, an actual multiplicity `k`, and its scalar factor from (3). All
remaining core-excess terms have the same cycle polynomial `y^k P(y)`.
Compute `G(y^k P)` **once before** distributing excess. This avoids repeating
the coefficient-cell multiplication and merging for every terminal excess.
It is particularly effective for signed pair-free-core enumeration. It does
not assert that the ordinary unit-tail identity itself is new.

Two independently checked integer realizations are provided.

1. **Global polynomial (`list_poly_count`).** Store `N! P` coefficientwise.
   Each raw history is integral: `g>=0`, `z` divides `M!`, and `M<=d<=N`.
   Reuse a multiplicity ladder. Between `k-1` and `k`, its numerator gains
   `2^{((ell-1)/2)(2k-1)+S_ell}` and its denominator is `ell*k<=N`.
   Every coefficient division is checked. Integrality follows for raw
   histories and then for their signed sums; it is not assumed merely from
   one observed aggregate divisibility.
2. **Local polynomial (`module_dp`).** Store `d! P`. Before excess, the
   integer multiplier is
   `T=(d+ell*k)!/(d!*ell^k*k!)`. It is an integer because it equals a binomial
   choice of the new labels times the number of permutations of those labels
   with `k` cycles of length `ell`. Excess adds the integral factorial ratio
   `d'!/(d+ell*k)!` and its signed binomial coefficient. These transfers depend
   on degree and multiplicities, not response or polynomial exponent.

Both scalars are sufficient, with no minimality claim. The earlier minimal
odd denominator belonged to the **unnormalized** raw weights; it is not
silently reused for (2). All final divisions use the matching normalization.

`hybrid_dp` implements the exact change of basis after a scalar degree-local
prefix. It uses the live fixed-point response `R_1=t`, converts
`O(d)w` to `d! w / 2^{binom(t,2)}`, and converts other responses by (1).
When no degree remains, it settles the scalar directly. This hybrid API
requires fixed points: it does not pretend an omitted `R_1` can be recovered
in the no-fixed observer. The full polynomial APIs support all four modes.

### Representation probes and loss boundaries

The calendar ladder skips forced zero-multiplicity stages, but its queue
bookkeeping can exceed the saving. Its shortcut count is per emitted
contribution, **not** the reduction in unique merged transitions. Directly
measured transition counts are the comparison.

A balanced-radix polynomial prototype is exact in bounded controls but slow:
it reserves large gaps between coefficients. The production polynomial
payload is an ordinary sparse dictionary. The prototype uses a conservative
coefficient budget: at most `2^N` relevant partitions, total absolute excess
coefficient at a fixed degree for each partition at most `2^N`, raw exponent
at most `binom(N,2)`, and common `N!`; its radix
has two extra guard bits. This is a representation control, not the selected
fast implementation. For the excess bound, a partition of mass `M` with `t`
cycles has absolute excess series `product_i(1-x^{a_i})^{-1}`; its degree
`d-M` coefficient is at most `binom(d-M+t-1,t-1)<=2^N`, since `t<=M`.
The empty history is handled separately.

There is also a sharper mathematical size bound before the unit stage. If
all actual parts are at least 3, then

\[
 g\le t(M-t)/2\le M^2/9.                       \tag{5}
\]

Indeed `gcd(a,b)-1 <= ((a-1)+(b-1))/2`; summing gives the first inequality.
Since `t<=M/3`, the second follows, and equal 3-cycles attain it. This bound
describes coefficient sizes, not the implemented prototype's radix choice.

## 4. Hostiles and connection contract

* Histories `(9)` and `(3,3,3)` have the same mass and zero future residual
  at the unit stage. Their normalized weights are `16/9` and `256/81`, but
  their actual weights are `16/9` and `2048/81`. Forgetting `y^t` fails.
* Histories `(5,13)` and `(9,9)` have the same mass and cycle count but future
  `S_3=0` and `4`. The cycle polynomial does not replace the gcd sidecar.
* Core `k=1,s=1` creates one actual cycle, not two; restoring a two-cycle
  Gaussian factor gives the wrong coefficient.
* The no-fixed ordinary coefficient at degree 3 is `2/3`. That observer is
  retained as a rational weight, not rounded or promoted to a class count.

**Source:** odd-cycle weighted histories and signed core substitutions.
**Target/map:** residual future-response keys with a normalized cycle-count
polynomial, then the exact linear readout (4). **Preserved:** every coefficient
in the declared modes, actual cycle count, all affordable future responses.
**Lost:** individual cycle partitions after summation, individual tournaments,
automorphisms, metagraph vertices, and Hamiltonian-path ports. **Sidecars:**
degree, stage, mode, response length, polynomial exponents, and denominator
normalization. A metagraph construction cannot be recovered from these counts.

## 5. Exact controls and measured scope

Run:

```text
python -B 05-knowledge/results/tournament_tail_kernels_20261021b/check.py
python -B -O 05-knowledge/results/tournament_tail_kernels_20261021b/check.py --verify
python -B 05-knowledge/results/tournament_tail_kernels_20261021b/check.py --benchmark
```

The deterministic verifier checks all four modes through bound 24 against
the previous scalar implementation and frozen incoming degree-local kernel;
unmerged rational partitions independently cover all modes through 18.
Raw histories through mass 30 test the normalized identity and integer
transfers. Large controls retain every ordinary coefficient through 200 and
all three other modes through 100. The stored independent A000568 table
through 76 and direct partitions through 60 provide additional routes.
Malformed literal bounds/modes, empty histories, omitted fixed points,
signed coefficients, and changes of basis at several tail thresholds are
explicitly covered. No assertion is load-bearing.

The fresh paired experiments include the incoming algorithm, not just the
older 7.94-second observation on a differently loaded session. Ordinary
bound 200 structural transitions fall from **2,994,814 to 1,515,929** and
peak profile keys from **96,836 to 63,299**. Coefficient cells remain, and the
global polynomial kernel performs **2,355,764 coefficient ladder divisions**;
there is no claim that the total number of scalar operations halves.
Core bound 100 structural transitions fall from **1,697,975 to 289,706**.
The global polynomial terminal fold is the principal measured core gain.
Ordinary runtime can still favor the incoming scalar degree-local method.

`benchmarks.json` records three alternating-order paired samples, all exact
output hashes, setup/output conversion, and separate small-bound traced
memory runs. Timings are shared-host observations, not isolated attribution
to one algebraic identity or a claim of universal fastest performance.
The inherited implementation also records more telemetry. These limitations
are retained rather than hiding a slower representation behind a key count.

The completed three-sample medians (seconds) are:

| Implementation | Ordinary through 200 | Core through 100 |
|---|---:|---:|
| Previous common odd scalar | 4.6453 | 1.0660 |
| Incoming degree-local scalar | 2.1967 | 0.8517 |
| Global normalized polynomial | 3.4817 | 0.2486 |
| Local normalized polynomial | 3.2262 | 0.3150 |
| Incoming scalar prefix, polynomial unit tail | 2.1777 | 0.3222 |

The ordinary hybrid is effectively tied with the incoming kernel in these
samples; no improvement over it is claimed. The selected core kernel is
about 3.43 times faster than that incoming implementation in this comparison.
For ordinary bound 100, traced peaks are 1,449,988 bytes for the incoming
scalar and 2,109,552 for the global polynomial: fewer keys do **not** give a
memory win there. For core bound 60, the corresponding peaks are 4,428,468 and
315,624 bytes; the previous common-odd kernel is smaller still at 278,020.
These are separate traced runs, not the untraced timing samples. Peak payload
coefficient cells at ordinary 200 are 96,842, compared with 96,836 previous
scalar keys; at core 100 they are 5,493 versus 5,480. Structural sharing and
the terminal fold, not deletion of those coefficient cells, explain the gain.

The next precise algorithmic obligation is to reduce polynomial coefficient
arithmetic without losing the final Gaussian functional, or to accelerate
signed nonterminal excess transforms. A polynomial bound on the number of
future-response states and an all-order efficient tournament classification
remain open.

## 6. Final reproduction and audit

Normal and independent optimized execution match the saved **79,565** active
checks and every exact coefficient. Receipt LF SHA-256:
`d45245f0b1ae5321beb3152c37ce8e5246a355eb65d2349cf31031df9c11c5d2`.
The separate benchmark receipt has LF SHA-256
`2bb79e980b5f11b825db6ab3c591cbc1b45cbd0163f46ab8eb2cc2b9784f1ca9`.

The parent independently reviewed the normalized transfer, global and local
denominator arguments, Gaussian terminal fold, coefficient-size bound, loss
boundaries and production module API: **PASS**. No timings were inferred from
that proof review. The polynomial prototypes and alternate scheduling
experiments remain labelled controls; they are not selected by hiding their
slower or larger results. All five owned files are frozen for integration.
