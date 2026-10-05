# Source-specific floors that survive refinement

2026-10-05. **PROVED:** a uniform signed kernel minorant, monotone atom
recovery, conditional finite ROOT extraction, and preservation under truthful
information refinement. **FINITE-EXACT:** the stated controls.
**OPEN:** an independent positive bound for every source. No new Collatz
source is certified by the numerical examples, which reuse completed words.

[Script](../../04-computation/experiments/collatz_localized_resolvent_floor_20261005.py)
and [exact output](collatz_localized_resolvent_floor_20261005.out).

## 1. Inheritance, portfolio, and what refinement means

The anchor is one source's injection weight, not a positive residue sum.
The niche is a signed dual from moment geometry. The wildcard is transporting
the old quadratic trace coordinate into a source-local kernel.
The board is **source identity / signed selector / observable access /
error budget / finite ROOT deadline / bounded-weight witness compactness**.

The closest current mechanism is
[Fourier atom positivity](collatz_fourier_atom_positivity_20261005.md):
`p_m=lambda(6m+3)` is a probability on nonnegative indices; an atom is positive
exactly when that particular source has a finite ROOT word. Positive residue
sums at every precision do not imply a positive atom. Its missing-atom
geometric law is the canonical hostile.

The incoming [Pascal classification](collatz_pascal_boundary_leaf_section_20261005.md)
classifies the **price** variable. Its source-floor audit corrects the claim
that all finite-dimensional encodings are excluded; the actual theorem only
excludes specified observers. Here the moments concern the **source**
distribution. A price prior and a source measurement have different types.

Two less-used routes supply the constructive move:
[THM-2237, Boolean moment atom duals](../../01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md)
uses a pointwise polynomial inequality to read out one atom;
[THM-2842, positive-cone multiplier observability](../../01-canon/theorems/THM-2842-ordered-positive-cone-vandermonde-multiplier-observability.md)
keeps the multiplier readout as an explicit extra input. Their mechanisms,
not their original support or Gaussian hypotheses, transfer here.

The [ordinary polynomial dual](collatz_atom_polynomial_dual_20261005.md)
already constructs `A_(m,d)<=p_m<=A_(m,d)+(32/9)4^(-d)` from ordinary
moments of `x_j=2^(-j)`. Its raw coefficient norm grows with the source index.
The construction below instead has coefficient norm17 for **direct kernel
measurements**. This change of observation primitive is a required sidecar,
not free access to better-conditioned information.

Three operations must remain distinct:

* Increasing measurement precision about one fixed source preserves a valid
  floor. Nested intervals make the lower bound for a fixed selector increase.
* Increasing the selector degree improves its exact lower readout. Retaining
  the maximum of all previously certified floors preserves it even if new
  measurement intervals are less precise.
* Prepending actual inverse edges changes the source. Its floor needs the
  exact guarded transport bound in the
  [deadline companion](collatz_floor_transport_deadlines_20261005.md).
  Arbitrarily long source changes have no uniform positive floor.

## 2. A two-term source-local selector

For a selected index `m>=0` and every `j>=0`, set

    t=2^(m-j),       h_m(j)=4t/(1+t)^2.

The kernel is symmetric in the index difference and satisfies

    h_m(m)=1;
    h_m(m-1)=h_m(m+1)=8/9       when the index exists;
    0<h_m(j)<=16/25             if |j-m|>=2.             (1)

Indeed h(t)=h(1/t), and h decreases for t>=1 since
`h'(t)=4(1-t)/(1+t)^3`. The nearest non-neighbor has t=4 or1/4.
There is no fictitious negative-index atom when m=0.

For integer d>=0 define

    q_(m,d)(j)=(9h_m(j)-8) h_m(j)^d,
    H_(m,d)=sum_j p_j h_m(j)^d,
    A_(m,d)=9H_(m,d+1)-8H_(m,d).                       (2)

Here `H_(m,0)=1`; all these moments are in[0,1]. Equations(1) imply
`q(m)=1`, q=0 at the available neighboring indices, and q<=0 everywhere else.
On the remaining support `|q|<=8(16/25)^d`. Integrating gives the theorem

    A_(m,d) <= p_m <= A_(m,d)+8(16/25)^d.              (3)

For every probability on the declared support, the exact lower readout
increases to p_m: a negative off-target term is multiplied by a number in(0,1)
when d increases, while the target stays1 and the neighbors stay0. Hence

    p_m>0  iff  A_(m,d)>0 for some finite d.            (4)

The convergence error in(3) is uniform in m. If an external argument already
gives `p_m>=epsilon`, it suffices to choose d with
`8(16/25)^d<epsilon` to force an exact positive readout. Equation(4) itself
does not prove that all source atoms are positive.

## 3. Exact interval receipt and preservation

Suppose independently justified measurements of the same actual measure give

    H_(m,d) in[a,b],       H_(m,d+1) in[c,e].

Then the signed interval computation gives

    L=9c-8b <= p_m <= 9e-8a+8(16/25)^d.                (5)

A strict L>0 is a source-specific floor. With absolute error at most eta in
each measured moment, the error of the readout is at most17eta, independently
of m and d. A positive floating-point central value is insufficient.

For nested valid intervals, c can only increase and b can only decrease, so
L can only increase. Across different d, retain `max(0,L_0,...,L_d)`; each
stored positive floor remains valid. The retained receipt consists of the
source index, kernel degree, exact signed coefficients, moment bounds, and
their proof provenance. Forgetting source identity or combining bounds from
different measures invalidates the inference.

The API checks exact numeric types and elementary interval consistency; it
does not certify arbitrary supplied intervals as facts about lambda. Its
extraction stage additionally verifies a literal first-hit ROOT word, so an
unsupported positive packet cannot become a false convergence certificate.

## 4. The quadratic trace becomes a localization map

Use the earlier quadratic trace coordinate on this new positive scalar:

    J=t+t^(-1),       h=4/(J+2).

Doubling the **source-index separation**, hence replacing t by t^2, gives

    J -> J^2-2,       h -> h^2/(2-h)^2.                (6)

This is the same algebraic trace identity as the earlier affine-word slope
construction, now used to classify separation from a selected source. The
source and target types differ: this t is a coordinate ratio, not a Collatz
word slope. No legal Collatz orbit or source is transported by (6).
The preserved predicate is the target t=1 and its ordered separation shells;
the lost datum under J is the sign of j-m. The source-index sidecar recovers
that orientation when needed. For the selector the loss is harmless because
both sides receive the same nonpositive value.

This supplies a precise use of a small recurring structure: trace doubling
organizes a bounded observable whose signed expectation isolates one atom.
The cheapest decisive test is (1), including the two neighbors, rather than
a numerical resemblance between unrelated constants.

## 5. Where the observation cost goes

Let `X=2^(-j)` and `s=2^(-m)`. Then

    h=4sX/(X+s)^2,
    R_k(s)=E[(s/(X+s))^k].

Expanding with `r=s/(X+s)` yields

    H_(m,d)=4^d sum_(a=0)^d (-1)^a C(d,a) R_(d+a)(s). (7)

The coefficient norm in(7) is8^d. Thus a constant17 error multiplier applies
to direct H-measurements, not to an implementation which first estimates all
the R-moments independently. Combining(2) and(7) has the safe error budget
`(9*8^(d+1)+8*8^d)eta=80*8^d eta` for equal errors at those primitive readouts.
The inverse powers are resolvent observables at the source-dependent scale s;
ordinary moments alone do not supply them at constant precision cost.

Changing base2 to b>1 changes the neighbor value to `c=4b/(1+b)^2`.
The same argument uses `(h-c)h^d/(1-c)` and the next-shell bound
`h<=4b^2/(1+b^2)^2`. Taking a very large base makes the measurement itself
nearly an atom query. Improved formal constants therefore do not establish
cheaper access to the missing source information. Base2 is kept fixed here.

The old missing-atom hostile also survives: if p_m=0, every exact A_(m,d)<=0,
even though the underlying law can have positive mass in every residue class
and strictly positive finite Hankel matrices. The signed readout retains the
inequality direction that a generic positivity statistic loses.

## 6. A positive floor pays for a finite literal route

The companion [floor/deadline theorem](collatz_floor_transport_deadlines_20261005.md)
proves, for a rooted nonroot source with counters N=L+K and odd ROOT time tau,

    tau<=N,        W(n)<=2/((N+1)(N+2)).               (8)

For rational epsilon=p/q>0 let

    M=floor(2q/p),       B=(isqrt(1+4M)-3)//2.

Any true `W(n)>=epsilon` forces tau<=B. For the selected leaf `n=6m+3`,
`p_m=W(n)`, so a positive lower receipt in(5) yields a literal U replay of at
most B odd steps that must reach1. For a unit source one may use its explicitly guarded
three-divisible predecessor from
[inverse predecessor sections](inverse_predecessor_sections_20261005.md),
then retain the actual edge to that unit. Positivity at every selected leaf
is still the missing universal premise.

The [weight-threshold compiler](collatz_weight_threshold_receipts_20261005.md)
provides a second extraction route: enumerate the entire finite superlevel
set. Its largest source is exactly `(4^(B+1)-1)/3`, so it can contain very
large integers while still being finite. Individual forward replay is usually
preferable when only one target is requested.

The script's known-head controls are:

| source | first tested positive degree | retained dyadic floor | finite deadline | actual odd steps |
|---:|---:|---:|---:|---:|
| 3 | 2 | 1/32 | 6 | 2 |
| 9 | 9 | 1/1024 | 43 | 6 |
| 27 | 42 | 1/68719476736 | 370726 | 41 |

These moment intervals use sixteen already completed leaf atoms and the
probability normalization to bound the remaining tail. They test the entire
representation-to-receipt interface but add no convergence coverage.
A fabricated packet claiming `W(27)>=1/4` fails actual bounded extraction.

## 7. Broader connections that change the next experiment

| Recovered source | Transferred operation and preserved predicate | Lost information / required sidecar | Decisive test |
|---|---|---|---|
| THM-2237 and THM-2842, linked above | Signed pointwise dual becomes a source-atom floor | External moment access and sign-aware error | Missing target must never give a positive valid floor |
| [THM-4027, Sun universal modular solubility](../../01-canon/theorems/THM-4027-sun-two-four-six-eight-universal-modular-solubility.md), with [THM-4026, its integer counterexample](../../01-canon/theorems/THM-4026-sun-two-four-six-eight-binomial-counterexample.md) | Keep ordinary height when passing through all congruences | Local witnesses may change or escape | A uniform **individual witness** weight, not a residue mass sum, confines witnesses to the finite superlevel bank |
| [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md) | Compose certificates only at a matching ordered boundary state | Positive totals can have disjoint state supports | Preserve actual source and endpoint in floor transport |
| [THM-4210, Rule30 lossless current tree](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md) | A lossless refinement retains its sidecar | Stored information does not imply physical admissibility or termination | Demand an actual replayed ROOT word after a positive floor |
| [THM-2253, online dyadic contrast extractor](../../01-canon/theorems/THM-2253-online-dyadic-contrast-tournament-extractor.md) | A whole-episode pathwise deadline limits accumulated cost | Fairness and per-stage stopping do not provide such a deadline | Use the total guarded counter budget, as in the floor companion |

The threshold companion makes the Sun connection quantitative: if every
refinement modulus has a rooted representative of a fixed source residue
whose **own** W-weight is at least the same epsilon, then the source itself
belongs to the finite bank. At epsilon=1/1000, source27 has changing witnesses
17 mod5,227 mod25,277 mod125, then none mod625. Positive coarse witnesses
have not supplied a persistent source floor.

The next productive target is now narrow: prove a source-indexed inequality
`9H_(m,d+1)>8H_(m,d)` with a rigorous error margin, using an arithmetic
identity or guarded decomposition that contains information beyond the same
finite ROOT bank. One need not estimate a whole distribution or require one
d for all m. Multiple constructions may cover different source families;
each successful inequality compiles into a finite, independently checked word.

## 8. Reproduction and audit

```text
python -B 04-computation/experiments/collatz_localized_resolvent_floor_20261005.py
python -B -O 04-computation/experiments/collatz_localized_resolvent_floor_20261005.py
```

There are29,483 explicit rational checks: source indices0..12, degrees0..8,
nodes0..40; independent finite probability controls with missing atoms through
degree20; trace and resolvent identities through distance12; exact interval
error; three actual pipeline controls and malformed input/false-premise tests.
The infinite-support inequality and convergence are proved in section2,
not extrapolated from this finite universe. Independent agent audits of the
kernel theorem and its oracle-conditioning distinction accompany the session.
