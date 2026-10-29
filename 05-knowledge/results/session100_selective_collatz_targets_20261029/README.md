# Session100: selective Collatz targets without circular exclusions

**Status:** exact bijections, arithmetic exclusions and conditional certificates,
independently audited. Collatz convergence remains an open obligation here.
No literature-priority claim is made for classical valuation conjugacy.

Starting snapshot: `92514b16e2f6a13b7a1f194a6e17d5e5c5fc9ee7`.
Core inherited results are pinned to Session99 commit
`1e742862a5a442c87a6d7271cbc5e8305e2c4482`; incoming pre-audit additions are
not used as proved dependencies.

## 1. A genuine set of bijections that excludes the named 5x+1 hostile

For a valuation word w=(a0,a1,...), with prefix sums S0=0,Sj=a0+...+a(j-1),
set

\[
 \Phi_b(w)=-\sum_{j\ge0}\frac{2^{S_j}}{b^{j+1}}\quad\text{in }\mathbb Z_2.
\]

This exactly decodes odd b*x+1 dynamics, excluding only the countable sources
that eventually encounter an infinite valuation. Positive integer sources
are never excluded.

For each finite prefix u, its exact source cylinder has least positive odd
representative

\[
 r_b(u)=(2^D-C_b(u))b^{-J}\pmod {2^{D+1}},
 \qquad C_b(u)=\sum_{i=0}^{J-1}b^{J-1-i}2^{S_i}.
\]

These representatives are nondecreasing as ordinary integers. The decoded
source is a positive integer **iff the representatives remain bounded**, or
equivalently eventually become constant. This is an arithmetic predicate,
not a test that assumes convergence.

Let W3+ be the rays with bounded3-source representatives. We now have exact,
clock-preserving bijections

    positive odd Collatz sources <-> W3+ <-> Y=Phi5(W3+).

They also connect to Session99's odd-index pole readout via n<->4^n.
Every positive Collatz source is retained. The target Y lies in the2-adic
5x+1 system, not in its positive-integer subsystem. Its intended fixed point
is -1, corresponding to the Collatz fixed odd source1.

The hostile5-cycle13->33->83->13 has word(1,1,5). Under3-decoding it becomes

    19/101 ->79/101 ->169/101 ->19/101.

It is positive real, so **sign is not the separator**. Its reduced denominator
101 persists through every inverse prehistory. Thus its entire eventual
itinerary basin is rigorously absent from W3+ and Y. This is a finite
arithmetic exclusion, not a circular definition discarding counterexamples.
[Exact maps, inverse and source sidecar](bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md).

**Remaining obligation:** prove every bounded-residue3-ray is eventually
constant2, equivalently every point of Y eventually reaches-1. The bijections
alone do not prove this. They remove the named5-counterexample while locating
the integer-source content that still has to be established.

## 2. A sharp limit on simpler target designs

Every finite valuation word has infinitely many positive integer sources.
Therefore no finite-prefix forbidden-pattern rule can keep ALL positive
Collatz sources while excluding even one other infinite word. A successful
selective target must retain unbounded source arithmetic, a global acceptance
condition, or an additional theorem. This is a safety-filter obstruction,
not a prohibition on all automata or all possible proof methods.

A separate rigidity result says a rational change of coordinates preserving
both exact branches a=1 and a=2 cannot change the multiplier at all: it must
be the identity. The2-adic parameter conjugacy is not such a rational map.

A small affine-chart test detects the lost parameter: for the three words
112,121,211, the ratio of successive carry contrasts is exactly b/2.
All two-slot order differences are independent of b. Thus a three-row guarded
atlas distinguishes3 from5 even after a common affine output normalization.
The rows have different native source cylinders and must not be silently
queried on one common integer.
[Atlas fingerprint and rational rigidity](atlas/AFFINE_ATLAS_FINGERPRINT.md).

## 3. The arithmetic guard distinguishes an actual controller trap

For a ternary queue branching when r divides the native valuation, use the
shifted potential (n+1/(b-2))(b/2)^Q. Its neutral boundary is

\[
 H=\frac{b^2+(b+2)2^r}{2^{r+2}-b^3},\qquad 2^{r+2}>b^3.
\]

This is exactly the fixed point of the balanced block(r,1,1).

- For3x+1 with r=3, H=49/5 lies below13, the first legal branch source.
  Every legal branch contracts; the auxiliary first-return procedure ends.
- For5x+1 with r=5, H=83 is an actual integer cycle entrance. The queue
  repeats Q->Q+2->Q+1->Q on83->13->33->83.
- We classify this auxiliary failure completely: the(5,5) controller runs
  forever iff it reaches83 while active. Input19 is an expanding exception,
  but its subsequent leaf-only orbit exhausts every finite queue.

The general proof reduces auxiliary termination to a finite set of legal
branch inputs at or below H. A completed-word slope bound cannot replace
this guard: it says nothing about paths that never complete.
[Finite exceptional-input theorem and exact trap](controllers/RETURN_CONTROLLER_TRAP.md).
This strengthens the return-controller framework, not repeated Collatz descent.

## 4. A smaller arithmetic target than full random genericity

For each positive odd3-source, with total valuation clock S_k,

\[
 \text{reaches1}\quad\Longleftrightarrow\quad
 \limsup_{k\to\infty}S_k/k>\log_2(10/3).
\]

The proof is a source-marked product bound: if1 is avoided, all odd states
are at least3. This identifies a sufficient clock target weaker than the
Haar mean2, while retaining the exact integer source. The random drift
threshold b=4 does distinguish3 from5, but a Haar-null exceptional set may
still contain a particular integer source.

A finite inequality2^S_k>(b+1/n0)^k certifies a descent below the ORIGINAL
source by time k. It does not necessarily certify a smaller endpoint:
the independent hostile starts at165 and ends at167 after17 steps, despite
an earlier descent. Source, endpoint and clock remain separate coordinates.
[Clock barriers, density bound and certified-basin refinement](drift/DRIFT_TARGET.md).

## Verification

Run `python3 verify_packet.py` in this directory. Five standard-library
mathematical checkers run normally and with optimization, compare outputs
byte-for-byte, and compare stored JSON receipt values (ignoring formatting). No checks depend on removable
Python assertions. The [independent audit](audit/INDEPENDENT_AUDIT.md) includes
exact cylinder precision, rational prehistories, arbitrary initial queues,
source/endpoint hostiles, and the rational-function proof review.

Finite checks validate implementations. The all-source and all-prehistory
quantifiers come from the proofs, not from a finite orbit census. Published
bytes, manifest hashes and relative links are checked separately. Empty CI
queries are not reported as a CI pass.
