# Independent audit of the reduced-scale LRC packet controls

**Verdict: PASS.** The complete declared bounded class consists of actual
primitive boxed 6+7 decoder equality entries with the full retained profiles.
The three first-denominator and incompatibility controls are exact. The
individual-interval count already proves the whole class strictly safe.
No new unresolved LRC class or universal closure follows, and the producer
does not claim one. No mathematical repair was requested.

This audit read the full producer proof and source, and the proved incoming
`05-knowledge/results/third_20260906_grid.md` and actual atlas in
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`.
It independently rebuilt the arithmetic without importing the producer or
using its transcript as a computational oracle.

## 1. Actual entry and all-parameter typing

The cofactor products, six-star V, seven-body U, primitive gcds and disjoint
prime supports all check exactly. The independent atlas is generated as
products of inert primes with exponents0,1,2 and sum bound356; it contains
the inherited5855 unordered coprime ratios. Its complete graph is the
claimed five-edge star and eight-edge seven-body graph.

For every1334<=t<=97096 with gcd(t,H)=1, the two crossing orientations are
excluded for all bounded coefficients, not just a selected relation:

- A pair of tV labels has gcd t*d. Its gcd with any U label is1, so the
  latter's nonzero coefficient is a multiple of t*d>Q.
- A U pair is D*p,D*q with gcd(D,t*v)=1. A nonzero outside coefficient is
  kD. Dividing a proposed relation by D forces
  |k|t*v<=Q(p+q), contradicted by t*min V>430Q.

The strict inequalities hold at1334 and persist with t. Both zero inner
coefficients are allowed, so mixed supports of size two are covered as
well. There are105 forward and126 reverse mixed triple supports. After
excluding them, every remaining relation is internal to a component's
weighted kernel. The actual bounded edge rows span those two kernels:
coordinate scaling turns them into incidence rows of connected graphs.
Thus both inclusions in W_(Q,3)=V_dec are justified and the rank is11.

The source's upper-endpoint box bound is valid. Disjoint prime support plus
gcd(t,H)=1 rules out equal labels across the two components. Every mixed
subset has gcd1. A subset of size at least7 contained in a single component
must be the whole U, also primitive. All retained complement words are
therefore exactly the all-one words, whose full bank membership was checked
against the pinned inherited JSON. This is a proof for the entire class,
not an inference from the three finite controls.

## 2. Exact packets and their limited obstruction

The referee uses unsafe-residue bitsets and enumerates every numerator at
every denominator below the claimed first one. It obtains:

|t|first denominator|complete first packet|one clearance|
|---:|---:|---|---:|
|1369|31|4,9,11,20,22,27|4/31|
|1373|29|2,4,5,24,25,27|3/29|
|1583|37|1,13,17,18,19,20,24,36|5/37|

For every denominator2 through28, a literal U label is divisible by it.
Both primitive packets at29 and31 are nonempty and all three scales are
units there. The full29 scale-residue profile matches the producer. The
located identity B_d(tV union U)={j in B_d(U): tj mod d in B_d(V)} is exact.
In particular1583 misses both of those prescribed useful U clocks.

The failure is a fixed-clock compatibility failure in actual entries. It
does not obstruct adaptive selection of another denominator, and the
physical31/29/37 witnesses already show safety.

## 3. Whole-class grid safety and the rounding distinction

For the U pair with primitive coordinates81,101, an independent enumeration
of all8181 pairs of danger-arc centers reconstructs25 open circle components
and their complete length multiset. Their sum is389/19089. The common
physical factor is coprime to t, so multiplication by it permutes any
translated t-grid; passing to the primitive pair clock is legitimate.

An open circle interval of length ell has at least ceil(t*ell)-1 grid
points. Summing these individual counts gives13 at1334 and a nondecreasing
function of t. The seven marginal dangers contain at most7ceil(t/7)
points in total because every gcd(t,u_i)=1. Subtracting this selected pair
intersection is a valid pointwise union correction. Hence there are at
least13-E>0 safe grid points, E=7ceil(t/7)-t<=6, throughout the class.

At1369 and1373 all21 primitive U-pair geometries were rebuilt by literal
center-pair enumeration. The sole positive pooled pair credit comes from
81,101 and is respectively55316/19089 and56872/19089. Even the sum of all
positive pooled credits fails against E=3 and6. Ceiling only the total
length still fails. The individual interval credit is15 in both cases.
Thus the difference comes from retained rounding information, as claimed.

All V labels are odd, so the half-phase lifts (2j+1)/(2t) preserve clearance
1/2 on V. Since7 divides H but is coprime to t, no U coordinate on these
lifts can have distance exactly1/14: that equality would imply
7u(2j+1)=t(14k+/-1), impossible modulo7. The weak-safe grid points are
therefore strictly safe. The independent integer engine obtains exactly
462,458,530 such lifts and verifies literal full-row strict clearance on
all1450. This proves the advertised stopping boundary.

The referee retains a noncoprime hostile: the unshifted7-grid is entirely
dangerous for speed7, contradicting a coprime marginal bound of1. The
individual clock gcd assumptions are used, not imported from global
primitivity.

## 4. Reproduction and frozen identities

The [independent source](../../04-computation/continuing4_20260906_lrc_packets_audit.py) and
[transcript](continuing4_20260906_lrc_packets_audit.out) pass **14720
always-active exact gates** in normal and optimized Python. Both outputs
are byte-identical actual LF, with no newline postprocessing.

The explicit finite universe includes all5855 atlas ratios; all693 mixed
support checks and12285 full profile instances across three actual rows;
the complete denominator banks; all21 literal center-pair overlap geometries;
and1450 literal safe half-phase lifts. Weighted rank uses fraction-free
integer elimination, distinct from the producer's rational row reduction.
Profiles are enumerated by omitted-label masks, and packets by unsafe bitsets.

From the repository root:

    python 04-computation/continuing4_20260906_lrc_packets_audit.py --profiles PATH_TO_PROFILES
    python -O 04-computation/continuing4_20260906_lrc_packets_audit.py --profiles PATH_TO_PROFILES

The optional default first looks beside the audit source (portable when
filed in04-computation), then in the current C:/w/s0905 computation directory.
The data pin is on actual raw LF bytes. The audit pins the producer source
and transcript but imports neither.

- Producer source: `22db0faeddf433aed8d1daeb664976813e4ed8a0b5a0fec57fdeaac027a23b08`.
- Producer output: `46c71d72aefab9d3e120c58317cc616c8ad37c1db3b7de9e2335c887a8fae151`.
- Producer report reviewed: `cf8aff56ad020e232597253549ffb87f3e98eb235aae93741b7fe4881f4debb5`.
- Inherited profiles: `935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
- Audit source: `c5008b8e176425da3ad5f63769efe1a958db15c42c7c968260ab6c65a872a88d`.
- Audit output: `54f832cb9edf4f0b8997359235a507bf19eb0b3c1d5b4181b28a5516871cb2a1`.

All work remained outside the repository; no producer, maintained file or
Git state was edited by this referee. The accepted conclusion is the stated
actual-entry construction and information-loss boundary, with inherited grid
safety of the full displayed class.
