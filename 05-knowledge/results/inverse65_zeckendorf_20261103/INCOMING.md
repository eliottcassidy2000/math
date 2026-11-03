# Concurrent inverse results and dual error/charge channels

The first checkpoint was rebased onto4a08185e87 and pushed as d2f73cbbe5.
Selected incoming claims are independently replayed here. No priority is
claimed for independently overlapping inverse results.

## N1. Incoming mechanisms and their boundaries

[golden_inverse_defects/THEORY I1-I5](../golden_inverse_defects_20261103/THEORY.md),
commit f6d5304849, independently proves the same65-bank realization inverse.
Its physical defects are `(d,e,f)=(x+y-h,y-x,h+x-y)`. It supplies finite fibre
polynomials and a stronger branch classification: no retained source supports
both3h+2 and3h+4 by same-end appendings. Its66 words in41 states with d=1 or2
independently agree with our maximum-H155 state census. Our norm9/5 improves
the per-symbol depth estimate over its two-step factor3; it is not a claim
of first discovering a realization inverse.

Its [TRANSPORT T1-T5](../golden_inverse_defects_20261103/TRANSPORT.md) replaces
a bounded source suffix. In our processing convention its seeds are

    U=r5 a a, U2=a a q1 a a, U4=q0 a a a a a,
    v(U)=(31,14,25), v(U2)=(95,41,75), v(U4)=(97,40,75).

The differences from3v(U) are tau=(2,-1,0) and2tau. Appending b reverses
tau; q2 preserves it. Even-b continuations give both positive offsets under
this explicit source rewrite. Its four proposed next states
`(287,122,225),(289,121,225),(293,119,225),(295,118,225)` are absent; our
inverse independently rejects all four. The construction is not closed.

The [adaptive packet](../adaptive_metagraph_reformulation_20261103/README.md)
supplies bounded-port carriers and distinguishes object reconstruction from
target-integer realization. Its repaired
[negative-offset descent review](../adaptive_metagraph_reformulation_20261103/INCOMING.md)
requires the correct threshold: for h'=a h+b, increasing h requires
`h>-b/(a-1)`. The negative-offset examples below respect this correction.

## N2. Dual error and charge coordinates have an integral gluing index

Let primitive integer vectors ell (column) and tau (row) satisfy
`P ell=sigma ell`, `tau P=sigma tau`, with nonzero `c=tau ell`. Then

    Pi=ell tau/c,  Pi^2=Pi,  P Pi=Pi P=sigma Pi.

A row v splits rationally into `(v ell/c)tau` and `v-(v ell/c)tau`, the
second having zero charge. Both parts are integral iff c divides v ell,
since tau is primitive. Moreover

    [Z^3 : (ker ell intersect Z^3) + Z tau] = |c|.     (N1)

Proof: the primitive charge map Z^3->Z is onto, its kernel is the first
summand, and Z tau maps to cZ. This elementary integral linear algebra does
not assert that either rational projection is an actual source state.

For the seven charge symbols in C1:

| Letter | Primitive error tau | c | Gluing index |
|---|---|---:|---:|
| b,q2 |(2,-1,0)|2|2|
| q1 |(3,-2,1)|5|5|
| q4 |(3,-2,0)|5|5|
| q5 |(2,-2,1)|-2|2|
| r21,r38 |(3,-1,-1)|13|13|

The concurrent b/q2 error is dual to our charge h-y; factor2 is its parity
gluing. For the norm13 family, J=2x+2y-3h gives charge component
`-J(3,-1,-1)/13` and first bulk coordinate `h+3J/13`. Clearing the
denominator gives precisely U=13h+3J from C3. The integral source lattice
explains that normalization independently of fitting terms. The set{2,5,13}
classifies these seven matrices; no Fibonacci law for larger atoms is inferred.
`integration_probe.py` verifies all projectors and their minimal denominators.

## N3. Actual negative-offset source replacements

For r21/r38, tau=(3,-1,-1) obeys tau P=tau; `(tau P_b)_h=-5`. Thus actual
sources with v(V)=3v(U)+m tau give, for every word w over{r21,r38},

    H(V w b)=3H(U w b)-5m.                            (N2)

Two actual processing-word pairs are:

| m | U / v(U) | V / v(V) | Initial H(Ub) -> H(Vb) |
|---|---|---|---|
|2|r4 a a / (39,16,31)|a a r27 a / (123,46,91)|227 ->671=3*227-10|
|4|a r12 / (33,33,39)|a a a r10 / (111,95,113)|321 ->943=3*321-20|

After k transport symbols and b, vertex clocks are10+5k and11+5k for the
first pair,9+5k and11+5k for the second. Words and ports remain explicit,
but V replaces U: this is not an induced extension or an asserted single
native metagraph edge. The retained spine closes into a directed Hamilton
cycle, proving strong connectivity.

The base sequence M_k=H(U r21^k b) obeys

    M_(k+2)=40M_(k+1)-13M_k-10J,
    (J,M0,M1)=(-23,227,9215) or (45,321,12081).         (N3)

Their generating functions have common denominator
`D(z)=(1-z)(1-40z+13z^2)` and respective numerators
`227-92z+95z^2` and `321-1080z+309z^2`. The target series are
`3M(z)-10/(1-z)` and `3M(z)-20/(1-z)`. This follows by summing(N3), with
the literal vertex clocks kept separately as specified above.

Literal subset DP checks both pairs at k=0,1, through16 vertices. Initial
targets671=11*61 and943=23*41 are composite. No prime-generator claim follows.

The offsets supply two ternary residue classes modulo6 in **different
families**. The complete H<=2000 census has20703 states, two separate m=2
lifts, two separate m=4 lifts, and none supporting both. This is finite
evidence with every participating state restricted to H<=2000, not an
all-order obstruction or an exclusion for all seeds H<=2000 with larger
outputs. No return rule closes these families.
Negative-offset descent itself is valid: H>10 implies3H-10 and3H-20 exceed H.

Next seek a common seed or source-compatible return, retaining joint state
and13-adic clock. Any finite scalar affine encoding must also pass the
[incoming return-cost spectral test](../session153_gadgets_20261103/rewriting/NOTE.md):
local residue coverage can still be defeated by an expensive return cycle.
