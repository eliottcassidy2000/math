# Paley geometry, projective sheets, and certified finite-observation aliases

2026-10-04. **PROVED** elementary transport, two-bit quotient and all-modulus
alias statements below. **FINITE-EXACT** for the declared group and word
censuses. **CITED** for classical switching-class context. Universal Collatz
convergence is **OPEN**. No literature-priority claim is made.

The positive construction is a finite-field tournament representation of
every legal Collatz word. Its sharp limitation can also be proved: all
quadratic Paley orientation probes retain only two exponent parities. A
stronger, explicit construction gives two already-convergent integers with
identical arbitrarily prescribed finite residue histories but opposite net
drifts at one common later time. Their route certificates are constructed,
not assumed from the conjecture.

## 1. Inheritance and the carrier choice

Closest mechanisms are the affine Paley action in
[the Fano/design note](collatz_mod6_20260922_paley_fano_octonion_design.md),
the guarded word algebra in [the blueprint audit](collatz_blueprint_20260921_affine.md),
and the real projective-sheet construction of
[THM-4464, checksum projective tournament](../../01-canon/theorems/THM-4464-checksum-projective-tournament-and-output-sheet-cocycle.md).
The latter already proves that pairwise determinant data can forget a common
lift sign. The finite-field determinant construction below is its appropriate
arithmetic analogue, not a newly discovered general switching principle.

The corrected near miss is the earlier identification of the Paley period3
code with a general Collatz symmetry; see
[the audited Paley bridge](collatz_paley_bridge_20261001.md). Its translation
symmetries are not global symmetries of the integer Collatz function. The
canonical hostile here is the pair of formal valuation words `(1)` and `(3)`:
their real multipliers have opposite drift, but every quadratic multiplier
character agrees. The least-used sidecars are the projective observer,
the retained orientation sheet, and the chosen word's exact domain.

Classical context: Babai--Cameron,
[Automorphisms and Enumeration of Switching Classes of Tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v7i1r38/pdf/),
Electronic Journal of Combinatorics7 (2000), R38, Theorems5.1--5.2,
classify the groups admitted by tournament switching classes. Those results
are not claimed here. All finite-field formulas used below are derived
directly, and the order-eight case is checked exhaustively.

Live board: legal affine word; intrinsic pair orientation; projective
observer; signed lift; real drift/axis; certified finite history. The
Anchor is the Collatz observation boundary, the Niche is the finite
projective action, and the Wildcard is the signed tournament-adjacent cover.

## 2. The intrinsic tournament and its projective switching law

Let p be a prime with `p=3 mod4`. On the projective line
`X=Fp union {infinity}`, choose

\[
 v_x=(x,1),\quad v_\infty=(1,0),\qquad
 S(x,y)=\chi_p(\det(v_x,v_y)).                       \tag{1}
\]

The quadratic character is nonzero for distinct projective points, and
`chi_p(-1)=-1`; hence this is a genuine tie-free tournament. Its vertices,
observable and orientation gauge are explicit. Infinity dominates the
finite points, which carry a Paley tournament in the convention
`x->y iff chi_p(x-y)=1`. Reversing the convention reverses all arcs.

For any invertible two-by-two matrix M over Fp, write

\[
 Mv_x=\lambda_x v_{g(x)},\qquad
 \epsilon_x=\chi_p(\lambda_x),\quad
 \delta=\chi_p(\det M).
\]

Taking determinants proves

\[
 \boxed{S(gx,gy)=\delta\epsilon_x\epsilon_yS(x,y).} \tag{2}
\]

Thus a square-determinant projectivity acts by vertex switching; a nonsquare
determinant adds global reversal. Scalar changes of the matrix representative
do not change delta and change every epsilon by a common sign, which cancels
in(2). The cyclic triple product `S(x,y)S(y,z)S(z,x)` removes vertex switching
but retains delta. It is an *oriented* triple invariant, not an unoriented
subset test.

At p=7, exhaustive permutation checks give:

| Retained object | Permutations |
|---|---:|
| The chosen normalized eight-vertex tournament |21|
| Its switching class |168|
| Permutations reversing its oriented triple invariant |168|

The168 switching permutations are exactly the projective square-determinant
matrices. With

\[
 S_0=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
 R_0=\begin{pmatrix}0&-1\\1&1\end{pmatrix},
\]

their projective orders are2 and3, their product has order7, and they
generate all168 permutations. This is an explicit finite quotient of the
`(2,3,7)` triangle-group presentation. The geometry comes from that specified
presentation and its surface realization;168 alone is not a genus.

The full group cannot preserve a single tournament faithfully: a nontrivial
involution would exchange some two vertices and reverse their arc. The
switching class retains a weaker predicate, and that distinction permits
the even-order projective group.

### Retaining the sheet gives a faithful signed cover

Use vertices `(x,e)` with `x in X,e=+-1`, and pair observable
`e*f*S(x,y)`. The two vertices over the same x are tied, so this is a
directed complete multipartite graph, not a tournament. At p=7 it has
16 vertices,112 oriented unordered pairs and8 antipodal ties.

Determinant-one matrices act by

\[
 (x,e)\longmapsto(gx,e\chi_p(\lambda_x)).            \tag{3}
\]

This preserves the complete observable. The336 matrices of SL2(F7) act
faithfully; the central matrix `-I` fixes all projective addresses and
reverses every sheet. Keeping ties is essential: forcing orientations on
the antipodal pairs would destroy that involution's invariance.

At p=5, `chi_p(-1)=+1`, so(1) is symmetric rather than antisymmetric.
It supplies a graph-sign construction, not a tournament. The spherical
five and hyperbolic seven cases therefore cannot be joined by silently
using the same directed carrier at both primes.

## 3. A Collatz word acts on Paley arcs, but the orientation sees two bits

For a positive valuation word `w=(a_1,...,a_j)`, let

\[
 A_i=\sum_{r=1}^i a_r,\quad P=3^j,\quad Q=2^{A_j},\quad
 B=\sum_{i=0}^{j-1}3^{j-1-i}2^{A_i}.
\]

On its legal exact cylinder the odd Collatz endpoint is

\[
 F_w(n)=\frac{Pn+B}{Q}.                              \tag{4}
\]

For p different from2,3, this rational affine map reduces to a permutation
of Fp, and for distinct x,y,

\[
 \boxed{\chi_p(F_w(x)-F_w(y))
 =\chi_p(3)^j\chi_p(2)^{A_j}\chi_p(x-y).}           \tag{5}
\]

The carry cancels in the difference. The formula holds for actual integers
sharing the exact valuation word, after reduction; it is not a claim that
the piecewise function U itself is one permutation modulo p.

For p congruent to3 modulo4 and p>3, quadratic reciprocity gives the table
(also directly checked at the displayed prime representatives):

| p mod24 | Orientation multiplier |
|---:|---|
|7|`(-1)^j`|
|11|`(-1)^A_j`|
|19|`(-1)^(j+A_j)`|
|23|`+1`|

Primes7 and11 recover both parities. Consequently the complete collection
of these quadratic orientation multipliers, even across *all* such primes,
is exactly a two-bit quotient of the word. It loses B, exponent magnitudes,
the real multiplier's comparison with1, and source position. This is a
statement about orientation *multipliers*, not about all labelled finite
field images or every possible arithmetic observable.

The words `(1)` and `(3)` have `j=1` and odd A. Their multipliers are
`3/2>1` and `3/8<1`, while(5) agrees at every prime. The words `(1,2)`
and `(2,1)` both have `(P,Q)=(9,8)` but carries5 and7 and fixed points
`-5` and `-7`. Even the signed variant `(3n-1)/2^a` has the same
determinant character as `(3n+1)/2^a`; the arithmetic carry sign must also
be retained.

Here the arithmetic carry sign `B -> -B` is distinct from the projective
representative sheet `M -> -M`. The latter leaves the rational affine map
unchanged. The matrices `[[P,+B],[0,Q]]` and `[[P,-B],[0,Q]]` have the same
sheet cocycle as well as the same determinant character, so retaining the
projective sheet alone does not recover the arithmetic carry sign.

Every matrix in(4) fixes infinity. These word maps remain in an affine
point-stabilizer; the extra projective transformations in section2 require
changing the observer and are not generated by Collatz words. This is the
same missing-operation boundary seen in
[THM-2626, the Paley--Borel frame](../../01-canon/theorems/THM-2626-paley-borel-projective-frame-torsor-and-physical-c13-boundary.md),
which provides a richer frame but no automatic physical transition.

## 4. An all-modulus hostile made entirely from completed routes

The following strengthens a same-character example into actual integer
Collatz trajectories. It does not assume unknown sources converge.

**PROVED.** For every integer `M>=1` and observation depth `D>=0`, there
exist distinct positive odd integers n,n', a common odd endpoint u and
an integer j>D such that:

* both reach u after exactly j accelerated steps, and u then reaches1;
* `n<u<n'`, so their net changes at time j have opposite signs;
* `U^i(n)=U^i(n') modM` for every `0<=i<=D`;
* their first D valuation exponents agree, all equal1;
* both have first-hit odd rank j+1.

The observation includes vertices0 through D and edges1 through D. It
does not include the outgoing valuation at vertex D, which can be the
first differing edge.

### Explicit construction

Factor `M=2^h*3^k*m`, where `gcd(m,6)=1`, and put

\[
 j=D+h+1,\quad
 u=\frac{2^{3^j+1}-1}{3},\quad
 L=2j\varphi(3^{j+k}m).
\]

Here `varphi` is Euler's totient, with `varphi(1)=1`; the expression
inside it is at least3. Equivalently,
`L=4j*3^(j+k-1)*varphi(m)`. Set

\[
 \boxed{n=\frac{2^j(u+1)}{3^j}-1,\qquad
 n'=\frac{2^j(2^L u+1)}{3^j}-1.}                    \tag{6}
\]

Since j>=1, lifting the exponent in `2^(3^j)+1` gives
`v3(u+1)=j`. Euler's theorem gives

\[
 2^L\equiv1\pmod{3^{j+k}m},\qquad 2^L>3^j.
\]

Both quotients in(6) are integers. The quotient `(u+1)/3^j` is positive
and even, so n>=3; the quotient `(2^L u+1)/3^j` is positive and odd.
These parity statements give the exact valuation words

\[
 n\xrightarrow{\;1^j\;}u,
 \qquad n'\xrightarrow{\;1^{j-1}(1+L)\;}u.          \tag{7}
\]

For example, iterate `x+1 -> 3(x+1)/2` while the prescribed valuation
is1. At the last edge of the second word, direct substitution gives
`3x+1=2^(1+L)u`, proving the exact final valuation. Also
`3u+1=2^(3^j+1)`, so U(u)=1. Every earlier value is odd and greater
than1; the first-hit assertions follow. The ordinary first-hit ranks are

\[
 2j+3^j+2,\qquad 2j+3^j+2+L.                       \tag{8}
\]

The first inequality n<u follows from
`n+1=(2/3)^j(u+1)`. For the second, writing `B=3^j-2^j`,

\[
 n'-u=\frac{(2^{j+L}-3^j)u-B}{3^j}>0
\]

because `2^L>3^j` and u>=1. Finally, while i<=D<j, both paths have
performed only valuation-one edges, so

\[
 U^i(n')-U^i(n)
 =\frac{2^{j-i}3^i u(2^L-1)}{3^j}.                 \tag{9}
\]

The expression is divisible by `3^k*m` by the choice of L, and by
`2^h` since `j-i>=h+1`. This proves every required congruence.

The smallest parameter control `M=1,D=0` gives `n=3,n'=53,u=5`, with
both sources mapping to5 and then1. At larger M and D the same proof
constructs the pair without searching for trajectories.

### The stored object and the boundary of this obstruction

The implementation builds the two certificates using
[the parent/row/block codec](inverse_ray_ternary_addresses_20261004.md):
construct u over root1 with exponent `3^j+1`, prepend either exponent1
or1+L, then prepend j-1 more ones. The full certificates retain exactly
the later differing exponent. Mixed residue queries use the binary and
ternary evaluators and an explicit coprime-to-three evaluator with CRT.

The huge control takes `M=2^64*3^12*35,D=10,j=75`. Both sources have
more than `10^35` binary digits; they are never expanded. Their finite
histories and certificates are computed symbolically. Small controls
independently expand and replay both integers.

The theorem rules out determining this later net-drift predicate from a
fixed modulus and fixed initial observation depth alone, even after the
common endpoint and common odd first-hit rank are supplied. It does not
rule out adaptive precision, an unbounded observation, source size, the
retained valuation word, or a Collatz proof using those data. Both examples
already converge, so this is not a divergence construction or a convergence
counterexample. It gives no new density claim beyond inherited inverse
completion families.

## 5. Connection ledger and reproducible universe

| Source / map | Preserved target | Lost data / repair | Decisive hostile |
|---|---|---|---|
| Projective points plus oriented representatives -> determinant tournament | Intrinsic pair orientation | Choice of sheet; retain signed cover | -I fixes addresses but reverses sheets |
| Projectivity -> switching action | Oriented triple relation when determinant square | A chosen normalized tournament |168 switching symmetries versus21 strict ones |
| Legal word -> Paley arc multiplier | Parities of j,A, exactly | Carry, magnitudes, source and sign sheet | Words1 versus3;12 versus21 |
| Completed route -> fixed modular initial history | Prescribed D-step observations | Later branch exponent and source height | Formula(6), both with known routes to1 |
| Full word and marked source -> affine geometric carrier | Exact endpoint and descent comparison | Nothing needed beyond declared arithmetic guards | Companion [axis note](geometry_collatz_drift_carriers_20261004.md) |

Reproduce from the repository root:

```text
python -X utf8 -B 04-computation/experiments/paley_geometry_observations_20261004.py
python -O -X utf8 -B 04-computation/experiments/paley_geometry_observations_20261004.py
```

The [output](paley_geometry_observations_20261004.out) records all336
projective matrices at p=7, every permutation of8 points, all336 SL2
signed-cover actions, and the p=5 symmetric hostile. The word universe
is all340 words of lengths1--4 with entries1--4, tested at seven named
primes and every finite pair. The alias universe is `M=1..48,D=0..4`
(240 pairs), with122 literal expanded controls under the declared4096-bit
cap and one unexpanded huge control. All checks use exact integer arithmetic
and explicit exceptions under normal and optimized Python. The all-height
theorems are proved above rather than inferred from this census.
