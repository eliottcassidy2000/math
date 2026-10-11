# Atomic source weights with distortion independent of composition depth

2026-10-11. **PROVED, scoped and self-audited:** C1–C3 below; C4 specializes
an inherited obstruction to the new priors. Native
receipts are **INHERITED** from the linked packets; finite replays are
**FINITE-EXACT**. Universal entry and independently positive ROOT scores
at every source remain **OPEN**.

## Inheritance and the precise extension

The [new bridge law](collatz_bridge_weight_transport_20261011.md) already
compares actual source weights with child weights on one bridge family,
uniformly under residue refinement. Its conditional distortion is below
1.078 at r=6. The [atomic-rank interface](collatz_atomic_rank_interface_20261011.md)
already pays arbitrary legal episodes and explicitly leaves their stopped
ROOT obligations open. We do not repeat either result as new coverage.

The present extension compares **arbitrarily long finite compositions**
of four kinds of contracting dependencies. Repeated multiplication of
one-step error factors is unnecessary: the affine offset of the composite
is controlled by one interval of fixed points.

Closest proved mechanism: bridge weight transport. Canonical hostile:
the authenticated neutral loop283->251->283. Corrected near miss:
conditioning on a guard without transporting the source law. Least-used
sidecar: the affine offset, regarded as a weighted average of fixed
points. Board: **native word / affine anchor / source atom / conditional
law / independent score / marked terminal**.

The map takes a native receipt word to its composite affine map. It
preserves the source-child equation and the weight estimate; it loses the
guards and the endpoint ROOT certificate. Both remain compulsory sidecars.
An arbitrary word of the affine maps is not automatically native.

## C1. A common interval contains all composite anchors

Represent an odd integer n by X=(n+1)/2. For a dependency n->h put
Y=(h+1)/2. The inherited contracting bank is:

| Dependency | Native source | Y=alpha X+beta | Fixed point gamma=beta/(1-alpha) |
|---|---|---|---|
| H | n=155 mod2048 | alpha=729/1024, beta=241/512 | 482/295 |
| R_r | n=3 mod4, r>=4, 3^r divides2n+1 | alpha=2^(r+2)/3^r, beta=-alpha/4 | -alpha/[4(1-alpha)] |
| bundled GR | n=251+288t, t>=0 | alpha=8/9, beta=0 | 0 |
| J_j, the incoming atomic G-to-run bridge | n=(2264K-23)/27+96Kt, K=2^(18j), j>=1,t>=0 | alpha=8/(3K), beta=16(K-1)/(81K) | 16(K-1)/[27(3K-8)] |

The name J_j here distinguishes the last **atomic bridge** from the
expanding letter G=(9n+5)/8. It does not rename its source document.
Receipts and guards are in the
[run credit interface](collatz_run_credit_interface_20261011.md) and
[atomic run bridge](collatz_atomic_run_bridge_20261011b.md).

Every alpha lies in (0,8/9]. The R_r anchor is at least -16/17, with
equality at r=4. The J_j anchor lies between0 and1/5, because K>136.
Consequently every anchor belongs to

\[
                  I=[-16/17, 482/295].                           \tag{1}
\]

For a nonempty finite composition of these affine maps,

\[
          Y=A X+(1-A)\gamma,qquad 0<A\le(8/9)^\ell,quad\gamma\in I,
                                                                    \tag{2}
\]

where ell is the number of maps. Indeed composing a map with slope a
and anchor g after a map with slope A and anchor gamma gives numerator
a(1-A)gamma+(1-a)g for the new anchor; its denominator is1-aA, the sum
of the two nonnegative coefficients. Thus the anchor is a convex average.
Induction proves (2), without any growth of the anchor interval.

For an actual Collatz dependency, retain every native guard and common-
future receipt before applying (2). The inherited two-R4 ray
13891+26244t ->10975+20736t ->8671+16384t is one genuine repeatable
example. The affine envelope also includes words whose guards are empty;
it makes no existence claim about those words.

## C2. One distortion bound, regardless of word length

Let nu_r(m)=(m+1)^(-r)-(m+2)^(-r), r>0. Suppose a native composed word
sends source index a to child index b, so X=a+1, Y=b+1 in (2), and
Y>=56 (the child odd integer is at least111). Define

\[
                 l=5346/5605,qquad u=121/119.
\]

Then

\[
 \boxed{u^{-(r+1)}A^{r+1}\nu_r(b)
          \le\nu_r(a)\le
        l^{-(r+1)}A^{r+1}\nu_r(b).}                              \tag{3}
\]

**Proof.** For z in [0,1], solve (2) for X to obtain

\[
 \frac{A(X+z)}{Y+z}
       =\frac{Y-(1-A)\gamma+Az}{Y+z}.                              \tag{4}
\]

Its derivative in z is (1-A)(gamma-Y)/(Y+z)²<0. At z=0 it is at most
1+(16/17)/56=121/119. At z=1 it is at least
(56-482/295)/57=5346/5605. Integrate the resulting reciprocal-power
inequalities in
nu_r(m)=r integral_0^1(m+1+z)^(-r-1) dz. This proves (3). QED.

These constants are independent of word length, the individual run/bridge
parameters, and every later restriction of the native domain. The bound
is weaker than the specialized1.078 bound for one J family; it adds a
different strength, uniformity under arbitrary allowed composition.

Fix one composed word, hence a fixed slope A, and any nonempty set B
of its admissible child indices >=55. An affine map with positive slope
has at most one source for each child. Normalize nu_r on B and normalize
its pulled-back source weights on B. For every subset E of B, or every
nonnegative score V on B, the induced normalized measures obey

\[
 D_r^{-1}\,\mathbb E_B[V]\le\mathbb E_{\mathrm{source}|B}[V]
                        \le D_r\,\mathbb E_B[V],qquad
 D_r=(u/l)^{r+1}.                                                \tag{5}
\]

This follows by dividing the lower numerator bound by the upper total
bound, and conversely. At r=6,

\[
                        D_6=1.564932104\ldots<1.565.             \tag{6}
\]

The set B may be a singleton, an arbitrary deep residue cell, or an
irregular guard refinement. No probability independence assumption is
used. For a union of different words, slopes differ and sources may
overlap: retain a disjoint address partition and the factors A^(r+1)
before summing. Equation (5) is not an unqualified union formula.

## C3. What this supplies to the independent ROOT-floor compiler

For any score V(m) in [0,1], the earlier source-floor theorem R6 in
[renewal/source weights](collatz_renewal_source_20261011.md) gives, for q>m,

\[
 W_q(m)=\sum_{k\ge0}\nu_r(m+kq)V(m+kq),qquad
 E_q(m)=\sum_{k\ge1}\nu_r(m+kq)\le(r+1)q^{-r-1}.
\]

Thus an independently established lower bound W_q(m)>=L_q implies

\[
             V(m)\ge\frac{L_q-(r+1)q^{-r-1}}{\nu_r(m)}.          \tag{7}
\]

The numerator must be strictly positive. For r=6 its error is7/q^7;
the sharper scaled error tends to6*zeta(7), an inherited series identity
that does not require irrationality of zeta(7).

Equations (3)–(5) can now transport a supplied child score through long
native words without losing a factor D_6 at every step. This removes one
quantitative obstruction to a refinement-stable floor. It does not
establish the child floor, and it does not identify a merely positive
component score with ROOT. The acceptance condition must remain a literal
checked terminal or an independently justified ROOT event.

Equation (5) transports the measure of a specified child function. For
the ROOT indicator, native common-future coalescence identifies the source
and child values exactly. For a discounted or randomized controller score,
one must additionally retain its routing probability and clock multiplier
before identifying that transported function with the actual source score.
Weight distortion alone does not preserve those different score semantics.

The mechanism also says exactly why two tempting shortcuts fail:

* More letters shrink A^(r+1). Bounded relative distortion does not create
  extra unnormalized source mass or a positive score from a zero score.
* An episode may stop at a state outside all current guards. No amount of
  transported mass authenticates that state as ROOT. Adding inverse G and
  expanding G without a new payment rule admits the neutral283/251 loop.

**Next independently testable questions.** Can a new guard attach each
stopped component to a child family with a proved positive score? Can an
estimate of W_q use native incidence and its exact carry to beat7/q^7
without first assuming ROOT(m)? For the current bank, can the required
disjoint word-address partition be maintained while adding genuinely new
source residues? These questions retain the target rather than replacing
it by the size of an unmarked equivalence class.

## C4. Keep the source law separate from a global flow weight

The stronger existing result P6–P8 of
[effective prefix mass](collatz_effective_prefix_mass_20261005.md) says
universal ROOT is equivalent to a strictly contracting summable incoming
flow, and already excludes bounded corrections of its original atomic
prior. The same hostile applies to the present nu_r, not just the older
prefix code. For every integer d>=1 the actual word of d valuations equal
to one sends

\[
             2^{d+1}-1\longrightarrow 2\,3^d-1.
\]

The corresponding source indices are 2^d-1 and 3^d-1, so the integral
bounds for nu_r give

\[
 \frac{\nu_r(2^d-1)}{\nu_r(3^d-1)}
       \ge\left(\frac{3^d}{2^d+1}\right)^{r+1}\longrightarrow\infty.
\]

If c nu_r<=w<=C nu_r for fixed positive c,C, then a discounted incoming
flow inequality would imply along that word
w(source)<=rho^d w(endpoint), rho<1, contradicting the displayed bound.
Thus no bounded positive correction of this prior can provide that
particular universal discounted-flow certificate.

This is an inherited mechanism specialized to the new weight, not a new
no-go against all positive scores. The promising state must keep the
source prior for measurement and a separate, potentially unbounded
climb/refuel account for flow payment. C2 concerns contracting authenticated
dependencies and is consistent with this expanding actual-orbit hostile.
The prior packet's9->7,41->31,169->127 refuel family remains a required
control for any proposed new account.

The [Vitali companion](vitali_atom_recovery_20261011.md) explains why the
full refinement sequence may recover a target erased by one finite
quotient. The [planar companion](mobius_affine_source_return_20261011.md)
provides a positive example of restoring the missing affine offset. The
analogy guides the choice of state; it does not transfer a planar theorem
into a ROOT proof.

## Reproduction

[Program](../../04-computation/experiments/mobius_source_anchor_weights_20261011.py),
[output](mobius_source_anchor_weights_20261011.json). The shared16,872 gates
include1,463 algebraic words of length1–3 (explicitly not all native),
250 actual or inherited source-child pairs, four exact integer r values,
and independent valuation replays for R, GR, J, and the two-run family.
H's native receipt is inherited; its affine index relation is rechecked.
The generic proof handles arbitrary length and real r>0. Normal and
optimized Python outputs agree byte for byte. Universal entry and
independently positive scores have not been proved by this experiment.
