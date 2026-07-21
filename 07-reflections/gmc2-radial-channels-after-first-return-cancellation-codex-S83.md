# GMC(2) after the first-return correction: radial channels, not isolated atoms

**codex-2026-07-21-S83.** Owner: synthesize past and incoming work, pull
often, and work creatively toward the two-dimensional Gaussian Moments
conjecture.

## Where the external problem now stands

The Gaussian Moments conjecture uses independent **real** standard Gaussians
and complex polynomial coefficients. Derksen–van den Essen–Zhao proved the
one-dimensional conjecture and the homogeneous two-dimensional case. Long's
2026 counterexamples make GMC false in every dimension at least three, leaving
two real dimensions as the unique open dimension.

In circular coordinates (Z=(X+iY)/\sqrt2), (W=\bar Z), the exact functional is

\[
\mathbb E[Z^aW^b]=\delta_{ab}a!.
\]

The repo's strong nullcone target NC2 says that all moments zero force strict
charge-one-sidedness. NC2 implies GMC(2), but it is stronger than the literal
GMC conclusion and is not known to be equivalent to it. Likewise the published
Jacobian indexing is the global/fixed-r implication GMC((2r)) to the
Special Image Conjecture in (r) variables; it does not currently give the
shortcut “GMC(2) implies JC(2).”

Rigorous closed regions include homogeneous polynomials, one-sided charge
support, the two-character/binomial stratum, the one-variable factorial/radial
case, and named bounded supports with complete resultant certificates. Incoming
commit `888c763ea` is one such valid certificate:

\[
P=aZ^3+bZ+cW+dW^3,quad M_2=M_4=0
\Longrightarrow M_6=466560(ad)^3
\]

on the all-nonzero branch. It is a local rung, not a universal cutoff theorem.

## What the audit changes

THM-1770 correctly observed that every balanced word at the least return level
is primitive. It incorrectly concluded that distinct primitive coefficient
monomials vanish separately. The moment equation is scalar, so they can cancel.
The exact star witness

\[
P=aZ^6+bW^2+cW^{18},qquad
M_4=4\cdot6!ab^3+4\cdot18!a^3c
\]

cancels with (abc\ne0). MISTAKE-211 and the resolved court case retract the
atomwise and star-closure consequences while preserving the primitive-word
lemma.

The correct formula for a monomial support is

\[
M_m=
\sum_{|r|=m,\ q\cdot r=0}
\frac{m!}{\prod r_i!}
\left(\sum a_ir_i\right)!\prod c_i^{r_i}.
\]

This exposes the missing state. Charge (q=a-b) determines survival under the
circle average, but radial height (h=a+b) determines the factorial. Two words
with the same charge can have different heights and phases. The faithful
Newton object is bivariate `(charge,height)`, not the Laurent charge polygon.

This also explains why three other attractive shortcuts do not finish the
problem:

1. A fixed-radius Duistermaat–van der Kallen theorem does not commute through
   the radial integral; cancellations across heights remain.
2. Holonomicity for each fixed polynomial does not provide an exact recurrence
   order or a uniform coefficient-family stopping bound.
3. A unique exposed factorial degree does not by itself dominate multinomial
   entropy. Near-face replacements can occupy a growing boundary layer.

## A new proved infinite slice

THM-2014 closes

\[
P=aZ+b(ZW)+cW
\]

for constants (a,c) and an **arbitrary-degree complex radial polynomial**
(b). The exact nullcone is

\[
b=0,\qquad ac=0.
\]

For (d=\deg b\ge2), leading coefficient \(\beta\), the proof establishes the
uniform asymptotic

\[
M_m=\beta^m(dm)!
\left[e^{b_{d-1}/(d\beta)}+o(1)\right].
\]

The (k)-pair charged correction loses (2d-1) factorial degrees while gaining
only two multinomial powers, so its normalized size is bounded by

\[
O\!\left(
\frac{(C m^{-(2d-3)})^k}{(k!)^2}
\right)
\]

uniformly for the entire range (0\le k\le m/2). This is the quantitative
separation the old “leading factorial” argument lacked. Degrees zero and one
close from exact exponential generating functions.

The proof identifies a useful general rule: a factorial face is safe only after
subtracting its entropy codimension. A two-power multinomial gain means degree
gaps at least three are directly separable; gaps zero, one, and two are the
resummation boundary.

## The corrected endgame

Exact computation now falsifies a universal (2R) cutoff but suggests the
rank-sensitive replacement

\[
M^*(S)\le(k-1)R(S)
\]

for a (k)-monomial support. A four-monomial radial-channel example has
(R=2), vanishes through four, and dies at six. Another has (R=3), vanishes
through six, and dies at nine. All tested trinomials through degree eight close
at (2R); all tested degree-at-most-four quads close at (3R); all tested
degree-at-most-three five-monomials close at (4R).

HYP-8765 therefore replaces “isolate each first-return atom” by:

1. localize at a chosen positive-negative coefficient pair;
2. order balanced relations by primitive return and radial-face defect;
3. pass from raw moments to cumulants or successive resultants, removing
   decomposable lower circuits;
4. prove the remaining radial-channel matrix is a nonzero confluent
   factorial-Hankel/Vandermonde determinant;
5. put every pair product in the radical, forcing the one-sided nullcone.

The factorial Hankel kernel (((r+s)!)) is strictly totally positive, so this is
a concrete determinant target rather than an appeal to termwise dominance.
The main risk is triangularity: one must prove that the chosen cumulant or
resultant ordering really removes every mixed circuit without introducing a
new same-defect collision.

## Assumption challenge and Tournament Analysis

I considered runners/charges, monomials, fixed radial shells, primitive
circuits, wall-crossing events, residues, Newton faces, Fourier modes, and proof
obligations as possible vertices. The charge quotient preserves exactly the
angular survival predicate but destroys height and phase; this session's
counterexample is the cost of that destruction. Radial channels and circuit
obligations are the least lossy current state.

There is no natural tournament on charges for this analytic problem: no
antisymmetric pair observable determines a many-term scalar cancellation.
Forcing one would hide the obstruction. A proof-obligation tournament is useful
only as navigation: orient (A\to B) when (A) must be discharged before (B).
Its tie Hamiltonian path is

`moment identity → circuit localization → channel determinant → pair radical
→ one-sided nullcone → GMC(2)`.

The incoming continuum-coordinate work reinforces the same methodological
lesson: a coarse coordinate can freeze while a large hidden shell still moves.
Here charge is the coarse coordinate; radial height and coefficient phase are
the hidden shell.

## Immediate next experiment

Build the localized cumulant matrices for the smallest unproved radial-degree-one
supports `{-2,-1,+1}` and `{-2,+1,+2}`. Compute their determinants symbolically
after quotienting earlier circuits, and compare them with minors of the
factorial Hankel matrix. A single uniform formula in the channel index would
turn the observed ((k-1)R) law into an induction and would be the first route
in this workspace that genuinely scales beyond bounded radial degree.
