# Finite moment feasibility and a signed localizer for one source

2026-10-05. **PROVED:** the conditional localizer floor, finite convex
optimization, rational dual verification, and bounded tail relaxation.
**FINITE-EXACT:** the declared synthetic measurement controls.
**OPEN:** independent measurement inequalities for every Collatz source.
No target ROOT word, orbit search, or completed-weight bank is used by this
package. Synthetic measurements demonstrate the certificate mechanism;
they are not asserted to be measurements of the Collatz distribution.

[Script](../../04-computation/experiments/collatz_moment_localizer_feasibility_20261005.py)
and [output](collatz_moment_localizer_feasibility_20261005.out).

## 1. Inheritance and the new finite obligation

The closest mechanisms are:

* [THM-2237, truncated Boolean moment interval](../../01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md):
  atom information is a feasible fibre, and a missing observable can leave
  both a zero-target and a positive-target law in that fibre.
* The proved dual-feasibility part of
  [THM-534, sector moment LP](../../01-canon/theorems/THM-534-lrc-sector-moment-lp-dual-certificate.md):
  an explicit pointwise polynomial inequality certifies the direction of
  an expectation bound. Its separate LRC extremality claim is not imported.
* [THM-2842, positive-cone multiplier observability](../../01-canon/theorems/THM-2842-ordered-positive-cone-vandermonde-multiplier-observability.md):
  the extra readout must be supplied; choosing a multiplier does not measure it.
* The [signed polynomial atom dual](collatz_atom_polynomial_dual_20261005.md)
  and [localized kernel](collatz_localized_resolvent_floor_20261005.md):
  ordinary Hankel positivity gives no target lower bound, whereas a signed
  minorant can. Their independent-oracle contract is retained.

The anchor is a positive floor for a named source. The niche is exact
finite-dimensional feasibility; the wildcard is a canonical optimizer.
The board is **source / moment provenance / signed localizer / tail budget /
finite dual / unique selection versus existence**. The canonical hostile is
an infinite-support law missing its selected atom. The corrected near miss
is an ordinary positive moment matrix, which by itself supplies upper atom
bounds. The least-used sidecar is a rational separating vector with its
signed measurement-error budget.

Use the inherited injection probability

    p_j=lambda(6j+3),   sum p_j=1,
    h_m(j)=4t/(1+t)^2,  t=2^(m-j),
    H_k=sum_j p_j h_m(j)^k.

Here m is fixed and retained. The target h_m(m)=1 is unique. Every other
value is at most c=8/9; the two nearest indices attain c when present.
These H_k are source-dependent kernel measurements, not the Pascal price
moments and not arbitrary finite-residue statistics.

| Source | Target and map | Preserved predicate | Loss and required sidecar |
|---|---|---|---|
| Actual injection law with selected m | H_0,...,H_(2d+1) | Every stated linear expectation | Source labels outside the kernel fibre; all observations must concern the same actual law |
| Moment packet | Signed matrix C_d | Expectation of (c-h)P(h)^2 | Finite order does not determine the law or target support |
| Rational vector a | Explicit polynomial minorant | Certified lower atom bound | Signed error and normalization P(1)=1 |
| Infinite law plus tail cap | Finite outer polytope | Every actual law remains feasible | Tail correlations are relaxed; feasibility can admit ghost laws |
| Canonical optimizer | One reproducible witness | The optimum of this finite problem | Uniqueness does not force that optimum to be negative |

## 2. The signed matrix and its exact lower bound

For d>=0 define the (d+1)-square rational matrix when the moments are rational,

    C_d[i,j]=c H_(i+j)-H_(i+j+1),   0<=i,j<=d.       (1)

For a polynomial P(h)=sum_(i=0)^d a_i h^i with P(1)=sum a_i=1,

    a^T C_d a = E[(c-h)P(h)^2]
              = -(1-c)p_m
                + sum_(j!=m) p_j(c-h_m(j))P(h_m(j))^2.

The second term is nonnegative. Equivalently,

    q_P(h)=(h-c)P(h)^2/(1-c) <= 1_(h=1)
              on {1} union [0,c].

Thus the finite certificate is

    p_m >= -a^T C_d a/(1-c).                         (2)

A rational vector with negative quadratic value supplies a positive floor.
The whole infinite off-target support is covered by the sign proof; no
finite-node extrapolation is used. If p_m=0 then C_d is PSD. A negative
quadratic direction therefore excludes the zero-target hypothesis. Under a
truthful packet such a direction has P(1)!=0, so it can be normalized to1.

There is at most one negative eigenvalue, counted with multiplicity. Indeed
C_d is the PSD off-target Gram matrix minus
(1-c)p_m times the rank-one matrix11^T. Alternatively, a two-dimensional
negative subspace would intersect the hyperplane sum a_i=0, where the form
is nonnegative. When a negative eigenvalue exists its eigenspace is thus
one-dimensional. This is uniqueness of a negative spectral mode, not a
proof that such a mode exists at the chosen degree.

Completeness holds only with increasing degree and adequate measurements.
Take P(h)=h^d. Then (2) is the previous monomial selector
9H_(2d+1)-8H_(2d). Its negative off-target error tends to zero, so

    p_m>0 iff some degree has a negative normalized localizer value.       (3)

For fixed d, failure to find a negative value does not imply p_m=0.
If truthful intervals of arbitrary precision are available, any strictly
negative exact rational-vector value eventually has a certified sign.
The existence of those measurements is a separate input obligation.

## 3. Finite convex optimization and a canonical witness

Although C_d can be indefinite, its restriction to sum a_i=0 is PSD:
P(1)=0 removes the negative target term. Write

    a=e_0 + D z,   columns(D)=e_i-e_0,  1<=i<=d,
    G=D^T C_d D,   b=D^T C_d e_0.

Then the normalized objective is

    C_d[0,0]+2 b^T z+z^T Gz,   G>=0.                 (4)

For every truthful packet, b lies in the range of G. To see this, regard
G as the Gram matrix on the off-target positive measure
(c-h) times the original law. If z lies in its kernel, the associated
polynomial vanishes in that weighted L2 space. Its inner product with the
constant polynomial is zero, so b^T z=0. A symmetric matrix has range equal
to the orthogonal complement of its kernel.

Consequently Gz=-b has a solution and every solution minimizes (4). For
rational moments, rational elimination supplies a rational minimizer. The
program selects one canonically by RREF with free variables zero. It is
not claiming that the mathematical minimizer is unique. If the tangent
matrix is not PSD or the equations are inconsistent, the packet violates a
necessary condition for a truthful law. Passing these tests is not a full
realizability check.

This is a decidable finite task for an exact rational packet: solve a
finite system and check the sign of its rational minimum. Choosing the
canonical solution has no effect on existence. For example, the law placing
half its mass at c and half at16/25 has a canonical optimum0 at degree1
and all larger tested degrees. It has no target mass. The same deterministic
selection rule therefore returns either a useful negative certificate or
an unhelpful nonnegative optimum.

**Actual-law uniqueness, without assuming the selected atom is positive.**
The actual injection law has a separately established infinite positive
subfamily. In [inverse predecessor sections, section4](inverse_predecessor_sections_20261005.md),
take k=9t+1, n_k=S^k(1), and z_t=(2n_k-1)/3. Each z_t is a leaf with the
explicit route z_t -> n_k ->1 and positive mass
4/[(k+1)(k+2)(k+3)]. Their indices j_t=(z_t-3)/6 are distinct and unbounded.
For every fixed selected m, infinitely many lie beyond m+1, giving infinitely
many distinct kernel values h_m(j_t) in(0,c) with positive off-target weight.
A nonzero polynomial cannot vanish at all of them. Thus the off-target
weighted Gram form is positive definite at every finite degree, and the
tangent matrix G is positive definite for d>=1. The normalized optimizer
is therefore unique for the actual law, even if its selected p_m is zero;
degree0 already has the sole vector(1).

These known positive atoms justify strictness of the objective, not a
positive floor at an omitted target. The actual unique optimizer can have
real, nonrational coefficients because its moments need not be rational.
The exact RREF implementation concerns rational packets; signed rational
witnesses with interval verification remain the finite proof interface.
If the actual selected atom were zero, its finite-degree minimum would be
strictly positive and every such optimized atom lower readout strictly
negative. Uniqueness alone therefore still leaves the sign obligation open.

An enumeration that asks increasing degrees and increasing certified
precision can similarly select the *first* positive floor uniquely under a
fixed enumeration rule. Equation(3) only proves its termination when the
target is positive and the moment access is available. It does not create
existence from uniqueness, nor add a Collatz coverage theorem.

## 4. Signed interval errors and a strict same-order improvement

Expand q_P(h)=sum q_k h^k. Independently justified intervals H_k in[l_k,u_k]
yield

    p_m >= sum_(q_k>=0) q_k l_k + sum_(q_k<0) q_k u_k. (5)

The script preserves exact rational signs. Uniform absolute measurement
error eta costs at most eta sum|q_k|. Its elementary packet validation does
not establish oracle provenance, joint realizability, or a Collatz identity.
An optimizer computed from a guessed midpoint is harmless only if its final
fixed vector is rechecked against truthful intervals using(5).

Here is a concrete improvement at the same maximum measured order3.
Consider the independently specified probability law

    mass1/100 at h=1;
    mass1/2   at h=8/9;
    mass49/100 at h=16/25.

The exact optimizer is P(h)=(25h-16)/9. Its factor vanishes at16/25, and
the signed factor h-c vanishes at8/9, so (2) recovers the exact floor1/100.
In contrast all three earlier monomial tests available from H_0,...,H_3 are

    9H_1-8H_0 = -2719/2500;
    9H_2-8H_1 = -43279/62500;
    9H_3-8H_2 = -686839/1562500.

The gain comes from using measured correlations to select the polynomial,
not from importing a completed Collatz target. This law is a method control,
not a claim about the actual injection weights. The optimized selector has
coefficient norm28577/81, compared with17 for each displayed monomial
readout. It is more selective here, but does not promise cheaper precision;
the interval calculation below pays its larger signed error budget.

## 5. A finite tail polytope and a rational Farkas-style receipt

Choose J>m, retain head masses x_0,...,x_(J-1), and one tail mass t. Put
u=h_m(J); every omitted index has0<h<=u. Require nonnegative variables,
sum x_j+t=1, and an independently proved tail cap t<=delta.
For moment intervals [l_k,u_k], k>=1, retain the finite inequalities

    sum_(j<J) h_m(j)^k x_j <= u_k,
    sum_(j<J) h_m(j)^k x_j + u^k t >= l_k.            (6)

They form an outer relaxation: every actual probability satisfying the
measurements and tail cap is feasible. Tail moments in different rows are
not forced to arise from one common tail law, so feasibility alone is weaker
than realizability. The case delta=1 remains valid when no sharper tail
bound is known. Setting delta=0 without proof would discard the hard part.

Write all inequalities as Ax<=b, with x>=0, and let e be the coordinate
vector selecting the target mass. An explicit rational vector y gives

    y>=0,       e+A^T y>=0
        ==>    e^T x >= -b^T y.                     (7)

Indeed e^T x>=-y^T Ax>=-y^T b. This is the entire dual-certificate proof;
no numerical solver status is used. A positive right side rules out x_m=0
for every feasible law. The exporter checks every rational multiplier and
column slack. An inconsistent packet must not be promoted into information
about the actual measure: truth of the measurement and tail premises remains
essential.

The bounded finite model also makes the zero-target obligation decidable.
Add x_m<=0 and enumerate rational vertices. Every nonempty bounded polytope
has a vertex: unless a feasible point has full-rank active constraints,
move along a nonzero common null direction until another boundary becomes
active, increasing rank; boundedness ensures a boundary is reached.
Thus complete full-rank active-set enumeration decides its emptiness. A
positive LP floor suffices for the actual infinite model; the converse need
not hold because this outer relaxation has forgotten tail correlations.
Within the fixed nonempty bounded relaxation itself, however, zero-target
infeasibility is equivalent to a strictly positive minimum target mass:
the continuous objective attains its minimum on that compact set.

The companion retains the law from section4 with unknown tail mass
1/100000, subtracting that amount from the16/25 head mass. The tail may be
any distribution on indices>=3, where h<=32/81. It derives moment intervals
from that information rather than selecting a convenient tail. The results are

    signed interval floor =63188909/6643012500 >0;
    checked finite LP floor=555875861/59787112500 >0;
    exact outer-LP minimum=3902080991/391937737500 >0.

The feasible bounded polytope has ten vertices; zero-target vertex
enumeration is empty. The output retains the exact dual multipliers. An
explicit feasible point demonstrates that the positive certificate is not
an artefact of contradictory measurement assumptions. These are synthetic
measurement controls with0 ROOT inputs.

## 6. Relaxation hostiles and the Collatz interface

Finite PSD is not target absence. Section4's positive-target law already
has a nonnegative degree0 localizer. A separate ghost illustrates loss of
discrete support: the point mass at h=7/10 gives PSD localizers at every
order because7/10<c. It also gives valid ordinary moment matrices. But7/10
is not any kernel value: it lies strictly between16/25 and8/9. The polynomial
(h-16/25)(h-8/9) is nonnegative on every actual kernel value and negative at
7/10. Keeping only the relaxed interval[0,c] admits that ghost. No positive
atom conclusion is invalidated; negative-certificate sufficiency is intact.

For the actual law, a truthful positive value epsilon from(2), (5), or(7)
at index m gives W(6m+3)>=epsilon. The inherited
[floor deadline](collatz_floor_transport_deadlines_20261005.md) then bounds
the actual odd ROOT time by

    B=(isqrt(1+4 floor(2/epsilon))-3)//2.

Literal replay or the [finite weight bank](collatz_weight_threshold_receipts_20261005.md)
extracts and checks the strict ROOT word. The present experiment intentionally
does not recycle such a word as its measurement input. What remains is a
proof of independently truthful source-indexed moment or tail constraints
that excludes the zero-target feasible set for a new desired family.

## 7. Reproduction and scope

    python -B 04-computation/experiments/collatz_moment_localizer_feasibility_20261005.py
    python -B -O 04-computation/experiments/collatz_moment_localizer_feasibility_20261005.py

The1,476 exact checks use four independently specified finite laws,
degrees0..4, support probes0..40, rational competitor coefficients in{-1,0,1},
an unknown-tail model tested at every single tail location3..44, an independent
full vertex enumeration, the ghost law at degrees0..4, and six malformed
input/false-dual controls. Infinite-support validity and finite optimization
are proved above. No assertion depends on Python assert, and no import runs
a Collatz computation. The program returns conditional mathematical receipts,
not verification of the independent measurement premises.
