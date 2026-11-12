# Compatible attachment fibers: five observables, dyadic rigidity, and symmetry release

**Status: PROVED (scoped algebraic statements), FINITE-EXACT (small-order
minimality controls).** This is a bridge through an explicit *compatibility
fiber*, not an identification of tournament attachments with plane stress.
The inherited quadratic attachment law and seven-coordinate fusion carrier
are not claimed as new. The new organization separates affine dimension,
real compatible directions, binary realizability, and scalar freedom.
Chan–Pak remains OPEN.

## Inheritance and connection contract

The fixed-root carrier is [Session264, F1](../../session264_constructive_20261111/fusion/PROOF.md);
its complete observation quotient and limits are
[Session281](../../session281_chanpak_20261112/wildcard/PROOF.md).
The actual41/43 central switch is inherited from
[Session282](../../session282_chanpak_20261112/central/PROOF.md).
The general quadratic attachment mechanism predates these sessions
(THM-4097/4114, routed by THM-4942). We rederive its coefficients only to
place it in the exact same coordinate system as the fusion carrier.

A distant warning is the [rooted-five Gold obstruction](../../rooted_five_core_gold_obstruction_20261019.md):
complete probabilities on a five-point core did not make core-only moves
sufficient; exterior operations changed feasibility. Another is
[the rooted-five motif/marked-phase study](../../metagraph_five_port_phase_frontier_20261023.md):
small observers may be complete on one universe and blind on another.
Here the analogous lost datum is **binary feasibility inside a real fiber**.

Connection contract: fixed geometry with compatibility equations -> fixed
deletion tournament with four preserved continuation responses. The map is
explicit below. Preserved predicate is equality of all noncentral fusion
coordinates. Lost information is whether a real root-orientation vector is a
binary vector, hence an actual tournament. Restoring that sidecar prevents a
continuous family from being mistaken for a constructive Chan–Pak operation.
No assertion that the mechanics theorem has this same matrix, dimension or
parameter interpretation is made.

## 1. The five-coordinate attachment surface

Fix a nonempty tournament D on n labelled vertices. A root z is attached by
bits x_i=[z->i]. Write h=H(D), H(empty)=1, and

    Q = sum_(S subset D) H(S) H(D-S).

Let s_i(S),t_i(S) count Hamilton paths of D[S] starting/ending at i, zero if
i is absent. Define the following fixed coefficients:

    u_i = sum_nonempty S s_i(S) H(D-S),
    v_i = sum_nonempty S t_i(S) H(D-S),
    a_i = sum_nonempty S (-1)^|S| s_i(S) H(D-S),
    b_i = sum_nonempty S (-1)^|S| t_i(S) H(D-S),
    B_ij = sum_(nonempty proper S) t_i(S) s_j(D-S),
    W_ij=B_ij+B_ji, i!=j,       W_ii=0.

All four rooted response coordinates are affine:

    U=h+u.x,       V=h+v.(1-x),
    d=h+a.x,       e=h+b.(1-x).                     (1)

Indeed remove the initial/final root from the weighted partial path defining
the coordinate. The only root arc tested is its first/last arc; the remainder
is a path of D and its complement contributes H(D-S).

Cut a Hamilton path at z. A root at an endpoint contributes s_i(D)x_i or
t_i(D)(1-x_i); an interior root contributes B_ij(1-x_i)x_j. Therefore

    H(x)=h + sum_i (u_i-t_i(D))x_i
              - sum_(i<j) W_ij x_i x_j.             (2)

Here sum_i B_ij=u_j-s_j(D). This proves the displayed simplification without
assuming that deleting an interior root leaves a Hamilton path: the two
remaining paths need not concatenate in D. Formula(2) is the multilinear
extension of the attachment count from the binary cube to R^n.

**Theorem1.** For a fixed D, the entire nontrivial fusion carrier has at most
five variable affine coordinates `(H,U,V,d,e)`; h,Q are fixed. This bound on
the affine hull is attained, already when D is the transitive order4.

For D with rows(14,12,8,0), use cuts0,1,2,3,4,8 in little-endian bit order.
Relative to cut0 their five-coordinate difference matrix is

    [ 8  8 -1  0 1 ]
    [ 4  4 -2  0 0 ]
    [10 12 -3  0 1 ]
    [ 2  2 -4  0 0 ]
    [ 0  1 -8 -1 0 ],       determinant24.          (3)

More precisely, let r be the rank of the four response rows(u,v,a,b).
For EVERY deletion tournament of order n>=2 the attachment affine hull has
exact dimension r+1. The response coordinates span r independent linear
functions on the binary cube. At least one W_ij is positive: split a
Hamilton path of D into two nonempty pieces. Thus(2) has a nonzero quadratic
multilinear coefficient and cannot be an affine combination of those linear
functions on the cube. This proves both inequalities. The n1 boundary has
constant H and dimension1. The values2,3,4,5 are attained respectively by
C3,TT2,TT3,TT4. Every TT_n for n>=4 has dimension5, by(4) and the same
argument.

No smaller deletion order attains dimension5: there are at most3 independent
linear response directions, plus the scalar quadratic coordinate. This is
an affine-span statement; the attachment set is finite and is not R^5.
Thus the first tournament order for this five-dimensional span is5.

A difference y between two real root vectors preserves all four responses
iff it lies in the kernel of the four rows u,v,a,b. The count change is the
exact linear-plus-quadratic expression from(2). For binary vectors this is a
physical central switch iff y is binary-compatible with the original vector.
Under arbitrary two-sided fusion its count difference is multiplied by the
left and right deletion counts, by the inherited fusion law.

## 2. Transitive spines: a local compatibility operator, but no binary motion

Order the transitive D as0->1->...->n-1. For n>=4, direct subset counting gives

    h=1, Q=2^n,
    U=1+sum_i 2^(n-1-i) x_i,
    V=2^n-sum_i 2^i x_i,
    d=1-x_(n-1),          e=x_0.                    (4)

The four coefficient rows are independent. Thus every nonempty real fiber
of (U,V,d,e) is an affine space of dimension n-4. At the midpoint x_i=1/2
it intersects the open cube in a relatively open neighborhood of that
dimension. However **every binary fiber contains at most one point**, since
U-1 is the binary code of x. No count-preserving or count-changing distinct
binary attachment lies in the same continuation fiber.

This is an all-order positive classification and an all-order obstruction:
real compatibility freedom grows without bound while physical binary
freedom stays zero. Five observables do not imply five physical parameters.

The compatible tangent y has endpoints zero and satisfies

    sum_i 2^i y_i = sum_i 2^(-i)y_i =0.

Equivalently its coefficient polynomial Y(t)=sum_i y_i t^i has the unique
factorization

    Y(t)=t(2-5t+2t^2) Z(t),       deg Z<=n-5.         (5)

The roots of the annihilator are2 and1/2, NOT golden-ratio roots. This local
second-difference filter converts the two global dyadic moments into a
finite-range operator.

## 3. Exact first curvature transition, and an unexpected sixfold factor

For the transitive spine the unordered interaction weights are

    W_ij = 2 if |i-j|=1,
           2^(|i-j|-1) if |i-j|>=2.                 (6)

For an interior split ending at i and starting at j, B_ij=1 when j=i+1;
when i>j, B_ij=2^(i-j-1); all remaining cases vanish. This proves(6).

Put y=t(2,-5,2)*z as in(5). The coefficient of epsilon^2 in
H(x+epsilon y) is

    q(z)=z^T K_m z,       m=n-4,
    (K_m)_(ij)=32,-25,11,-2 for |i-j|=0,1,2,3,
                  0 otherwise.                    (7)

To verify for every m, convolve the symmetric weights(6) with(2,-5,2)
on both sides. The resulting diagonals are(-64,50,-22,4,0,0,...);
multiplication by−1/2 gives(7). Exponential tails vanish exactly.

**Theorem2.** The restricted quadratic form is positive definite for n=5,6,7,8
and indefinite for EVERY n>=9. It is zero-dimensional at n=4.

Proof: leading determinants of K_1,...,K_4 are32,399,2646,8988, all positive,
so Sylvester proves the first four claims. For m=5, z=(10,25,32,25,10)
has q(z)=−2, whereas q(e_1)=32>0. Zero-padding these two vectors embeds both
signs in every larger m. Thus the first indefinite real compatibility fiber
has deletion order9, i.e. tournament order10. No numeric eigenvalue
approximation is used. ∎

The finite-range Fourier symbol factors exactly:

    32-50 cos(theta)+22 cos(2theta)-4 cos(3theta)
      =2(1-c)(5-4c)(1-2c),   c=cos(theta).           (8)

The new zero at cos(theta)=1/2 is a sixth-root-of-unity resonance. The factor
5−4c comes from the reciprocal dyadic filter, not a pentagonal symmetry.
These are infinite-line symbol zeros, not assertions that finite K_m has
those eigenvectors or zero eigenvalues. They explain why small fibers can
look convex before long fibers admit signed cancellation.

At the symmetric point x=1/2 the linear variation along the displayed
symmetric negative witness direction is zero:
for transitive D, (u-t)-W1/2 is antisymmetric under index reversal and the
negative witness y from the symmetric z above is symmetric. More generally
only the quadratic coefficient is asserted; no universal stationarity claim
is needed for Theorem2. Real epsilon can be chosen small enough to stay in
the unit cube. None of this manufactures a binary tournament.

## 4. Symmetry releases actual compatible moves, then too much symmetry kills them

Let D be vertex-transitive of order n>=3. All four response coefficient
vectors are constant by automorphism invariance. Their sums are Q-h,Q-h,
−h,−h respectively. The alternating sums follow from the inherited
K_D(−1)=0 identity. Thus every k-subset root cut has

    U=h+k(Q-h)/n,       V=h+(n-k)(Q-h)/n,
    d=h(1-k/n),         e=hk/n.                     (9)

Consequently all k-subset attachments agree in EVERY fusion coordinate
except H, with no relaxation of binary feasibility. Also s_i=t_i=h/n and
every row sum of B is r=(Q-2h)/n, giving

    H(C)=h+r k-sum_{{i,j} subset C} W_ij.           (10)

This converts the root-attachment problem into a weighted induced-subset
energy problem *inside an actual full-central physical fiber*. The weight
matrix is intrinsic path-cover data, not the original adjacency matrix.

**Theorem3.** A vertex-transitive D has nonconstant count on its two-element
root cuts iff W has at least two distinct off-diagonal entries. A unit
odd-count switch occurs on two such cuts iff two entries differ by2.
The same-cardinality scalar output is constant for every k when W is constant.

If Aut(D) is transitive on unordered pairs, W is constant and

    H(C)=h + (Q-2h)/(n(n-1)) * k(n-k).              (11)

This applies in particular to Paley tournaments of prime order p=3 mod4:
use the square-affine unordered-pair action proved in
[THM-309](../../../../01-canon/theorems/THM-309-five-cycle-design-paley.md)
.
The same elementary group argument applies over finite fields of odd prime
power order q=3 mod4. No global extremal claim about Paley tournaments is
made. More symmetry can create large compatible fibers but erase ALL scalar
freedom on them.

### Exact illustrations, inherited versus new explanatory content

* C3: h3,Q12,W=2, H(k)=3+k(3-k); no same-cardinality freedom.
* Carousel C5 with rows(6,12,24,17,3): h15,Q120, W pair-orbit weights8,10.
  Its k2 fiber has counts41,43. This is the INHERITED central pair, now
  explained by the split into two unordered-pair energies.
* Paley7: h189,Q2016,W=78, H(k)=189+39k(7-k). Its entire cardinality3
  fiber is count657, despite many distinct attachment choices.
* Carousel7: h175,Q1988,W values68,82,84; k2 outputs559,561,575 and k3
  outputs629,641,643,659. Equal-cardinality gives full-central compatibility,
  not arbitrary interval coverage.
* Paley11: h95095,Q1776500,W=28842; the coefficient in(11) is14421.

These values are independently checked against literal Hamilton-path counts
for the small examples. The general proof is by path decomposition and
symmetry, not an extrapolation from the table.

## 5. What five explains, and what remains to be constructed

Three different first scales must not be conflated:

1. Tournament order5: full affine dimension5 of a fixed-deletion attachment
   carrier first appears, even for a completely rigid transitive deletion.
2. Tournament order6: the first fixed-deletion binary central count switch
   appears (the inherited C5 example); exhaustive cores through order4 give
   the stated lower-order boundary, not a new all-order classification.
3. Tournament order10: transitive real compatibility curvature first becomes
   indefinite, despite physical fibers remaining singleton forever.

The creative positive target is therefore **anisotropic symmetric cores**:
vertex symmetry supplies a large binary-compatible bank, while pair-orbit
splitting supplies count variation. One can enumerate k-subset energies using
association schemes or recursive orbital algebras. That is a concrete search
program; neither multiple orbits nor a coefficient gcd proves a complete
interval. The exact41/43 repeated-gate mesh obstruction from Session282 still
applies, and no universal supply or absolute target selector is asserted.

**OPEN conjecture candidate (deliberately stronger than the evidence):**
there exist vertex-transitive cores D_j and cardinalities k_j whose normalized
energy spectra contain intervals of lengths tending to infinity, while their
central deletion contents admit sufficiently diverse coprime choices for
an absolute-target construction. These are two separate obligations; the
small C5/C7 evidence supports only tiny intervals, not the conjunction or
Chan–Pak.

## 6. Turning the orbital geometry yields a new physical pair

The cyclic order9 core with steps `{1,7,6,5}` has outgoing rows

    (226,452,393,275,39,78,156,312,113).

Attach the root beating cut `{0,2,4}` (mask21) or `{0,1,3}` (mask11).
Literal pointed subset-path counting, independently of the energy compiler,
gives the two complete packets

    (14207,3357,50004,18906,34455,2238,1119),
    (14209,3357,50004,18906,34455,2238,1119).         (12)

Thus these are actual ten-vertex full-central unit switches. The shared
deletion count3357=9*373 is coprime to the carousel7 deletion175=25*7.
Search of available theorem/results Markdown for the exact pair found no
prior pair; this is a scoped provenance search, not a claim that every past
binary output was checked. The C5 pair is explicitly inherited, not renamed.

The simple carousel order9 has NO unit pair at any fixed cardinality,
even though all its nonconstant fibers have difference gcd2. Rotating the
*choice of core step set*, rather than iterating one successful order5
shape, restores a unit switch. This is a concrete example of changing the
background geometry after a fixed compatibility fiber fails.

However coprime deletion types only remove an arithmetic obstruction.
Combining one unit switch of deletion175 and one of3357 supplies normalized
digit offsets `{0,175,3357,3512}` relative to its actual baseline. Their gcd
is1, but these four integers are not an interval. The absolute-carry problem
from Session282 is untouched.

### Controlled hypothesis generation

The experiment `scan_orbits.py` has two distinct finite universes:

* every fixed-cardinality cut of the carousel cores of odd orders5 through17;
* every cyclic orientation on Z_n for n=9,11,13, with every cardinality2
  through(n−1)/2. This includes duplicate isomorphism classes deliberately;
  it is not a census of all vertex-transitive tournaments.

In the second universe there are18 consecutive-count pairs counting those
duplicate presentations at n9, and none at n11 or13. The maximum consecutive
run is2 at n9. In the carousel-only universe the unit pairs at orders5 and7
disappear at every tested order9 through17. These observations motivate:

* **H-A, OPEN:** carousel cores of every odd order n>=9 have no same-cardinality
  unit count switch. This is a precise proposed rigidity family, not a theorem.
* **H-B, OPEN:** other vertex-transitive cores admit same-cardinality central
  unit switches at unbounded orders. The new n9 pair supports continued
  flexibility beyond the two smallest examples, while the complete cyclic
  n11/13 negatives show there is no simple 'every larger order' extrapolation.
* **H-C, OPEN, more ambitious:** some recursive orbital-algebra family admits
  growing consecutive energy intervals and enough deletion-content diversity
  to feed an absolute-carry selector. This has two independent obligations;
  H-B would not imply it.

The desired unification is not 'five causes everything'. It is the
**compatible-fiber/physical-lattice split**: preserve a declared future
interface, identify the exact allowed lattice inside that fiber, then study
the scalar energy on it. Symmetry changes both the lattice and the energy;
maximal symmetry may flatten the latter completely. Real curvature and
five-dimensional affine spans provide useful geometry but no substitute for
that lattice step.
