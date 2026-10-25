# Power maps, retained arithmetic, and the scope of a unification

2026-10-25. **PROVED elementary constructions and obstructions** below;
standard toric geometry and the lambda interpretation are **CITED**.
The proposed cross-problem research program is **OPEN**, not an equivalence
between RH, Collatz, P versus NP, or tournament counting.

## Inheritance and the board

Closest mechanism: the literal product/involution model in
[microcosm.md](../zero_heavy_ports_20261020g/microcosm.md), where a pointed
Boolean set of size q has q^2 paths and q blue diagonal paths. Its primitive
triple is (q,(q^2-1)/2,(q^2+1)/2), with consecutive Euclid parameters.
Canonical hostile: an identical capped weighted core can have blue counts
17 and 0. Corrected near miss: preserving a count or a conic does not
preserve path incidence or its special arithmetic subfamily.

The least-used sidecar is the **integral model at the prime 2**. The recent
[absolute-geometry audit](../session58_hourglass_complement_geometry_20261025/absolute-geometry/ABSOLUTE_GEOMETRY_BRIDGES.md)
distinguishes a genuine all-prime carrier from a counting polynomial, while
the [projective-support bridge](../session58_hourglass_complement_geometry_20261025/projective-support-bridge/PROJECTIVE_SUPPORT_BRIDGE.md)
preserves divisibility but explicitly loses addition. Neither is an RH proof.

Board: exact tournament ports; canonical owners across rank; integer gap
roles; toric power maps; global arithmetic correlations; uniform proof
obligations. The useful common question is whether the chosen map preserves
the predicate after the next operation. Shared vocabulary alone does not.

## 1. A genuine commuting Frobenius carrier for Pythagorean triples

Let C be the projective scheme UV=W^2 over Z. The Veronese map

    [m:n] -> [m^2:n^2:mn]

identifies C with P^1_Z. On U!=0 use W/U as coordinate; the relation gives
V/U=(W/U)^2. The analogous V-chart gives the other affine line, and these
two charts cover C. Their transition is inversion. This proves the stated
identification, including characteristic two.

For every positive integer r define

    psi_r[U:V:W] = [U^r:V^r:W^r].                         (1)

The coordinates never vanish simultaneously at a geometric point, the
equation is preserved, and psi_r psi_s=psi_(rs). Under the parametrization,
this is [m:n] -> [m^r:n^r]. For prime p, reduction modulo p is the absolute
Frobenius. Thus this is an explicit instance of the standard toric
commuting-Frobenius structure, in the lambda framework of Borger, sections
1 and 2.4: [primary paper](https://arxiv.org/html/0906.3146v1).
No novelty is claimed for the toric construction.

Over Z[1/2], the coordinate change

    [U:V:W] -> [x:y:z]=[U-V:2W:U+V]

identifies C with x^2+y^2=z^2. In the latter coordinates (1) is

    Phi_r(x,y,z) = (((z+x)^r-(z-x)^r)/2,
                    y^r,
                   ((z+x)^r+(z-x)^r)/2).                 (2)

The halves in (2) are integer polynomials: in the sum/difference all
remaining binomial coefficients occur twice. Also

    z'^2-x'^2 = ((z+x)(z-x))^r = y^(2r),
    Phi_r Phi_s = Phi_(rs).                              (3)

For odd prime p, (2) reduces to coordinatewise p-th powers modulo p. The
original conic model is different at 2: x^2+y^2-z^2 becomes
(x+y+z)^2, a nonreduced double line. The coordinate change has determinant
divisible by 2. For example Phi_2=(2zx,y^2,z^2+x^2), which is not the
coordinate Frobenius of that bad model. At [1:0:1] modulo 2 its three
coordinates are all zero, so these formulas do not even define a global
projective map on that model. C repairs the model; this is not a proof
that every integral model admits the same lift.

If m>n>0 are coprime of opposite parity, (2) applied to their primitive
Euclid triple is

    2^(r-1) (m^(2r)-n^(2r), 2m^r n^r, m^(2r)+n^(2r)).   (4)

The inner triple is primitive. Indeed the new parameters remain coprime
and of opposite parity. Formula (4) therefore specifies the exact scale
removed by primitive normalization; that normalization is an arithmetic
operation, not an unrecorded projective identification.

## 2. The first hostile: the blue-count subfamily is not stable

The tournament triple uses m=n+1. For every r>1 and n>=1,

    (n+1)^r-n^r > 1.                                    (5)

Consequently toric powering leaves the consecutive-parameter subfamily.
This is stronger than failing to match its powers-of-two subsequence.
Explicitly psi_3 sends the primitive triple (3,4,5) to (63,16,65).
Its actual Phi_3 value before primitive normalization is (252,64,260).
The blue-diagonal recipe with first leg 63 would instead give
(63,1984,1985). Both satisfy a conic equation, but they are different
arithmetic constructions.

Likewise the literal tournament deletion/growth q' = 2q-1 is not the
prime-power map. Starting at (3,4,5), growth gives (5,12,13). There is
no commuting square between those two operations in the present model.

The original coordinatewise-power guess fails even earlier: for a
nondegenerate Pythagorean triple and r>1,

    x^(2r)+y^(2r) < (x^2+y^2)^r = z^(2r).

Repairing that guess required the exact multiplicative coordinates U,V,W.
Repairing the tournament connection would additionally require transport
of the consecutive-parameter condition and of actual path incidences.

### A positive repair: allow the exceptional set to grow

There is an exact finite-set construction behind the missing parameter.
Let B be a subset of A, with |A|=m>|B|=n>0. Form

    X=A disjoint-union B,

and let kappa exchange the two copies of each b in B, fixing A minus B.
Write q=m+n and d=m-n. On X squared, J=kappa times kappa has d^2 fixed
points. Therefore the three actual cardinalities

    |X times Fix(kappa)| = qd,
    nontrivial J-orbits = (q^2-d^2)/2,
    all J-orbits        = (q^2+d^2)/2                     (9)

are the Euclid triple (m^2-n^2,2mn,m^2+n^2). Coprime, opposite-parity
m,n give a primitive triple. This generalizes the one-exception orbit
count without asserting that the first set in (9) consists of blue paths.

The operation (A,B) -> (A^r,B^r), using the actual subset inclusion,
realizes toric powering. The two copies of B^r are exchanged, while
A^r minus B^r is fixed. Hence

    q_r=m^r+n^r,                 d_r=m^r-n^r,
    d_r=d sum_(j=0)^(r-1) m^(r-1-j)n^j.                  (10)

Flattening tuples makes successive powers compose as r*s. Equations
(9)--(10) give an explicit commuting diagram between nested finite sets,
their integer parameters, and the toric conic over Q. For primitive integral
homogeneous coordinates first divide gcd(m,n), or restrict to coprime m,n.
The growing exceptional
set is the exact coordinate lost by insisting on d=1. Starting with
m=2,n=1, the third power has q_3=9 and d_3=7; its orbit counts are
63,16,65. No fixed bound on d survives all powers of a nondegenerate pair.

At d=1 and n=2^(k-2), k>=2, this pointed involution has the same cycle type
as the complement involution of the old Boolean choice set of size
q=2^(k-1)+1, so an isomorphism of involutive sets exists. For other
consecutive m,n the orbit construction still applies, but that specific
Boolean-family identification is not asserted. The actual path interpretation additionally singles out the
swap R on X squared. Its fixed set still has size q, whereas the first
set in (9) has size qd. For d>1 these are different. Thus the construction
repairs arithmetic and composition, while marking exactly where the
blue-path predicate stops being transported. It is a usable generalized
counting carrier, not a claimed enlargement of the tournament theorem.

## 3. Why an honest toric model still does not prove RH

For a fixed finite field F_Q, C has Q+1 points. Hence

    Z(C/F_Q,t) = 1/((1-t)(1-Qt)).                         (6)

Its good integral model has, in Re(s)>2, the global finite-prime product

    product_p Z(C/F_p,p^(-s)) = zeta(s) zeta(s-1).         (7)

This identity follows directly from the two Euler factors. More generally,
if an actual all-prime model has #X(F_(p^m))=sum_j a_j p^(jm) for one
fixed finite integer polynomial (finitely many j>=0), then

    product_p Z(X/F_p,p^(-s)) = product_j zeta(s-j)^a_j    (8)

in their common absolute-convergence half-plane. Integer coefficients of
either sign are allowed in this algebraic identity when the stated point
counts actually hold. Equation (8) does not construct the model.

These formulas place the same unresolved Riemann zeta factors in another
geometry. The elementary eigenvalues 1,Q of (6) do not locate the
nontrivial zeros of the global factors in (7). A geometric carrier and
commuting Frobenius maps are substantial data, but the requisite global
positivity has not followed from them.

The same circularity test applies to A000568: taking A000568(n) disjoint
copies of Spec Z gives zeta(s)^A000568(n). This is a valid realization
once the copies are specified, and provides no method to compute the
unknown exponent with less work.

## 4. A precise RH-facing target: support and normalization together

**CITED:** Weil's criterion reformulates RH as a sign condition on the
explicit-formula functional for autocorrelation test functions satisfying
two moment constraints. In the symmetric multiplicative-Haar normalization
these are integral h(x) x^(1/2) d*x = integral h(x) x^(-1/2) d*x = 0.
Connes--Consani discuss this in section 4 of
[The Scaling Hamiltonian](https://alainconnes.org/wp-content/uploads/scalingH.pdf).
Their section 3 also shows why a tempting local positivity argument fails.
We import no unproved positivity statement from that paper.

An elementary support calculation sharpens the proposed interface. Let h
have support in [exp(-L),exp(L)] and use multiplicative Haar measure. Its
autocorrelation has support in [exp(-2L),exp(2L)]: any nonzero term requires
both y and xy to lie in the support. Thus prime-power evaluations in the
finite-prime explicit formula can involve only p^j<=exp(2L). The global
normalization, archimedean term, and two moment conditions remain required.

**OPEN construction target:** for every L, construct a positive
representation of the correctly normalized Weil form on all admissible
test functions with this support, compatible under enlargement of L.
If it proves the full cited inequality on those spaces, every compactly
supported admissible test function is covered, and the cited criterion
would imply RH. The unsolved step is exactly that representation/inequality.
No arbitrary tournament Gram matrix is identified with this form.

The useful inspiration from owner languages is retaining boundary data
through enlargement. Its transfer here would have to preserve the actual
form, moments and normalization. This is a proposed proof architecture,
not a reduction from tournament connectivity to RH.

## 5. Classification by proof duty, not asserted equivalence

| Problem | Exact object that must survive | Remaining universal duty |
|---|---|---|
|Unrestricted A000568 recurrence|Canonical class owner and attachment incidences|Control the growing rank-transition state and its exact quotient|
|Black-core connectivity|Actual directed paths, positional blackness, labelled target|Cover every surviving endpoint/port configuration|
|Collatz|One fixed positive integer, exact affine carries and chronology|Certify eventual entry/descent for every starting integer|
|RH through Weil's criterion|The normalized global quadratic form and moment constraints|Prove the sign for every admissible compactly supported test function|
|P versus NP|An input language, verifier and input-length resource bound|One deterministic polynomial-time algorithm for every NP language, or a separation|

This is a classification of bottlenecks. It is not evidence of mutual
logical determination. In particular, NP-completeness is an explicit
polynomial-time many-one equivalence between languages: a map must preserve
yes/no instances in both directions. See Cook's
[official problem description](https://www.claymath.org/wp-content/uploads/2022/06/pvsnp.pdf),
section 2. A000568 is a counting function, and RH and Collatz are universal
mathematical assertions; none has been supplied with such a reduction here.

P=NP would make proof search with a unary-encoded length cap efficient in
the formula length and that cap, for a polynomial-time proof verifier.
A binary-encoded cap does not give the needed polynomial witness bound.
P=NP supplies neither a bound on the length of
a proof of RH nor an assertion that RH is true. Likewise a finite-state
machine for every fixed parameter need not give one machine of polynomial
size when the parameter becomes part of the input.

To claim a stronger unification, provide: source instances; target
instances; a computable map; iff preservation of the target predicate;
resource/parameter distortion; and a way to lift a witness or obstruction.
The present positive maps satisfy this contract only for their stated
counting, power-map, support, or ownership predicates.

## 6. Exact checks and stopping boundary

`f1_check.py` verifies the toric chart on projective points over the primes
2,3,5,7,11,13; all small primitive Euclid parameters m<=20; composition,
exact scale and the first hostile; explicit nested-set powers and their
involution orbits for 1<=n<m<=4 and r<=3; odd-prime coefficient congruences; and
local zeta coefficients by an independent fixed-count recursion.
All gates are explicit exceptions and survive Python -O. General proofs
are above; no finite experiment is an RH test or a global descent proof.

The proved positive bridge is an actual all-prime toric carrier. The
proved negative boundary is failure to preserve the blue-count subfamily.
The proposed RH route preserves a named target form but currently has no
positivity certificate. This distinction is the research outcome.
