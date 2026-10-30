# The arithmetic of p^k+1: clocks, lattices, and paired groups

**Status: PROVED elementary deductions; CITED finite-field cyclicity and
the specified Zsigmondy input; FINITE-EXACT implementations.**
Classical norm-one and trace mechanisms are credited to the recovered
sources in RECOVERY.md. No Selmer/Sha realization or finiteness theorem
is asserted.

## A1. Integer factorization keeps the exponent's parity

For an odd prime p and k>=1 put N_k=p^k+1. Elementary factorization gives

    N_k = product_(d|k, k/d odd) Phi_(2d)(p).             (1)

For positive a,b, with g=gcd(a,b),

    gcd(p^a+1,p^b+1) =
       p^g+1,  if a/g and b/g are both odd;
       2,      otherwise.                               (2)

For a complete proof let D be the gcd and choose integers u,v with
ua+vb=g. The residue p is a unit modulo D, so p^g=(-1)^(u+v).
If a/g or b/g is even, raising this identity to that exponent forces
1=-1 modulo D; since D is even, D=2. If both are odd, the equation
u(a/g)+v(b/g)=1 forces u+v odd, giving D | p^g+1. The reverse
divisibility follows from odd-power factorization.
The two-adic valuation is also explicit:

    v_2(p^k+1)=v_2(p+1), k odd;  =1, k even.              (3)

For odd k the quotient by p+1 is an odd alternating sum. For even k,
p^k=1 mod8. These also settle the two-adic part of (2).

For every odd p and k>=2, N_k has a prime ell with

    ord_ell(p)=2k, hence ell=1 mod 2k.                    (4)

For k>=4 this follows from the **CITED** classical primitive-divisor
theorem for p^(2k)-1, as stated in the introduction of
[Flatters--Ward](https://arxiv.org/pdf/1002.4829), with original
[Zsigmondy source](https://doi.org/10.1007/BF01692444).
A primitive divisor cannot divide p^k-1, so lies in N_k.
No new proof of Zsigmondy's theorem is claimed.

The two small cases need no deep input. For k=2, the odd part of p^2+1
is greater than1 and each of its primes has order4. For k=3 use
p^2-p+1. Every prime divisor other than3 has order6. If3 divides this
factor then its valuation is exactly one; for odd p>=3 the factor is
larger than3. Thus a prime of order6 exists. At k=1 the missing case
is p+1 a power of two. Base2 is a genuine hostile: 2^3+1=9 has no
prime of order6, as the older Bang/Zsigmondy packet already established.

The verifier factors all56 integers for p=3,5,7,11,13,17,19 and k=1..8,
checks every prime's order, and checks (2) for a,b<=12. None of these
statements says that N_k itself is prime: for odd p it is even.

## A2. A projective line is a complete norm-one clock

Let q=p^k be odd. In F_(q^2)=F_q(alpha), choose alpha^2=delta with
delta a nonsquare in F_q. Frobenius over F_q sends alpha to -alpha.
The norm-one group

    U={z: z^(q+1)=1}

is cyclic of order q+1. A complete parametrization by P1(F_q) is

    z_t=(t+alpha)/(t-alpha),  t in F_q;  z_infinity=1.     (5)

The denominator never vanishes. The formula is injective by cross
multiplication; with q+1 values it exhausts U. Multiplication becomes

    t star s = (ts+delta)/(t+s),                         (6)

with infinity as identity and a zero denominator giving infinity.
Thus p^k+1 is a literal group size with a closed operation and a retained
coordinate chart, not merely a numerical resemblance to a tournament size.

If z generates U, its degree over F_p is exactly2k. A shorter
Frobenius period d would be a proper divisor of2k, hence d<=k, whereas
p^d-1<p^k+1 cannot be divisible by ord(z). The reciprocal trace
z+z^-1 lies in F_q and has degree exactly k: equality after p^e,
e<k, would force z^(p^e)=z or z^-1, requiring p^k+1 to divide
p^e-1 or p^e+1.

These are different operations:

- Squaring z changes its order from q+1 to (q+1)/2.
- Reciprocal trace identifies z and z^-1, leaving (q+3)/2 values.
- The p-Frobenius orbit has length2k at z and k at its trace.

Trace is not closed under arbitrary future multiplication. For a
named h not equal to +/-1,

    tr(hz)-tr(hz^-1)=(h-h^-1)(z-z^-1),                   (7)

so the two branches are distinguishable unless z=+/-1. This is the
inherited reciprocal-branch mechanism, specialized to the complete
q+1 clock. It does not say that a traced time sequence has half its
period. The finite-field hypothesis is essential; composite quadratic
rings can have more branches and ramified fibers.

## A3. The determinant phase separates switching from literal symmetry

Represent projective points by v_t=(t,1) and v_infinity=(1,0), and put

    K_xy=chi(det(v_x,v_y)).

This is skew-symmetric exactly when q=3 mod4. It then defines a
tournament with infinity a source and a Paley tournament on finite
points (up to the harmless global convention for the finite arcs).
For q=1 mod4 it is symmetric, and is not a tournament relation.

For g in GL2(F_q), write g v_x=lambda_x v_(pi(x)). Then

    K_(pi(x),pi(y)) =
       chi(det g) chi(lambda_x) chi(lambda_y) K_xy.       (8)

This proves the classical PSL2(q) switching action directly. We do
not claim a classification of every symmetry of a Hadamard matrix,
where independent row/column operations and exceptional small cases
are different actions.

The star multiplication by s in (6) is represented by

    M_s=[[s,delta],[1,s]], det M_s=s^2-delta.

The determinant character is the quadratic character in U:

    chi(s^2-delta)=z_s^((q+1)/2).                        (9)

For example, apply Euler's criterion to s+alpha in F_(q^2):
z_s=(s+alpha)^(1-q), so both sides equal the same sign.
Consequently a primitive (q+1)-clock has determinant character -1.
Its square lies in the determinant-square half and gives a clock of
order(q+1)/2. At q=3 mod4 the full clock reverses orientation up to
switching; the half-clock preserves it up to switching.

Neither assertion makes this even-order clock an ordinary tournament
automorphism. Tournament automorphism groups have odd order, and
the earlier regular-parent theorem fixes the adjoined source in
the literal orientation. The signed scales lambda_x and the determinant
phase in (8) are the needed sidecars.

In particular every 5^k field has q=1 mod4. Recursive carousel
tournaments on5^k vertices still exist, but the quadratic determinant
relation on that field is symmetric. The same issue occurs at even
exponents for p=3 or7.

## A4. Two different27-vertex parents and a Frobenius attachment residue

The iterated triangle tournament C3^[3] and the Paley tournament on
F_27 both have27 vertices. Their full automorphism group orders are
3^13=1594323 and1053, respectively; their transitive maxima are8 and5.
The former values follow from the inherited wreath and substitution
laws; the latter are checked by native full-group and exact
transitive-set searches.

The Paley construction has actual automorphisms x->a*x^(3^j)+b,
where a is a nonzero square and j=0,1,2. These give27*13*3=1053
different maps; native enumeration confirms completeness.
The subset {0,1} retains the order-three Frobenius group. Therefore
the unmarked one-vertex attachment beating those two points has
automorphism group C3, by regular-parent R1. With the documented
irreducible polynomial [1,0,2,1], the three-set {0,1,alpha}
(encoded as {0,1,3}) has trivial stabilizer and gives a rigid28-child.

Thus the two-point recipe from prime fields does not remove field
Frobenius automatically. The literal sources and the field basis
are retained in the certificate. For general Paley F_(p^k), k>1,
Frobenius fixes {0,1} pointwise, so that cut cannot be distinguishing.

## L1. The integral root-lattice cokernel is completely determined

Let H be any doubly regular tournament of order q>=3 and S its skew
matrix. On the root lattice

    L={x in Z^q: sum x_i=0},

set theta=(I+S)/2, C=(q+1)/4, r=(q-1)/2.
The integral action and relation

    theta^2-theta+C I=0, or theta(1-theta)=C I,           (10)

are **inherited** from
[THM-4814, section E](../../../01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md).
No freeness over a nonmaximal quadratic order is needed below.

**Cokernel theorem.** In integral bases,

    Smith(theta)=diag(1 repeated r, C repeated r),
    L/theta L is isomorphic to (Z/C)^r.                 (11)

The same holds for 1-theta.

**Proof.** The characteristic polynomial of theta on L tensor Q is
(X^2-X+C)^r, since S^2=-qI on the sum-zero space and its conjugate
eigenvalues have equal multiplicity. Hence det(theta)=C^r.
Equation (10) makes every Smith invariant divide C.

For a prime ell dividing C, reduction gives an idempotent theta and
characteristic polynomial X^r(X-1)^r. Its two roots are distinct even
when ell=2. Thus exactly r Smith invariants are divisible by ell.
If a=v_ell(C), their valuations are at most a and sum to ra, so all
r equal a. In the ordered Smith list these occupy its last r positions
for every ell. Combining all primes gives (11). QED.

The controlling quantity is exactly (q+1)/4, rather than an unspecified
quadratic determinant. Thus every odd prime in q+1 appears with known
exponent and multiplicity in this cokernel.

| q | C | cokernel theta |
|---|---:|---|
| 7 | 2 | (Z/2)^3 |
| 11 | 3 | (Z/3)^5 |
| 15 | 4 | (Z/4)^7 |
| 19 | 5 | (Z/5)^9 |
| 23 | 6 | (Z/6)^11 |
| 27 | 7 | (Z/7)^13 |
| 59 | 15 | (Z/15)^29 |

Full Smith forms are computed exactly, together with independent checks
of integrality, (10), and the Gram-adjoint identity. The q=27 action
is over the nonmaximal order of discriminant-27, already qualified in
the source; this does not alter (11). The q=15 example is the doubled
skew-Hadamard construction in GROWTH.md; it tests a non-field source
and a non-squarefree norm C=4.

## L2. Complementary cokernels carry a perfect pairing

Put X=L/theta L, Y=L/(1-theta)L. For representatives x,y in L define

    lambda([x],[y]) = <(1-theta)x,y>/C mod Z.             (12)

The ordinary Euclidean pairing restricts to L. The adjoint of theta
is1-theta. Replacing x by x+theta a changes (12) by <a,y>; replacing
y by y+(1-theta)b changes it by <x,b>. Thus lambda is well defined.

In the basis e_i-e_q the Gram matrix of L is I+J, with determinant q.
Since gcd(q,C)=1, that Gram matrix is invertible modulo C.
If lambda([x],Y)=0, then (1-theta)x lies in CL. Therefore
x=theta x+(1-theta)x belongs to theta L, because CL is contained in
theta L by (10). This proves left nondegeneracy; the adjoint argument
proves right nondegeneracy.

On X direct-sum Y the form

    omega((x,y),(x',y'))=lambda(x,y')-lambda(x',y)         (13)

is a perfect alternating pairing. Its group order is C^(q-1), a square.
The construction works at2 as well: alternation follows from the
displayed subtraction, not from dividing a skew form by2.

This is an actual paired finite group arising from a tournament
operator. Identifying it with a Selmer group or Tate--Shafarevich group
requires an arithmetic realization and compatible local conditions
that have not been supplied.

For prime C=ell, (13) is a2r-dimensional symplectic F_ell space.
Its number of Lagrangians is

    product_(i=1)^r (ell^i+1).                          (14)

Count pairs consisting of a Lagrangian and a line in it, then quotient
by that line to obtain the recurrence L_r=(ell^r+1)L_(r-1).
The first nontrivial instance here, q=7, has r=3 and135 Lagrangians.
The familiar fifteen in dimension four is not its full-space count.
Choosing an isotropic line and reducing to its orthogonal complement
modulo that line DOES give dimension four and fifteen Lagrangians.
That choice is an additional marked datum, not a canonical identification
with the separate fifteen Ogg primes.

This also explains the three requested families directly. Write L_j(p)
for the number of Lagrangians in a2j-dimensional symplectic F_p space.
Then the classical recurrence is exactly

    L_0(p)=1,  L_j(p)=(p^j+1)L_(j-1)(p).               (15)

The paired tournament lattices at q=11,19,27 have respectively
(p,r)=(3,5),(5,9),(7,13). Their Lagrangian counts therefore contain
the products of3^j+1,5^j+1,7^j+1 through j=5,9,13. Symplectic
reduction by an isotropic line supplies the lower-dimensional spaces.
This is a map between actual objects, not just equal group orders.
Independent RREF enumeration in dimension four gives15,40,156,400
Lagrangians for p=2,3,5,7. These dimension-four counts are not the
full Lagrangian counts of the larger tournament-lattice examples.

The formula here requires a FIELD: a Z/4 lattice quotient is not an
F_4 vector space merely because both rings have four elements.

## L3. New clock primes force new paired torsion

Fix p=3 mod4 prime and let k>=3 be odd, so q=p^k supports a Paley
tournament. A1 gives an odd prime ell with ord_ell(p)=2k. Therefore
ell divides C=(q+1)/4, and divides none of p^j+1 for 1<=j<k.
If a=v_ell(C), L1--L2 give the actual primary groups

    X[ell-infinity] = (Z/ell^a)^((q-1)/2),
    (X direct-sum Y)[ell-infinity] = (Z/ell^a)^(q-1),

with the restricted perfect alternating pairing. This deduction uses
the CITED primitive-divisor input specified in A1, not a new proof of it.

The standard affine--Frobenius maps x -> u*x^(p^j)+v, with u a
nonzero square, form a literal tournament subgroup of order
q(q-1)k/2. The new ell divides none of its factors: ell is not p,
gcd(q-1,q+1)=2, and ell>=2k+1>k. Thus new norm-clock primes become
paired torsion while being absent from this supplied symmetry group.
No general full-automorphism classification is needed for that claim.

| p | k | q | C | a primitive ell |
|---|---:|---:|---:|---:|
| 3 | 3 | 27 | 7 | 7 |
| 3 | 5 | 243 | 61 | 61 |
| 7 | 3 | 343 | 86 | 43 |
| 11 | 3 | 1331 | 333 | 37 |

At27 the full-group calculation in A4 makes the absence a statement
about the entire literal group. At the other displayed larger fields
we assert absence only from the explicit affine--Frobenius subgroup.
The finite factorization bank checks these entries; their lattice
groups follow from L1, not from unperformed large Smith computations.

## What these maps do and do not identify

At q=27, the literal Paley group has order1053 with prime support
{3,13}, the norm clock has order28, and the integral cokernel is
(Z/7)^13. These are three different objects with explicit maps.
A prime can be absent from the literal automorphism group while being
forced in the integral arithmetic carrier. The paired construction
supplies algebraic structure, not a finiteness proof for an unrelated
infinite arithmetic family.
