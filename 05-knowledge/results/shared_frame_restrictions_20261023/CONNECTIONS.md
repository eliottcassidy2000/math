# Connections that retain their missing coordinate

23 October 2026. Elementary deductions are marked **PROVED scoped**;
recovered results are attributed with their actual scope. This is a map of
usable interfaces, not a claim that shared numbers identify different objects.

## 1. Incoming work changes the portfolio

While this session was computing restrictions, `origin/main` advanced from
8b71b746f7a9 to **1eb9c41c5192**, the concurrent Session43 packet. We read its
proofs and both audit reports. Its
[divisor-height construction and least-realization criterion](../session43_observers_phase_guards_20261023/cycle-geometry/PROOF.md)
independently obtain the same literal cycle and census we had just computed,
and go further by solving arbitrary finite height prescriptions. Credit this
as concurrent/inherited rather than claiming our cycle as a new discovery.
Our additional control computes the complete N/2 census as well.

Its [finite-frame Gaussian results](../session43_observers_phase_guards_20261023/gaussian-jacobian/THEOREMS.md)
prove that finite Reynolds averaging cannot imitate infinite Gaussian
power-nullness in a reduced Q-algebra: the first |G| vanishing orbit power
sums already force the element to vanish, by Newton identities. Its
[finite-degree theorem](../session43_observers_phase_guards_20261023/FINITE_GUARDS_AND_DEGREE.md)
proves that exact C_L equivariance of a planar polynomial map of degree D
lifts to full weight-(1,-1) torus equivariance when L>D+1. A constant nonzero
Jacobian then forces diagonal linearity. Exact C8 and C9 actions combine to
C72, giving the stated conclusion through degree 70; the degree-71 shear
(x+y^71,y) is the sharp nonlinear boundary. These are useful authenticated
bandwidth rules, not arbitrary mod-8 tournament guards.

The incoming monodromy construction also realizes the optional (5)(2)
permutation completion by a double-star polynomial, with a separate joining
branch generator. We retain its distinction between a sphere as underlying
surface and a hyperbolic marked orbifold. There is no need to reproduce a
weaker branch-cover analogy in this packet.

## 2. Gaussian order moments and planar brackets: contract after pairing

The signed order moment in RESTRICTIONS.md has a Gaussian realization:
take independent continuous identically distributed variables Z_i, e.g.
independent real standard Gaussians. Their rank order is uniform, so

    mu(F)=E product_(i<j, ij in F) sign(Z_j-Z_i).          (1)

This follows by partition into the n! order chambers, whose probabilities
are equal by exchangeability and whose boundaries have probability zero.
For a subgroup fixing a port, condition its label to occupy that rank;
the remaining uniform rank assignments give the hypergeometric sign sums
in RESTRICTIONS.md. Thus the fixed-port Krawtchouk coefficients describe a
real conditional order experiment. The signs are not polynomials, however.

The actual Gaussian Moments Conjecture concerns complex polynomials P with
E(P^m)=0 for every positive m, and eventual vanishing of E(QP^m) for each Q.
Derksen–van den Essen–Zhao give an all-dimensional implication to the
Jacobian conjecture; this does not identify (1) with that hypothesis or
give a fixed-dimensional shortcut. See their
[definition and reduction](https://arxiv.org/html/1506.05192v1).
We do not certify the repo's inherited GMC2 formalization claims in this
session; the incoming provenance explicitly separates those from its audit.

Here is a precise small bridge to planar Jacobian work. Let a finite abelian
group act by Poisson automorphisms on C[x,p], and decompose f=sum_chi f_chi,
g=sum_psi g_psi into character sectors. For the Reynolds projection R,

    R{f,g}=sum_chi {f_chi,g_(chi^-1)}.                    (2)

**Proof.** Equivariance and bilinearity imply that {f_chi,g_psi} has
character chi*psi. Averaging kills exactly the nontrivial characters.
Equivalently, R{f,g}={Rf,Rg}+R{f-Rf,g-Rg}; the two mixed terms average to
zero. QED. This is an application of the inherited weight-additive bracket
in [THM-1345](../../../01-canon/theorems/THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow.md),
not a new general Poisson theory. The related
[THM-2375 charge tomography](../../../01-canon/theorems/THM-2375-gaussian-angular-complete-line-charge-tomography.md)
already retains labelled sectors, and explicitly distinguishes Hermitian
energies from scalar moments and missing relative phases.

With {f,g}=f_x g_p-f_p g_x and (x,p)->(zeta*x,zeta^-1*p), x and p are
opposite nonzero charges. Their individual averages vanish, while

    R{x,p}=1 != {Rx,Rp}=0.                               (3)

Retaining the paired sectors recovers exactly the missing constant. Yet
even the averaged bracket is insufficient to test the Keller condition:
P=x(1+p), Q=p has {P,Q}=1+p and R{P,Q}=1 for every L>=2, although the
actual determinant vanishes at p=-1. The correct target is the entire
bracket, not only its charge-zero component. This is the same logical
boundary as preserving endpoint incidence rather than an averaged degree.

An all-size invariant-carrier control makes the lost branch explicit.
For C_m with m>=2, put u=x^m,v=p^m,w=xp. Every invariant monomial factors
through these three generators and uv=w^m. Directly,

    {u,v}=m^2 w^(m-1),   {w,u}=-m u,   {w,v}=m v.        (4)

Every invariant polynomial has zero differential at the origin: invariant
monomials have no degree-one terms. Therefore two invariant polynomial
coordinates cannot have a nonzero constant Jacobian. Opposite-charge
coordinates such as x,p avoid this obstruction. This excludes one naive
quotient carrier, not an arbitrary polynomial Keller pair.

**Sign repair recovered during this check.** Under the bracket convention
above, {xp,x}=-x and {xp,p}=p. The derivation x*d_x-p*d_p is {-,xp},
not {xp,-}. THM-1345's local notation is repaired accordingly. Its
top-weight cancellation also requires a positive total top weight;
the pair (x,p) with top weights 1,-1 is the minimal zero-weight exception.
Neither correction changes the elementary exact-equivariance rigidity
argument. Other claims in that historical theorem are not re-audited here.

## 3. Branch sheets, polynomial sources, and geometry are different layers

The spiral picture is exact for the logarithm's covering, not for the
complex plane itself. On exp:C->C*, one winding adds 2*pi*i to the lifted
coordinate. The negative-real cut selects a branch; it is not an intrinsic
singular line. This convention and its two boundary values are stated in
[DLMF 4.2](https://dlmf.nist.gov/4.2).

For m>=2 the cotangent lift

    (x,p) -> (x^m, p/(m x^(m-1)))                         (5)

has determinant 1 on x!=0 and has m sheets. Its deck action is exactly
(x,p)->(zeta*x,zeta^-1*p); this ties the finite invariant carrier in (4)
to a literal missing branch label. The denominator and punctured domain
prevent (5) from being a polynomial Keller map on the original plane.
The repo already has stronger source-sensitive exclusions, including the
[two-component invariant-divisor gate](../planar_carrier_factor_rigidity_20261019.md)
and [rational-carrier saturation](../planar_rational_carrier_saturation_20261011.md).
We do not reopen those closed carrier routes. A useful new candidate must
explain how changing carriers preserve source regularity and the complete
bracket, including the terms suppressed by (2).

The recovered [5/6/7 geometry work](../geometry_567_collatz_observers_20261004.md)
already supplies exact maps: for triangle type (2,3,q), the sign of
1/2+1/3+1/q-1 gives spherical/Euclidean/hyperbolic geometry at q=5/6/7.
The same Heawood graph has genus-one or genus-three embeddings depending
on its retained local rotation signs. A graph, a rotation system, an
orbifold and an affine dynamical multiplier are distinct carriers. Their
common lesson is to preserve the datum that makes the requested operation
well-defined, not to identify curvature from a shared integer label.

## 4. A genuine 13/4 family: projective incidence as another exact coupling

The least speculative 13/4 connection is the projective plane over F3:
13 points, 13 lines, four points on each line. The family is

    v=q^2+q+1,   k=q+1,                                  (6)

so q=2,3,4 gives (7,3),(13,4),(21,5). This is already represented in the
[Session35 Singer/line-flip work](../session35_dimensions_dualities_20261023/singer-audit/AUDIT.md).
Its prior odd-q theorem distinguishes only two regular states in a given
line-flip family; it does not classify all regular tournaments. Its
full projective-group action cannot orient a tournament invariantly,
since a projectivity can exchange any given pair. A Singer odd cyclic
subgroup is a different, admissible symmetry.

For the point-line incidence matrix M of any finite projective plane of
order q, the axioms give MM^t=q I+J: a point lies on q+1 lines, and two
different points share exactly one. The walk point->uniform incident
line->uniform point has

    K=MM^t/(q+1)^2,
    eigenvalues 1 and q/(q+1)^2,
    spectral gap (q^2+q+1)/(q+1)^2.                       (7)

This is an exact incidence-coupling parallel to the graph/tournament
operator. The spectra come from different Gram identities: at q=3 the
gap is 13/16, not 8/9. The source, target, incidence and uniform weights
are all specified; no isomorphism with tournament classes is asserted.
Our independent coordinate construction checks (7) at q=2,3,5, using
normalized nonzero triples and orthogonality over F_q. The general proof
uses the plane axioms, not those three tests.

## 5. The valuation cycle is constructive; its normalization still matters

For N=2^5*3^7*5^3*7^4*11*13^2, define D_N(a)=max{k:a^k divides N}
for a>=2. This is a composite-base divisibility height, not a p-adic
valuation. For a=product p^r_p,

    D_N(a)=min_(r_p>0) floor(v_p(N)/r_p),
    #{a>=2:D_N(a)>=k}=product_p(floor(v_p(N)/k)+1)-1.      (8)

If a includes a prime not dividing N, its height is zero. The values 0
and 1 are terminal; `check.py` adds bookkeeping self-loops by explicit
convention, not by applying the maximum formula at 1.

Equation (8) makes 4->2->5->3->7->4 a literal five-cycle and 13->2 an
entrance. Our complete divisor calculation agrees with the incoming
stronger proof: 5,616 nonunit divisors enter terminal 1 and 143 enter
the nontrivial cycle. For N/2, the exact height polynomial is

    4656 z + 120 z^2 + 16 z^3 + 6 z^4 + z^7.             (9)

The nontrivial cycle becomes 2<->4, but its divisor basin still has size
143. The terminal-1 basin falls to 4,656. Changing from the Baby Monster
order to the actual involution-centralizer order 2.B is therefore visible
in the dynamics. The group-order origin and centralizer correction remain
inherited from the prior Monster packet; no group-theoretic new theorem
is inferred from this arithmetic census.

**A simple infinite-family extension.** Multiply N by a squarefree M coprime
to N. For every k>=2, the product in (8) is unchanged because each new
prime has exponent one. The actual bases of height at least two are
unchanged as well: none can use a new prime. No new prime is <=7 because
N already contains 2,3,5,7. The recurrent core and its 143-base basin
therefore remain identical. Only the height-one basin grows. The same
argument applies to N/2. This is stability under a specified operation,
not a claim that all multiples preserve the dynamics.

The numerical equality 5616=|PSL_3(F3)| is tempting alongside the 13-point
plane, but the counting mechanisms differ. The group order follows from
[(27-1)(27-3)(27-9)]/2, whereas 5616 here is 5760-144 from exponent boxes.
No canonical bijection intertwining their operations is supplied. We keep
this as a **candidate question**, not a transfer: does a natural group
action exist on those height-one divisors that respects any arithmetic
observable? The cheapest hostile is to demand preservation of height
alone, which is vacuous on this fibre; a meaningful test must also retain
divisibility or the prime-exponent coordinates.

## 6. What this changes for the long-term targets

For Collatz, the reusable target is a refinement rule that retains the
actual source, terminal, integer valuation clock and conditional weights.
An auxiliary uniformly contracting chain does not prove that every source
enters its carrier. The fixed-port modes give an explicit warning against
transporting a spectral bound through arbitrary conditioning.

For native wiring, cuts give an all-size positive spectral theorem, while
kernel-conditioned guard histograms give an exact finite consumer. The
remaining comparison must retain the SAME path ports on every transition.
For planar Jacobian work, paired character sectors and polynomial source
regularity are the corresponding retained information. These are precise
reusable research moves; no inter-problem reduction has been proved.
