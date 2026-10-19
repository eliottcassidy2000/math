# A precise common mechanism: first lifting obstructions and integral normalization

2026-10-19. Status: elementary proved linear-algebra synthesis of this session's
results. The lifting lemma is standard in substance; no claim of novelty for
it or for Bockstein methods is made. The applications identify actual maps
and lost data, rather than equating the underlying conjectures.

## 1. The first obstruction for an integer matrix

Let M:Z^a->Z^b be an integer matrix and let p be prime. Write Mbar for its
reduction modulo p. For v in ker Mbar choose an integer lift vtilde and put

    beta(v) = [M vtilde / p] in F_p^b / im Mbar.

This is a well-defined F_p-linear map. Changing vtilde by p w changes the
bracket representative by Mbar w, and the same observation proves linearity.
Furthermore,

    beta(v)=0 iff v has a lift vtilde+p w with M(vtilde+p w)=0 mod p^2.

Indeed the equation is exactly Mbar w=−M vtilde/p modulo p. Thus an augmented
rank computation is a complete first-lift test, not merely a heuristic.

In Smith coordinates, a nonzero diagonal entry p^e u with p not dividing u
contributes a kernel direction modulo p exactly when e>=1. Its beta image
is nonzero exactly when e=1. A zero Smith entry contributes a direction
with zero obstruction. Therefore beta distinguishes torsion killed at the
first p-adic level from classes surviving to higher precision. It is
stronger than either mod-p nullity or a characteristic-zero rank alone.

The defect Mv/p must be retained before reducing. Reducing Mv to zero first
and discarding the quotient destroys the entire obstruction.

## 2. Ternary norms: two different reconstruction losses

In the period-paper extension, M is a literal ternary norm matrix restricted
to cyclic homogeneous inputs. The independent cyclicity restriction is
essential; the order-three norm alone has an averaging projector when p!=3.

For odd degree n and two-variable-supported inputs, the session proves

    ker(M mod p^m) = p^(m-1) ker(M mod p), p>=5, m>=2,
    dim ker(M mod p) = floor(log_p n).

Consequently beta is injective on this kernel: every supported odd modular
class is obstructed at the first lift. The corresponding integral p-primary
torsion consists of floor(log_p n) copies of Z/p, not a higher cyclic group.
This all-degree statement concerns support-preserving corrections.

The full cyclic carrier is different. A modular class may have genuinely
three-variable support yet still be visible on a coordinate plane. Another
class may be invisible on every coordinate plane. These are distinct losses:

- failure to reconstruct from boundary values is a kernel of restriction;
- failure to lift a modular relation is a beta obstruction in a norm cokernel.

A class may exhibit both, but the predicates are not interchangeable.
See the period packet for the separately proved universal families and
finite-exact exceptional prime computations.

## 3. Tournament deletion: carry count is the entire obstruction spectrum

The integral source-free basis gives a unimodular diagonal presentation of
D_m:C_n->C_(n+m), with entries

    a_r=binom(n-r+m,m), multiplicity t_r−t_(r−1), 0<=r<=n,
    with t_(-1)=0 and t_r the number of tournament isomorphism classes.

Therefore the first obstruction beta is nonzero precisely on primitive
string coordinates with v_p(a_r)=1. Coordinates with v_p(a_r)>=2 have a
first lift. More precisely, a nonzero mod-p coordinate of valuation e has
a kernel lift modulo p^k exactly for k<=e, and fails at p^(e+1). For a
mixed mod-p class, take the minimum valuation among its nonzero coordinates. By Kummer's
binomial valuation formula, these cases mean exactly one carry versus two
or more carries in adding n-r and m in base p. Thus carry information is
not a visual resemblance to Pascal triangles: it classifies the actual
integral matrix's lifting obstructions and cokernel torsion.

Ordinary repeated deletion differs:

    d^m=m! D_m.

Its valuations acquire v_p(m!), an additional normalization contribution.
The adjoint statement in the divided-power/polynomial pairing is

    ordinary derivative of order m = m! times Hasse derivative of order m.

This is the same integral factor that changes the Wronskian paper's jet
valuation budget. It cannot be canceled after reduction modulo a prime.
The underlying tournament, arithmetic sampling nodes, and prime-sum target
are still different objects; only this module operation and lattice factor
are identified.

## 4. Scalar Frobenius compatibility: a nonlinear companion, not the same matrix

For A in Z[X,X^-1], let E(A)=A(X)^p−A(X^p). At first precision, an alteration
A->A+pB changes the divided defect by

    E(A+pB)/p = E(A)/p−B(X^p) mod p.

Hence its transverse class in the additive quotient by p-grid polynomials
is intrinsic to f=A mod p. This parallels the matrix beta calculation, but
E is nonlinear and the correction image is spatial substitution, not a
ternary norm or vertex deletion.

The session goes further: at precision p^m a compatible lift exists iff
support(f) lies in one residue class modulo p^(m−1), and its coefficient
class is unique. The lift is multiplicative on its domain but not additive.
For example at p=2,m=2, the lift of 1+X^2 is 1+2X+X^2, rather than 1+X^2.

Evaluation at X=1 may erase a nonzero polynomial obstruction: for A=1+X and odd p,
E(A)/p at 1 is twice the base-two Fermat quotient. At a Wieferich prime that
number vanishes modulo p, while the coefficient of X in E(A)/p remains1.
Thus a zero scalar measurement does not imply a compatible polynomial lift.
This is a concrete hostile inside the same carrier, not an analogy between
unrelated exceptional primes.

## 5. Connection ledger and next proof questions

The common workflow is exact: construct an integral relation, divide its
known prime factor before reduction, identify the permitted correction
image, and test the remaining class. The target, correction space, and
normalization must be retained at every step.

This suggests three precise directions:

1. Full cyclic ternary lifts: find a uniform dual obstruction for the
   degree-p supported class when genuinely interior corrections are allowed.
   The support-preserving theorem does not discharge that larger problem.
2. Modular reconstruction exceptions: characterize which integral norm
   torsion classes produce zero-boundary polynomial solutions, rather than
   extrapolating the smallest observed degree from universal constructions.
3. Nonlinear automata: isolate an integral correction space compatible with
   the physical local rule before applying any linear lifting instrument.
   A formal linear kernel alone does not establish admissibility.

These are proof targets, not consequences already established. Neither
Gilbreath's prime conjecture, Rule30, Collatz, the Jacobian conjecture, nor
finiteness of Wilson/Wieferich primes follows from this synthesis.

## 6. A literal prime-power and zeta carrier inside the norm lattice

For odd n, let T(n) be the order of the prime-to-six torsion of the supported
cyclic norm cokernel in Section 2. The proved Smith classification gives

    T(n)=product_(p>=5 prime) p^(floor(log_p n)).

Set T(1)=1. Because every positive power of an odd prime is odd, the ratio
at successive odd degrees has the exact law

    T(n)/T(n-2) = p if n=p^r for a prime p>=5 and r>=1,
                 1 otherwise.                          (n>=3 odd)

Thus the elementary Z/p torsion multiplicity increases by one at each
p-power degree. No comparison map between the degree cokernels is claimed.
It does not become Z/p^2 at degree p^2: the new event adds a second
independent Z/p class. This is a concrete distinction between multiplicity
of prime-power arrivals and depth of p-adic lifting.

Writing psi(x)=sum_(p^r<=x) log p, one obtains exactly

    log T(n)=psi(n)-floor(log_2 n)log2-floor(log_3 n)log3.

Equivalently, for Re(s)>1, absolute convergence permits the Euler-product
calculation

    sum_(n>=3 odd) log(T(n)/T(n-2))/n^s
      = -zeta'(s)/zeta(s)
        -log2/(2^s-1)-log3/(3^s-1).

This is an explicit map from a family of integral relation matrices to the
prime-power von Mangoldt measure and a logarithmic zeta derivative. It
connects the user's earlier prime-power theme to the supplied ternary-period
paper through a proved intermediate lattice, without assigning significance
to a matching small number.

The new matrix carrier supplies no new error term for psi and no assertion
about Riemann zeros or prime races. The elementary Smith proof itself uses
the prime-power structure of binomial coefficients, so the identity does not
create independent analytic information. The exceptional even-degree
interior 29-torsion class is outside this supported odd family and must not
be forced into this formula.

## Artifact links

- [Integral tournament lattice](../tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md)
- [Ternary norm reconstruction packet](../ternary_frobenius_reconstruction_20261019/)
- [Wronskian and scalar Frobenius packet](../wronskian_frobenius_20261019/)
- [Finite Gilbreath entry](../gilbreath_binary_entry_20261019/BINARY_ENTRY_AND_EXACT_MIXING.md)

The prime-power Dirichlet series is the classical Euler-product identity for
−ζ′/ζ, applied to the independently stated integral torsion formula. Kummer's
carry theorem and first-lifting/Bockstein linear algebra are classical inputs.
