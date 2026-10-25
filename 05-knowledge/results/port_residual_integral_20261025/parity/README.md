# Odd symmetry gives a universal parity law without new operations

2026-10-25. **PROVED:** the user's general modulo-two identity, the parity
of pair-free tournament counts, and a two-term parity law for the positive
clone interpolation. Integral counting, preserved operations and positivity
are kept as three separate requirements throughout.

**Concurrent inheritance:** incoming commit f093b00695 independently proves
the general identity (P1) in Theorem1 of
[the all-odd-action interface](../../session63_integral_observers_root_arrangements_20261025/cancellation/RESULT.md).
We recover that identity here as the user's stated input. The canonical
inflation, pair-free parity and positive-interpolation consequences below
are the extensions in this packet; no historical-priority claim is made.

## 1. The general identity and its exact reason

Let a finite odd-order group A act on an n-element gap set G, faithfully
or otherwise. Let P_A(t) count its orbits on subsets, graded by subset size.
Every subset orbit has odd cardinality, by orbit-stabilizer. Therefore,
grade by grade,

    P_A(t) = (1+t)^n modulo 2.                       (P1)

This argument neither divides a group average in a bad characteristic nor
replaces a full automorphism group by a chosen subgroup. Both orbit sets
individually satisfy the congruence; their integer counts can differ.

For the positive cancellation polynomial

    Q_A(u)=(-1)^n P_A(-1-u),

the same statement is

    Q_A(u)=u^n modulo 2.                            (P2)

Thus its characteristic-two depth is exactly n, including n=0 with depth0.
The [rational minimum-cycle formula](../../unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md)
can have much smaller depth. Formula (P2) does not change the available
operations: constant XOR toggles that descend to subset orbits are still
exactly the unions of A-orbits on G. Their rank is the number of such
orbits, independent of the coefficient field. In the elementary abelian
regular family it stays one while n grows without bound.

Oddness is load-bearing. A transposition on two gaps has P=1+t+t^2,
which differs from (1+t)^2 modulo two. Its middle coefficient is the
smallest failure of the unqualified congruence.

## 2. Positivity supplies a correction, not an equivariant quotient

There is an exact nonnegative integral correction before reduction:

    B_A(t)=((1+t)^n-P_A(t))/2
          =sum_(subset orbits O) (|O|-1)/2 * t^rank(O). (P3)

Its coefficients are nonnegative integers because each |O| is odd.
This is a genuine positive weighted enumeration of the excess labelled
objects in each orbit. It does not construct an equivariant object of
each prescribed smaller cardinality.

For example let C5 act regularly on five gaps. At rank one, (P3) has
coefficient two. There is no equivariant surjection from that transitive
five-element singleton orbit onto a two-element C5-set: its image would
be transitive, with cardinality dividing five. An arbitrary two-element
set realizes the coefficient, but not the requested quotient operation.
Likewise pairing all but one point in that orbit requires an extra choice
of the unpaired point; it cannot be invariant under the transitive action.

This is the three-part distinction in a minimal example. Integral
coordinates are the coefficients of P and B. The operation requirement
is descent through the actual action. Positivity makes abstract counts
possible, but cannot by itself supply that descent or a geometric model.

## 3. Canonical clone inflation makes pair-free parity explicit

Write a_N=A000568(N), and let c_N count pair-free tournament classes,
with a_0=c_0=1. The canonical maximal-transitive-module decomposition
gives a unique pair-free quotient. This uses
[THM-4463, homogeneous tournament contacts](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
and the exact [clone-cut carrier](../../ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md).
For a fixed quotient R on q>0 vertices, inflations of total order N are
Aut(R)-orbits of positive integer weight vectors summing to N. Every
such orbit is odd. Their count is therefore

    binom(N-1,q-1) modulo 2.                        (P4)

No new clone merger is introduced: R is pair-free, the inflated blocks
are exactly the maximal transitive modules, and their actual ordered
weights are retained until taking the automorphism orbits.

Let A(x)=sum a_N*x^N and C(x)=sum c_N*x^N. In F2[[x]], (P4) gives

    A(x)=C(x/(1-x)),
    C(x)=A(x/(1+x)).                                (P5)

The inherited parity theorem, [THM-4750 section6](../../../../01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md),
gives A(x)=1+x+x^2 modulo two. Its elementary reason is that a_N has the
same parity as the groupoid mass 2^binom(N,2)/N!: every Aut(T) is odd,
and its reciprocal is one in F2. The latter mass is even for N>=3.
Substituting in (P5) gives

    C(x)=1+x/(1+x^2) modulo 2.                      (P6)

Hence the pair-free count is odd exactly at positive odd orders; its
order-zero value is one. This is an all-order consequence of canonical
inflation, not extrapolation from the computed prefix.

The integer version of (P5) is false. At N=4 its putative right-hand
coefficient is 1+3*c_3+c_4=6, whereas a_4=4. The C3 quotient has three
positive weight vectors of type (2,1,1), which form just one orbit.
The missing two objects vanish modulo two, exactly as (P3) predicts.
This is the first integer failure, retained as a hostile control.

## 4. The whole positive interpolation has only two surviving coefficients

Use the positive interpolation already proved in
[POSITIVE_JET.md section6](../../unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md):

    R_N(u)=sum_M b_(M,N-M)*(1+u)^M
          =sum_([T], |T|=N) (1+u)^q(T)*Q_T(u),       (P7)

where q(T) is the number of maximal transitive modules and the b's are
the signed clone-excess coefficients. Its coefficients over Z are
nonnegative; its constant term is c_N and its leading coefficient is a_N.

Apply (P2) to every T, then group by its pair-free quotient order q and
use (P4),(P6). For N>0,

    R_N(u)=sum_(q odd) binom(N-1,q-1)
                        (1+u)^q*u^(N-q) modulo 2.   (P8)

If N=2m+1, write q=2j+1 and use
binom(2m,2j)=binom(m,j) modulo2, obtained by comparing coefficients of
(1+x)^(2m)=(1+x^2)^m in F2[x]. Then (P8) becomes

    (1+u)*((1+u)^2+u^2)^m=1+u.

For N=2m, the analogous identity
binom(2m-1,2j)=binom(m-1,j) modulo2 gives u*(1+u).
Consequently the complete parity pattern is

    R_0(u)=1;
    R_N(u)=1+u modulo2,       N odd;
    R_N(u)=u+u^2 modulo2,     N positive even.        (P9)

In particular R_(N+2)=R_N modulo2 for every N>=1. A low-dimensional
modular observable has emerged from the exact clone semantics. It
does not imply that the integer polynomials, the action, or the literal
tournament-move graph have only two coordinates. The higher integral
coefficients remain nonnegative even integers, not absent objects.

## 5. Inheritance, checks and next use

The closest mechanism is the inherited positive-jet and canonical-cut
theorem. The smallest hostile is the even transposition action; the
integer clone-transform hostile is C3 inflated to order4. The previously
unused operation is reduction grade by grade before the binomial
inflation transform. The concept board is odd orbit cardinality /
integer positive correction / canonical weighted quotient / field
reduction / actual descended toggles.

The checker constructs small permutation actions, enumerates their subset
orbits literally, and compares P, Q and B. A second path enumerates every
tournament class through order6, computes its actual homogeneous-pair
components and full automorphisms, and sums (P7). It compares with a
direct signed odd-partition/substitution calculation through20, checking
integrality, positivity and (P9). These are bounded controls of the
stated proofs, with explicit gates under normal and `-O` execution.

The useful next application is a control for the new integral recurrence:
its final pair-free parity must satisfy (P6), and its bivariate extension
must satisfy (P9). Neither identity licenses replacing the exact recurrence
by a Boolean toggle model. The [integral gauge](../counting/README.md)
preserves actual transitions; this parity calculation specifies a smaller
observable of their final result.
