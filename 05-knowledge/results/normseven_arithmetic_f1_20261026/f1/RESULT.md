# The exact determinant obstruction to a marked Laurent repair

2026-10-26. **PROVED by the self-contained arguments below; independently
audited in [AUDIT.md](AUDIT.md).** The finite controls are **FINITE-EXACT**. Congruence
lifting is classical; the contribution here is its exact application to
the repository's marked root-fibre operation, its minimal compensating
character, and the affine-boundary separation. No historical priority is
claimed for the underlying linear algebra.

## Inheritance and scope

The closest mechanism is [K1 and K3--K4](../../normseven_f1_20261025/KUMMER_AND_F1.md):
an invertible root-character matrix modulo two reconstructs a marked root
fibre over every coefficient ring, and a congruent unimodular integer lift
exists. Its corrected near miss was treating the determinant three of one
chosen Boolean lift as intrinsic. The canonical hostile for general root
order is exponent two modulo five. The least-used sidecar is the determinant
class **modulo the root order**, with the individual character rows fixed.

The anchor here is the exact Laurent-lift criterion; the niche is minimal
stabilization; the wildcard is the distinction between a torus and its
positive affine monoid. The concept board is:

| Object / operation | Preserved target | Sidecar / decisive test |
|---|---|---|
| Root algebra / character matrix C | Entire marked root fibre | Root order m and base units; invertibility mod m |
| Laurent torus / integral lift Q | Individual character rows | det(C) mod m; require det(Q)=+/-1 |
| Added root coordinate | Original rows plus compensating character | Its exponent and radicand; block determinant |
| Affine positive monoid | Boundary as well as open torus | Positive cone; atom images must be a permutation |

This extends the inherited rank-one six/seven test to all ranks. It does not
identify root order seven with the Eisenstein norm-seven ideal: the latter
still selects **quadratic** character packets, so its root order is two.
The cyclic conductor [C1--C3](../../normseven_boundary_20261025/CYCLIC.md)
remains a separate boundary invariant; a Laurent lift does not remove it.

## L1. Complete criterion, with each character fixed

Let m>=2, n>=1, and C be an n-by-n integer matrix invertible modulo m.
There exists Q in GL_n(Z) with Q=C mod m if and only if

    det(C) = +1 or -1 mod m.                              (L1)

For a nonzero commutative ring R and units b_1,...,b_n, put

    A_b = R[u_i^+/-1]/(u_i^m-b_i),
    w_j = product_i u_i^C_ji.

The inherited universal root-algebra theorem already reconstructs A_b from
the w_j with their base radicands whenever C mod m is invertible. L1 answers
the stronger, different question: can multiplying the w_j by Laurent
monomials in the b_i make their **ambient torus exponent matrix** unimodular,
while retaining each marked character? Writing Q=C+mD, the required change is

    w'_j = w_j product_i b_i^D_ji.                       (L2)

It exists exactly under L1. It works over every R, including residue
characteristic dividing m; it is not then an assertion of an etale cover.

**Proof of necessity.** A unimodular Q has determinant +/-1, and congruence
of matrices implies congruence of determinants.

**Proof of sufficiency, without a field assumption.** The case n=1 is
immediate. Write S=Z/m. Every determinant-one matrix over S is a product of
elementary matrices. Here is a direct elimination proof, including the
composite-modulus issue. Its first column is unimodular: for each prime p
dividing m, at least one entry is nonzero modulo p. Choose coefficients t_i
such that a_1+sum_(i>1)t_i a_i is nonzero modulo every such p. For a given p,
take all t_i=0 if a_1 is already nonzero; otherwise take one t_i=1 at a
nonzero entry. The Chinese remainder theorem combines these choices. Row
additions therefore make the top entry a unit v in S. Clear the entries
below it by elementary row additions.

The diagonal operation diag(u,u^-1), acting on two rows, is elementary,
because, with E_ij(t)=I+t e_ij,

    diag(u,u^-1)
      = E12(u) E21(-u^-1) E12(u) E12(-1) E21(1) E12(-1).

Use u=v^-1 to make the first pivot one, and elementary column additions
to clear the rest of its row. The remaining block has determinant one.
Repeat. At size one the last entry is one. Thus the original matrix is a
product of elementary matrices and their inverses. Lifting every elementary
coefficient to an integer gives a determinant-one lift. If det(C)=-1 mod m,
first multiply its first row by -1, lift, and then multiply that row by -1
again. This proves L1. QED.

The elimination uses operations across rows to construct Q, but the final
congruence Q=C mod m retains **every** input character row. Intermediate
operations are not being declared permissible relabelings of observations.
The marked norm-seven U/D quotient subgroups are therefore retained when
m=2. The coefficient corrections in L2 remain essential.

## L2. Exactly which root orders allow every marked repair

For each fixed n>=1, all matrices in GL_n(Z/m) admit a unimodular lift iff

    m is one of 2,3,4,6.                                (L3)

Indeed the determinant map onto (Z/m)^* is surjective, via diag(a,1,...,1).
So L1 applies to all matrices exactly when every unit is +/-1. This is
equivalent to phi(m)<=2. The elementary formula for phi shows that any
prime divisor at least five is excluded; writing m=2^a 3^b then forces
a<=2, b<=1, and excludes a=2,b=1. The resulting list is L3.

More precisely the fraction of GL_n(Z/m) admitting a lift is

    |{+1,-1} mod m| / phi(m).                           (L4)

All determinant fibres have equal size, by multiplication by diag(a,1,...,1).
Thus **every** sixth-root character automorphism lifts, while exactly
**one-third** of seventh-root character automorphisms lift, at every rank.
Modulo five the fraction is one-half, so seven is not the first failure.

An actual sixth-root repair is u^5/b=u^-1 for u^6=b. By contrast, for
u^7=b the root-fibre coordinate w=u^2 is invertible with u=w^4/b, yet no
integer exponent congruent to two modulo seven equals +/-1. This is a
failure of a marked ambient torus lift, not a failure of fibre recovery.

## L3. One compensating character is sufficient and, when needed, minimal

Let delta=det(C) in (Z/m)^*. Adjoining identity coordinates cannot repair a
failure of L1: det(diag(C,I_r))=delta at every r. This remains true even if
the eventual lift may mix all coordinates by multiples of m.

But adjoining **one specified nonidentity character** always suffices:

    C^+ = diag(C,delta^-1) mod m                       (L5)

has determinant one and hence admits an integral unimodular lift. If delta
is not +/-1, zero auxiliary coordinates are impossible and one is possible,
so the auxiliary dimension one is minimal in this precisely stated operation
class. If delta is +/-1, the minimum is zero.

Arithmetic realization: add a unit c and a root v^m=c. On the enlarged root
algebra the character operation is (w_1,...,w_n,v^(delta^-1)). A lift Q^+
from L5 changes these coordinates only by monomials in b_1,...,b_n,c.
In particular each of the first n lifted coordinates has v-exponent
divisible by m, so any dependence on v reduces to the known base unit c.
Each original marked character is preserved; inside the enlarged algebra,
each subalgebra generated by a specified set of original character monomials
is unchanged after multiplying those generators by the known base units.
In the field case the same statement applies to the generated subfields.
The enlarged algebra is A_b tensor_R R[v^+/-1]/(v^m-c);
we are not claiming it is A_b or that it has a rational point when A_b does.

An exact seventh-root example is

    Q = [2  7] ,   det Q=1,   Q mod7=diag(2,4),
        [7 25]
    Q^-1 = [25 -7].
           [-7  2]

For u^7=b and v^7=c the two monomials are

    W=u^2 v^7=u^2 c,
    V=u^7 v^25=b v^4 c^3,
    u=W^25/V^7,    v=V^2/W^7.                          (L6)

Their seventh-power equations are W^7=b^2 c^7 and V^7=b^7 c^25.
The first character remains exponent two on u, and the compensating
character is exponent four on v. The analogous earliest hostile m=5 uses
[[2,5],[5,13]], congruent to diag(2,3), again with determinant one.

Allowing a free re-basing of target characters changes the question:
multiplying C by a matrix of determinant delta^-1 removes the obstruction.
This may preserve the generated quotient field while changing individual
observables. L1 concerns **fixed** character rows, exactly as the quadratic
fixed-character repair did. Neither a re-basing nor stabilization may be
silently called the same marked map.

## L4. No such criterion fills in the affine boundary

Adjoin an absorbing zero to the exponent groups. The resulting split-torus
monoid has automorphisms GL_n(Z). By contrast, the free positive monoid
N^n (also with absorbing zero adjoined) has only permutation automorphisms:
the nonzero additive atoms are exactly e_1,...,e_n, so a monoid isomorphism
must permute them. Equivalently an integral matrix and its inverse both
preserve N^n iff the matrix is a permutation matrix.

Consequently there exists Q=C mod m whose exponent map is an automorphism
of this **fixed affine positive monoid**, with no added units, iff C mod m
is a permutation matrix. This is much stronger than L1. It is an existence
claim about a congruent exponent map, not an assertion that the coordinate
corrections L2 extend regularly across the old zero locus. Already C=[m+1]
is congruent to [1], but replacing u^(m+1) by u requires division by b=u^m;
the Laurent decoder is not a regular decoder from that old observation at
b=0.
For example [[1,1],[0,1]] is already unimodular, but its affine monomial map

    (u,v) -> (uv,v)

collapses every point with v=0 to (0,0). The inverse divides by v.
This is the smallest dimension in which a unimodular nonnegative shear
exhibits that collapse. In dimension one the inversion u->u^-1 already
fails to extend as a positive monoid map.

The result is about a specified monoid and boundary. It excludes neither
different fans/toric modifications nor other F1 formalisms. In particular
it supplies no torification of the additive branch arrangement, K3 surface,
or Campedelli quotient. The inherited additive equations and normalization
remain separate and indispensable.

## Connection contracts, controls, and stopping reasons

The source is the marked m-th-root algebra with its units; the target is a
split Laurent torus coordinate presentation. The map is Q=C+mD with the
base-unit corrections L2. It preserves individual root characters and the
full arithmetic fibre. Forgetting m or the base units destroys the inverse;
forgetting the positive cone destroys boundary extension. The decisive
tests are det(C) modulo m, followed by the atom test for the affine monoid.

Every result updates all four board rows: the determinant class separates
root recovery from torus repair; the compensating character changes the
ambient torus but preserves the old characters; the affine atom test blocks
the false next implication in both dimensions. The inherited cyclic carry
conductor and residual square classes are not determined by any of these
determinants. A boundary atlas remains a separate next target.

Reproduce from this directory:

```sh
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
```

The two receipts must agree byte for byte. The checker constructs actual
unimodular lifts and independently computes the generated congruence subgroup
by breadth-first search for every rank-two matrix at m=2,...,12. It also
exhausts rank three at m=2,3 and runs deterministic seeded invertible inputs
at ranks1,...,6 and m=2,...,30, applying the one-character repair exactly
when needed. All failures use explicit exceptions, not Python assertions.
Its finite universe does not prove the all-rank or all-modulus statement.

**Source boundary.** The congruence mechanism is classical: Conrad,
[*SL2(Z)*, Theorem3.2](https://kconrad.math.uconn.edu/blurbs/grouptheory/SL%282%2CZ%29.pdf),
proves the rank-two surjectivity and explains the GL determinant obstruction.
The all-rank elimination above is included so its exact use is reviewable.
The monoid/base-extension convention follows Deitmar,
[*Schemes over F1*, section1](https://arxiv.org/html/math/0404185v7).
His monoids are unpointed; adjoining zero and contracting that zero in the
monoid algebra gives the pointed convention used here. Both sources checked
2026-10-26. The automorphisms of this torus are not being identified with
Deitmar's separately defined GL_n(F1) of a finite free module.
