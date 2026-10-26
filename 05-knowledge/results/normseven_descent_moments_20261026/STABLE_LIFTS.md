# Root automorphisms, determinant classes, and the exact price of a torus lift

2026-10-26. **PROVED** L1-L3, with exhaustive bounded controls. This extends
the earlier quadratic Laurent-lift observation to general m and identifies
its minimal stabilization. The linear-algebra ingredients are classical;
the point here is the marked root/torus/section contract, not historical novelty.

**Concurrent integration:** incoming `2daf5e00af`,
[F1 L1-L4](../normseven_arithmetic_f1_20261026/f1/RESULT.md), proves the same
criterion, minimal compensation and explicit 2-by-2 block. These are shared
results. This packet adds an explicit v=1 section equation and independent
controls over noncyclic unit groups; it does not claim a second new theorem.

## L1. Fixed-rank criterion

Let m>=2, n>=1 and C an integer n-by-n matrix invertible modulo m. As in
`normseven_f1_20261025/KUMMER_AND_F1.md` K1, for any nonzero commutative ring
R and units b_i, the monomial map on the m-th-root algebra is an isomorphism.
A root coordinate may be multiplied by a Laurent monomial in the known b_i.
This replaces its exponent row by that row plus m times an integer row,
without changing its individual character modulo m.

There exists Q in GL_n(Z) with Q=C modulo m **if and only if**

    det C = 1 or -1 modulo m.                              (L1)

Necessity is immediate from det Q=+/-1. For n=1 it is also sufficient.
For n>=2, SL_n(Z)->SL_n(Z/m) is surjective by the elementary construction
below. Multiply the first row of C by epsilon=+/-1 to make its determinant
one modulo m; lift that matrix to determinant one over Z and multiply the
first row back by epsilon. This gives Q.

For completeness, an invertible first column modulo m generates the unit
ideal. For each prime p dividing m, add either no row (if its first entry
is nonzero mod p), or one row with nonzero first entry, to row one.
Choose the row-addition coefficients simultaneously by CRT; the new first
entry is a unit modulo m. Clear the remaining first-column entries.
Diagonal scaling diag(u,u^-1) on two rows is elementary, since

    w(u)=E12(u) E21(-u^-1) E12(u),
    w(u) w(-1)=diag(u,u^-1).

Thus normalize the pivot using determinant-one operations, repeat on the
remaining block, and clear the last upper entries. The last pivot is one
because the determinant is one. Every elementary matrix lifts to the same
integer elementary matrix; reverse the operations to get an integer lift.
No division by a nonunit or assumption that m is prime occurs.

When Q=C+mD, changing w_j to w_j product_i b_i^D_ji gives exponent Q.
The source units are essential. In characteristic dividing m this is still
a root-algebra isomorphism, but it is not an etale cover of distinct roots.

## L2. One auxiliary root always suffices, and sometimes is necessary

Choose a modulo m with a det C=1. Then diag(C,a) has determinant one
modulo m, so it has a lift Q in SL_(n+1)(Z). Add a root v^m=c for a specified
unit c; modulo m this has no cross terms with the old characters. As above,
the unimodular lift differs only by explicit base-unit monomials.

If det C is not +/-1 modulo m, no same-rank lift exists by L1, whereas this
one-coordinate stabilization exists. The minimum number of extra coordinates
under these operations is therefore exactly zero or one. Allowing an arbitrary
replacement of the original characters is a different problem and can erase
the question by selecting the identity basis; the marking is load-bearing.

Taking c=1 adds the group scheme mu_m, **not a geometrically invisible point**.
To retain the original root fibre, select the rational section v=1 and carry
its image under Q. If V_i are target torus coordinates, that image is the
primitive monomial equation

    product_i V_i^((Q^-1)_(n+1,i))=1.                     (L2)

Discarding it enlarges the object. Over fields not containing all m-th roots
of unity the added group scheme also has nontrivial Frobenius action, so the
constant-group averaging theorem in MOMENTS does not automatically apply.

## L3. An explicit six/seven map and its boundaries

For root order seven, C=diag(2,1,1,1,1,1) is an automorphism of the six-root
algebra but has no same-rank torus lift. Add one coordinate and use the block

    Q_(u,v) = [2  7],      Q^-1 = [25 -7].
              [7 25]                 [-7  2]

Its determinant is one and it reduces to diag(2,4) modulo seven. Identity
on the other five coordinates gives an explicit seven-dimensional torus
automorphism. On u^7=b,v^7=1 it gives

    W=u^2 v^7=u^2,       Z=u^7 v^25=b v^4,
    u=W^25 Z^-7,         v=W^-7 Z^2.                       (L3)

On v=1 the target must retain Z^2=W^7; with the fixed source b it also has
Z=b. This is an actual Laurent monomial map with an inverse and a marked
subtorus. Over the positive free monoid it is not an automorphism: its
inverse has negative entries. Adjoining this coordinate does not compactify
the map over branch zeros or absorb the additions E+x, E+z.

For root order six every unit residue is +/-1, so every invertible character
matrix has a same-rank lift. The full list of m>=2 with this universal
property is 2,3,4,6. Proof: the property is equivalent to phi(m)<=2. The
formula phi(m)=product_(p^e || m) p^(e-1)(p-1) restricts prime factors to
2,3 and immediately leaves those four m. For m=7, only two of the six
determinant classes lift. The determinant map on GL_n(F_7) is surjective
with equal-size fibres, hence exactly one third lift at fixed rank.

This six/seven distinction is about **root order**. The explicit six-to-seven
stabilization is about **torus dimension**. Neither is the Eisenstein norm
seven, seven branch forms, nor a K3/general-type transition. The original
norm-seven construction has quadratic roots and is already in the always
liftable m=2 case. This is a precise operation contrast, not an identification
of the two geometric families.

## F1 and arithmetic contract

In the pointed-monoid convention, the split torus is Z^n with an absorbing
zero. Its monomial automorphisms are exactly GL_n(Z), and base extension
gives the corresponding Laurent torus over Z. Thus Q and its marked primitive
subtorus are actual multiplicative F1 data. The convention is inherited from
the earlier packet and checked against [Deitmar, sections 1-3](https://arxiv.org/html/math/0404185v7)
on 2026-10-26 (his original monoids are unpointed).

Source: marked m-th-root algebra and source units. Target: a possibly
stabilized Laurent torus with its image root equations and marked section.
Operation: determinant compensation, elementary integral lift, and base-unit
correction. Preserved: each old root character, the full root algebra and
the chosen original fibre. Destroyed if omitted: arithmetic constants,
section, branch regularity and additive specialization. Cheapest hostile:
m=7,C=(2), compared with m=6,C=(5). No F1 torification of the surface,
new rational point, or open number-theory theorem follows.
