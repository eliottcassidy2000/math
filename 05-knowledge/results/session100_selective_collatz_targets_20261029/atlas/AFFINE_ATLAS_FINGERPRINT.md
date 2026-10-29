# A multiplier-selective affine atlas: ternary order and rational rigidity

**PROVED:** branch-preserving rational rigidity and the minimal-length
permutation-contrast fingerprint below. Finite checks are separate. This is
an obstruction to losing the multiplier in a proposed reframe, not a
convergence proof. No global novelty claim is made for rational-function
rigidity or elementary affine composition.

## 1. The target must carry the branch charts, not just their symbolic graph

For odd positive multiplier b>=3 define the rational affine chart

    f_(b,a)(x)=(b x+1)/2^a, a>=1.

Its native domain is the exact valuation cylinder v2(bx+1)=a. These domains
are part of the structure. The valuation-word conjugacy in ../bijections/
forgets b at the symbolic level; its inverse source-residue map restores it.
Here we give a second way to detect whether a candidate target has lost b:
retain the affine charts and their order contrasts.

A word w=(a1,...,aJ) has D=sum a_i and

    F_(b,w)(x)=(b^J x+C_b(w))/2^D,
    C_b(w)=sum_(i=0)^(J-1) b^(J-1-i)2^S_i,
    S_i=a1+...+ai, S_0=0.                                (1)

A positive odd source follows the word exactly iff

    x ≡ (2^D-C_b(w)) b^(-J) mod 2^(D+1).                  (2)

The extra bit in the modulus enforces an odd final quotient and hence exact
valuations, rather than only divisibility. The inverse exists because b is
odd. The bijection lane proves this guard and its nested-residue property.

## 2. Two order contrasts forget b; three recover it

Specify the observer precisely: compare affine outputs of permutations of a
fixed valuation multiset, retaining their differences but discarding a
common additive offset. Since J and D agree, the input-dependent slope in
(1) cancels. This is a counterfactual chart observer, not the assertion that
one native integer source can legally follow all the words.

For two slots,

    C_b(a,c)=b+2^a,
    C_b(a,c)-C_b(c,a)=2^a-2^c,

independent of b. Absolute two-slot carries still see b; the blindness claim
is ONLY for this difference observer.

At three slots, use the smallest unequal positive multiset {1,1,2}. In the
orders112,121,211 the carries are

    b^2+2b+4,  b^2+2b+8,  b^2+4b+8.

Their differences are4 and2b. Therefore

    [F_(b,211)(x)-F_(b,121)(x)] /
    [F_(b,121)(x)-F_(b,112)(x)] = b/2.                    (3)

The common denominator16 and the common slope cancel. This ratio is unchanged
by any common nondegenerate affine re-expression of the three outputs. It is
3/2 for Collatz and5/2 for the hostile multiplier. Among word-length order
contrasts of this specified kind, three is the first length that can recover
b: length one has no nontrivial contrast, and length two was proved blind.
No minimum claim is made for arbitrary possible observers.

The exact source cylinders differ as well. For b=3 the residues of112,121,211
are7,27,9 modulo32; for b=5 they are29,9,23. They are pairwise disjoint within
each parameter. A target may store this three-row chart atlas, but may not
pretend to run all rows on the same source without carrying their guards.

This makes a selective reframe possible: require the chart ratio3/2, in
addition to the source-residue sidecar. The5x+1 atlas fails that predicate
for an algebraic reason independent of whether any orbit converges. It does
not prove that all systems passing the predicate converge.

## 3. Rational changes of coordinates cannot erase the multiplier

**Rigidity theorem.** Let b,c>=3 be odd integers and R(x) a nonconstant
rational function over a characteristic-zero field. Suppose, as rational
identities,

    R((b x+1)/2^a)=(c R(x)+1)/2^a, a=1,2.                (4)

Then b=c and R(x)=x. The conclusion applies if (4) is known only on each
native branch cylinder wherever R is defined: infinitely many points on
that cylinder force the rational identity after denominators are cleared.
The two cylinders need not intersect.

Proof: f_(b,2)(x)=f_(b,1)(x)/2, and similarly on the right of (4). Since
f_(b,1) is invertible affine, comparison gives

    R(y/2)=R(y)/2.                                      (5)

Expand R as a formal Laurent series at zero. Equation(5) forces every
coefficient except degree1 to vanish, since2^(-j)!=2^(-1) for j!=1.
Thus R(y)=A y with A nonzero. Substituting into (4) forces A=1 by the
constant terms and then b=c by the linear terms. QED.

A useful affine-only derivation writes R(x)=Ax+B: the x coefficient forces
b=c; the constants at a=1 and2 force B=0,A=1. The checker exhausts small
integer A,B as hostile controls, while the proof above handles arbitrary
rational functions, not just that finite bank.

Thus the genuine2-adic valuation conjugacy between3x+1 and5x+1 cannot be a
rational coordinate change preserving these branch charts. There is no
conflict with the conjugacy theorem: it is a much less rigid homeomorphism
and does not preserve ordinary integer sources or affine output differences.
This prevents a false transfer of real/rational geometry through an abstract
symbolic bijection.

## 4. Relation to the repository's joint-order theme

The structural parallel with earlier tournament port work is scoped:
a lower-order response can lose relative order that a later compound
operation detects. Here we have an explicit multiplier recovered from three
order rows, and explicit two-row difference blindness. No isomorphism with
the 81/27 tournament pair or a tournament Hamilton-path count is asserted.
The reusable lesson is to test joint order before declaring a representation
faithful. The chart guards are the arithmetic analogue of named port data.

## 5. Reproduction and proof boundary

Run `python3 check_affine_atlas.py` and `python3 -O check_affine_atlas.py`.
The checker tests odd multipliers3..31, exact three-row source cylinders and
native trajectories at25 lifts each, every two-slot contrast with digits1..7,
and the finite affine rigidity controls. No symbolic algebra dependency is
needed. Numerical values of an allegedly convergent orbit are not inputs.
The infinite rational-function theorem follows from the Laurent proof, not
from finite sampling.
