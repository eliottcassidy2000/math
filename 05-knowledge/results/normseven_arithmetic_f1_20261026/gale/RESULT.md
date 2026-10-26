# The outer Gale sign is an arithmetic isomorphism invariant

2026-10-26. **PROVED from the inherited exact geometric rank and the
standard elliptic-surface divisor/section theorem.** The computations below
are **FINITE-EXACT** controls. This is an application of classical quadratic
twisting and Galois descent, not a historical-priority claim.

## 1. Source, correction boundary, and exact statement

The starting result is [BG2](../../six_seven_gale_cube_20261025/BUDGET_AND_GALE.md):
the outer six-line Gale discriminant double plane retains a minus-one twist.
That note deliberately did not decide whether another map could absorb the
sign. We now decide this for its actual outer specialization `u=v=1`.

Put

    F(E,x,z)=E(E^2-x^2)(E^2-z^2)(E+x+z).

For a number field k and delta in k*, let Z_delta be the smooth minimal
resolution of the normal double plane `w^2=delta F`. Define
`rho_k(Z)=dim_Q(NS(Z_kbar) tensor Q)^(Gal(kbar/k))`; this is the
arithmetic Neron-Severi rank, distinguished from geometric Picard rank.
Let chi_delta be the quadratic Galois character of k(sqrt(delta))/k,
including the trivial character when delta is a square.

**Theorem GA.** For every such k and delta:

    NS((Z_delta)_kbar) tensor Q = Q^18 direct-sum Q(chi_delta)

as Galois representations. In particular:

1. `rho_k(Z_delta)=19` for square delta and `18` otherwise.
2. `Z_delta` and `Z_epsilon` are isomorphic over k **if and only if**
   delta/epsilon is a square in k. Necessity holds for arbitrary abstract
   surface isomorphisms, without preserving the double-plane map or fibration.
3. In the marked elliptic fibration over k(t), the generic curve has rank
   one for square delta and rank zero otherwise. Its torsion is `(C2)^2`.

Thus the actual outer Gale double plane `Z_-1` is not Q-isomorphic to
`Z_1`; it becomes isomorphic over Q(i). All these surfaces are isomorphic
geometrically and have geometric Picard rank nineteen. No statement here
identifies them with the degree-eight six-root K3, whose map to the
double-plane quotient has degree sixteen.

Inheritance: [exact outer rank](../../session60_six_seven_geometry_20261025/outer-picard/EXACT_RANK_COROLLARY.md)
proves geometric rank19, geometric Mordell-Weil rank1, and the rational
nontorsion section P used below. [The independent upper-bound proof,
Sections 1-2](../../session60_six_seven_geometry_20261025/picard-upper/INDEPENDENT_UPPER_BOUND.md)
supplies the rational branch-resolution lattice. These are load-bearing,
not inferred from the new finite counts.

The closest proved mechanism is that rational nontorsion section; the
canonical hostile is dropping a constant square class from a geometric
isomorphism. The corrected near miss is mistaking the determinant of a
chosen monomial lift for an intrinsic cover invariant, repaired in
[KUMMER_AND_F1](../../normseven_f1_20261025/KUMMER_AND_F1.md).
Here the retained invariant is an actual Galois representation. The
least-used sidecar is the **anti-invariant divisor class**, rather than
the geometric rank alone.

## 2. One explicit anti-invariant class completes the rational divisor basis

On E=1 set `t=x`, `d=1-t^2`, `X=-dz`, `Y=dw`. For Z_1 this gives

    Y^2=(X-d)(X+d)(X-d(t+1)),
    P=(dt,d^2).

The involution iota:w->-w is exactly elliptic negation `(X,Y)->(X,-Y)`.
It fixes the zero section. Let Pbar denote the section divisor, and put

    A=[Pbar]-[(-P)bar] in NS((Z_1)_kbar) tensor Q.

The standard divisor/section isomorphism `NS/Triv = MW` sends A to 2P,
which is nontorsion. Therefore A is nonzero, and `iota(A)=-A`.
Both section divisors are defined over Q.

The double plane has nine ordinary double points and two ordinary triple
points. Its resolution contributes nine A1 curves and two D4 trees, all
defined over Q; the pullback of a plane line gives one more class. Their
intersection lattice is

    <2> direct-sum A1(-1)^9 direct-sum D4(-1)^2,

of rank18 and determinant -16384. Every one of these classes is fixed by
iota. For a node, its unique exceptional curve is invariant. For a triple,
the central branch exceptional curve is invariant, and each of the three
arms lies over its own rational crossing, hence is invariant individually.
These statements concern divisor classes; iota need not fix every point
of every exceptional curve.

The invariant rank18 lattice and the nonzero anti-invariant class A are
independent. Since the inherited geometric rank is exactly19, they span
the entire rational Neron-Severi space. Thus Galois acts trivially on
`NS((Z_1)_kbar) tensor Q`, and iota has eigenspaces of dimensions18 and1.

Over k(sqrt(delta)) the explicit isomorphism from Z_delta to Z_1 is
`w -> w/sqrt(delta)` on the normal models. It lifts uniquely to their
minimal resolutions. For sigma in Gal(kbar/k), its descent cocycle is
identity or iota according to chi_delta(sigma). Twisting the preceding
19-dimensional representation therefore gives exactly the formula in GA.

Any k-isomorphism of surfaces induces a Galois-equivariant isomorphism of
their Neron-Severi spaces. Their characters have traces
`18+chi_delta(sigma)`, so an isomorphism forces chi_delta=chi_epsilon,
equivalently delta/epsilon square. Conversely, if delta/epsilon=c^2,
`w -> w/c` gives the isomorphism. This proves the abstract classification.
It does not compute the full automorphism group or an integral NS lattice.

## 3. Exact Mordell-Weil consequence and its limits

The twisted generic fiber is

    E_delta: Y^2=delta (X-d)(X+d)(X-d(t+1)).

Its point at infinity is the zero section, and the three roots give the
same rational two-torsion points for every delta. Over kbar(t) it is the
untwisted curve. That curve has free rank1, spanned over Q by P. Since P
is rational, Galois acts trivially on this one-dimensional free part.
The twist changes this action to chi_delta. A rational point on a
nonsquare twist would give an invariant free class, so its free class is
zero. Hence every rational section is torsion when delta is nonsquare.
For square delta, the section P transports and the rank is one.

The geometric fibration has an I0* additive fiber at t=1. In characteristic
zero its torsion injects into that fiber's component group `(C2)^2`.
The three displayed rational two-torsion points already exhaust this bound.
Thus the geometric torsion, and the torsion of every constant twist over
any number field, are exactly `(C2)^2`. For nonsquare delta this proves the
complete group statement `E_delta(k(t))=(C2)^2`.

The section's loss is visible before any rank theory: on the old line
`z=-t`, the equation becomes `w^2=delta(1-t^2)^2`. For nonsquare delta it
has no solution in k(t), while after adjoining sqrt(delta) the two sections
are `w=+/-sqrt(delta)(1-t^2)`. The rank argument proves that replacing this
line by another rational section cannot evade the obstruction.

This does **not** say Z_delta has no rational points, no rational curves,
or rank zero in every specialized elliptic fiber. Branch points and the
zero/two-torsion sections remain, and specialization can gain rank. The
rank assertion keeps this marked fibration and the original constant field.

The classical inputs in this section and Section2 are **CITED** from
[Schuett-Shioda, Elliptic Surfaces, arXiv:0907.0298v3](https://arxiv.org/pdf/0907.0298):
Theorem6.3 (printed p.23) identifies the divisor quotient with sections;
Lemma7.3 (p.32) gives the I0* component group, and Lemma7.8 (p.35) gives
the prime-to-characteristic torsion injection for an additive fiber.
The primary text was checked 2026-10-26. The inherited fibration has
singular fibers and characteristic zero, as those inputs require.

## 4. The Gale map and an exact arithmetic residue witness

For the actual outer forms `(E,E+x,E-x,E+z,E-z,E+x+z)`, the net relation
matrix is

    C=[-2,1,1,0,0,0; -2,0,0,1,1,0; 1,-1,0,-1,0,1].

Put `lambda=(E+x,E+z,2E)`. Directly, its six diagonal entries are

    C^t lambda=(-2(E+x+z), -(E-x), E+x, -(E-z), E+z, 2E).

Their product is `-4F`. Dividing the double-cover coordinate by2 therefore
gives Z_-1 over Q. This restates the inherited explicit map in the
coordinates to which GA applies; the new conclusion is its rational
nonisomorphism to Z_1 under **any** possible alternative isomorphism.

For a finite field F_q of odd characteristic greater than3 and nonzero
delta in F_q, write

    S_q=sum_(x,z in F_q) chi((1-x^2)(1-z^2)(1+x+z)).

Then the same branch resolution proves

    #Z_delta(F_q)=q^2+18q+1+chi(delta) S_q.             (GA-count)

There are nine A1 corrections q and two D4 corrections4q. This remains
true for nonsquare delta: an exceptional node conic is a double cover of
P1 branched at two rational points, so it has a rational ramification point
and is P1; the D4 central component and its three arms are rational too.
The affine raw double cover contributes `q^2+chi(delta) S_q`, and the
line E=0 contributes q+1. This proves GA-count without changing the twist
on extension fields.

At7, direct independent counting gives `S_7=4`, hence

    #Z_1(F_7)=180,             #Z_-1(F_7)=172.

At13, -1 is square and both counts are408. At49 the fixed -1 twist becomes
square, so both counts are3384. These are useful hostile/positive controls:
equality at one split prime, or over an even extension, does not erase the
arithmetic character. The all-twist average and even-extension behavior
are inherited more generally in
[Session63](../../session63_integral_observers_root_arrangements_20261025/twists/ZETA_AND_EFFECTIVITY.md);
the present advance is the exact NS character and twist classification.

## 5. Connection contract, board update, and reproduction

Source -> target: the outer Gale determinant double plane -> constant
quadratic twists of the outer K3 -> Galois action on its divisor space.
Map: the displayed lambda substitution and `w/2`, followed by
`w/sqrt(delta)` and its descent cocycle. Preserved: the geometric surface,
branch labels, and elliptic generic fiber after constant extension.
Lost by geometric identification: the one-dimensional anti-divisor
character and rational free section. Necessary sidecar: delta in k*/k*^2.
Decisive test: the class `[P]-[-P]`, with finite witness180 versus172 at7.

The lane's board is (i) geometric rank19, (ii) the rational branch lattice,
(iii) P and its anti-class, (iv) the Gale sign, (v) finite Frobenius counts.
The new anti-class connects every item: it fills rank19, lies outside the
branch lattice, transforms by the Gale sign, and explains why equal
geometric invariants cannot justify equal arithmetic counts. It does not
produce an F1 torification: the proof uses an additive branch map and
Galois descent, not just the monomial carrier.

Anchor outcome: the minus-one ambiguity is closed over Q for the actual
outer specialization. Niche outcome: all constant twists are classified
by square classes, with exact arithmetic Picard/Mordell-Weil ranks.
Wildcard hostile: split-prime/even-extension count equality is recorded
as insufficient. Integral NS discriminants, primitive generators and
nonconstant twists are not determined. No new META-PATTERNS card is proposed;
this uses "Inventory retained power classes" and "Controlled forgetting".

From the repository root:

    python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/check.py --output 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/normal.json
    python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/check.py --output 05-knowledge/results/normseven_arithmetic_f1_20261026/gale/optimized.json

The standard-library checker counts the whole declared projective base,
uses both root multiplicities and Euler characters, tests every nonzero
twist over the declared fields, verifies the branch incidence and Gale
identity, and checks rational-section identities. Computation corroborates
the formula and witnesses; the inherited geometric rank and the proof
above carry the unbounded classification. Normal/optimized receipts must
agree byte for byte. Parent-session audit and hashes are recorded separately.
