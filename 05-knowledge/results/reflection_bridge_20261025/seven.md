# Seven objects: boundary contraction, projective supports, and the historical 3-1-3

**PROVED:** the finite-set quotient, reflection, incidence obstructions,
and coordinate formulas below. **FINITE-EXACT:** independent controls in
`seven_check.py`. The chosen pointed-monoid projective convention is
inherited from the audited Session59 source identified below. The earlier
Keller fragment is recovered as provenance; no Keller, magic-square,
tournament-count, or F1-descent theorem is inferred from matching sevens.

## 1. Inheritance and the actual earlier fragment

The closest operative mechanism is the full-coordinate divisor chart in
[arithmetic.md, sections 3--4](../trijection_reset_20261020/arithmetic.md).
It represents the proper divisors of `p^2qr` by ten planar points and
keeps complement distinct from exchange of the two unrepeated primes.
The canonical hostile is complement sending a nonempty full support to
the excluded empty support. The corrected near miss is the historical
`mu_3` interpretation below. The least-used sidecar is the distinction
between a finite-set quotient, a specialization poset, and Fano incidence.

The working board is divisor supports / boundary orbit / reflection /
projective point convention / incidence / iterated-map covariance. The
anchor is the 3-1-3 bridge, the niche is the precise lost incidence, and
the wildcard is checking whether the historical reflection iterates.

The exact old reference is
[FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md](../FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md).
Its correction banner and **MISTAKE-319, root audit of klein-S691** in
[MISTAKES.md](../../../01-canon/MISTAKES.md) identify the fragment as three
successive inverse fibers of the map in
[THM-2473, sporadic Keller branch tower](../../../01-canon/theorems/THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md).
The current theorem file labels its exact algebraic claims (1)--(6)
`PROVED`, finite tower checks `FINITE-EXACT`, and its off-spine tower law
`CONDITIONAL`. This note does not use those global claims as dependencies.

The fragment records first-layer targets

    P-=(-1,3/2,13/2), P0=(0,0,-1/4), P+=(1,-3/2,13/2)

and respective second-layer fiber sizes `3,1,3`. It corrects a former
finite-collision/Arithmetic-Kakeya interpretation to escape of sheets
at infinity and withdraws the `mu_3` label. Session58's
[hourglass theorem, section 1](../session58_hourglass_complement_geometry_20261025/hourglasses/THEOREM.md)
already recovered this provenance and separated it from a seven-cell
magic hourglass. The Boolean seven considered here is a third object.

There is a further elementary scope check on the displayed reflection.
Write the fragment's map as

    u=1+xy,
    F=(u^3 z+y^2 u(4+3xy),
       y+3xu^2z+3xy^2(4+3xy),
       2x-3x^2y-x^3z).

Let `sigma(x,y,z)=(-x,-y,z)` and `tau(a,b,c)=(a,-b,-c)`.
Substitution, term by term, gives

    F sigma = tau F,              not F sigma = sigma F.       (1)

Since `tau` fixes `(-1/4,0,0)`, `sigma` acts on that first fiber,
exchanging `P-,P+` and fixing `P0`. It does not thereby act on the
second-level seven leaves. Indeed

    w=(-1/8,0,0),       F(w)=P0,
    F(sigma(w))=(0,0,1/4) != sigma(P0)=P0.                     (2)

Thus the first-layer fixing symmetry is valid, but that same coordinate
reflection fails already on the central second-layer leaf. A reflection
of all seven leaves requires another specified map. Counts `3,1,3`
alone do not provide it. This observation refines the connection scope;
it does not contradict the fragment's explicitly first-layer equation.

## 2. A precise Boolean 3-1-3 and the boundary it contracts

Let `E={p,q,r}`, `B=2^E`, and `C(A)=E\A`. On all eight subsets, `C`
is a fixed-point-free involution with four two-element orbits. It does
not restrict to the seven nonempty subsets because `C(E)=empty`.

Instead form the finite-set quotient

    Y=B/(empty identified with E),

and call the contracted class `star`. This relation is `C`-invariant,
so complement descends to an involution `J` of `Y`. Its explicit form is

    J(star)=star,
    J({p})={q,r}, J({q})={p,r}, J({r})={p,q},

with the reverse swaps on pairs. Thus

    Y = three singletons + one boundary class + three pairs,  (3)
    cycle type(J)=1^1 2^3.

Choosing `E` as representative of `star` identifies this set with the
seven nonempty subsets. Under this identification `J` is **modified
complement**, fixing `E`; it is not the ordinary complement map. The
choice of the two cardinality-one endpoint states to contract is
canonical under permutations of the three prime labels. Every such
permutation commutes with `J`.

The mechanism generalizes: identifying empty and full subsets of an
`m`-element set turns the complement action on `2^m` points into an
action on `2^m-1` points with one fixed point and `2^(m-1)-1` exchanged
pairs. The distinct `3,1,3` cardinality layers occur specifically at
`m=3`; the general rank layers are the binomial layers with the two
endpoints represented by one distinguished class.

This is a quotient of **sets with an involution**, not of ordered
lattices. If an order-preserving quotient of the Boolean lattice
identifies its bottom and top, every intermediate element must also
be identified with them by antisymmetry. Equivalently, the nonempty
subset poset has one maximum and three minima, so admits no order
reversing bijection. The repaired `J` cannot be such a bijection.

## 3. Actual arithmetic and graph maps into this seven-set

The seven squarefree proper divisors of `N=p^2qr` are exactly

    p,q,r; pqr; pq,pr,qr.

The support map `A -> prod_{ell in A} ell` is a bijection from the
nonempty subsets of `E` to these seven divisors. It preserves inclusion
as divisibility. The singleton/full/pair order in (3) gives their
`3-1-3` decomposition and refines the earlier `7+3` divisor-balance split.

The transported modified reflection swaps each prime with the product
of the other two and fixes `pqr`. It arises from complement relative to
the radical `pqr`, after contracting its empty/full boundary orbit.
It is not complement relative to `p^2qr`: the latter sends `q` to
`p^2r`, outside the squarefree seven. Thus the ten-point divisor
reflection and the seven-point reflection have different domains and
different retained data. The numerical prime labels restore magnitudes;
the full exponents are needed to restore the ten-point lattice.

There is also a literal graph model. Choose a bijection from `E` to
the three edges of a labelled `K3`. Then `B` is the set of all eight
simple graphs on those vertices, and `C` is ordinary graph complement.
The quotient `Y` identifies the empty graph and the triangle. Its
three one-edge graphs and three two-edge paths are exchanged; its
boundary class is fixed. Choosing the triangle representative gives
the nonempty-graph seven, with a modified complement at the triangle.
Edge incidence is retained by the full mask and the chosen edge labels.
The contraction loses the difference between the empty graph and
triangle, so cannot preserve connectivity or cycle existence there.

## 4. The same seven labels are not the same projective incidence

Identify a nonempty subset with its indicator vector in `F2^3`.
These are the seven nonzero vectors, hence the points of the Fano
plane `P^2(F2)`. A line is exactly a triple of distinct nonzero vectors
whose sum is zero, or in subset notation

    {A,B,A symmetric-difference B}.                          (4)

The repaired reflection `J` is not a collineation. The line
`{{p},{q},{p,q}}` is sent to `{{q,r},{p,r},{r}}`, whose vector sum is
`(1,1,1)`, not zero. More precisely, `J` preserves exactly the three
lines through the full-support point. Each such line is
`{E,A,E\A}` and is fixed setwise. For any of the other four lines,
all three vectors are changed by adding `(1,1,1)`, so their new sum
is `(1,1,1)` and the image is not a line. This proves both the positive
scope and the obstruction without identifying structures by counts.

The chosen coordinate permutations do preserve Fano lines and commute
with `J`. They also preserve subset inclusion. The extra Fano lines
come from addition of vectors; multiplication and support in a pointed
monoid alone do not supply that additive incidence.

## 5. The chosen F1 convention: seven Hom points, three rank-zero points

We use exactly the pointed-monoid projective model in Session59,
[Cyclotomic support frames, sections 1--2](../session59_f1_geometry_and_information_20261025/cyclotomic-support/CYCLOTOMIC_SUPPORT_DYNAMICS.md).
It states its primary Thas/Deitmar references and distinguishes relative
base-preserving morphisms from unrestricted coefficient morphisms.

Let `A_m=C_m union {0}`, with `0` absorbing and `C_m` cyclic of order
`m>=1`. For projective dimension two, its point set is

    P_m(3)=(A_m^3 minus {(0,0,0)})/C_m.

An exact support of size `s` has `m^(s-1)` points: normalize one
nonzero coordinate to the identity unit, leaving `s-1` independent
units. Thus its support layers and total are

    singletons: 3,    pairs: 3m,    full support: m^2,
    |P_m(3)|=3+3m+m^2.                                      (5)

At `m=1`, there is exactly one point per nonempty support: seven
literal base-preserving Hom points. The closed/rank-zero points are
the three singleton supports. The ordinary projective counting
polynomial is

    |P^2(F_q)|=1+q+q^2 = 3+3(q-1)+(q-1)^2.

Its formal value at `q=1` is three. The literal `A_1` count substitutes
`m=1`, which corresponds to `q=2` in this polynomial, and gives seven.
These are three explicitly different questions: literal Hom points,
rank-zero points, and a counting-polynomial specialization. There is
no cyclic group `C_0` in the construction, and no ordinary field with
one element is being introduced.

The support bijection to the Boolean seven preserves the specialization
poset only when that poset is retained separately: `S subset T` means
the point with support `S` lies in the closure of the point with support
`T`. Thus larger supports are more generic. This is not an
identification with Fano additive incidence. The reflected set `Y`
does not become a monoid-scheme automorphism or order duality merely
because it has the same seven underlying points: Section2's extremal
count obstruction already rules out the order-duality reading.

## 6. Connection contract and reproduction

The genuine shared operation is boundary-sensitive involution descent:
start with a total complement map, identify an invariant endpoint orbit,
and keep track of the fixed point thereby created. The exact maps above
connect Boolean supports, squarefree divisor labels, and labelled K3
edge masks. They preserve the specified finite-set action. The missing
data are order at the contracted endpoints, graph connectivity/cycles
there, and Fano additive incidence. The historical inverse-fiber seven
has no supplied intertwining map into this action; equation (2) blocks
the immediately suggested coordinate reflection.

Run from the repository root:

    python -B 05-knowledge/results/reflection_bridge_20261025/seven_check.py
    python -B -O 05-knowledge/results/reflection_bridge_20261025/seven_check.py

The checker independently enumerates all Boolean supports at sizes 1
through 8, all seven Fano lines and their images, all six coordinate
permutations, numeric squarefree divisors of 60, and projective unit
frames for `1<=m<=10` by normalization of actual vectors. It checks
the two displayed rational map witnesses and covariance on a finite
grid using exact fractions. Explicit exceptions remain active under
optimization. Normal and optimized stdout agree byte for byte; the saved
receipt is [seven_check.json](seven_check.json). These finite checks
supplement the all-size elementary arguments; they do not re-audit
THM-2473's global claims.
