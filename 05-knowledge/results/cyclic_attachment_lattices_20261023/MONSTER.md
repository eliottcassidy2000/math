# Four Monster/Baby Monster views and the shared-quotient research move

**CITED** group/modular inputs; **PROVED** elementary consequences and the
abstract transfers stated below. No new moonshine theorem or direct map
from the Monster to Collatz is claimed. This note recovers audited repo
work before choosing the analogy.

## 1. Prime support, valuation, action, and level are different resolutions

The fifteen prime divisors of the Monster order are

    2,3,5,7,11,13,17,19,23,29,31,41,47,59,71.

The [ATLAS Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) and
[ATLAS Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/) give the
orders used in the audited [level-filtration note](../monster_levels_square_sum_7_21_20261022.md)
and [Session32's centralizer note](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md).
Their three strata are:

| Primes | Group-order view | Action on the Monster's2A involutions |
|---|---|---|
|17,19,23,31,47|same exponent in M and B|every Sylow p-subgroup fixes some point|
|2,3,5,7,11,13|larger exponent in M|some order-p subgroup fixes a point, no full Sylow does|
|29,41,59,71|absent from B|no order-p subgroup fixes a point|

The action statement is the inherited [Sylow fixed-point lemma](../session32_square_monster_20261023/SYLOW_FIXED_POINT_PARTITION.md):
for G/H, a subgroup fixes a point iff it lies in a conjugate of H. Here
H=C_M(t)=2.B. The subgroup index is |M|/(2|B|), not |M|/|B|. Its deficits
at the six small primes are(4,7,3,4,1,2). This preserves MISTAKE-668's
factor-two repair; automatic nesting of rationality lists is not claimed
as a new group-theoretic theorem.

The external geometric input is [Aricheta, Theorems1.1 and1.3](https://arxiv.org/html/1809.07421):
away from the level, supersingular nonrational points measure a genus
difference, and rationality on a monstrous curve is equivalent to its
associated centralizer containing a prime-order Fricke element. The repo's
audited1|2|6 filtration recovers the prime supports of M, B, and Fi22.
Those are existing results, not a new classification here.

The especially useful retained example is the repo's independently checked
level-two computation at29,41,59,71: rational j-vertices acquire two
parallel coarse2-isogeny edges exchanged by Frobenius. Forgetting the edge
lift erases this obstruction even though both endpoints remain rational.
This is a multigraph with Frobenius action, not a tournament to be oriented
artificially. Its role here is a hostile to endpoint-only reconstruction.

## 2. The central sign is visible in an actual representation

An additional primary-source view comes from Breuer, Magaard and Wilson's
[character-table verification, section2, p.3](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf).
The196883-dimensional Monster character restricts to2.B with constituent
degrees1,4371,96255,96256. The central involution acts as+1 on the first
three and-1 on the last. Thus its two eigenspaces have dimensions100627
and96256, and its trace is4371.

The elementary deduction is important for the present task. Passing from
2.B to B is legitimate for an action only when the central involution acts
trivially. The96256-dimensional sector fails that test. Replacing it by a
B representation would erase a measured central charge, not just relabel
the same object. The projectors(1+t)/2 and(1-t)/2 explicitly retain the
two characters of the kernel C2.

This is the same OPERATION used in CHARGE.md: separate the characters of
a forgotten kernel before adding contributions. It is not an isomorphism
between a coloured Fibonacci tree and a Monster representation. The source
groups, dimensions, and target predicates differ. The transfer is a proof
move: test whether the kernel acts trivially on the actual observable before
descending it to a quotient. The quaternion and ternary central-extension
controls are cheap examples of failure of that test.

## 3. A Monster local group really is glued over one common quotient

The same [primary supplement, sections4.13 and5](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf)
describes the7A cyclic-subgroup normalizer as the subdirect product of
7:6 and He.2 over C2. This is the source of the existing
[Paley7 bridge](../session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md).

Write epsilon for the square/nonsquare multiplier sign of AGL(1,7)=7:6,
and delta for the outer coset map He.2->C2. The subdirect product is

    {(a,h):epsilon(a)=delta(h)}.                        (1)

The exact shared-coordinate requirement is the displayed equality.
Each projection is surjective, but independently taking both projections
admits the two mismatched sign combinations. The compatible product has
half the size of the unrestricted product; each pair's admissibility
depends on its actual common sign, not just either projection's order.

Map to this session: replace(7:6,He.2,C2) by(P_i,P_j,P_i/K_i) in the
block-stabilizer theorem. Preserved predicate: one pair has a common
quotient image. Destroyed data: its independent kernel lifts. Required
sidecar: the quotient value and the identification of its two appearances.
The cyclic lattice record is an arithmetic implementation of this fibre
product, including the independent lifts. This is an exact common logical
structure, not a numerical relationship between the orders of the groups.

The He kernel is not the Baby Monster, and this7A statement is not a new
description of2.B. It is a different local view inside M. It provides a
particularly clear model of why two plausible local choices cannot be
combined independently at a shared port.

## 4. Connections that survive the type checks

| Source -> target | Map and preserved predicate | Information that must remain |
|---|---|---|
|actual block pair -> compact cyclic record|translation stabilizer coset -> two congruences; preserves all chart comparisons|same cell variable on every edge|
|level-two moduli -> j-vertices|forget isogeny structure; preserves endpoints, not Frobenius-fixed edges|which parallel edge/lift|
|2.B representation -> potential B action|test triviality of central C2; descent iff its action is trivial|central-sign sectors|
|7A local group -> two quotient actions|common-C2 fibre product; preserves compatible pairs|one shared sign|
|zero-slot tree -> integer length|first-register projection; preserves length, loses z-charge|split/recolour debt and occurrence owners|

The all-size positive result here is the compact cyclic compatibility and
restricted-character certificate. The Monster themes guide WHICH quotient
loss to test; they do not certify any omitted tournament block, establish
a uniform source-specific Collatz bound, or turn the fifteen-prime list
into a universal cutoff principle.

Primary sources were reopened this session. The larger supersingular field
enumerations and Monster character-table computations are credited as
existing input; they were not rerun here. Only the explicit small-group,
dimension-arithmetic, lattice, and character controls in this packet were
executed anew.
