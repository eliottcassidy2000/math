# Nine modular symbols suffice when converse reflection is authenticated

October 23, 2026. **PROVED:** the reflection-marginal forest implication, by
the elementary criterion in `../CONVERSE_MARGINAL_REPAIR.md`. **FINITE-EXACT:**
the complete 6,880-class order-eight modular menu below. This is a separate
information model from the one-direction class-separation theorem in
`README.md`; that theorem's twelve-symbol minimum is unchanged.

## Information retained

Write

    E_pq(T) = { { (c_w mod p,c_vw mod q) : v->w } : v in V(T) }.

Three questions must remain distinct:

1. Does E_pq(T) identify an ordinary class?
2. Does the explicitly supplied ORDERED pair
   (E_pq(T),E_pq(T^op)) identify an ordinary class?
3. Do the separate transported E_pq margins, their counterpart margins
   forced by authenticated converse reflection, and the actual endpoint/
   color data determine the port table?

The second does not imply the third merely by injectivity. A four-cycle of
signature pairs retains distinct joint labels while admitting an alternating
switch that preserves both separate marginal families.

For the third question, form the bipartite multigraph Gamma: its left and
right vertices are signature values, and each ordinary class T contributes
one distinct edge joining E_pq(T) to E_pq(T^op). If Gamma is a forest,
leaf elimination uniquely recovers the mass on every class edge from its
two marginal families. Parallel edges count as a two-cycle and are never
suppressed. Restricting target classes by endpoints/color preserves a forest.
Native reflection supplies the second margins for the same source-category
unknowns, not invented joint bucket counts. The separate root proof gives
all hypotheses and the paired-injectivity hostile.

## Full order-eight results

The independent checker computes triangle-edge incidence by two-step arc
intersections and directly verifies every incoming-row signature against the
supplied converse-class involution. It tests p=1..8, q=1..7, plus exact-head/
parity-edge, parity-head/exact-edge, head-mod3/exact-edge, and exact/exact.

| Signature | Ordinary labels | Ordered paired labels | Gamma cycle rank |
|---|---:|---:|---:|
| E_22 | 6,668 | 6,800 | 81 |
| E_23 | 6,863 | 6,878 | 2 |
| E_32 | 6,865 | 6,874 | 6 |
| E_33 | 6,876 | 6,880 | 0 |
| E_42 | 6,871 | 6,875 | 5 |
| E_43 | 6,880 | 6,880 | 0 |
| exact head, parity edge | 6,874 | 6,876 | 4 |
| parity head, exact edge | 6,863 | 6,878 | 2 |
| head mod3, exact edge | 6,876 | 6,880 | 0 |

The complete grid and paired collision witnesses are in
`converse_receipt.json`. E_33 is the unique coordinatewise-minimal success
in the tested grid for paired class separation and for the global forest
criterion. Its unordered paired labels separate all 3,528 converse-merged
nodes.

### The repair consists of eight tiny stars

E_33 has exactly the same class-collision partition as head-mod3/exact-edge.
Its only collisions are:

| Classes, canonical words | Their converse canonical words |
|---|---|
| 4328706, 4330626 | 25337872, 25727236 |
| 25694272, 25694466 | 8397924, 8393858 |
| 27467792, 27467812 | 2297090, 2303121 |
| 50676048, 50676504 | 4589843, 4399241 |

Each of these eight converse classes has a globally singleton E_33 label.
Thus every target class is identified from its original marginal or its
reflected marginal. Gamma has 13,752 vertices, 6,880 class edges, and 6,872
components: 6,864 single-edge components and eight two-edge stars. Four
collisions occur on each shore, with no edge joining two ambiguous centers.
The repair is an exact reconstruction statement, not failure to find a
counterfeit.

Consequently E_33's nine-symbol alphabet determines the complete ordinary
port table, hence its colored converse-merged multigraph, when the actual
reflection equality is imposed. Without reflection, its ordinary bucket
system still has the inherited 1,351 units of unresolved mass and actual
alternative tables; those are a different model.

## Minimality within the modular family, not arbitrary alphabets

The exact-head/parity-edge and parity-head/exact-edge models admit explicit
nonnegative reflected switches in `../residual_counterfeits.json`, checked
by `../check_residual_counterfeits.py` and the independent rule audit.
Each changes eight ordinary port cells and four merged black entries by
+/-2 while preserving the specified source signature buckets, target
endpoint/color margins, class reflection, blue data, black evenness and
recorded cyclic-triangle currents.

Every E_pq with q<=2 is a coarsening of exact-head/parity-edge or its
constant-edge projection. Every E_pq with p<=2 is similarly coarser than
parity-head/exact-edge or its constant-head projection. The same actual
switches remain feasible after coarsening. Every p*q<9 has p<=2 or q<=2,
so it fails. E_33 attains nine. This proves a nine-symbol minimum for the
specified reflection-constrained wiring model within the positive-integer
product-residue family, not an absolute information-theoretic minimum.

The separate optional quotient-alphabet test deliberately goes beyond this
product-residue family by merging arbitrary E_33 symbols. It must not be
folded into the modular minimality claim.

## All-order and computation boundaries

The forest theorem is all-order, but the specific E_33 forest is finite at
order eight. The strongly connected order-fifteen circulant witnesses in
`README.md` defeat exact E and all its modular quotients; converse does not
separate that pair either. Every vertex has the same incoming and outgoing
cyclic-edge profile, so both shores' signatures agree. No universal class
or wiring claim follows from the small-order forest.

Run from the published synthesis root:

    python -B incidence/converse_check.py --classes wiring/classes_input_n8.tsv
    python -O -B incidence/converse_check.py --classes wiring/classes_input_n8.tsv

Normal and optimized runs match the frozen receipt. The checker also writes
exact `converse_signature_ids.tsv`; no numerical-hash equality is used for
class decisions. No original incidence packet file is changed by this
addendum. The original full-port replay supplies the authenticated feasible
table; this checker tests the complete class-level forest needed to force
that table from the two marginal families.
