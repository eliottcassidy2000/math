# 33, 11, and 7: values, symmetries, and metagraph fibres

Status: recovered proved mechanisms plus exact refinements. These numbers do
not have one interchangeable meaning. The 33 missing coefficients and the
eleven-vertex roots from the preceding session are separate roles until a
map is supplied; the minimum-root proof in PROOFS.md supplies one such map.

## 1. An actual 33 -> 11 quotient and an actual metagraph line

At order six the H33 classes have Hamilton-frame representatives

| Frame code | Outgoing bit rows | Aut | tc=H/Aut |
|---:|---|---:|---:|
|20|18,12,9,17,38,15|1|33|
|30|58,12,9,16,38,14|1|33|
|12|50,4,9,19,38,14|1|33|
|163|14,20,24,18,33,15|3|11|
|105|38,52,8,19,37,12|3|11|

Converse pairs are (20,12) and (163,105); 30 is self-converse. These are
ordinary classes, not five merged nodes. The 33 ordered Hamiltonian paths
of either symmetric class fall into eleven free three-element orbits.

The literal merged-line census has a black line pair joining merged
representatives20 (H33) and79 (H11) at order six. Its multiplicity is two.
Complementing all nonpath edges, with the path itself retained, is the map;
it does not preserve H. This line relation differs from dividing H by Aut.
`number_line_witnesses.json` records the complete 11/33 endpoint-filtered
line rows through eight, including H11 self-lines at five.

## 2. Recover the existing all-order H33 sequence

The [Session24 multiplicity proof](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
already gives, for n>=6,

    # {ordinary T of order n: H(T)=33}=n^2-5n-1,
    # {converse-merged such T}=n(n-5)/2,
    # {self-converse such T}=1.

Session224 and [the eternal-colour recovery](../eternal_colours_self_converse_values_20261108.md)
also recover this result; no new priority is claimed here. Our independent
small census confirms the ordinary terms5,13,23 at n=6,7,8.

Mechanism: there are five strong H33 cores at order six and three at seven;
H11 has one rigid strong core at five; C3 has H3. The ordered strong-component
grammar is therefore

    F_33(x)=(5x^6+3x^7)/(1-x)^2 + 2x^8/(1-x)^3,
    F_11(x)=x^5/(1-x)^2.

The second H33 term is the two orders of C3 and the H11 core, with three
singleton gaps. Omitting it would falsely turn a quadratic sequence into a
linear one. The finite strong inventory is globally complete by the inherited
Busch bound, pinned in THM-4744,
[`THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).

## 3. A refinement that tracks symmetry and path value together

THM-4812(g),
[`THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md),
classifies the complete strong tc11 fibre: four types, all visible here:

    (order,H,Aut)=(5,11,1), (6,33,3), (6,33,3), (7,99,9).

We use only its audited v<=21 part, not its later v23 extension. The only
strong tc1 blocks are K1 and C3. Since tc multiplies across ordered strong
components and 11 is prime, every tc11 tournament contains exactly one of
these four exceptional blocks, with arbitrary K1/C3 words on either side.
If u records the exponent of three in Aut, the exact bivariate series is

    Q_11(x,u)=x^5(1+u x)^2 / (1-x-u x^3)^2.          (2)

Thus every tc11 tournament has Aut=3^k and H=11*3^k. Setting u=1 gives the
all-order count sequence beginning1,4,8,14,26,46 at n=5,6,7,8,9,10.
This is an elementary refinement of the inherited finite strong-fibre
classification, not a newly classified inventory.

At u^1, (2) counts precisely the H33/Aut3 subfamily:

    2x^6/(1-x)^2 + 2x^8/(1-x)^3.

For n>=6 its count is n^2-11n+32. The remaining rigid H33 family has6n-33
classes. These add to the inherited n^2-5n-1 law. This refinement distinguishes
how a fixed H33 value occupies tc11 and tc33 nodes at every size.

## 4. Seven requires its boundary and its acting group

Inherited THM-4987,
[`THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`](../../../01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md),
excludes7 and21 from H(T)/|A| for every A<=Aut(T). Its inputs are the full-group
gap theorem and the free path action. This does not exclude an automorphism
of order7 or a selected endpoint fibre of size7.

THM-4988,
[`THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md`](../../../01-canon/theorems/THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md),
gives the endpoint-orbit denominator Aut(T)_v and first boundary values7 at
order six and21 at seven. Our unquotiented start/end-count sequence records
0,0,0,0,0,22,142,513 classes with an endpoint count7 at n=1..8. It is a
different statistic from the stabilizer-divided endpoint-orbit count.

Paley7 gives H189 and full group order21, hence tc9; quotienting only by
translations C7 gives27 path orbits. THM-5094,
[`THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md`](../../../01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md),
proves the sharp C7 lower bound H/7>=25, attained by the carousel, whose
H is175. A numeral alone does not specify which quotient is being counted.

Finally, c=7 is an allowed odd ray coefficient: it asks for H=1+7*2^s,
not H=7. `number_roles.json` keeps H, tc, Aut, endpoint entries, and seed
coefficient occurrences in separate columns at each order.
