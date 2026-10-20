# Incoming connections: defect templates, cubic junctions and non-leaf bridges

**2026-10-20. PROVED scoped identities and corollaries; FINITE-EXACT word
checks. Actual tournament realization and a universal snark reduction
remain OPEN.**

Two commits arrived during close-out: `602286dc6f`, the
[finite-defect grammar](defect_boundary_20261020/README.md), and
`b57c87b03d`, the
[pure-blue / strong-core packet](merged_metagraph_purity_and_strong_core_20261020/README.md).
The root and the matching/metagraph agents read and audited their relevant
proofs. This addendum records the resulting connections, rather than
treating either packet's status labels as a substitute for that review.

## 1. The same 105 objects now enter tournament enumeration

The incoming [sparse theorem](defect_boundary_20261020/sparse_theory.md),
section 5, classifies the top-degree sector at feedback budget k:
a perfect matching on 2k ordered defect endpoints, a unit weight at
each endpoint, and a heavy coordinate in every gap. Here the cap is k+1.

At k=4 there are eight endpoints and nine heavy gaps. Thus each matching
gives a distinct capped weighted core of order 17 and total weight 53.
Each contributes `x^53/(1-x)^9`, so this sector contributes

    105*x^53/(1-x)^9,

with leading coefficient `105/8!=1/384` in the eventual order-n count.
The incoming theorem proves these weighted cores are distinct: an
isomorphism fixes the heavy vertices in their intrinsic transitive order,
then fixes each intervening light endpoint, recovering the matching.

Fix the endpoint address `0,...,7 -> F2^3`. The matching is exactly the
same object used by our eight-pole kernel B_M and by the H8+I2 chart atlas.
Its pair differences give the same 7+42+56 split. This retains the full
marked matching, not just a coincidence of cardinalities. The endpoint
order cannot be quotiented by arbitrary affine permutations.

Our rank-91 kernel map loses 14 linear directions among these 105
tournament templates. It does not reduce their class count to 91.
All kernel relations have coefficient sum zero, by evaluation at a
constant boundary word. Consequently total coefficient mass survives
this linear quotient, but individual templates need not.

For a word s, write n_c for its number of occurrences of colour c.
The combined response has the particularly simple exact formula

    sum_M B_M(s) = product_c (n_c-1)!!   if all n_c are even,
                  0                    otherwise,

where `(-1)!!=1`. Proof: each colour class must be independently
perfect-matched. The responses on eight ports are 105, 15, 9 or 3
on the parity-allowed multiplicity patterns, and zero otherwise.

Multiplying this response by `x^53/(1-x)^9` gives the common-weight
generating function of the marked leading-template responses. This is
an explicitly defined joint observable, not a theorem that tournament
Hamiltonian-path counts equal cubic colouring counts.

## 2. A failed operation and a concrete cubic repair

The incoming lower-degree corrections shorten a gap, omit a gap, or
share a defect endpoint. The first two preserve the endpoint matching,
so B_M cannot see their changed tournament weights or clone mergers.

Sharing an endpoint instead gives `P3 + 2K2` on seven ports. If every
defect edge is naively translated into equality, the resulting kernel is

    E(s)=[s0=s1=s2] [s3=s4] [s5=s6].

It has 27 words. Exactly one colour has odd multiplicity in every word
of E. But the cubic seven-pole parity lemma requires all three colour
multiplicities odd. Thus `<f_P,E>=0` for EVERY cubic seven-pole P.
This is a type obstruction, not evidence of a universal snark family.

A genuine cubic closure uses a cubic vertex at the three exceptional ports:

    J(s)=[s0,s1,s2 all distinct] [s3=s4] [s5=s6].

This kernel has 54 words and every one satisfies seven-pole parity.
Attach the first three ports to one new vertex and pair the remaining
ports. The closure count is exactly `<f_P,J>`, with multiplicity one for
each permitted boundary word. Loops or parallel edges require their
usual graph checks; a loop-producing pairing automatically contributes
zero in a proper colouring.

This supplies a valid candidate operation on the cubic side. It does
not establish that the repaired junction preserves the tournament
defect count, the template normalization, or an existing obstruction.
The next decisive test would retain an actual seven-pole extension
vector and compare a proposed local replacement in every exterior.

## 3. Strong connectivity supplies a forbidden wiring region

The incoming [cycle-reduction proof](merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md),
section 4, is elementary: a Hamiltonian path of a nonstrong tournament
starts and ends in distinct ordered strong components. Its first
vertex beats its last. Complementing all nonpath arcs reverses that
endpoint arc and closes the retained path to a directed Hamiltonian cycle.

Therefore, for every n>=3, all complement-line neighbors of a nonstrong
node are strong. The nonstrong locus has no internal lines or loops.
This is additional necessary wiring information, beyond degrees,
parity and black twins.

The same packet's
[colour-resolved recurrence](merged_metagraph_purity_and_strong_core_20261020/STRONG_CORE_COUNT_RECURSIONS.md)
turns the forbidden region into exact cut masses. If T_n and B_n count all
tilings and blue tilings, and s_n,b_n their strong subsets, the crossing
line counts are `T_n-s_n` in total and `B_n-b_n` in blue. No factor of two
is lost: each crossing line has exactly one nonstrong endpoint.

## 4. All even-order bridges split into two exact kinds

The incoming pure-blue proof establishes, for even n>=4:

    pure-blue node <=> degree-one node,

and each is attached to a mixed strong node. This uses the inherited
strong-component center law, cycle-cut anti-automorphisms and
[Lu's published crossing theorem](https://www.sciencedirect.com/science/article/pii/S0012365X20300959).
The root independently checked the primary indexed theorem statement;
no fresh full-PDF or independent proof of that cited theorem is claimed.

Combine this with our black-component bridge quotient Q. Remove its
pure-blue leaf vertices and their edges, obtaining Q*. Then

    Bridges(G) = {pure-blue leaf edges}
                 disjoint_union lifts(Bridges(Q*)).

The second part consists precisely of non-leaf bridges of G. Every
remaining quotient vertex represents a black-active component; it can
be a single vertex with black loops. Each such component contains a
strong node, by the preceding endpoint-arc argument. The original
endpoints of a lifted non-leaf bridge are mixed self-converse nodes:
the bridge is blue, but neither endpoint is pure-blue.

Proof of the bridge split: deleting leaf vertices changes no alternative
path between retained vertices. Apply the proved black-contraction
criterion to every retained blue edge. All deleted edges were bridges
by the degree-one condition. This requires no black-core connectivity
conjecture.

Thus the remaining all-order issue is the bridge structure of Q*, not
whether pure-blue even-order nodes are leaves. At odd orders the degree-three
pure-blue example prevents the same leaf deletion rule.

## 5. The old hostile now identifies a useful missing coordinate

Our order-five surgery changes blue edges `(0,11),(2,6)` to
`(2,11),(0,6)`, retaining every fibre size and coloured degree.
The involved classes have cyclic-triangle counts

    c3(0)=3, c3(2)=4, c3(6)=5, c3(11)=0.

It preserves the strong/periphery cut: nodes 0,2,6 are strong and node11
is transitive. That new constraint alone therefore does not reject it.
But the inherited blue triangle-current rule at n=5 requires an odd
jump of absolute value at most three. The original jumps 3 and 1
become 4 and 2. This rejects both altered edges locally.

The incoming [foundation note](merged_metagraph_purity_and_strong_core_20261020/FOUNDATIONS_AND_MIXED_LANDING.md)
routes that rule to
[THM-649, affine flip trilogy](../../01-canon/theorems/THM-649-flip-parity-trilogy-affine-D-binomial-law.md)
and [THM-790, blue parity and leg law](../../01-canon/theorems/THM-790-blue-parity-law-proved.md).
This is an additional necessary screen for complete wiring.
The full presentation certificate in our metagraph note remains
sufficient; the listed scalar screens are not promoted to sufficiency.

## 6. A final incoming proof closes the pure-blue alphabet gap

Commit `fcf72e76dd` supplied the
[uniform prefix-rerouting proof](pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).
The root and metagraph agent checked its symbolic cases: F1/F2/F3 path
increments and violated arcs, the residual sign patterns, the m=4
boundary, and the unit reindexings by plus or minus two at composite
orders. Unlike the earlier affine-only argument, it tests the actual
positional reversal on an explicit arc. The proof passes this audit.

Together with the audited odd-circulant reduction, the nontrivial
strong pure-blue tournaments are exactly the directed triangle and the
regular five-vertex tournament. Call the latter R5 here, to distinguish
it from an undirected pentagon.

The incoming
[branch-count note](pure_blue_atom_classification_20261020/PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md)
also refers to companion rigidity and endpoint-transport documents.
The following deductions give self-contained replacements for those
dependencies using the proved alphabet and the already-read center law.

### A short rigid-strong lemma

If a strong tournament T has `H(T)/|Aut(T)|=1`, its automorphism group
acts transitively on ordered Hamiltonian paths: the action is free,
and its orbit count is that ratio.

Choose a Hamiltonian cycle `v0,...,v(n-1)`. Its two successive cuts
`(v0,...,v(n-1))` and `(v1,...,v(n-1),v0)` must be related by an
automorphism. That automorphism is precisely `vi -> v(i+1)`, so T
is circulant in these coordinates. The reflection
`r(vi)=v(n-1-i)` is therefore an anti-automorphism.

Every other Hamiltonian path is `alpha(P0)` for an automorphism alpha.
Its positional reversal is `alpha r alpha^(-1)`, also an
anti-automorphism. Hence **every rigid strong tournament is pure-blue**.
The new alphabet theorem now leaves only C3 or R5; the latter has
`H/Aut=15/5=3`. Including the singleton, the rigid strong alphabet
is exactly `{1,C3}`. This argument does not use a Hamiltonian-bypass
theorem.

### Normal forms and an all-order coefficient of D_n

By [THM-2454, the proved center law in section 2](../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md),
every pure-blue class has the unique form

    A -> Z -> A^op,

where A is an ordered word in singleton and directed-triangle
components, and Z is empty, singleton, C3 or R5.
Let `a_0=1`, `a_j=0` for j<0, and `a_j=a_(j-1)+a_(j-3)`.
The pure-blue node counting series, including the empty object, is

    P(z)=(1+z+z^3+z^5)/(1-z^2-z^6).

This P is not our degree polynomial D_n. Pure-blue fibre sizes are
exactly one or three. Their respective counts are

    L_(2m)=a_m,             J_(2m)=0,
    L_(2m+1)=a_m+a_(m-1),  J_(2m+1)=a_(m-2).

Conversely any degree-one node is pure-blue: black line multiplicities
are even, and each loop consumes two incidences. Therefore, for n>=3,

    [z]D_n(z)=L_n,

and the entire pure-blue contribution to D_n is
`L_n*z+J_n*z^3`. The other nodes are not asserted to have distinct
degrees; in particular `[z^3]D_n >= J_n` is all that follows here.

### Mixed landing and exact branch geometry

At n>5 every pure-blue tournament is nonstrong, because its possible
strong atoms have orders at most five. Its complement partners are
strong by section 3, and cannot be pure-blue at that order. The same
reason works at n=4. At n=5 the three pure-blue types are TT5,
`1 -> C3 -> 1`, and R5.

The root independently enumerated every Hamiltonian path and every
vertex permutation for these three controls. Their complement
partners are all strong, with `(H,|Aut|)` profiles respectively

    TT5:             (9,1),
    1 -> C3 -> 1:    (13,1),
    R5:              (11,1), (13,1), (15,3).

None is R5, whose profile is `(15,5)`. Thus **every pure-blue neighbor
is mixed for every n>=4**. This avoids an additional endpoint-transport
dependency. At n=3 the two pure-blue leaves are adjacent, as before.

The incoming score-shell argument further proves the three neighbors
of any degree-three pure-blue node are distinct, even after converse
merging. For `A -> R5 -> A^op` with `|A|=r>=1`, its complement partner's
central five vertices are exactly those of score r+2; the other scores
are at most r or at least r+4. The induced central tournament is one
of the three order-five partners with H=11,13,15. Hence no isomorphism,
or converse isomorphism, can identify these three full partners.

Consequently, at every n>=4, there are **exactly L_n leaf bridges**.
All remaining bridges are classified by the black-component quotient
after deleting these leaves. Its remaining pure-blue vertices are
degree-three vertices with three distinct mixed neighbors in the
original metagraph; contraction may identify those neighbors.
Every other quotient vertex is black-active. This sharpens the
all-even statement of section 4 to a structural description at all
orders, but it does not prove that the remaining quotient has no bridges.

The six-node abstract hostile has a bridge joining its two pure-blue
nodes. The mixed-landing theorem now excludes it from the actual
tournament metagraph at every n>=4; it also is not the two-node graph
at n=3. Its abstract involution realization remains valid. Thus degree
parity alone still cannot prove this exclusion, while the retained
normal forms and partner constraints can.

The metagraph agent additionally replayed the incoming independent
prefix checker read-only: all 32,760 normalized odd-circulant
presentations from order 9 through 31, all 591 residual-pattern cases
through m=200, the four order-seven words, and the small positive
controls passed. The symbolic proof supplies the all-order quantifier.

## Reproduction and scope

The [independent script](../../04-computation/boundary_support_20261020_incoming.py)
enumerates all 105 matchings, all 6,561 eight-port words and all 2,187
seven-port words. The [receipt](../../04-computation/boundary_support_20261020_incoming.out)
records the matching-sum identity, the budget-four normalization,
the hostile/repair parity counts, and the all-path order-five
pure-blue partner audit.

    python3 04-computation/boundary_support_20261020_incoming.py
    python3 -O 04-computation/boundary_support_20261020_incoming.py

The incoming packets retain their own universes and reproduction scripts.
No fixed-budget counting theorem is asserted to determine maximum
metagraph degrees, and no cubic-junction construction here forces a
Petersen minor.
