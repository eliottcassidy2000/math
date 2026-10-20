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

## Reproduction and scope

The [independent script](../../04-computation/boundary_support_20261020_incoming.py)
enumerates all 105 matchings, all 6,561 eight-port words and all 2,187
seven-port words. The [receipt](../../04-computation/boundary_support_20261020_incoming.out)
records the matching-sum identity, the budget-four normalization,
and the hostile/repair parity counts.

    python3 04-computation/boundary_support_20261020_incoming.py
    python3 -O 04-computation/boundary_support_20261020_incoming.py

The incoming packets retain their own universes and reproduction scripts.
No fixed-budget counting theorem is asserted to determine maximum
metagraph degrees, and no cubic-junction construction here forces a
Petersen minor.
