# Session38 independent audit

23 October 2026. **PASS / FINITE-EXACT** for the complete order-eight class,
port, colored multigraph and dyadic computations listed here. **PROOF REVIEW
PASS** for the scoped incidence and electrical arguments. No Git files were
modified by this audit. The inherited input is pinned to repository commit
`2130e7633805521d4e24355c37596bd2cd79f044`.

## 1. Independence and exhaustive coverage

This is not a sample and does not use the producer's Hamiltonian-path map
as a classification oracle. The audit independently:

1. Decodes the 6,880 inherited big-endian upper-triangle representatives.
   It minimizes over every ascending-score-block permutation, checks each
   automorphism count, checks distinct canonical words, and verifies the
   orbit-size sum is exactly 268,435,456 = 2^28. There are 403,193 checked
   score-block permutations. Independent subset DP verifies all 6,880
   Hamiltonian-path counts. Converse labels are separately canonicalized.
2. Computes exact E signatures directly from cyclic triples and proves them
   distinct on that independently certified complete class universe.
3. Builds each of the 2,097,152 fixed-path masks directly and classifies it by
   E. This differs from the producer's enumeration and normalization of
   Hamiltonian paths of representatives. Complete class coverage plus
   invariant injectivity justifies every classification, not a hash guess.
   C++ hash-table equality compares the complete exact signatures.
4. Rebuilds every colored ordinary endpoint-port and merged edge entry,
   including loops. The independently produced complete u16 class map is
   byte-for-byte identical to the Hamiltonian-path producer:

   `636f0ae8f69455d4d39143a7c612260cfbb2a9736c91ca67fe2ac3ea65449be1`.

The C++ audit checks signed-root-curl reconstruction on every one of the
2^21 masks. It checks every endpoint transport and every reflection/converse
relation. Reflection is position reversal plus converse, not free-arc
complement. Blue means a fixed actual presentation under that reflection,
not merely a self-converse class label. All black twin multiplicities are
checked even, for directed and undirected supports, including loops.

## 2. Exact wiring results and entrywise comparison

- 6,880 ordinary classes; 176 self-converse classes; 3,528 merged nodes.
- 4,096 blue presentations; 2,048 blue lines.
- 1,046,528 literal black lines; 523,264 paired-black lines.
- 956,339 ordinary endpoint/color/target entries, all compared exactly.
- 98,580 endpoint-capacity entries, all compared exactly.
- 460,699 directed merged endpoint pairs, or 460,816 color-separated cells.
- 438,216 undirected endpoint pairs, or 438,374 color-separated cells.
- 618 literal black loops on 232 supports, and four blue loops.

Every class-field and colored merged multiplicity comparison passes. A further
independent check verifies all four remaining class metadata fields and all
twelve node fields, including cyclic triples, scores, cycle counts from
closed-path orbits, and every integer triangle current.
Independent recursive graph traversal/bridge computation agrees with the
producer's iterative implementation on all five inventories: literal full,
blue, literal black, paired black, and fundamental blue-plus-paired-black.
The full support is connected. Black support has one 3,525-node component
and three isolated nodes; blue support has one 176-node component and 3,352
isolated nodes. The full graph has eight support bridges but only three
literal multigraph bridges. All degrees are checked against actual port
fibers; loops contribute twice to degree. These are finite n8 results,
not all-order connectivity or bridge classifications.

## 3. Incidence claims

A second implementation calculates edge incidences using directed two-step
intersections rather than cyclic-triple enumeration. It independently
verifies the complete representative census and a 13-by-7 modular grid:

- E and E_(4,3) each distinguish all 6,880 classes.
- Dropping the head count gives 6,749 signatures; dropping the edge count
  gives 6,535; unpairing their rootwise marginal lists gives 6,871.
- Global arc triples give 6,832, J4 gives 6,878, and second-round J4 gives
  6,880.
- Head modulus 1, 2 or 3 fails even with the exact edge count. Edge modulus
  1 or 2 fails even with the exact head count. Thus (4,3) is the minimum
  positive-integer modular pair in the coordinatewise sense and its
  twelve-symbol alphabet is minimal **within this modular class-separation
  family**. Moduli are not in general ordered by refinement; the lower-bound
  proof uses exact opposite coordinates, not that false monotonicity.
- Independent exact Q values are (2881,2887) and (3008,3000) on the two J4
  collision pairs, whose complete unassociated edge histograms agree.
  Thus (J4,Q) is another complete finite class invariant.

The all-size Eulerian triangle-balance identity and path-local incidence
transport proof pass review. Independent local transport checks cover every
fixed-path mask through n7 (33,866 masks), beyond the producer's n6 controls.
They also verify the n4 hostile showing source E alone does not provide its
actual-path partner E.

Class-injective transported invariants plus actual endpoint/color capacities
force the complete n8 aggregated table. This uses the actual complement's
signature. Source-class signatures or node margins alone do not supply it.
The mirrored refinement gives a smaller, separately typed result. E_(3,3)
has 6,876 direct signatures, but every one of the 6,880 classes is a singleton
either directly or under converse. Full reflection therefore forces every
ordinary table cell. An independent modular/converse scan and explicit
verification of both mirrored lower-bound witnesses prove that nine is the
minimum alphabet size in the positive-integer modular family for this
reflection-compatible wiring model. This is distinct from the twelve-symbol
minimum for ordinary class separation. The exact-head/edge-parity counterfeit
excludes every q<=2; the head-parity/exact-edge counterfeit excludes every
p<=2, even when the other coordinate is exact. Their changes preserve all
actual source signature capacities and target endpoint capacities, reflection,
and black evenness, but change the merged graph. Every affected native cell
and removed literal mask was checked against this audit's own reconstruction.

The general residual-cycle/potential theorem was also reviewed and is valid
for its stated integer transportation model, with labelled channels and
explicit cell bounds. Peeling exhaustion proves uniqueness; a nonempty
peeling residual alone does not establish an alternative table or graph.

## 4. Explicit all-order boundary

The two recorded 15-vertex circulants are independently verified to have
identical exact E, uniform tags (score,c_v,strong4_v)=(7,28,184), and outgoing
edge incidences (3,3,3,4,4,5,6). Fixed-start recursive completion DP, multiplied
by cyclic transitivity, gives H=196,264,905 and 196,277,235. All induced-subset
completion counts separately recover these H values and give ordered
F2 totals 4,857,720,210 and 4,857,619,470.

The difference in H proves nonisomorphism. Uniform weighted neighbor data
remain indistinguishable at every stated one-dimensional refinement round.
Adding an ordered chain of universal sinks preserves the E collision and
unequal H at every n>=15; explicit controls cover n16, n17 and n18. This is
not a claim that fifteen is the first failure, not a same-H collision, and
not itself a counterfeit metagraph. The root strong-lift continuation was independently reviewed: cyclic blocks
T -> {a} -> P_r -> T give H_lift=(2^r+1)H+(2^r-1)F2 and preserve exact E.
The difference 113070-88410*2^r is nonzero for every r>=1, proving strong
collisions at order15 and all orders n>=17. The H and F2 differences have
opposite signs; the explicit formula, not a same-sign argument, supplies this conclusion.

## 5. Electrical proof and independent numeric audit

The switch-detection proof is sound after its explicit distinction between
comparison on the common unreduced blue kernel and on a common component
quotient. The latter requires the stated preserved-component assumption.
The perfect-pairing proof of rad(q_k)=V_(k+1), including qualifying integral
lifts and lower-bit corrections, is sound. The mismatch-energy theorem is
correctly conditional on actual authenticated target signatures.

For n8 the audit independently builds the blue kernel and first half-form
from its own undirected edge ledger, grounding vertex 0 rather than the
producer's vertex 3527 and using a different elimination convention:

- dim V1=3363, rank q1=3358, dim V2=5; blue half-form rank=4.
- All five mod-four lift vectors, every coordinate of A(lift)/4, and their
  complete q2 matrix pass direct integer-edge verification. Their binary
  classes have rank five, hence span V2. Rank q2=4, so dim V3=1.
- All six further complete lift/image certificates through modulus 256 pass
  direct integer-edge verification. The nonzero binary class is preserved;
  q3,...,q7 vanish, while q8=1.

Consequently the complete result is

    K(G8)[2-primary] = (Z/2)^3358 + (Z/4)^4 + Z/256,
    v2(number of spanning trees) = 3374.

No producer higher-lifting elimination code is imported. The certificates
and the independently established initial dimensions suffice for completeness.
This concerns the actual full literal network, not merely its paired-black
subgraph. The odd-primary part is not computed.

The weighted support-bridge unimodular splitting proof also passes review.
Independently removing the eight verified leaves gives a connected 3,520-node
core with no support bridge and minimum support degree three. The bridge
weights are 1,1,1,2,2,2,2,2. Direct restriction of all eleven mod-four and
higher actual lift certificates preserves the five-dimensional V2 and every
higher form. Therefore the core has

    K(core)[2-primary] = (Z/2)^3353 + (Z/4)^4 + Z/256,
    v2(number of core spanning trees) = 3369.

The high-order factor is genuinely in the bridge-free core, not a collection
of pendant doubled-black edges. This is a finite localization, not a uniform
explanation of high torsion at every order.

The n5 invisible-switch control is also independently checked: every binary
harmonic vector and all half-pairings agree, the form is nondegenerate, and
fraction-free Bareiss determinants give 41,096 versus 36,680. Both 2-primary
groups are (Z/2)^3, so even the complete dyadic tower does not authenticate
arbitrary wiring. The control is synthetic, not another tournament realization.

## 6. Reproduction and publication scope

Requirements: C++17 compiler and Python 3 standard library. No external solver,
network, Python package or producer source import is required. See README.md
for commands. Full regenerated audit ledgers are intentionally not duplicated
in the publication: they are generated by independent_wiring.cpp, compared
entrywise, and pinned by hashes; the primary wiring packet already preserves
the same complete table. FREEZE_RECEIPT.json lists the explicit publication
whitelist and source/output hashes.
