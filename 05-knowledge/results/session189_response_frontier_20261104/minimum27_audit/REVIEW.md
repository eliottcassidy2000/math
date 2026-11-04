# Independent audit of the 27-name original-module boundary

**PASS for the stated finite C3-child theorem, the root-containing twenty-query family battery, and the source-relative lower bound of 28 old vertices.** All independent physical recounts, metadata checks and normal/optimized verifiers pass. The reviewed producer PROOF.md SHA-256 is `7e0612183fd4c546d94874a7bd3d042bd1800ef68c1137f1c3934e346ecc08a0`. Its pending-review status line is superseded by this separately frozen approval. The producer source and mathematical data are frozen by `PRODUCER_PINS.json`.

## Scope

The finite claim concerns a prime nine-name quotient Q with exactly nine directed-triangle children on fixed physical names. A coarse query is constant on each actual child. The source Q[C3,...,C3] and quotient-converse target use the same child sets and ports x->y. It is not an arbitrary-child or fixed-ambient order-nine theorem.

The global consequence, after combining the separately reviewed highest-prime and nontransitive-child restrictions, is that a nonidentity simultaneous XOR of modules of an ORIGINAL source with the same full named two-port response needs at least 28 old vertices. This is a source-relative lower bound, not an exhibited collision, unrestricted reconstruction theorem, comparison of arbitrary family outputs, or permission to use later-created modules.

## 1. Independent census domain

The separate `c3_census_audit` packet uses a degree-bucket canonicalizer, not the producer's color-refinement canonicalizer. It enumerates all degree-respecting permutations of each supplied representative and its converse, counts automorphisms as the multiplicity of the minimum, checks ordinary canonical codes are pairwise distinct, and sums their exact orbit sizes. At order nine:

- 97,144 converse-merged representatives expand to 191,536 ordinary classes;
- their disjoint labeled-orbit mass is exactly 68,719,476,736 = 2^36;
- 33,867,002 degree-respecting permutations are checked;
- 69,445 merged representatives are prime, equivalent to 136,655 ordinary prime classes;
- 67,097 prime merged classes have nonzero endpoint parity and 2,348 have zero endpoint parity.

Full labeled-orbit mass of pairwise disjoint isomorphism classes proves completeness without trusting a published count or the incoming generator. The independent pair-seeded distinguisher-closure module test is valid: every module containing the seed must contain every added distinguisher, and a fixed proper closure itself is a module. Literal directed-path DFS computes every prime representative's exact start/end arrays. All 2,348 residual physical row arrays, source codes, and 42,264 residual endpoint entries match the producer exactly. This audit additionally cross-checks the producer filter against every independent prime row: all 1,250,010 start/end coordinates and all 67,097 expanded physical endpoint-witness records agree. This audit read and checked the separate implementation; its final receipt and review are pinned dependencies.

## 2. Producer residue logic

For C3 the complete ordered spanning-path-cover profile is (q1,q2,q3)=(3,6,6). Splitting a physical Hamilton path into maximal child visits yields a directed word on the nine child names and the two singleton ports, with no consecutive equal child names. If child i occurs r times, its ordered r-cover contributes q_r. The word and ordered internal covers reconstruct the path uniquely.

Divide each child profile by 3 to obtain (1,2,2). Each child that appears twice or three times therefore contributes a factor 2, so any term with at least three repeated children vanishes modulo eight. The producer's state truncation is exactly this rule. Its visit multipliers are 1,2,1, producing normalized cumulative weights 1,2,2, and its final multiplication by 3^9 restores the physical count. No division by a run factorial, count of distinguishable repeat copies, or omission of a three-visit term occurs.

The source and target quotient rows differ by exact edge converse. Both port rows retain x->y, with a_i=[child i -> x] and c_i=[child i -> y]. The comparison is on the same physical names, not independently relabeled response polynomials.

## 3. Independent physical-subset recurrence

`orbit_physical.cpp` does not count occupied-run words, normalize child profiles, or truncate repeated children. It starts from the ordinary physical-vertex Hamilton endpoint recurrence and uses only the independent cyclic-rotation automorphisms of the nine C3 cells.

A physical subset in one C3 is classified by its size d=0,1,2,3. For d=2 choose its canonical directed edge 0->1. Let T_i(s) be the total Hamilton endpoint count in child i for the canonical physical subset represented by state s. Let B_i(s), stored when d_i=2, be the count ending at vertex 1 of that directed edge. If E is the external-predecessor sum for the state with the last vertex removed, then:

- d_i=1: T_i=E;
- d_i=2: T_i=2E+T_i(previous), B_i=E+T_i(previous);
- d_i=3: T_i=3(E+B_i(previous)).

For d_i=3 the three endpoint counts agree by cyclic rotation. Removing the canonical last vertex leaves the directed two-vertex edge, and only its terminal vertex points to the removed vertex. A singleton port uses the ordinary predecessor sum. The singleton-state base count is one. All predecessor terms are evaluated on smaller states.

There are 4^9 times 4 = 1,048,576 physical subset orbits per side, with both port-presence bits retained. Summing the full-state endpoint counts gives the actual 29-vertex augmented Hamilton count modulo eight. Every state, including arbitrarily many interrupted child visits, is present. This is a materially different exact counting derivation from the producer's normalized run-word recurrence.

The complete independent recount retains one certified separator for each of the 2,348 residual representatives. Every source and target residue agrees with the producer, and every source-minus-target residue is 4 modulo eight. The audit does not claim to independently recount all entries of the producer's larger twelve-query bank: the independently verified 2,348 witnesses suffice for the exclusion.

The new recurrence is additionally checked against literal uncompressed physical-vertex subset counting on 580 query pairs: every quotient and every coarse query for one, two and three C3 children (548 pairs), plus 32 four-child cases. These 1,160 exact physical graphs have order at most 14. Independent child converses and common arbitrary physical permutations are included. Every exact integer count reduces to the independent orbit residue. Normal and optimized Python modes preserve all checks and yield identical receipts.

## 4. Battery, endpoint branch and physical transport

The twelve residual coarse queries are:

    (0,1), (1,0), (1,2), (1,3), (3,1), (3,3),
    (3,5), (0,7), (0,18), (40,256), (4,3), (130,128).

All masks name nine quotient positions. Child i uses physical names 3i,3i+1,3i+2; x=27 and y=28. The lift of mask m is the sum of 7*2^(3i) over its set bits. `physical_records.json` retains every actual 29-row source and target, both physical port masks, the shared child sets, and both independently recounted residues. The verifier checks every edge, both masks, the fixed port edge, and that exactly the between-child edges are conversed.

For the 67,097 nonzero-endpoint classes, the separately reviewed all-order endpoint theorem in the pinned `prime8/PROOF.md`, Section 2, applies at arbitrary quotient order: if p_i=S_i+E_i is odd, the actual query is (e_i,0), and its gap is 2 modulo four. It is not (0,e_i). No prime-eight census is extrapolated here; only the explicitly general theorem is used.

The union of all nine endpoint anchors and the twelve residual queries has twenty distinct members, since (1,0) occurs in both. Thus twenty fixed coarse queries suffice in the retained representative coordinates. An arbitrary actual quotient is handled by one common isomorphism, possibly after quotient converse, transporting each entire child set and both masks together. Converse merely swaps the compared source and target residues. This does not assert that the same twenty raw masks work without transport under every arbitrary labeling.

Within a C3 child, any orientation change preserves the complete profile (3,6,6). A uniform query therefore has exactly the same count for any such target child interiors by occupied-run substitution. Equivalently, an internal isomorphism fixes the coarse query because all three physical names have the same port bits. This is equality of counts at the same physical query, not equality after separately renaming a nonuniform response.

## 5. All original-module families and the order bound

Assume T and a different original-module XOR target U have equal complete named responses. The pinned, separately reviewed highest-prime integration says every highest changed quotient in T's ORIGINAL strong-module tree is prime of order at least nine. The pinned crossed-TT4 review, together with its reviewed common-end and distance dependencies, gives the branch-complete restriction that every child of such a node is nontransitive. Its old conditional prose is superseded by the explicitly completed independent review; there is no remaining TT4 exception.

Choose a highest changed node M with m children. Every nontransitive tournament has at least three vertices, so |T| >= |M| >= 3m >= 27. If |T| <= 27, equality is forced throughout. Hence |T|=|M|=27, m=9, M is the whole source, and all nine children are directed triangles. The quotient is prime and its changed target is its converse. The child target profiles agree by hereditary equality; alternatively, proper original modules of this prime substitution are confined to a triangle, where the only nonsingleton module is the whole triangle. Therefore every child target is one of the two orientations of C3, regardless of which original-module family produced it.

The finite coarse separator contradicts full response equality. The argument includes arbitrary simultaneous original-module families, not merely a single root reversal. Families that do not change a nine-child highest quotient were already excluded by the reviewed structural restrictions.

There is no exterior ambient at the equality boundary and no second highest changed node: the selected M occupies all 27 names. Consequently this new argument needs neither an unjustified ambient lift of a scalar modulo-eight gap nor a cancellation assumption among multiple local changes. Same-query restoration is part of the inherited lower-quotient and transitive-child restrictions; those previously reviewed results are the exact inputs, not replaced by finite C3 experiments.

The twenty-query coarse battery is family-uniform only over root-containing families. It does not separate a family that merely converses one or more C3 interiors while retaining Q: all complete child profiles then agree, so every coarse query agrees. Such child-only changes are excluded from a full-response collision by the inherited lower highest-node theorem, not by this coarse battery.

The conclusion is a lower bound of 28, not a claim that order 28 admits a collision or that all order-nine quotients with arbitrary nontransitive children are separated.


## 6. Reproduction and final approval boundary

`verify.py --recount` compiles the independent physical counter and reproduces every large witness. Normal and optimized Python produce identical audit receipts, using explicit exceptions for every check. The final verifier cross-checks all original physical source/target rows and masks against the producer, all prime endpoint arrays against the independent census, every endpoint anchor record, every independently recounted chosen bank entry, and range/shape/first-witness consistency for all 56,352 producer bank entries. It does not mislabel those latter consistency checks as independent recounts of the entire bank.

`PUBLICATION_FILES.json` is the exact release allowlist. Its checksums exclude binaries and runtime logs. A separate address/undefined-behavior sanitizer run of four large physical cases also agreed; leak scanning is unavailable under the execution environment's process tracing and was not used as evidence.

The only mathematical extension beyond the displayed quotient-only comparison is exact complete-profile substitution for the 512 cyclic child-orientation vectors. The global size bound then uses the pinned, independently reviewed structural theorems exactly as stated. No new arbitrary-ambient lifting statement, unrestricted child theorem, later-created-module action, or unrestricted named-response injectivity claim is approved.
