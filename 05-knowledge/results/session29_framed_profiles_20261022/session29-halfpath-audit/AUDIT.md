# Independent audit: mirror half-path algorithm and census

Status: ACCEPTED, 2026-10-22. No correction required. This audit modifies no repository or producer files.

## Proof review

1. **Anti-involution existence and conjugacy.** A tournament automorphism cannot exchange a pair in an involution, so its automorphism group has odd order. An anti-automorphism has order twice an odd number, and the corresponding odd power is an involution. It fixes at most one vertex. If a,b are anti-involutions and g=ab has odd order d, then g^((d-1)/2) a g^(-(d-1)/2)=b. This verifies the stated conjugacy without assuming the automorphism group is abelian. Order one is handled separately.
2. **Half-path and fixed-cycle bijections.** The unrecorded half is forced by reversal and the anti-involution. The only additional condition is the last fold arc, to the unique fixed vertex in odd order or to the partner in even order. A fixed directed odd cycle has the unique fixed vertex as its canonical anchor. The first half chooses at most one vertex from each exchanged pair. There is no double counting by rotation after anchoring. Fixed cycles and fixed full paths remain distinct.
3. **Dynamic programming.** Deleting the final chosen vertex gives the stated recurrence uniquely. There are 2m3^(m-1) valid state/endpoint pairs, each with at most m-1 predecessors. The O(m²3^m) arithmetic and O(m3^m) integer storage claims are correct; the report explicitly distinguishes bit complexity and its rectangular implementation. The anti-involution is required input, not obtained for free.
4. **Blue tile quotient.** Aut(T) acts freely on ordered Hamiltonian paths. A path determines at most one reversing anti-involution. A mirror-fixed tile orbit thus contributes exactly |Aut(T)| incidences (P,tau). Counting instead by the conjugate anti-involutions yields b=f/|C_Aut(T)(sigma)|. This is not f/|Aut(T)|. The parity and pure-blue criterion follow correctly.
5. **Bound and expectation.** The fold condition permits exactly one sign of the last pair, giving the stated bound. In the fixed-sigma labelled ensemble each signed half-ordering tests m distinct independent edge-orbit bits. Thus its probability is 2^(-m), and the exact ensemble mean is m!, not an unlabelled-class average.
6. **Signed-pair census.** Conjugacy of anti-involutions adjusts any ordinary isomorphism to commute with sigma. The acting signed-pair group has size 2^m m!, with stabilizer C_Aut(T)(sigma). Orbit-stabilizer and the labelled mean prove both weighted checksums. They are not formulas for the unweighted number of self-converse classes.

## Independent exact computation

`independent.py` imports no producer code. It uses a vertex-bitmask recursion with partner exclusion, rather than the producer's ternary array. Separate directed DFS enumerates full Hamiltonian paths and all cycles. Full vertex permutations compute automorphisms, centralizers and direct spine-normalized tile sets. Signed-pair canonical keys are independently compared with the full spine class keys.

All 599 fixed-reversal tournaments of orders one through six pass. Their labelled means, per-class orbit sizes, signed-pair class partitions and both weighted census checksums pass. Recovered self-converse class counts are 1,1,2,2,8,12. Normal and optimized outputs are byte-identical. A fresh producer replay also matches its frozen results exactly.

Named-case tuples (H,Aut,centralizer,f,blue,total tiles) are independently reproduced:

- C3: (3,3,1,1,1,1).
- C5{1,2}: (15,5,1,3,3,3).
- P7: (189,21,3,9,3,9).
- Seed A: (15,3,1,1,1,5).
- Seed B: (13,1,1,3,3,13).

For P7 the independent fixed-cycle count is 12, whereas fixed full paths number 9. In particular P7 is mixed, not pure blue. No primitive-exception characterization was assumed or audited.

## Scope

The general claims follow from the reviewed bijections and group actions, not extrapolation from the finite checks. The result provides faster symmetry-refined counts and exact census checks, not full H(T), an unrestricted A000568 algorithm, or complete metagraph classification. Historical novelty claims are appropriately limited to the stated extension; the prior Paley half-path idea is credited.
