# Audit scope and result

Final mathematical review: PASS, 2026-10-20 15:44 UTC. No outstanding corrections were identified in the scoped packet.

## Proofs reviewed

* All-order mixed landing of the unique transitive complement line, including its explicit normalized black mask and the n3 exception.
* Cycle-cut anti-automorphisms and their rotation2 product; odd-order regular circulant reduction; exclusion of orders divisible4; Lu-based exclusion of the remaining even strong orders.
* Even purity iff self-converse rigidity, with the complement-line leaf consequence restricted to n>=4.
* Injection of directed Hamiltonian cycles, counted up to cyclic rotation, into n-cycle automorphisms.
* Nonstrong-periphery independence via the endpoint/apex Hamiltonian-cycle closer.
* First-SCC recurrence, mirrored-outer-component blue recurrence, and all entries of the colour-resolved strong-core cut table.
* The explicit count-preserving n5 rewiring obstruction and the correction that nonrigid pure-blue nodes need not be leaves.

## Computation scope

The portable programs have no repository imports and use only Python's standard library. Exact saved-output replay passes:

1. n5: all64 tilings, all120 vertex permutations in canonicalization, both converse classes, complete blue-line list and hostile two-edge surgery.
2. n3..100: direct adjacency/path/mask/reflection/score checks for the transitive-pipe family.
3. n3..7: every literal fixed-path tiling, independent directed reachability from/to vertex0, literal coordinate reflection, and complement-pair classification. All six strong/core/cut count columns match their recurrences. No periphery--periphery line occurs.

A separate independently implemented finite audit checked pure-blue fibre masses, strongness, successor automorphisms, and cycle counts through n6. A separate sorted-score equality criterion applied to the census verified the strong/core/cut table through n8. These corroborating audits do not extend the portable scripts' claimed direct enumeration range; n8..12 output in the portable count script is recurrence-derived.

## External and inherited dependencies

Lu's Theorem7 was checked in primary indexed HTML, including its precise exception and the theorem statement; source access is documented in LU_CROSSING_SOURCE.md. Full-PDF access is not claimed. Camion's classical Hamiltonian-cycle theorem is cited input.

The first-SCC recurrence is already implied by THM-349 together with THM-354; the apex-closing mechanism is inherited THM-333. The through62 numerical pure-blue consequence is conditional on inherited THM-2477 and THM-2453, whose lengthy proofs were not reaudited here. It is not part of the newly audited all-order structural proof.

Emerging asymptotic-finiteness and unrestricted rigid-atom classification arguments are outside this packet and this approval. No claim that they have been proved is inferred from the present audit.
