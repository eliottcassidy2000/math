# Even the complete path and cycle length spectra do not determine the metagraph neighborhood

Take the order-seven tournaments with the following zero-based
outgoing row masks:

    A: 34,68,33,87,71,90,5
    B: 70,64,98,39,47,3,56

Both are rigid, neither is self-converse, and they are not
isomorphic even after identifying converse. The score multisets
are (2,2,2,2,4,4,5) and (1,2,3,3,3,4,5), and the converse of
the former has (1,2,2,4,4,4,4). All three are different.

The complete counts of ordered directed simple paths by vertex
length 1 through 7 agree: (7,21,53,107,161,151,65).
The counts of unrooted directed simple cycles by length 3 through 7
also agree: (9,15,16,11,3). This is more agreement than the two
Hamiltonian counts and stabilizer size alone.

However the neighborhoods in the multiplicity-retaining merged
metagraph are different even if only the Hamilton counts of mates
are observed. Both have 65 line incidences, including 21 closed
incidences, but the closed incidences whose open mate has no
Hamilton cycle are **three for A and six for B**. The lists of
mate path counts on these nonstrong incidences are respectively

    A: 33,37,37
    B: 23,25,33,33,37,37

For the calculation, for every Hamilton path p flip its path arcs
and read the resulting tournament. This is the converse of the
nonpath-arc complement in [THM-4729](
../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md).
Path and cycle tallies, and strong connectivity, are unchanged by
converse, so it gives the same observable mate. A p is closed exactly
when its final vertex beats its first vertex in the original graph.
The absence of automorphisms and anti-automorphisms on either
given node means no line incidences need to be divided by a further
symmetry factor at that node. Distinct lines can nevertheless have
the same unmarked mate; the claim is about the metagraph with
line multiplicities, not any simple graph obtained by removing them.

The [strong-bridge fan formula](CYCLE_NECKLACES.md) describes the
same boundary from the cycle-reversed side: it depends on bridges
in each actual reversed-cycle source. The complete marginal
length spectra do not retain that information. Thus an exact
recurrence on those spectra alone cannot return the correct
labeled-by-(H,C) neighborhood for every tournament class.

Run `python metagraph_collision.py` for literal checks using
vertex permutations only. No minimal-order claim is made.
