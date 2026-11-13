# Independent audit of the arithmetic-tree containment obstruction

**PASS: PROVED all-order reduction plus inherited complete census and
independently reproduced finite subset counts.** The underlying isomorphism
census was not regenerated in this audit.

The audited producer is `../containment_edges.py`. Its complete target
universe consists of all odd H<=73, at every tournament order. Busch's
bound for a strong tournament S is `H(S)>=5^((|S|-1)/3)`; the exact integer
inequality `5^8>73^3` excludes strong order9 and above. Thus the inherited
complete through-eight bank contains every possible nontrivial strong
component of every target in this universe.

Some inherited small-core representatives identify converse pairs. This
loses no subset-count values: reversal bijects Hamiltonian paths on each
induced subset, so a tournament and its converse have identical labelled
subset-count tables. Labels within a representative also do not change
the union of values. The order-eight bank supplies all ordinary classes.

For each of148 retained strong representatives the independent checker
enumerates individual simple directed paths recursively. It increments
the cell indexed by the path's visited set, rather than using either of
the producer's forward/backward endpoint dynamic programs. It enumerated
87,388 literal paths, recovered all23,512 subset counts, and checked every
full count and every strongness flag.

The all-order weak-target step is exact. Every tournament is an ordinal
sum of its strong components, and its H is their product. Every induced
subset similarly decomposes into the restrictions to those components,
whose path counts multiply. Singleton components and empty restrictions
contribute1. Conversely any selected component representatives and
subsets can be joined in an ordinal sum to realize the corresponding
product. Therefore arbitrary weak padding adds no unaccounted values.

The producer closes these sets by ascending divisor-pair convolution.
The audit instead enumerates nondecreasing multisets of component counts
and takes their products, including repeated factors. All37 target
unions agree, including the empty unions at7 and21. No count value is
mistaken for a particular core or a chosen root.

The independent binary-heap parent calculation also agrees. With
`m=(h-1)/2`, the ordinary odd-label parent is `2*floor(m/2)+1` for h>1.
Contracting the forbidden nodes7 and21 produces exactly the splice edges
`3->13,15` and `11->41,43`; all other retained edges are unchanged.
Every parent is smaller, so this is a rooted spanning tree of the allowed
odd labels. Branching is ternary at3 and11 after contraction, rather than
binary everywhere.

Exactly these five prescribed edges fail ordinary induced containment in
the checked target universe:

    13 ->27, 25 ->51, 29 ->57, 31 ->63, 35 ->69.

The three prime-source failures13->27,29->57,31->63 are inherited from
[the November2 containment packet](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md)
and independently reproduced here. The full target-union method also
handles composite sources, yielding25->51 and35->69 in this table.
No novelty is claimed for the inherited prime rows or their mechanism.

Thus no pure-birth lift of this arithmetic tree can work, even with an
arbitrary source representative of the required parent count and no
unit-ready constraint. Every sequence of births retains the entire old
tournament as an induced subset. Changing old arcs or deleting vertices
is outside this obstruction; so is changing the arithmetic tree.

Reproduce from the repository root:

    python -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/containment_audit.py
    python -O -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/containment_audit.py

Both receipts agree. They pin the inherited data/proof, the producer and
its output, and the arithmetic-tree code. The numerical agreement is an
independent check of the finite stage, not a second isomorphism census.
