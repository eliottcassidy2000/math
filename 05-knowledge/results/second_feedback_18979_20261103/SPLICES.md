# Exact connector recurrences and the necessary joint boundary

**PROVED general formula; FINITE-EXACT application.** The complete square-sum
classification being used is the fixed-anchor minimum-connector catalogue
in [SC1--SC2](../return18979_ports_20261103/SQUARES.md), not a classification
of all square-sum graphs. No nonedge is assigned an artificial orientation.

## SP1. Turn the distinct-name constraint into an exact recurrence

Fix an ordering of c old components and an oriented Hamilton path inside
each component. At minimum connector count r=c-1, the inherited splice
bijection says the path must alternate these old blocks and new vertices.
Let C_i be the finite set of legal new names at gap i, excluding all old
names. For square sums, each C_i is the complete factor-pair list

    uv=a-b, u<=v, u=v mod2, z=((v-u)/2)^2-b>0.

For nonempty B of gap indices write m(B)=|intersection_(i in B) C_i|.
The number J(S) of legal distinct connector assignments on slots S obeys

    J(empty)=1,
    J(S)=sum_(B subset S, min(S) in B)
         (-1)^(|B|-1) (|B|-1)! m(B) J(S\B).              (SP1)

This is an exact recurrence on all slot counts. It computes the matching
count of the slot/name incidence matrix without replacing it by endpoint
marginals. For two and three slots it becomes

    J12=m1*m2-m12,
    J123=m1*m2*m3-m12*m3-m13*m2-m23*m1+2*m123.             (SP2)

**Proof.** Sum over set partitions pi of S, weighting a block B by
(-1)^(|B|-1)(|B|-1)! m(B). A chosen function f from slots to their legal
names contributes the sum of these coefficients over partitions refining
its equality partition. On a fibre of size t the coefficient sum equals
the signed sum of permutations of t letters, grouped by their cycles:
each block has (|B|-1)! cyclic orders. It is one for t=1 and zero for
t>1 (pair permutations by a fixed transposition). Thus only injective f
survive. Isolate the block containing min(S) to obtain SP1. This is the
standard partition-lattice inclusion-exclusion, derived here for the
actual attachment menus; no priority claim is made.

Summing J over component orders and internal oriented paths, and dividing
by two for whole-path reversal, gives exactly the number of minimum
square-splice paths. `splice_mobius.py` obtains9 and17 for the original
eight anchors and those anchors plus8. A separate subset matching DP
checks each guard. The old catalogue also establishes that each of these
paths has its own connector set; counting paths and counting connector
sets are not interchangeable in general.

## SP2. Every depth of joint data can be necessary

For r>=2, take the even-weight and odd-weight binary words of length r as
two separate connector universes. Let C_i consist of words whose ith bit
is one. Every intersection of a proper nonempty set of slots has size
2^(r-|B|-1), in both universes. The full intersection has size one in
exactly one universe, according to the parity of r.

SP1 therefore gives a difference of absolute size (r-1)! in the numbers
of injective assignments: only the one-block term changes. Hence **all
proper intersections together do not determine depth-r composition**.
At r=3 the two counts are2 and4, despite identical single-slot sizes and
pair intersections. Exact independent controls cover r=2..5; the proof
has no finite-r cutoff.

## The useful transfer to tournaments

The source is an endpoint-resolved ordered path-cover skeleton. The map
replaces each old block by its actual path and each connector slot by a
distinct named vertex. It preserves legal concatenation and vertex use.
It forgets unused arcs, internal automorphisms, and alternative skeletons;
those must remain as sidecars when the target is a whole tournament or
an isomorphism-class metagraph rather than this selected path family.

For a tournament, deleting added vertices from a Hamilton path again
produces an ordered old path cover. But SCC blocks need not remain whole
after attachments, and other covers contribute. Therefore SP1 is an
exact ingredient on a declared skeleton, not an equality with the total
tournament H count. The inherited
[THM-3121 path-cover substitution kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)
supplies the wider sum. The actual81/27 control shows why equal H and
equal labelled cyclic triples cannot replace the underlying source.

The connection to IC1 is organizational rather than numerical: both
calculations aggregate structured partial assignments, retaining the
weights needed for later composition. IC1 uses permutation-cycle
placement integers; SP1 uses signed cycle weights to enforce distinct
names. Their targets and symmetries differ. Neither scalar formula
reconstructs the full tournament metagraph, and neither proves a Boolean
inverse for every allowed Hamilton count.
