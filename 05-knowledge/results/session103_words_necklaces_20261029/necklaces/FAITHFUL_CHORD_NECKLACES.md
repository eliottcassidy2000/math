# A faithful necklace model for unique-Hamilton-cycle tournaments

**PROVED:** the encoding and isomorphism/converse statements below. **FINITE-EXACT:** complete fixed-cycle enumeration through order7 and the displayed insertion/repair examples. This is a newly specified model, distinct from the recovered five-letter insertion alphabet. No literature-first claim is made for encoding a cyclically marked graph by its chords.

## 1. Use the intrinsic cycle as the coordinate frame

Let T have a directed Hamilton cycle C, and write its vertices in cycle order as Z/n, n≥3. For each i define a binary row-letter

    W_i(d)=A_(i,i+d),   d=1,...,n−1,

where A is the directed adjacency matrix and indices are modulo n. The cyclic word of row-letters has length n; each letter has n−1 bits. It satisfies

    W_i(1)=1, W_i(n−1)=0,
    W_i(d)+W_(i+d)(n−d)=1.

Conversely these conditions specify a unique tournament with the distinguished directed cycle0→1→...→n−1→0. The first pair records the cycle; the reciprocity equation says that every unordered vertex pair has exactly one direction. Thus the model is **lossless**. The alphabet grows with n and its letters have overlap constraints; it is not a free fixed-alphabet language.

**Necklace theorem.** Isomorphism classes of pairs(T,C) are exactly rotation necklaces of these constrained words. If T has exactly one Hamilton cycle, the same necklaces classify the ordinary unmarked tournaments in this family.

Proof: any isomorphism of directed cycle-marked pairs sends i to i+r for one r, hence rotates the row-letters. Conversely every such rotation is a pair isomorphism. When the Hamilton cycle is unique, an ordinary tournament isomorphism must preserve it, so the mark imposes no extra restriction. This is an exact instance where an intrinsic unique structure supplies a coordinate frame without additional marking data. No automorphism-triviality assumption is needed.

The ordinary group of a unique-cycle tournament is consequently its necklace rotation stabilizer. At orders4–7 the computed necklaces are aperiodic; C3 is the small periodic exception. We make no all-order rigidity claim here.

## 2. The uniqueness guard is essential

The cyclic5-tournament R5, with positive differences1 and2, has two Hamilton cycles:

    (0,1,2,3,4),   (0,2,4,1,3).

Their constant row-letters are respectively1100 and1010. These are distinct necklaces of marked objects but describe the same unmarked tournament. Thus arbitrary cycle-marked necklaces cannot be counted as ordinary tournament classes when uniqueness is dropped.

The exact checker verifies both cycles and the distinct necklace codes. This is a concrete marking-loss witness, not merely a caution about possible automorphisms.

## 3. Converse is a twisted reflection, not bare word reversal

To express T^op while restoring the forward cycle convention, relabel i↦−i. Its word is

    W_i^dagger(d)=1−W_(-i)(n−d).

Thus converse combines reversal of vertex order with reverse-complement of each row-letter. It preserves the reciprocity conditions and squares to the identity. Rotation together with this twisted reflection gives the appropriate dihedral action on the constrained words.

A plain reversal of the row-letter sequence can fail to encode any tournament. For the first parent below the row word is

    (1000,1000,1010,1110,1110).

Reversing only vertex indices violates reciprocity at some pairs. The checker verifies this directly. Therefore ordinary free-word bracelet formulas must not be applied without specifying the internal letter action.

For the native ordered two-run observer M2, reversal of both directed paths gives

    M2(T^op)=M2(T)^t

in the same labelled ordered-pair coordinates. Relabelling then simultaneously permutes rows and columns. Consequently determinant, rational rank and integral Smith data are converse-invariant. This symmetry is real; it does not prove that ordinary reversal of an unspecified insertion history induces converse.

If c_n is the number of unique-cycle necklaces and s_n the self-converse subset, the quotient up to converse has(c_n+s_n)/2 classes. This is Burnside for the proved twisted involution, not a guessed word quotient.

## 4. Complete small-family controls

We fix the directed cycle, enumerate all2^(binom(n,2)−n) chord orientations, test Hamilton-cycle uniqueness exactly, and quotient only by cyclic rotations. Every unique-cycle tournament has a representative in this set, by the theorem above. The resulting census is complete in its stated orders:

|n|fixed-cycle candidates|accepted fixed-cycle tournaments|isotypes|unimodular M2|singular M2|self-converse isotypes|
|---|---:|---:|---:|---:|---:|---:|
|3|1|1|1|1|0|1|
|4|4|4|1|1|0|1|
|5|32|15|3|3|0|3|
|6|512|48|8|7|1|2|
|7|16384|147|21|16|5|7|

These recover the applicable older HYP-9445 data by a different canonicalization route. They do not independently verify the user's new orders10–12. The procedure avoids a full n!-permutation canonicalization once uniqueness is known, but still has exponentially many chord candidates; it is not a closed all-order enumeration algorithm for the special word grammar.

## 5. A completely specified same-insertion hostile

The following matrices are the two order-five parents, rows in the marked cycle order0,1,2,3,4:

    A: 01000       B: 01000
       00100          00100
       10010          10011
       11001          11001
       11100          11000.

They are nonisomorphic, have a unique Hamilton cycle, share sorted score sequence(1,1,2,3,3), and both have det M2=1. They differ by the chord2↔4.

Apply the identical marked operation to each: add vertex5 with5→0 and j→5 for every j=1,2,3,4. The new distinguished cycle is0,1,2,3,4,5. Both children retain a unique Hamilton cycle and share sorted scores(1,1,2,3,4,4), but their observer determinants are0 and1.

This independently realizes the determinant/coarse-data obstruction described by the user. It is **not asserted to be the exact original grammar-type or b/e example**, because this explicitly marked operation uses a different convention. The operation here is fully specified and independently reproducible.

The singular child's complete primitive integer right-kernel is certified separately in [the kernel report](../singularity/MARKED_HOSTILE_KERNEL_CERTIFICATE.md). With q_ij=e_(i,j)−e_(j,i) in end-pair column coordinates, a basis is

    q02−q12,  q03−q13,  q02+q25,  q03+q23+q35.

It has Smith form1^26 0^4. Thus the observer is singular but its cokernel has no finite torsion. A blanket replacement of “singular” by “has torsion” would fail already here.

## 6. The kernel can be repaired by another unique-cycle insertion

Keep the singular six-vertex child and adjoin vertex6 with

    6→{0,1,5},   {2,3,4}→6.

The resulting tournament has unique Hamilton cycle(0,1,2,3,4,6,5), and det M2=−1. Thus singularity is not inherited by arbitrary vertex extensions even within the unique-cycle family.

This is not a counterexample to the now-proved be/ee criterion: this extension gives the new vertex indegree three, while every actual grammar insertion gives it indegree one. It is a counterexample to proving that criterion by an unqualified claim that an existing singular observer remains singular after every unique-cycle insertion. The [bridge–sink proof](../bridge/BRIDGE_SINK_STAR_LEMMA.md) supplies exactly that guarded kernel transport, choosing its witness at the last forbidden pair.

Structurally, the six-vertex child consists of two C3 blocks, all edges from the high block{2,3,4} to the low block{5,0,1} except the reversed bridge1→2. The new vertex is inserted between these blocks. The [independent repair report](../singularity/REPAIR_AND_RELAY_SCOPE.md) includes the exact certificate and a bounded relay-family probe; no all-order relay theorem is claimed.

## Reproduction

`python -B necklaces/check_chord_necklaces.py` and its optimized version produce the same saved `results.json`. The observer is built by permutation/path-cut enumeration, not imported from the subset-DP producer in the singularity lane. The companion lane checks the hostile and repair by both constructions and the universal insertion update. All test failures raise exceptions under optimization.
