# The same product-and-surplus mechanism explains 7, 21 and the role of 3

**CITED exclusions, with an independently checked local certificate.**
Neither exclusion is new: see [THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
and the [original bof–Royle discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament).
The purpose here is to expose the same ingredients that a constructive
coverage proof would need to control.

## The count is a product plus two nonnegative defects

For strong A and B sharing exactly v, keep every cross edge in T. Each
Hamilton path of A and B splits before/after v. Merging the two prefixes
and the two suffixes produces at least one Hamilton path of T. Define

* S: the sum, over pairs of input paths, of the number of merges minus one;
* O: the number of output paths whose restriction to A or B is not a path.

Then the inherited exact identity is

    H(T)=H(A)H(B)+S+O.                              (E1)

Both defects are nonnegative. Their sum is even because all three counts
H(T), H(A), H(B) are odd. For a triangle B, this is

    H(T)=3H(A)+2e,         e=(S+O)/2≥0.              (E2)

The 3 counts the three ways to cut the directed triangle into a Hamiltonian
path. Equality requires both unique merging and no broken restrictions.
Merely seeing a triangle does not imply an exact factor three.

## A short common proof framework for the two exclusions

Use the [Busch decomposition and merge lemmas](GRAMMAR.md). Every nontrivial
strong tournament with count below9 must be a single nearly transitive atom:
two nontrivial strong pieces would contribute at least3·3. The only atom
counts below9 are3 and5, on three and four vertices respectively.
Therefore no strong tournament has count7. SCC multiplication excludes7
for arbitrary tournaments as well.

If a strong tournament had count21, it could not be a nearly transitive
atom. In its two-piece decomposition each factor has count at most7, since
the other factor is at least3. Thus each piece is C3 or R4. Two R4 pieces
already give at least25; the only candidates are C3 with C3, or C3 with R4.

All16 cross matrices for C3 with C3 give counts {9,11,13,15}. For C3 with
R4, the following two 8×8 tables exhaust the possibilities up to converse.
Neither contains21. This closes the strong case. Finally21=3·7 cannot be
rescued by separate SCCs, since a component with count7 or21 would be
required. This proof uses the structural decomposition, not an extrapolated
minimum sequence or a large-order census.

## The entire small local obstruction

Use R4 with outgoing rows (6,12,8,1), vertices0,1,2,3. For a chosen root r,
add a,b with r→a→b→r. Let i,j∈{0,…,7} encode a's and b's outgoing edges
to the other old vertices, in their increasing order. Rows below are i,
columns are j. Roots2 and3 are the converses of roots1 and0 respectively.

Root0:

    19 33 29 37 25 37 31 37
    37 43 45 43 43 45 43 37
    29 43 31 37 33 43 29 31
    37 43 37 33 41 43 33 25
    23 41 33 41 23 37 29 33
    33 43 43 43 33 37 37 31
    23 41 25 33 23 37 19 23
    23 33 25 25 23 29 19 15

Root1:

    23 33 29 37 31 37 33 37
    37 41 45 45 43 43 41 37
    31 43 29 37 37 43 29 31
    37 43 37 37 43 45 33 29
    33 41 37 41 37 37 37 33
    37 41 45 45 37 33 37 29
    33 41 29 33 37 37 27 23
    29 33 29 29 31 29 23 17

The local obstruction is precise: above the product15, surplus0,2,4 occur
but surplus6 does not. Larger surplus8 does occur. The implementation
counts every row independently by literal vertex permutations, endpoint
dynamic programming and odd-cycle collections. It also records S and O
separately for all256 cases, verifying E1 rather than just a scalar table.

## Multiplying by 3 has several distinct meanings

|Construction or observer|Exact meaning|Information that must survive|
|---|---|---|
|Ordinal C3 component|H(C3→T)=3H(T)|The one-way SCC separation|
|Triangle sharing a vertex|3H+S+O|Cross matrix, merge multiplicity and broken restrictions|
|Two indegree-one births|H+2F=3H+2(F−H)|The actual ending two-cover count F at the root|
|Three path bins in a transitive block|3^m assignments before nonempty-bin correction|The run count and its boundary conditions|

The third law and its forced three path classes are inherited from
[the two-root normal form](../mutation_collatz_chanpak_20261112/mutation/TWO_ROOTS.md).
The fourth concerns ordered covers: exactly j nonempty bins in TT_m give
j!S(m,j). It is an exponential mode, not divisibility of the whole H count.
Conflating these four mechanisms would erase the data needed for recursion.

The reverse side of the exclusion proof is therefore not 'divide every
multiple of3 by a triangle'. A strong count63 witness has no ordinal
triangle factor and does not require a count21 parent. Our certificate uses
a count15 parent and surplus18. For coverage, the task is to choose a
smaller realizable core and a realizable surplus together.

## Complete same-count preparation barriers for triangle attachment

The exhaustive source bank and the published minimum f(8)=45 make the
following statements complete at every source order. For each listed H,
even arbitrary H-neutral preparation to a different **strong** core and
arbitrary re-marking cannot enable the indicated triangle attachment:

|Source H|Missing offsets d in 3H+d|
|---:|---|
|11,13|0|
|15|2,4|
|17|4|
|25|2|
|27|2,4|
|33|4|
|35|2,4|
|41,43|2,4|

The other attained strong counts below45—3,5,9,19,23,29,31,37,39—admit
all three offsets0,2,4 after a suitable choice of core and root. The table
excludes an operation on a source fibre, not the target count itself.
Weak cores and operations that remove part of the old core are outside
this particular source universe. The separate all-tournament H35 readiness
obstruction is [O1–O3](../ready_return_pentagon_20261103/OBSTRUCTION.md).
