# The30/66 compensation targets cannot restore unit deletion count

**PROVED by an elementary all-order reduction and a FINITE-EXACT exhaustive
certificate:** no tournament T of any order with H(T)=1409 or1427 has a
vertex w satisfying H(T-w)=1. This is stronger than the requested fixed-name,
fixed-order, root-end-zero repair exclusion. It does not exclude1409 or1427
from the Hamiltonian spectrum: the actual source children already attain
them with deletion count3.

No strong-order literature bound is needed. The entire remaining finite
universe consists of2047 one-apex bit words. All are counted independently
by literal Hamilton-path subset DP as well as two combinatorial formulas.

**Counting mechanism inherited:** the earlier
[`feedback_word_hamilton_spectrum_20261020/proofs.md`, Sections1--4](../../feedback_word_hamilton_spectrum_20261020/proofs.md)
already proves the weighted binary-inversion formula, SCC reduction,
four-coordinate positive transfer, sharp extrema, and the all-order missing
value13. Our bits are its complemented convention: here b_i means v_i->w,
whereas that packet's bit means w->v_i. The split calculation below is an
independent verification and specialization of that mechanism, not a new
counting formula or a novelty claim. The present new target is the specific
empty fibres1409/1427 and their implication for the actual30/66 debt repairs.

## Inherited source and the precise failed requirement

The closest actual mechanism is
[`entry_selector/PROOF.md`](../../c4_port_returns_20261103/entry_selector/PROOF.md):
two prepared q=1 targets of counts469 and475 undergo exact C2 births and
give strong children1409 and1427, retaining root99 and the seed mark0.
The deletion coordinate at0 changes1->3. The single triangle reversal
that retains strongness restores q=1 but lowers both H and E99(2):

| Original child | Strong q=1 repair | Required increase in(H,E99(2)) |
|---:|---:|---|
|1409|1379|(+30,+30)|
|1427|1361|(+66,+66)|

The original finite statement only rejected the three individual triangle
reversals. The present theorem rejects **every final tournament** at either
desired count with q=1, regardless of what intermediate operations, births,
deletions, relabelings, or changed deletion orders are allowed. Thus adding
an arbitrary compensation while retaining q=1 cannot solve these targets.

The canonical hostile to overstatement is the actual q=3 child, which has
the desired H and root response. The corrected near miss is treating a
failed fixed physical star as a general impossibility without first reducing
all possible deletion orders. The least-used coordinate here is the exact
transitivity forced by q=1. Our concept board is named source / deletion
resource / transitive order / apex split / endpoint response / native guard.

## 1. Unit deletion count is exactly a one-apex condition

For any tournament U, H(U)=1 iff U is transitive. One direction is immediate.
For the other, a nontransitive tournament contains a directed triangle,
which has three Hamilton paths. Insert every remaining vertex into each
path by a deterministic valid tournament insertion. Deleting those inserted
vertices recovers the original triangle path, so the three resulting paths
are distinct. Thus H(U)>=3.

Consequently H(T-w)=1 exactly when the other vertices have a unique
transitive order

```
v0 -> v1 -> ... -> v_(m-1).
```

Write b_i=1 iff v_i->w. The pair(T,w) is now specified completely by that
binary word, up to renaming the ordered vertices. For m>=2, T is strong iff

```
b0=0 and b_(m-1)=1.                              (1)
```

Necessity: otherwise v0 is a source or v_(m-1) is a sink. Sufficiency:
`w->v0->...->v_(m-1)->w` is a directed spanning cycle. This tests all
possible deletion orders; it is not a choice of one fixed physical order
inside the source graph.

## 2. Exact split formula and an elementary size cutoff

Every Hamilton path has a unique split `A,w,B`, where A and B are increasing
words partitioning the m ordered vertices. If both pieces are nonempty,
write i=max A and j=min B. The boundary arcs are legal precisely when
b_i=1,b_j=0.

If i<j, the split is possible only when j=i+1, and it is unique. If j<i,
all indices before j must lie in A and all indices after i in B; the
i-j-1 interior indices are arbitrary. There are exactly2^(i-j-1) splits.
The two empty-piece cases give the formula

```
H(b) = (1-b0) + b_(m-1)
       + sum_(i=0)^(m-2) b_i(1-b_(i+1))
       + sum_(0<=j<i<m) (1-b_j)b_i 2^(i-j-1).       (2)
```

It holds for every nonempty binary word, including weak graphs. For m=0,
the singleton count is1 and is handled separately.

Under strongness condition(1), the two empty-piece paths both exist, and
the term j=0,i=m-1 always contributes2^(m-2). All other terms are
nonnegative. Hence for strong one-apex tournaments of order n=m+1>=3,

```
H(T) >= 2 + 2^(n-3).                              (3)
```

At n>=14 the bound is at least2050, already greater than both targets.
Thus **every strong q=1 target of H1409 or1427 must have order3,...,13**.
This bound supplies the all-order reduction independently of the more
general Moon--Busch bound. No numerical asymptotic estimate or assumed
finite-template closure is used.

The inherited sharp one-apex minimum is stronger: for n>=5 it is
`(2^floor((n-3)/2)+1)(2^ceil((n-3)/2)+1)`, giving2145 at n14. The script
also verifies this minimum throughout its complete census. Either bound
gives the same sufficient cutoff here.

## 3. The complete finite result

For order n, condition(1) leaves exactly2^(n-3) middle-bit assignments.
Their total over n=3,...,13 is2047. The exact census is:

| Order | Strong words | Distinct H values | Minimum H | Maximum H |
|---:|---:|---:|---:|---:|
|3|1|1|3|3|
|4|2|1|5|5|
|5|4|2|9|11|
|6|8|4|15|23|
|7|16|9|25|51|
|8|32|14|45|107|
|9|64|31|81|227|
|10|128|57|153|467|
|11|256|116|289|963|
|12|512|211|561|1955|
|13|1024|435|1089|3971|

Neither1409 nor1427 occurs at any of these orders. These are marked
transitive-deletion words, not unlabeled tournament isomorphism classes;
duplicates of an unmarked tournament do not harm completeness.

`one_apex.py` verifies every count three ways:

1. Formula(2).
2. Literal enumeration of every subset A, checking its maximum incoming
   edge to w and the minimum outgoing edge from w to its complement.
3. A separately compiled literal directed-path subset DP on the actual
   adjacency matrix, independent of the transitive-split formulas.

All2047 matrices are checked strong. Another62 arbitrary words at orders2
through6 include weak boundary controls. The bit subfamily starting01 and
ending1 has1024 words through13; its empty target set follows already from
the broader strong census. No root-end condition is needed for the theorem.

## 4. Why arbitrary weak padding cannot evade the result

In a tournament whose deletion of w is transitive, every directed cycle
contains w. Hence at most one strongly connected component can have more
than one vertex, and that component contains w. Hamilton paths concatenate
in the total order of SCCs, so H is the product of their Hamilton counts.
If H>1, the sole nontrivial component C therefore satisfies

```
H(C)=H(T), C strong, and H(C-w)=1.
```

The strong no-hit result applies to C. Consequently no weak tournament of
any larger order can attain either target with a unit deletion count either.
This completes the stated all-order exclusion, including arbitrary choices
of the deleted vertex and every physical labeling.

By inherited tournament Hamilton-count oddness, every deletion count at
either target is therefore at least3. The actual C2 children attain3, so
the unit reset is asking for a resource value below the attainable minimum.
The current theorem leaves preparation with q=3, different endpoint guards,
and larger actual operations open.

## 5. Exact retained-name debt sources and fresh observations

Both strong q=1 repairs use the same unique deletion order after removing
physical name0:

```
99,10,11,2,1,3,4,12,200,401,400,201.
```

Their words in this order are

```
H1379: 011010001001;
H1361: 011011001001.
```

The debts have an explicit inherited-word explanation. The prepared targets
have deletion words `0110100011` and `0110110011`, both ending11. For such
a word b of length m put

```
L = sum_(j<m-1, b_j=0) 2^(m-j-2).
```

Inserting one0 immediately before the final1 changes H to `H+L+2`.
Inserting two0s there changes it to `H+3L+4`: the old contributions to the
final1 are multiplied by4, the two new zeros contribute3, and the preceding
1 acquires one additional adjacent cut. The actual first C2 birth is the
one-zero insertion after deleting0; the strong triangle repair is the
two-zero insertion. Since the C2 child itself has count `3H+2`, its loss
under that strong q1 repair is exactly

```
2H-3L-2.
```

For H469, L302 gives first count773, repaired count1379, and loss30.
For H475, L294 gives first count771, repaired count1361, and loss66.
The script checks these identities against the actual named stage matrices,
not just abstract words.

The independent counter reconstructs full H, P2, E1/S1/E2/S2, every ordered
J pair, and every single-vertex deletion count directly from the saved
actual matrices. For each debt source it verifies q0=1 and E99(2)=H.
For each original child it verifies the desired H, q0=3, and E99(2)=H.
The retained matrices, physical names and full profiles are in the receipt.
No endpoint profile is transported from an equal-count graph, and no native
legality is inferred from the numerical equalities.

## Reproduction and counter interface

```
python 05-knowledge/results/compensation_returns_20261103/recovery/one_apex.py
python -O 05-knowledge/results/compensation_returns_20261103/recovery/one_apex.py
```

The normal/optimized receipts `one_apex.json` and `one_apex_optimized.json`
are identical. The script compiles `audit/exact_counter.cpp` in a temporary
directory. That counter accepts one line `id n row0 ... row_(n-1)` per graph.
Default output is its full endpoint/joint/deletion profile; argument `h`
requests only H. It checks actual tournament rows and uses signed64-bit
arithmetic only for n<=19, where every ordered two-slot count is at most
20!<2^63. The present exhaustive matrix census stops at13.

The finite negative is part of the proof, not a guessed pattern. The result
closes the requested q=1 compensation target; it does not resolve the global
`{7,21}` conjecture, repeated positive-digit closure, or A000568.
