# Quinary paths, exact carriers, and regular-parent growth

2026-10-30. **Exact computation with proved orbit quotients.**
The Hamilton-path number below is recovered prior work; the present
contribution is the verified carrier and its precisely scoped minimization.

## 1. Recover the old count before claiming a new one

Let T5 be the carousel on Z/5 with connection set{1,2}, and let T=T5[T5]
be its25-vertex lexicographic square. The exact count is

    H(T)=2398512090439546875.

It already appears in
[the October21 circulant experiment](https://github.com/eliottcassidy2000/math/blob/1bd2d9b9b913ac1f5e092a1fd3edadb25b7ad9f9/04-computation/experiments/lattice_bridges_20261021/circ25_H.txt),
at connection set{1,2,5,6,7,10,11,12,16,17,21,22}.
The map x=child+5*inner identifies that circulant with T.

The full automorphism group has order5^6=15625. Its action on ordered
Hamilton paths is free, and

    H(T)/15625=153504773788131.

This quotient counts path orbits, not labelled paths.

## 2. The exact subset-and-endpoint quotient

The ordinary Hamilton-path DP has states(visited subset,last vertex).
For T, normalize the current child to0 by the outer C5 rotation, and
normalize its endpoint to0 by the inner rotation.

- The current child's visited mask contains0:16 possibilities.
- Each of the four other children has an unpointed visited subset modulo
  its C5 rotations:8 possibilities.
- Retain all five child states jointly.

Thus there are16*8^4=65536 orbit states. Each transition enumerates every
actual available target vertex, normalizes the successor, and sums
multiplicities when successors agree. Full wreath automorphisms make
this an equitable partition of the raw transition graph.

Initialize the singleton orbit with mass25, not1. Then the recurrence
counts total labelled prefixes and therefore the stated H(T). Forgetting
this initial orbit size or dividing transitions by stabilizers would
change the answer.

The [portable checker](portable_checks/quinary_paths.js) independently
compares transitions on all2304 raw pointed states of the nine-vertex
C3[C3] control and on1000 quinary raw states. It also performs a separate
raw subset DP at order9, giving3159. Forward and backward counts agree.

The current orbit states have

| Graph | Full orbit states | Reachable | Full weighted partition | Reachable weighted partition |
|---|---:|---:|---:|---:|
| C3[C3] |64|48|57|42|
| C5[C5] |65536|63488|25437|24297|

The last two columns are the coarsest partitions preserving the level
and the sums of outgoing weights into each partition class. Since all
edges increase level, reverse-level refinement proves minimality in
that stated category.

An initial validation-only branch had a dangling else, oversampling the
ternary control and skipping the quinary samples. It was repaired in
commit a039df2bacadc46f3aaece7c5188d0474681c3f5, with exact coverage
assertions. Hamilton counts and partition counts were unchanged.

### Why cardinality summaries are insufficient

The two visited child subsets{0,1} and{0,2} both have size2 and trivial
subset stabilizer. Their remaining three vertices induce, respectively,
a transitive tournament and a directed triangle. A later single run
through the remaining child therefore has1 versus3 Hamilton orders.
The eight actual subset types cannot simply be replaced by six sizes.

Likewise the earlier folded64->24->16 experiment uses paired selections
under an involution and a fixed terminal root. It is NOT this64-state
unrestricted subset/endpoint carrier. Equal state counts do not identify
representations or their permitted observations.

## 3. Partition size is not linear observation rank

Let A contain outgoing transition multiplicities and f indicate the
full-visited level. A^j f is supported only on level n-j. For0<=j<n it
is nonzero: a Hamilton path supplies an extendable prefix. Hence these
vectors have disjoint nonzero supports, while A^n f=0.

The acceptance-only linear observability rank is therefore n:25 here,
9 in the ternary control. The ternary support sizes, checked directly,
are1,2,4,7,10,11,8,4,1.

This does not contradict the partition minima. Acceptance-only linear
observations discard explicit level information on dead states; the
weighted partitions preserve it. Preserving all level indicators
linearly has the elementary upper bound n(n+1)/2, not a proved exact
rank here. Moreover the completion-count embedding itself requires
solving the DP. Rank25 is not a free25-state counting algorithm.

Named-vertex questions require another refinement: normalizing the last
vertex to0 does not retain whether it was literally vertex3.

## 4. Subset recurrences also count actual children

For the cyclic-prime wreath action let B_0=2 count unpointed subsets and
A_0=1 count included-pointed subsets. Then

    B_k=(B_(k-1)^p+(p-1)B_(k-1))/p,
    A_k=A_(k-1) B_(k-1)^(p-1).

These count subset orbit states, not paths. At p=5 they give B_1=8 and
A_2=65536.

The stronger stabilizer inventory is the classical wreath/necklace
recurrence. If P_k(z) counts binary subset orbits by the base-p
logarithm of their stabilizer order,

    P_0(z)=2,
    P_k(z)=(P_(k-1)(z)^p-P_(k-1)(z^p))/p
             +z P_(k-1)(z^p).

Nonconstant words have no nontrivial outer rotation because p is prime;
constant words acquire one extra factor p.

This recurrence and its regular-parent use are shared with the
concurrent [prime-power-plus-one work](https://github.com/eliottcassidy2000/math/blob/51cd593a3ef755f65d0d0014ff52922e183a34de/05-knowledge/results/prime_power_plus_one_20261029/GROWTH.md)
and [the full inventory](https://github.com/eliottcassidy2000/math/blob/b8d24efc9a9d30ab327f45001d5eab60abd227af/05-knowledge/results/prime_power_plus_one_symmetry_20261030/PROOFS.md).
The1554 count is a concurrent recovered result, not an exclusive
discovery claimed here.

### Forgetting the added root correctly

For a regular parent T, every automorphism of a one-vertex extension
fixes the new vertex. A child has at most two regular-deletion roots:
two are twins outside their mutual edge, and three are impossible.
Root forgetting merges exactly the open/closed clone-mask pair for each
vertex orbit. Thus unmarked child classes equal subset orbits minus
vertex orbits, refined by the corresponding stabilizer.

For the vertex-transitive tower, with
g_k=(p^k-1)/(p-1), this gives

    E_k(z)=P_k(z)-z^(g_k-k).

This is the [regular-parent theorem](https://github.com/eliottcassidy2000/math/blob/b57ed62de7477c9a98ce0bdb98b8f181451cf03b/05-knowledge/results/regular_parent_symmetry_20261029/PROOFS.md),
not an unrestricted child-enumeration shortcut for nonregular parents.

At p=5,k=2,

    P_2(z)=1554+2598z+1728z²+576z³+96z⁴+6z⁵+2z⁶.

Subtracting z^4 gives6559 ordinary26-vertex extension classes, including
1554 rigid children. These are children of THIS parent, not all
26-vertex tournaments. The [independent checker](portable_checks/regular_augmentation.js)
enumerates all32768 five-child type words and all32 base extensions
under all6! relabelings. At order6 it obtains seven classes, with
automorphism orders1,1,1,1,1,5,5.

The number of rigid subset types follows R_k=(R_(k-1)^p-R_(k-1))/p.
The binary seed2 is fixed exactly at the odd prime p=3; for p=5 it grows
2,6,1554,... . This is a structural explanation of the ternary/quinary
difference, not an identification of every previous occurrence of3 or5.
