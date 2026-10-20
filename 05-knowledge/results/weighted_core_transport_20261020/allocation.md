# Run-shuffle transport on weighted tournament cores

**Status:** PROVED identities and explicit counterexamples below;
FINITE-EXACT checks in `allocation_check.py/.out`. No bounded universal
state or new closed formula for A000568 is claimed.

## 1. What is inherited, and what this adds

The previous [growth packet](../metagraph_growth_20261020/growth.md) gave:
canonical transitive-chain roots, a unique scale from chain-length gcd,
the exact uniform-inflation complement mask, and the `ABBCCA` witness
showing that a block word alone loses the complement target class.

The canonical chain decomposition is inherited from
[THM-4463, homogeneous tournament contacts, section 2](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).
The Hamiltonian path-count mechanism of splitting into maximal block runs
is already present in
[THM-3121, path-cover walk-content substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
and uniform substitution responses in
[THM-3235, transitive blowup saddle and observer twins](../../../01-canon/theorems/THM-3235-transitive-tournament-blowup-saddle-and-hypergeometric-decimation.md).
Local single-clone growth also connects to
[THM-4145, rooted homogeneous-pair expansion](../../../01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md).

The additions here are an exact allocation carrier, its operation rules,
a finer *path* dilation scale, a compressed blue test, and two sharply
scoped failures of smaller states.

The board is: **weighted core; ordered visits; allocation shuffle;
path scale; sparse complement defect; canonical normalization**.

## 2. A lossless allocation carrier

Let `Q` be the quotient by maximal transitive chains, with positive chain
lengths `l_i`. Its expansion `T=Q[K_(l_1),...,K_(l_b)]` has vertices
`(i,a)`, `1<=a<=l_i`, with `(i,a)->(i,c)` iff `a<c`. Different blocks
inherit their arcs from `Q`. The quotient has no homogeneous pair.

Fix a Hamiltonian path of `T`. Split it into maximal consecutive runs
inside one chain. Store the ordered list

```text
R=((i_1,s_1),...,(i_r,s_r)),
```

where adjacent block indices differ and `Q` directs `i_t -> i_(t+1)`.
Number the visits to each chain `i` chronologically as `1,...,r_i`.

For every chain, store a word

```text
a_i=(a_i(1),...,a_i(l_i)),
```

read in **increasing clone-label order**: `a_i(c)=j` means clone `c`
is allocated to the `j`th visit to chain `i`. Its letter multiplicities
are exactly that chain's visit sizes.

**Bijection.** Given `Q`, lengths, `R`, and all `a_i`, the path is
reconstructed by replacing each visit by its allocated clone labels in
increasing order. Conversely, every Hamiltonian path gives exactly this
data. Increasing order is required within a consecutive visit, but not
between two separated visits to the same chain.

For a fixed run list, all choices of these multiset words are valid.
Their number is

```text
M(R)=product_i l_i! / product_(visits j to i) s_(ij)! .    (1)
```

Thus a lexicographic multiset rank stores each `a_i` as an integer in
`[0,M_i-1]`. Conditional on the run list, a joint mixed-radix rank in
`[0,M(R)-1]` is optimal among fixed-length binary encodings, requiring
`ceil(log2 M(R))` bits. This is a conditional information claim, not a
claim that the chosen core or run-list representation is globally minimal.
The companion implements exact ranking and unranking.

For the earlier `ABBCCA` example with two clones in each chain, the two
possible assignments to the two singleton visits to `A` are `a_A=12`
and `a_A=21`. The missing information is precisely one bit.

## 3. Reflection, dilation, and a refined scale

### Reflection

The grid reflection `sigma` takes the converse tournament and reverses
the observer path. On this carrier:

```text
Q -> Q^op,
R -> reverse(R),
a_i(c) -> r_i+1-a_i(l_i+1-c).                            (2)
```

The last rule both reverses clone order in the transitive chain and
reverses the numbering of its visits. Formula (2) is an involution.

### Uniform dilation

Inflate **both** tournament and marked path by replacing each path
vertex with `k` consecutive ordered clones. The exact rule is

```text
l_i -> k*l_i,
(i,s) in R -> (i,k*s),
each letter of a_i -> k consecutive copies of that letter. (3)
```

The core `Q` and the sequence of run block indices stay fixed. Reflection
commutes with this operation. In a run-length encoding of the allocation
words, dilation multiplies every length by `k` without adding any new
constant runs.

### The decorated-path scale

Let `delta` be the gcd of the lengths of all maximal constant runs inside
all allocation words `a_i`. Then

```text
delta divides gcd_i(l_i);
delta(D_k(T,path)) = k*delta(T,path).                     (4)
```

**Unique path-root criterion.** A marked tournament is the uniform
`k`-dilation of another marked tournament exactly when `k` divides
`delta`. In that case divide all constant-run lengths of the `a_i` by
`k`, and all visit sizes and chain lengths by `k`.

**Proof.** In a dilated path, the `k` consecutive clone labels replacing
one old clone belong to the same visit. Hence each allocation word is
constant on its consecutive length-`k` label groups, which is equivalent
to divisibility of every maximal constant-run length by `k`. Conversely,
that divisibility groups clone labels into unambiguous consecutive
`k`-tuples. Each visit receives whole tuples, so its size is divisible by
`k`; collapsing them gives a valid path and reconstructs the original
under (3). The canonical chain-root theorem from the previous packet
makes this unique up to marked isomorphism. ∎

This is a refinement of the node scale: a tournament can admit a uniform
root while a particular interleaved path does not. The `ABBCCA` paths in
`D_2(C3)` have node scale 2 but path scale 1.

## 4. Blue can be tested without expanding the tournament

Mirrored entries of the run list force a possible map `h` of core
vertices: if run `t` lies in `i`, and mirrored run `r+1-t` lies in `j`,
then `h(i)=j`. The presentation is blue if and only if all of the
following hold:

1. Mirrored run sizes agree and the forced `h` is consistent.
2. `h` is a weight-preserving involution and an anti-automorphism of `Q`.
3. For every `i`,

   ```text
   a_(h(i)) = reverse(a_i) with letters j -> r_i+1-j.     (5)
   ```

**Proof.** A blue presentation has the reversing-path anti-automorphism
`v_t -> v_(N+1-t)`. It must send canonical transitive chains to canonical
transitive chains, reverse their internal orders, and induce the displayed
`h`; this proves necessity. Conversely, conditions 1--3 build precisely
that anti-automorphism on every clone, so the path-coordinate matrix is
fixed by `sigma`. ∎

Every core vertex occurs in the path, so the run list forces at most one
candidate `h`; there is no search over all core permutations. After the
run records are indexed, this test costs `O(b^2+N)` on expanded words and
can compare run-length encoded allocations directly. It needs no
`N`-by-`N` adjacency matrix. Black is the failure of this exact blue test.

### Endpoint closure is another directly readable bit

The incoming
[THM-4729, closed and open tilings, law L1](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md)
records the closed/open partition of complement-line endpoints. Only its
elementary endpoint law is used here; no pure-blue classification is a
dependency.

Let `c(P)=1` when the final path vertex beats the first. If their chain
indices differ, read `c` directly from the corresponding arc of `Q`.
If both lie in chain `i`, find the smallest clone label allocated to the
first visit and the largest clone label allocated to the last visit.
Then `c=1` exactly when the latter label is smaller than the former.

For `N>=3`, that endpoint edge is a free path-coordinate edge, so ordinary
`tau` flips `c`. Reflection preserves `c`, because it both reverses the
path and takes the converse. Uniform marked-path dilation preserves `c`
for `N>=2`, since the two extreme clones belong to the two distinct old
endpoint vertices and inherit their arc.

Thus each ordinary line, blue or black, has one open and one closed
presentation. This is a bit of a **presentation**, not generally a bit
of an entire tournament node. It separates the original two-clone
`ABBCCA` example, but does not repair the inversion obstruction below.

## 5. Ordinary complement: exact sparse defect, then normalization

For any tournament `T` with Hamiltonian path `P`, all-free-bit complement
is

```text
tau_P(T) = T^op XOR E(P),                                (6)
```

where `E(P)` is the set of its `N-1` consecutive unordered path edges.
Converse reverses every arc; the XOR restores precisely the mandatory
path arcs. This formula holds independently of modular decomposition.

To use a weighted core, replace `Q` by `Q^op`, reverse the clone labels
inside each chain, transport `P` by those label reversals, and store
the transported path-edge mask. The mask has only `N-1` edges. Its
incidence is already specified by the run-shuffle carrier, so a full
dense target matrix need not be stored.

There are two honest ways to continue:

- Keep this weighted reference tournament plus a structured defect mask.
  Its weights describe a reference expansion, not necessarily the actual
  maximal transitive modules after the defects.
- Recover the actual homogeneous-pair components of the output, form its
  new canonical weighted core, and encode the still-valid observer path
  again by its visits and allocation shuffles.

The second route is closed under ordinary `tau`, and applying it twice
returns the original fixed-path tournament. Reflection commutes with it.
The companion checks these identities across every fixed-path cube
through order six.

This is exact semantic closure. It is not a proof that canonical
normalization has bounded state size, or that canonical class matching
can be performed using only a few scalar counters.

## 6. A local additive rule, with its exact scope

To append a new vertex that is another clone in an existing chain `i`,
choose its insertion rank `p` among that chain's ordered labels.

- If the current path ends in another chain `j`, require `j->i` in `Q`.
  Append a new singleton visit to `i`, and insert its new visit label at
  position `p` in the allocation word `a_i`.
- If the path already ends in chain `i`, require the new clone rank to
  exceed the rank of the path's final clone. Increase the final visit
  size by one and insert its existing visit label at position `p`.

Increase `l_i` by one. All old path vertices retain their relative
relations. These rules give every legal endpoint append within this
single-clone growth class and are reversed by deleting that appended
clone. The weighted core stays fixed.

Arbitrary new-vertex adjacency patterns are broader: they can split old
chains and change `Q`. For those, exact reconstruction plus canonical
normalization is required, or the nonuniform new adjacencies must remain
as defects. No local fixed-core rule is asserted for arbitrary growth.

## 7. Pairwise inversion counts are still insufficient

Take the cyclic core `A->B->C->A` with chain lengths `(4,1,2)`, and the
same run word in both cases:

```text
AA B CC AA.
```

Give the four ordered `A` clones the visit-allocation words

```text
a_A = 1221,      a'_A = 2112.
```

Both have the same visit sizes `(2,2)` and the same inversion count 2.
All other chains have one visit. Therefore their entire vectors of
pairwise visit-inversion counts agree. Their multinomial ranks are,
however, 2 and 3. Both presentations are open, so their endpoint closure
bits agree as well.

The ordinary complement targets have degree multisets

```text
(1,2,3,3,3,4,5),          (1,2,2,4,4,4,4).              (7)
```

These are unequal even after allowing `d -> 6-d`, so the targets are
different converse-merged nodes.

The first failed implication is **same pairwise run-inversion counts
and endpoint closure imply the same complement target**. The surviving exact compression is
the multiset rank or the complete run-length encoded shuffle. Inversion
counts can remain useful statistics; they cannot replace this coordinate
without a new sufficient-state proof.

## 8. What cannot remain bounded

### Canonical core size

Let `K_N` be the transitive tournament with its unique Hamiltonian path.
It has one canonical chain, so its weighted core has one vertex.
After ordinary complement, the arc from `i` to `j>i` is forward exactly
when `j=i+1`.

For every `N>=5`, this output has **no homogeneous pair**. To prove it,
take `i<j`.

- If `j=i+1`, use `i-1` when `i>1`, and vertex 3 when `i=1`, to
  distinguish their incidences.
- If `j=i+2`, vertex `i+1` distinguishes them.
- If `j>=i+3` and `i>1`, vertex `i-1` distinguishes them.
- If `i=1`, `j>=4`, and `j<N`, vertex `j+1` distinguishes them.
- For the remaining pair `(1,N)`, vertex 3 distinguishes them when `N>=5`.

Thus one ordinary complement sends a one-vertex weighted core to an
`N`-vertex core with all weights one. This proves that bounded canonical
core order is not closed under `tau`. It does **not** prove that every
symbolic description must be large: this particular family still has
the short formula above, or the sparse defect description (6).

### Allocation information

Even with the cyclic three-vertex core fixed, weights `(2m,1,1)` and run
word `A^m B C A^m` permit `binom(2m,m)` different allocation shuffles.
The weighted core has trivial automorphism group, so these are distinct
fixed-path presentations. Any lossless representation of all of them
needs at least

```text
log2 binom(2m,m) >= 2m-log2(2m+1)
```

bits in the worst case. The inequality follows because the middle
binomial coefficient is maximal among the `2m+1` coefficients summing
to `4^m`. This lower bound concerns literal presentations; it is not a
claim that all their next complement targets are distinct.

## 9. Reproduction and the next exact target

```text
python 05-knowledge/results/weighted_core_transport_20261020/allocation_check.py
python -O 05-knowledge/results/weighted_core_transport_20261020/allocation_check.py
```

The script imports the previous packet's elementary tournament utilities.
It has no assertions. It checks 2,785 shuffle rank/unrank words, all 1,100
fixed-path presentations through order six, 2,200 dilation transports,
1,109 decorated-root controls, and 5,164 admissible clone-appends; it
also checks 1,098 endpoint-closure flip controls, the inversion hostile
with its common open bit, and the core-growth family through `N=100`.
Normal and optimized Python pass with identical output.

The constructive survivor is a **weighted core plus run list plus exact
allocation ranks**, with sparse path defects when ordinary complement
destroys the module structure. The next reduction should compare states
by the actual continuation operations they must support, rather than
discarding allocation information because one statistic agrees. A
compact quotient of these allocation ranks that preserves the desired
operation language remains OPEN.
