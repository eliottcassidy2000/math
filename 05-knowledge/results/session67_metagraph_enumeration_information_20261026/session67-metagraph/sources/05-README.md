# Weighted cores, allocation shuffles, and an exact structural recurrence

**20 October 2026 — PROVED scoped constructions / FINITE-EXACT controls /
OPEN smaller scalar recurrence for A000568.**

There is now a well-founded structural recurrence: generate new pair-free
cores by alternating attachments to smaller weighted cores, and generate
weights by coordinate increments modulo the core automorphisms. It
reconstructs A000568 through order eight without constructing metagraphs
or using the Davis/Burnside formula in generation. A run-shuffle carrier
also preserves the information needed for doubling, reflection and
ordinary blue/black line transport.

The compression has an exact boundary. Complete core cycle indices suffice
for counting weight allocations but can fail to determine the next growth
response. Ordinary line complement can turn a one-vertex weighted core
into a core with arbitrarily many vertices.

## 1. Inheritance and terminology

The closest proved mechanism is
[THM-4463, homogeneous tournament contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md):
the homogeneous-pair graph is a union of directed paths, whose components
are the canonical maximal transitive modules. Contracting those chains
gives a tournament Q with no homogeneous pair. Every tournament is uniquely

    T = Q[TT_(l_1), ..., TT_(l_b)]

up to a core isomorphism preserving the positive weights. The chain orders
are rigid; Aut(T) is the weight stabilizer in Aut(Q).

Here **pair-free core** means this Q. **Scale-primitive** means
gcd(l_1,...,l_b)=1, as in the
[previous packet](../metagraph_growth_20261020/README.md).
Neither term means modular indecomposability.

Canonical hostile: the two ABBCCA paths have the same visits but different
complement targets. Corrected near miss: a present-day count or even an
entire weight response need not support the next operation. Least-used
sidecar: the allocation of ordered clone ranks among repeated visits.

Anchor: structural generation of A000568. Niche: exact allocation transport.
Wildcard: cycle-index differential equations and their missing boundary.
The live board is **weighted core / attachment runs / allocation shuffles /
automorphism action / line defects / fixed-point-free cycle terms**.

## 2. The growth dichotomy and the alternating recurrence

Let T have n vertices and b canonical chains. Add x, and read its adjacency
bits along each ordered chain. Let s be the total number of changes between
successive bits.

**PROVED.** Exactly one of the following happens:

1. x has a homogeneous neighbor. Then x is inserted into one existing
   chain, all other chains stay unchanged, and the core remains Q.
2. x has no homogeneous neighbor. Each old chain splits at its adjacency
   changes, x is a singleton, and the new core has b+1+s vertices.

Old vertices cannot acquire a new homogeneous-pair edge: any old
distinguishing vertex remains present. If x has a homogeneous neighbor u,
its entire external adjacency agrees with u; this forces precisely a
chain insertion. In particular, x cannot glue different old chains.
The full proof and boundary cases are in [inheritance.md](inheritance.md).

Thus a child on n+1 vertices is itself pair-free exactly when:

- x's adjacency strictly alternates along every parent chain;
- x has no homogeneous neighbor.

There are only 2^b alternating candidates, one initial bit per chain.
Of these, the numbers accepted are

    2^b - 2b          if n=b,
    2^b - 1           if n=b+1,
    2^b               if n>=b+2.

The exceptional rejected candidates are precisely clone insertions.

Let C_q be the set of pair-free core classes of order q and W_Q(n) the
positive weight vectors of sum n modulo Aut(Q). Then the following is an
**all-order structural set recurrence**:

    C_1 = {one vertex}.

    For N>=2:
      form all parents Q[l], with Q in C_q, q<=N-1, [l] in W_Q(N-1);
      try the 2^q alternating attachment phases;
      reject children with a homogeneous pair;
      identify isomorphic surviving children to obtain C_N.

    a_N = sum_(q<=N) sum_(Q in C_q) |W_Q(N)|.

Every core is obtained by deleting any vertex and then reversing that
deletion. Its parent has smaller total order, so the recurrence is
well-founded. Weight orbits themselves obey a smaller constructive
recurrence: start at (1,...,1), increment one coordinate, and identify
vectors under Aut(Q). See [core_count.md](core_count.md).

This is not a scalar recurrence using just previous a_N. It retains actual
cores for new-core generation and uses canonical isomorphism tests.

| n | pair-free core classes | A000568 |
|---:|---:|---:|
| 1 | 1 | 1 |
| 2 | 0 | 1 |
| 3 | 1 | 2 |
| 4 | 2 | 4 |
| 5 | 5 | 12 |
| 6 | 26 | 56 |
| 7 | 255 | 456 |
| 8 | 4602 | 6880 |

At order eight the generation tests 42,890 alternating attachments;
39,188 survive before isomorphism identification. They give 4,602 new
cores, and weighted smaller cores supply the other 2,278 classes.
No metagraph or order-eight path cube is constructed.

### Root-forgetting requires a symmetry weight

Counting each valid attachment as one class is incorrect. With a_N denoting
ordinary classes, the exact general identity is

    a_(n+1) = 1/(n+1) sum_[T], |T|=n
                         sum_(all attachments z) |Aut(T+x_z)| / |Aut(T)|.

Restricting to alternating attachments and pair-free children gives the
same identity for |C_(n+1)|. For each individual child the accumulated
weight is exactly one; the script checks this, not merely the grand total.
For example, one valid alternating attachment to TT2 gives C3.
Naively dividing by three gives 1/3; its automorphism ratio is three.

## 3. An unexpectedly small exact growth kernel

For **labelled** extension counts, the two integers (n,b) do close. If
K_(n,b)(y) records the number of attachment patterns by resulting core size,
then

    K_(n,b)(y) =
      (n+b)y^b
      + y^(b+1) [2^b(1+y)^(n-b) - 2b - (n-b)y].

There are 2^b binom(n-b,s) patterns with s switches. Of these, 2b
endpoint clone insertions have s=0, and n-b interior clone insertions
have s=1. Moving those contributions to core size b proves the formula.

This supplies a closed recurrence for labelled core-size distributions.
It agrees with the independent formal-series identity

    T(z) = C(z/(1-z)),       C(z) = T(z/(1+z)),
    T(z) = sum_(n>=0) 2^binom(n,2) z^n/n!,

where C counts labelled pair-free cores exponentially. Both computations
agree through order twenty. This counts automorphism-weighted mass,
not ordinary isomorphism classes.

A basic compression limit is also explicit. In a uniformly chosen labelled
tournament, a specified pair is homogeneous with probability 2^(-(n-2)).
The union bound therefore gives

    Pr(core size < n) <= binom(n,2) * 2^(-(n-2)).

Almost all labelled mass eventually has full core size. Weighted-chain
compression remains especially useful in structured families, rather than
eliminating the large-core part of the general census.

## 4. The count-only differential recurrence and its precise missing data

Let Z_C(p_1,p_2,...) be the full cycle-index series of pair-free cores,
including the empty constant 1:

    Z_C = sum_[Q] 1/|Aut Q| sum_(g in Aut Q) product_(cycles c of g) p_|c|.

Define Z_T likewise for all tournaments. Canonical ordered-chain substitution
and attachment of one marked vertex give

    Z_T(p) = Z_C(p_i/(1-p_i)),
    partial_(p_1) Z_T(p) = Z_T(2p_1,2p_2,...).

The second identity holds because an automorphism cycle sees a newly marked
vertex with one common adjacency bit, independently for each cycle.
Combining the two gives the exact core recurrence

    (1+p_1)^2 partial_(p_1) Z_C(p)
        = Z_C(2p_i/(1-p_i) : i>=1).                     (*)

A000568 is recovered by

    sum_(n>=0) a_n x^n = Z_C(x^i/(1-x^i) : i>=1).

This is a genuine recurrence on smaller algebraic data. However, (*) leaves
the boundary Z_C(0,p_2,p_3,...) undetermined. That boundary records
automorphism terms with **no fixed vertices**, rather than a subclass of
cores. Its beginning is

    1 + (2/3)p_3 + (4/5)p_5 + (10/9)p_3^2 + (8/7)p_7 + ...

Thus deleting or marking a single vertex cannot supply all symmetry data.
This is the actual-core version of the obstruction already isolated by
[THM-514, the A000568 unit-tail core kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md).
No new priority is claimed for rooting or cycle-index methods.

The full hierarchy obtained by adjoining odd automorphism cycles recovers
the known gcd orbit calculus; it computes the boundary but does not evade
Burnside's mathematical content. The useful new target is a structural
description of this fixed-point-free boundary in the weighted-core language.

### Two hostile tests locate the compression boundary

**Same core order and same automorphism-group order are insufficient even
for weights.** Two order-six cores both have cyclic automorphism group of
order three. One action has a 3-cycle and three fixed vertices; the other
has two 3-cycles. Their numbers of weight orbits at total weight seven
are four and two. The actual permutation cycle types are indispensable.

**Even the entire cycle index is insufficient for individual-core growth.**
Two rigid order-six cores both have cycle index p_1^6:

    rows A = (0,5,25,35,11,23)
    rows B = (8,17,3,38,13,23).

Among their 52 pair-free extensions, A has 50 children with automorphism
order one and two with order three; B has 51 and one. These count
attachment instances. Thus equal complete cycle indices do not determine
the next-core automorphism profile. This does not contradict the global
identity (*), which aggregates all cores.

## 5. The smaller allocation encoding that survives the operations

Store the weighted core, the path's ordered run list, and for each chain a
**shuffle word read in increasing clone-rank order**. A letter names the
visit receiving that clone. With visit sizes s_(ij), the number of possible
allocation words for chain i is

    l_i! / product_j s_(ij)!.

A multiset rank stores exactly that choice. For a fixed run list, a joint
mixed-radix rank uses ceil(log2 product_i(l_i!/product_j s_(ij)!)) bits,
which is optimal for fixed-length lossless encoding of those choices.
The ABBCCA ambiguity is repaired by one bit.

The operation laws are exact:

- **Doubling:** double weights and run sizes; repeat each shuffle letter
  twice. General k-fold dilation repeats it k times.
- **Reflection:** reverse the run list, reverse each shuffle and reverse
  its visit numbering; take the converse core.
- **Blue detection:** mirrored runs force a possible core anti-automorphism;
  test it and the corresponding shuffle reversals directly.
- **Ordinary line complement:** use T^op together with the N-1 path edges
  toggled back. Keep this sparse defect description, or recompute the
  actual canonical core and its allocation data.
- **Growth:** a clone insertion has a local shuffle update. Arbitrary
  attachments use the chain-splitting rule of section 2 and re-encode.

The gcd of all constant-run lengths inside the shuffle words is a unique
**marked-path dilation scale**, finer than the gcd of the chain lengths.
A node can have a doubling root while one of its interleaved path
presentations does not.

Two hostile examples prevent overcompression:

1. For cyclic core C3, weights (4,1,2), and word AABCCAA, the A-shuffles
   1221 and 2112 have identical pairwise inversion counts but complement
   targets with different degree multisets, even after converse merging.
2. Complementing the unique path presentation of TT_N gives a pair-free
   core of order N for every N>=5. Hence no bounded canonical-core size
   can survive ordinary line transport. A symbolic sparse-defect carrier
   still describes this family economically.

Full proofs, encoders and controls are in [allocation.md](allocation.md).
Exact semantic closure is established; a smaller quotient of allocations
preserving all requested continuation operations remains OPEN.

## 6. Connections, incoming work, and next targets

| Source -> target | Preserved | Lost or required data | Decisive test |
|---|---|---|---|
| Tournament -> weighted canonical core | full unmarked tournament class | path incidence needs runs and shuffles | ABBCCA |
| Core -> automorphism cycle index | every positive-weight orbit count | actual next-core response | rigid order-six pair above |
| Attachment -> (n,b) transition | labelled core-size distribution | target class and symmetry | source+C3 versus C3+sink |
| Pointed cycle index -> unpointed index | terms containing p_1 | fixed-point-free boundary | C3 coefficient 2p_3/3 |
| Path presentation -> sparse converse defect | exact complement tournament | reference chains may cease to be modules | TT_N -> full N-core |

During the session, origin/main acquired the pre-audit
[THM-4729 closed/open tiling note](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md)
and
[THM-4730 pure-blue alphabet claim](../../../01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md).
The elementary closed/open endpoint bit is compatible with this carrier.
The new alphabet claim is not a dependency of any result here. Its odd-order
argument used an invalid affine-only palindromicity test; the exact nonaffine
witness and surviving finite tests are retained in
[incoming_audit.py](incoming_audit.py). The original checker and theorem
status have been repaired: the classification remains unrefuted, but the
all-order odd-order proof and dependent consequences remain provisional.

The second incoming packet,
[THM-4731, the tribonacci partner](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md),
also connects directly to our strongest hostile. Its transitive-node
partner is exactly tau(TT_N), whose canonical core we prove has size N
for N>=5. Nevertheless, its Hamiltonian-path count has the inherited
three-term recurrence h(N)=h(N-1)+h(N-2)+h(N-3), with h(1)=h(2)=1,
h(3)=3. The March run-decomposition proof is recovered there. Thus
unbounded canonical-core order can coexist with a small symbolic counting
grammar. This gives a concrete positive control for the sparse-defect
direction, rather than a reason to abandon compression altogether.

The strongest remaining targets are now specific:

1. Find a structural description of fixed-point-free core automorphism
   terms, with repeated cycle lengths treated explicitly. Rooting alone
   provably leaves precisely this boundary.
2. Minimize allocation states relative to a specified operation language:
   mere inversion counts fail, while full multiset ranks work.
3. Develop a normalized sparse-defect syntax whose size is controlled by
   defect complexity rather than by the order of the expanded core.

The comparison with earlier arithmetic-certification work is methodological:
state reduction must preserve the next operation, and resetting a canonical
representation can reveal previously hidden complexity. No Collatz
conclusion is inferred from these tournament results.

Meta-patterns used: search the statement before the method; retain local
incidence through counting; correct the object before sharpening technique.
No new cross-domain meta-pattern is promoted from this packet.

## 7. Verification

The packet scripts use explicit exceptions, so checks survive python -O.
The growth script independently generates all classes through eight,
compares complete path cubes through seven, checks labelled orbit masses,
and tests every attachment to every labelled parent through five.
The allocation and weight scripts use separately organized controls.

See [VALIDATION.md](VALIDATION.md), [audit.py](audit.py), and
[manifest.json](manifest.json) for finite universes, reproduction commands,
retained output comparisons and normalized hashes. No theorem ID is reserved.
