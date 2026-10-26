# A complete tournament code, an exact necklace family, and the 81/27 boundary

**Status: PROVED elementary classification and counterexample; CITED canonical
decomposition and realization theorem; FINITE-EXACT independent controls.**
The necklace formulas and modular decomposition are classical mechanisms.
This note records an independent recovery of their exact interface for opposite
blocks, the full automorphism factor, and the cyclic-triple invariant repair.
Concurrent packets now supply the same emblem and necklace classification,
as well as stronger sharp-boundary results; their priority and scope are
recorded below.

## Inheritance and scope

**Concurrent integration, 2026-10-26.** Incoming commits `d761e4cf0f`
and `c9a8823436`, read after this packet's first proof/checker freeze, contain
[tournament composition classification](../../tournament_composition_classification_20261026/README.md)
and [cyclic-triple composition](../../cyclictriple_composition_20261026/README.md).
They already establish the decorated-tree interface, the 81/27 repair, and
the complete carousel fibres, period strata, groups and orbit masses for
`A=C3->1`. The former also proves the sharp `4p` threshold for the p-part of
automorphism order and uniqueness of the bad order-12 fibre; the latter
adds the multiplicity-retaining quotient spectrum and pointed growth.
Those stronger statements are inherited from their own proofs and audits,
not claimed as consequences of this packet's small finite census. The
present packet remains an independent recovery, with the explicit arbitrary
non-self-converse block formulation and full adjacency/isomorphism controls.
No external novelty claim is made for the shared classifications. The
incoming maintained repair of THM-4803 is the current correction lineage.

Closest proved mechanism: [Session56's strong-module and quotient interface](../../session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md),
especially its intrinsic recovery proof and cycle-index warning. The current
[A000568 family packet](../../a000568_intransitivity_families_20261026/FAMILIES.md)
supplies prime cycle-index data through order nine and family counts through ten.
[Session67](../../session67_metagraph_enumeration_information_20261026/README.md)
provides another complete mathematical enumeration carrier, with implementation
through eight; [the HN packet](../../hn_metagraph_information_20261026/README.md)
keeps ordinary and converse-merged classes separate.

The canonical hostile is Paley7: a cyclic tournament whose full automorphism
group has order 21, rather than seven. The corrected near miss is that
preserving a quotient's abstract automorphism group does not preserve the
stabilizer of its child-isomorphism labels. The least-used relevant sidecar
is that actual child-colouring, before passing to the cyclic-triple observable.
Our board is strong modules / prime quotient adjacency / child types /
cyclic periods / converse / information fibres. The anchor is exact taxonomy,
the niche is its symmetry strata, and the wildcard is the 81/27 failure.

“Prime tournament” means indecomposable with at least three vertices. It
does not mean that its number of vertices is a prime integer. A strongly
connected tournament need not itself be prime: its canonical top quotient
is prime, but its child modules can be nontrivial.

## 1. The full finite-tournament classification interface

A module M is a nonempty vertex set seen uniformly by every outside vertex.
A strong module overlaps no module except by containment. The strong modules
form a canonical rooted inclusion tree. Its children partition each node;
the quotient on those children is transitive or prime. This structural input
is inherited from Session56 and is stated as Theorem 16 in
[Boussaïri–Chergui–Ille–Zaidi](https://arxiv.org/html/1805.04607).

Define a recursive code C(T):

1. A singleton has one leaf code; the empty tournament has a separate code.
2. At a transitive node, retain the child codes in their unique source-to-sink
   order. These are the ordered strong components of that node.
3. At a prime node, retain the child-indexed quotient tournament together
   with its child codes, modulo **quotient isomorphisms carrying each child
   code to the matching code**. Equivalently, minimize the pair consisting of
   quotient adjacency and the ordered child-code list over all index permutations.

**Theorem 1 (complete code).** For all finite tournaments T and U,

    C(T)=C(U)  if and only if  T is isomorphic to U.

**Proof.** An isomorphism preserves modules, strongness, inclusion, and the
quotient arcs. Induction therefore carries each child to an equal child code;
it preserves the linear order or induces the permitted prime-quotient
isomorphism. Conversely, equal codes supply a quotient isomorphism and,
by induction, an isomorphism of each matched pair of children. Their union
preserves internal arcs by construction and cross arcs by the quotient.
Thus it is an isomorphism of the parent tournaments. The base cases finish
both directions. QED.

This is a complete algorithm on every finite input: brute force can compute
all modules and all quotient permutations. It is not a bounded catalogue of
prime tournaments or a new closed-form A000568 recurrence. An arbitrary
prime tournament appears as one root with n leaf children; its full adjacency
remains in the code. Bare tree shape, prime arities, or prime counts are not
the code. The unbounded prime registry and its coloured isomorphism actions
remain expensive. The necklace family below is a closed-form subfamily inside
this full classification, not all tournaments.

### Exact recovery condition and the coloured automorphism group

For a specified substitution Q[B_i], require that its displayed blocks be
exactly its maximal proper strong modules. Under that intrinsic condition,
two such substitutions are isomorphic **iff** there is a quotient
isomorphism pi with B_i isomorphic to B'_(pi(i)) for every i. Without it,
matching decorated quotients still suffices, but necessity is unjustified.
For example the transitive quotients with transitive block sizes (1,2,1)
and (2,1,1) both give the same four-vertex transitive tournament, although
those ordered block-size lists differ. Those cuts were not canonical.

For a prime Q of order at least three and arbitrary nonempty B_i, recovery
holds automatically. Here is the inherited elementary proof. If a module
M meets at least two blocks, its touched index set is a module of Q.
Primeness forces it to meet every block. If M omitted a vertex from block i,
that vertex's uniform relation to M would make i a source or sink of Q,
impossible for a prime of order at least three. Thus M is the whole
tournament. Every proper module lies in one block, so the blocks are exactly
the maximal proper modules and are strong. The B_i themselves need not be
strongly connected.

Consequently, putting K=Stab_(Aut Q)(([B_i])_i), one has

    Aut(Q[B_i]) ≅ (product_i Aut(B_i)) semidirect K,
    |Aut(Q[B_i])| = (product_i |Aut(B_i)|) |K|.       (1)

The kernel is the product of the internal groups; the image is exactly K.
Choosing a common labelled model for each child type supplies a splitting.
At a linear node K is trivial because the child order is intrinsic. At a
prime node it is generally smaller than Aut(Q). Actual isomorphism types,
not merely child sizes or child cyclic-triple structures, colour this action.

## 2. The odd carousel has exactly cyclic symmetry

For n=2m+1>=3 define R_n on Z/n by

    i -> j  iff  1 <= j-i mod n <= m.

**Lemma 2.** R_n is prime and Aut(R_n) is precisely its n rotations.

**Proof of primeness.** If two vertices a,b belong to a module and the
forward distance from a to b is at most m, every vertex strictly between
them on that short forward arc sees a and b oppositely, so must also lie
in the module. Every pair has one such orientation. Thus any module of
size at least two contains consecutive vertices; rotate them to 0,1.
Vertex m+1 sees 0 and 1 oppositely and must join the module. The short arcs
1 to m+1 and m+1 to 0 then force every remaining vertex into it. QED.

**Proof of the group assertion.** The out-neighbourhood of i is the transitive
tournament on i+1,...,i+m, whose unique source is i+1. Hence the successor
map i->i+1 is intrinsically recoverable from the adjacency relation.
Every automorphism commutes with it and is a rotation; every rotation is
indeed an automorphism. This also works for m=1. QED.

Thus R_9 and R_15 are prime tournaments even though 9 and 15 are composite
integers. The assertion is specific to this carousel. Paley7, oriented by
differences {1,2,4}, also has a cyclic regular action, but multiplication
by two is an additional automorphism. The full small oracle finds order21.
Its markings with ones at {0,1} and {0,2} yield isomorphic block tournaments
but are not rotationally equivalent binary words. For a general prime cyclic
base Q, replace necklaces by the full Aut(Q)-orbits of child labels.

## 3. Opposite blocks are classified exactly by necklaces

Fix a tournament B on b vertices such that B is **not** isomorphic to its
converse B*. Put a=|Aut(B)|=|Aut(B*)|, and on a common labelled vertex set
set B_0=B, B_1=B*. For a binary word w of odd length n>=3 define

    T(w)=R_n[B_(w_0),...,B_(w_(n-1))].

**Theorem 3 (necessary and sufficient classification).**

    T(w) ≅ T(v)  iff  v is a cyclic rotation of w.       (2)

If d is the least positive period of w, then d divides n and

    Aut(T(w)) ≅ Aut(B)^n semidirect C_(n/d),
    |Aut(T(w))|=a^n n/d.                               (3)

**Proof.** Lemma 2 and the recovery lemma make the n blocks intrinsic.
Every isomorphism therefore induces a rotation of R_n. The two block types
are nonisomorphic, so this rotation must carry the bit word to the other
word. Conversely such a rotation lifts using the identical block models.
The word stabilizer consists of rotations by multiples of d, so has order
n/d. Formula (1) gives (3); the cyclic factor acts by permuting the n
internal group coordinates. QED.

It is the *quotient* automorphism group that has order n/d. Forgetting the
factor a^n is wrong unless B is asymmetric. Conversely, omitting B!=B*
is wrong even with an intrinsic prime quotient: for the self-converse
two-vertex block all eight words over R_3 give isomorphic tournaments.

### Divisors, Möbius inversion, and exact symmetry strata

Let P(d) be the number of binary words of length d with least period d.
Every word decomposes by its unique least period, so

    2^d = sum_(e|d) P(e),
    P(d) = sum_(e|d) mu(e) 2^(d/e).                    (4)

The words of length n and least period d are repetitions of these words.
Their rotation orbits have size d. Thus the number of tournament classes
with least period d, equivalently automorphism order a^n n/d, is exactly

    P(d)/d,                    d|n.                   (5)

Summing gives the binary-necklace formula

    N(n)=sum_(d|n) P(d)/d
        =(1/n) sum_(e|n) phi(e) 2^(n/e).             (6)

The second expression also follows by counting words fixed by each rotation.
For prime integer n=p, two constant classes have order a^p p; the remaining
(2^p-2)/p classes have order a^p. Composite odd n has intermediate strata:

| n | Least period d : number of classes | Total N(n) |
|---:|---|---:|
| 3 | 1:2, 3:2 | 4 |
| 5 | 1:2, 5:6 | 8 |
| 7 | 1:2, 7:18 | 20 |
| 9 | 1:2, 3:2, 9:56 | 60 |
| 15 | 1:2, 3:2, 5:6, 15:2182 | 2192 |

### Converse: complemented reflection, not ordinary bracelets

The map i->-i identifies R_n* with R_n. Therefore T(w)* has necklace
represented by w'_i=1-w_(-i). The appropriate converse-merged action is
generated by rotations and **reflection combined with bit complementation**.
It is not ordinary reflection of the word. No member is self-converse:
isomorphism would equate the number of 1 bits with n minus that number,
which is impossible for odd n. Thus the converse-merged class count is
exactly N(n)/2. Individual block converse and total tournament converse
are different operations, even though both enter this formula.

## The 81 versus 27 boundary

Take B=A=C3 -> 1, with vertices 0,1,2 forming 0->1->2->0 and all beating
vertex3. Its out-neighbour bit masks are

    A=(10,12,9,0),       A*=(4,1,2,7).

Their degree multisets are {0,2,2,2} and {1,1,1,3}, so they are nonisomorphic;
both have automorphism group C3. Order four is the smallest possible size of
a non-self-converse block, since all tournaments on at most three vertices
are transitive or C3 and are self-converse.

Set U=R_3[A,A,A] and V=R_3[A,A,A*]. The blocks are maximal proper strong
modules even though A itself is not strongly connected. Formula (3) gives

    |Aut(U)|=3^3·3=81,       |Aut(V)|=3^3=27.          (7)

Nevertheless **C3(U)=C3(V) literally on the same twelve labelled vertices**.
Inside a block, converse preserves the cyclic triples. A triple meeting
exactly two blocks is never cyclic, because its singleton vertex sees the
other block uniformly. A transversal triple depends only on the unchanged
outer R_3. There are 3 internal and 4^3 transversal cyclic triples, total67.
More generally all T(w) for a fixed n and opposite block B share the same
labelled cyclic-triple hypergraph, while (3) realizes every divisor-period
symmetry stratum in that one hypergraph fibre.

Literal twelve-vertex rows, with bit j meaning an arc to vertex j, are:

    U=(250,252,249,240,4000,4032,3984,3840,2575,3087,2319,15),
    V=(250,252,249,240,4000,4032,3984,3840,1039,271,527,1807).

The independent adjacency-only backtracker finds exactly81 and27
automorphisms and no isomorphism. Exhausting all nonempty vertex subsets
finds maximal proper modules with masks15,240,3840 in both. The subset-DP
Hamiltonian-path counts are both462915. The counterexample refutes the
automorphism invariant, not the separate Hamiltonian-path invariant.
No global minimality claim for an arbitrary counterexample order is made.

### Failed implication and exact repair

The pre-repair claim in
[THM-4803](../../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md)
and [INTRANSITIVITY, section3](../../a000568_intransitivity_families_20261026/INTRANSITIVITY.md)
was that |Aut(T)| is a function of C3(T). Its displayed substitution formula
was correct, but the inference from Aut(Q*)=Aut(Q) omitted ancestor effects.
Changing a descendant can change whether two children are isomorphic, so
the coloured stabilizer K at an ancestor changes. In (7), an internal
linear-node swap A->A* destroys the outer threefold rotation. Local moves
are not confined to the automorphism factor of the node where they occur.

The strongest survivor is (1), with every child isomorphism label retained.
Global converse preserves |Aut|; independent module moves need not.
Finite agreement through order nine remains a finite observation.
The reported cyclic-triple realization theorem, finite class counts and
Hamiltonian-path argument are not refuted by this example.

Likewise C3 determines strong-component vertex sets and their cyclic-triple
structures, not their tournament isomorphism types. The strong tournament
R_3[A,1,1] and its converse have the same labelled C3 and are nonisomorphic;
the oracle checks this six-vertex control. Thus even the only strong
component can change its type.

### A complete realization fibre in this explicit family

For A=C3->1, the necklace family is stronger than a sampled subfamily.
The [cited realization theorem](https://arxiv.org/html/1805.04607),
Theorems10 and24, Proposition51 and Corollary57, says realizations share
strong modules and are obtained by independently choosing each quotient
realization. At prime nodes the choice is quotient versus converse; at
linear nodes it is a total order of the children.

Here the root is R_n, each of its n children has one linear node joining
C3 and a singleton, and each C3 has its own prime node. The root can be
R_n or R_n*, each linear node gives A or A*, and each internal C3 direction
gives the same child isomorphism type. Root reflection normalizes R_n*
back to R_n. Hence **every realization of this fixed C3 hypergraph is
isomorphic to exactly one of the N(n) necklace classes above**. There are
2^(n+1)·2^n=2^(2n+1) labelled realizations; that is not the unlabelled
necklace count. This deduction uses the cited realization theorem; the
elementary 81/27 refutation itself does not depend on its hard direction.

## Exact verification and boundaries

From this directory:

```sh
python3 -B check_taxonomy.py > /tmp/taxonomy.normal.json
python3 -B -O check_taxonomy.py > /tmp/taxonomy.optimized.json
cmp /tmp/taxonomy.normal.json /tmp/taxonomy.optimized.json
cmp /tmp/taxonomy.normal.json receipt.json
```

The standard-library checker imports no producer implementation. It checks:

- all1,099 labelled tournaments of orders1..5: the recursive strong-module
  code is compared in both directions against full n! adjacency canonization;
  the class counts are1,1,2,4,12;
- all8 binary words at n=3 and all32 at n=5, inflated to12 and20 vertices:
  full adjacency-only isomorphism for all36 and528 unordered word pairs
  including the diagonals, and exact automorphism counts for every word;
- the literal81/27 witness, its67 triples, both module universes and both
  Hamiltonian-path counts; the strong-six converse hostile;
- all necklaces at odd n=3,5,7,9,11,13,15 versus the divisor/Möbius formulas;
  carousel modules through13 and full automorphism groups through9;
- self-converse block collapse, arbitrary noncanonical linear cuts, and
  the Paley7 multiplier identifying two nonrotation words.

The proof handles every odd n and every finite non-self-converse B; the
finite controls neither substitute for that proof nor extend their own
enumeration scope. The full tree interface handles every finite tournament,
but prime quotient adjacency remains essential. No black-core connectivity,
new all-order A000568 numerical range, or inexpensive prime-registry formula
is asserted.
