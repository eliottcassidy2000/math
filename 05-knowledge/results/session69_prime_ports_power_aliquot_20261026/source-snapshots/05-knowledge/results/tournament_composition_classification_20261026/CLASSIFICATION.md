# Classify the decorated object and its operations

2026-10-26. **PROVED** elementary structural consequences below, using **CITED**
canonical modular decomposition and cyclic-triple realization results. This is
an exact recursive classification framework, not a completed atlas of all prime
tournaments or metagraphs. The packet's independent audits are in `AUDIT.md`.

## 1. The objects and the equivalences

An ordinary tournament is a finite complete oriented graph, with isomorphisms
preserving arcs. A path presentation is a pair `(T,P)` where `P` is a directed
Hamiltonian path, written in increasing position order. The free arcs are the
nonconsecutive pairs. The involution `tau` reverses exactly these free arcs;
it retains the path. It is not global converse and not a single-arc flip.

We keep four levels separate:

1. the labelled tournament;
2. the ordinary isomorphism class;
3. a cyclic-triple hypergraph, labelled or up to hypergraph isomorphism;
4. the path presentations and the operations between them.

An invariant at one level need not determine an operation at another. In
particular `H(T)/|Aut(T)|` counts path presentations up to isomorphism: the
automorphism action on directed Hamiltonian paths is free, since fixing an
ordered path fixes every vertex. This count is part of the metagraph fibre,
not a consequence of the cyclic-triple hypergraph alone.

## 2. The exact recursive classifier

The canonical tree of strong modules has a leaf at each vertex. At a linear
node retain the **ordered** list of its strong components. At a prime node
retain its prime quotient `Q` and the isomorphism type of each child tournament
as a vertex colour of `Q`. Quotient relabellings must preserve this colouring.
Use maximal proper **strong** modules; arbitrary modules need not give the
canonical partition.

**Decorated-tree lemma (PROVED, relative to canonical decomposition).** Two
tournaments are isomorphic if and only if their rooted decomposition trees are
isomorphic preserving linear orders and the child-coloured prime quotients.

Indeed, any tournament isomorphism preserves the defining module predicate
and hence the canonical partition at each node. It induces the requisite
coloured quotient isomorphism and child isomorphisms. Conversely these maps
assemble to an arc-preserving bijection: arcs within children are handled by
the child maps, and arcs between children by the quotient. Induct on size.

At a prime node put

    K = {g in Aut(Q) : T_g(i) is isomorphic to T_i for every i}.

Then

    Aut(Q[T_1,...,T_k]) ≅ (product_i Aut(T_i)) ⋊ K.              (2.1)

The kernel consists of independent internal automorphisms. To see splitting,
choose a single model for each child isomorphism type and coherent
identifications with its copies; each element of `K` then permutes those
copies through the chosen identifications. The split is not canonical. At a
linear node the quotient is rigid, so the group is the product of the child
groups. For an arbitrary, noncanonical substitution, (2.1) describes the
partition-preserving subgroup; one must first prove that the partition is
intrinsic before identifying it with the entire automorphism group.

This is the useful part of the simple-group analogy: distinguish irreducible
objects, composition data, and the induced actions. Prime tournaments play
the role of indecomposable quotient objects. Their automorphism groups are
not thereby simple. An involution in a tournament automorphism group would
interchange some pair and reverse its arc, so the group has odd order. No
classification of finite simple groups is needed for this packet.

**Failure boundary.** The multiset of prime factors, the undecorated tree,
the list of child group orders, and the total group order are each coarser
than the recursive classifier. The witness below changes the colouring
stabilizer `K` while leaving the internal group orders unchanged.

## 3. The observer fibre is an action, not only a cardinality

Write `C(T)` for the labelled cyclic-triple hypergraph and `R(C)` for its
tournament realizations on the same vertex set. We use the following cited
inputs from [Boussaïri–Chergui–Ille–Zaidi](https://arxiv.org/abs/1805.04607):
the shared strong modules (Theorem 10), prime realization up to converse
(Theorem 24, citing Boussaïri–Ille–Lopez–Thomassé), and the nodewise realization
bijection and count (Proposition 51, Corollary 57). They give

    |R(C)| = 2^(number of prime nodes) × product_(linear nodes v) arity(v)!.

This is a labelled count. It does not say the realizations have equal
automorphism groups, equal ordinary orbit sizes, or equal metagraph fibres.

**Lift-action lemma (PROVED).** Set `G=Aut(C)`, acting by relabelling on
`R(C)`. Its orbits are exactly the ordinary tournament isomorphism classes
realizing the hypergraph isomorphism type of `C`. For each lift `T`,
`Stab_G(T)=Aut(T)`. Therefore

    sum_[T] 1/|Aut(T)| = |R(C)|/|Aut(C)|,                     (3.1)
    number of ordinary lift classes = (1/|G|) sum_g |Fix_R(C)(g)|. (3.2)

Proof: any isomorphism of two lifts of this labelled `C` preserves `C`, so
belongs to `G`; conversely a `G`-orbit is an ordinary isomorphism class. The
stabilizer identity follows from the definition. Summing the orbit sizes
`|G|/|Aut(T)|` proves (3.1); counting pairs `(g,T)` with `gT=T` proves (3.2).
Every tournament realizing an isomorphic hypergraph can first be transported
to the fixed labelled `C`, so none are omitted. No equal-stabilizer premise
has been used.

For the 12-vertex example let `A=C3→1` and let `B=A*`, on the same four
labels. Compare `C3[A,A,A]` and `C3[A,A,B]`. The former permits the cyclic
permutation of its three children and the latter does not. The full observer
fibre has the following four ordinary classes:

| Type | Aut order | Orbit size inside the labelled fibre | H | H/Aut |
|---|---:|---:|---:|---:|
| AAA | 81 | 16 | 462915 | 5715 |
| AAB | 27 | 48 | 462915 | 17145 |
| ABB | 27 | 48 | 462915 | 17145 |
| BBB | 81 | 16 | 462915 | 5715 |

Here `|Aut(C)|=1296`, `|R(C)|=128`, and the mass is
`2/81+2/27=8/81=128/1296`. The common number of cyclic triples is 67.
The [fibre proof and exact checker](fibre/RESULT.md) establish the witness,
the sharp minimum 12, and the general sharp threshold `4p` for the `p`-part
of the automorphism order.

## 4. Size change requires operation data

An atlas should distinguish at least these operations:

| Operation | Carrier retained | Exact content | Lost if reduced to scalars |
|---|---|---|---|
| Modular substitution | canonical quotient with child colours | classifier and (2.1) | repeated-child stabilizer |
| Passage to cyclic triples | strong-module tree and quotient duality classes | labelled lift set and action (3.1) | orientation choices and lift stabilizers |
| Ordered join | two path carriers with their boundary positions | Cartesian product of legal-relocation graphs | endpoint and half-edge identity |
| Free-arc complement `tau` | fixed directed path | conjugates join to an associative hinge product | cannot descend from ordinary edge endpoints alone |
| `p`-adic symmetry audit | digit blocks and Sylow action | exact valuation maxima and carry defects | full isomorphism type and cycle index |

The [metagraph lane](metagraph/RESULT.md) gives the join/product and hinge
laws and a smallest order-five failure of composition from ordinary edge
endpoint classes. The [arithmetic lane](arithmetic/RESULT.md) classifies
tournaments containing a Sylow `p`-subgroup of `S_n`, proves the all-order
valuation maxima, and tracks composition deficits by base-`p` carries.

Neither scalar automorphism order nor its prime factorization is an adequate
substitute for the actual permutation action. The randomly recovered
[THM-4764](../../../01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md)
already distinguishes two order-25 actions by their subset orbit polynomials.
That is an exact downstream consumer of the cycle-index data retained here.

## 5. Research boundary and next decisive objects

**OPEN:** a concise structural atlas of all prime tournament species; full
metagraph composition under general substitution; a smallest finite rooted
port object sufficient for that composition; and any new implication for
LRC(14). Existing exact enumeration through a specified order remains finite
evidence, not a proof of an all-order structural regularity.

The immediate targets are concrete: classify the rooted port types that
separate the order-five hostile pair; determine which quotient actions
appear at the sharp `4p` boundary beyond the cyclic examples; and ask which
cycle-index/port enhancements are closed under general substitution.
The canonical decorated tree already classifies ordinary tournaments, but
an effective prime atlas and a minimal operation-compatible quotient are
separate problems.
