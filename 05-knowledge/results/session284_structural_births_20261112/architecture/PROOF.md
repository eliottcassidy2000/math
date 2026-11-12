# What a comprehensive recursive grammar can and cannot mean

**Status: PROVED elementary synthesis and corollaries; inherited algorithms are credited, not claimed new.**

## 1. The all-size exact rule already exists

Let G_n be the set of tournament isomorphism classes of order n. For each representative Q, adjoin x with binary attachment z in {0,1}^n, where x→i iff z_i=1. Then

    G_(n+1) = {can(Q+x_z): [Q] in G_n, z in {0,1}^n}.

Here can means the lexicographically least adjacency string over all vertex permutations. This is a complete, terminating, finite-syntax rule valid at every order, with no empirical extrapolation. Every child has a deletion card and therefore occurs; canonicalization identifies precisely its isomorphism class. It is not an efficient enumeration promise and does not by itself explain a child's properties. The rule is classical extension-and-canonicalization, and THM-4169 gives the prime-parent specialization with additional exact port information.

If the new vertex remains marked, the correct quotient is only by Aut(Q). The number of marked child types from Q is

    E(Q) = (1/|Aut(Q)|) sum_(g in Aut(Q)) 2^(c(g)),

by Burnside: an attachment fixed by g is constant on its cycles. Forgetting the mark introduces additional identifications between distinct attachment orbits, and between distinct parents. For a child T, the number of attachment orbits above Q equals the number of Aut(T)-orbits of vertices v with T-v isomorphic to Q. This identity is a bijection of rooted isotypes, not a claim that all vertices in such a deletion fibre form one orbit.

Consequently, if D(Q) is the number of distinct unmarked children,

    E(Q)/(n+1) <= D(Q) <= E(Q).

These formulas retain precisely the sidecar needed to turn a marked grammar into the actual extension metagraph. They distinguish the metagraph from its vertex count A000568(n).

## 2. Exponential branching inside a concise explicit family

The births lane proves Aut(W_n)=1 for every odd n>=5. The inherited prime augmentation classification, [THM-4169](../../../../01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md), says exactly2+2n of its2^n attachments are nonprime: two source/sink patterns and two clones for each old vertex. All other patterns are prime.

Thus the number of distinct prime unmarked children of W_n is at least

    ceil((2^n-2n-2)/(n+1)).

Proof: since the parent is rigid, its attachment patterns are distinct marked isotypes. Each unmarked child accounts for at most n+1 such types. No probabilistic or asymptotic assumption is involved. This corollary combines the inherited prime-parent classifier with the explicit rigidity proof; neither ingredient is advertised as a newly discovered universal generator.

There is therefore no eventual bounded number of isomorphism-class successors even along this one easily described family. Yet all these successors have the concise exact description 'every binary word except the explicitly listed source/sink/clone words, modulo child isomorphism'. Exponential branching and finite parametric rules are compatible.

## 3. Compactness must name what is compressed

The labelled universe has2^(n(n-1)/2) tournaments. A class contains at most n! labellings, hence

    A000568(n) >= 2^(n(n-1)/2)/n!.

Any injective fixed-length binary address system for all classes requires at least ceil(log2 A000568(n)) bits. In particular it requires at least ceil(n(n-1)/2-log2(n!)) bits. This standard counting bound does not rule out a short generator, symbolic family library, fast individual invariant, or unexpected theorem. It rules out uniformly short *individual lossless addresses* for all large tournaments. A scalar count, a generating function, a parametric grammar, a list of classes, and a list of their transitions have different information obligations.

We now have three compatible levels rather than a false choice between exhaustive tables and vague patterns:

1. Universal exact extension-and-canonicalization: complete at every order, with potentially huge state space.
2. Signature-run and feedback-one family grammars: compressed exact rules over unbounded size, retaining required lineage/isomorphism sidecars.
3. Observable recurrences: exact only for specified profiles or laws; not automatically sufficient for child isotypes or path-count realization.

## 4. Relation to Chan–Pak

A complete generation grammar guarantees that every tournament can be generated. Chan–Pak instead requires that every prescribed allowed odd integer is the Hamiltonian path count of some generated object. No surjectivity follows merely from completeness of the object grammar. The new recursive endpoint protection and family inventories identify structural constraints that any target selector must respect; they do not yet produce such a selector.
