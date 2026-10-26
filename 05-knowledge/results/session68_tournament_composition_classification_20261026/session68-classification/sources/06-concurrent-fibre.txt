# Cyclic triples, decorated symmetries, and the sharp prime thresholds

2026-10-26. **PROVED**, conditional only on the explicitly **CITED** classical
modular-decomposition and C3-realization theorems below; **FINITE-EXACT**
controls in `check.py`. Independent audit is recorded in `AUDIT.md`.
No historical-priority claim is made for modular decomposition or orbit-stabilizer.

## Inheritance and the error this repairs

The inherited [C3-fibre theorem](../../a000568_intransitivity_families_20261026/INTRANSITIVITY.md)
and [modular interface](../../session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md)
retain a canonical rooted tree, its child ports, prime quotients up to
duality, and orders at linear nodes. The repository claim that `|Aut(T)|`
is constant on every C3-fibre is false. Its induction held the child
isomorphism colours fixed while allowing internal tree moves that change
those colours. The exact subgroup of quotient automorphisms preserving
child colours can therefore shrink.

The hard inherited input is precise: realizations of a fixed labelled
cyclic-triple hypergraph have the same strong modules; at a prime node its
quotient is determined up to duality, and at a linear node the child order
is free. This is **CITED** from Boussaïri–Chergui–Ille–Zaidi,
[*3-uniform hypergraphs: modular decomposition and realization by tournaments*](https://arxiv.org/pdf/1805.04607),
Theorem10, Theorem24 and Corollary57 (printed pages4,7,32 in the retrieved
33-page version), primary text checked 2026-10-26. The all-order proofs below
use that structural input; the numerical tests are not substitutes for it.

Anchor: recover the 81/27 witness and determine its first possible size.
Niche: prime-adic symmetry thresholds under composition. Wildcard: distinguish
the two different meanings of chirality. The concept board is the marked
modular tree / child isomorphism colours / quotient stabilizer / path-cover
profile / prime divisibility. Every result below names which coordinates
survive the C3 map and which do not.

## F1. The exact automorphism recursion

At a prime modular node write `T=Q[T_i : i in V(Q)]`, using the canonical
child modules, and give quotient vertex i the colour `[T_i]`, its ordinary
tournament isomorphism class. Put

    Gamma(Q;[T_i])={sigma in Aut(Q): T_sigma(i) is isomorphic to T_i for every i}.

Then there is an exact sequence

    1 -> product_i Aut(T_i) -> Aut(T) -> Gamma(Q;[T_i]) -> 1,

and in particular

    |Aut(T)| = |Gamma(Q;[T_i])| product_i |Aut(T_i)|.       (F1)

**Proof.** Automorphisms preserve canonical child modules. The kernel acts
independently inside each child, so is the stated product. The image preserves
the quotient and the child isomorphism colours. Conversely, any such quotient
permutation lifts by independently choosing an isomorphism between every
source child and its target child; the uniform cross arcs make the combined
map an automorphism. This proves exactness and the order formula. At a linear
node the transitive quotient has no nontrivial automorphism, so the formula
reduces to the product of the child automorphism orders. QED.

Dualizing the quotient itself does not change `Aut(Q)`. But moving inside
one child can change its colour without changing its C3-hypergraph. Thus
the base group and even the factors' orders can stay fixed while Gamma
changes. The missing coordinate is the equality pattern of the **oriented
child isomorphism classes**, not just their sizes or C3-types.

## F2. The sharp p-adic onset is 4p

For every odd prime p and tournaments T,T' on n vertices with the same
labelled cyclic triples,

    n < 4p  =>  v_p(|Aut(T)|)=v_p(|Aut(T')|).               (F2)

The bound is sharp for every odd prime p: there is a pair on exactly 4p
vertices with identical labelled cyclic triples and automorphism-order
ratio p. Consequently every prime divisor of the reduced rational ratio
`|Aut(T)|/|Aut(T')|` is odd and at most n/4.

**Proof of stability.** A tournament automorphism group has odd order:
an involution would transpose some pair and reverse its arc, so no involution
exists; Cauchy's theorem excludes even order. For fixed p, induct on n<4p
using the common canonical modular tree. Linear nodes use the product rule
and induction. At a prime node, the quotient may be replaced by its dual,
which has the identical automorphism group on the same child ports.

A C3-fibre on at most three vertices has just one ordinary isomorphism type.
Therefore only children of size at least four can change their isomorphism
colour. Let L be these large-child ports; there are fewer than p of them.
Let B be the subgroup of Aut(Q) preserving every child size and the oriented
isomorphism colours of all small children. The group B is unchanged by any
C3-fibre moves below the node. Both refined colour stabilizers Gamma and
Gamma' lie inside B. Every p-subgroup of B fixes each port of L: its orbits
on a set of size less than p have size one. Thus every such subgroup lies
inside both Gamma and Gamma'. In particular a Sylow p-subgroup of B lies
in each, so both have the same p-part as B. Induction on the child factors
in F1 proves F2. The p=2 valuation is zero at every order. QED.

**Sharp construction.** Let `A=C3 -> 1` and `B=1 -> C3`; these are converse,
nonisomorphic four-tournaments with the same labelled C3-hypergraph and
automorphism order three. Choose any circulant tournament Q on a prime
p vertices; for example i beats j exactly when
`j-i mod p` belongs to `{1,...,(p-1)/2}`. Its translations act transitively.
Such Q is modular-prime: vertex-transitivity forces strong connectivity,
and its canonical maximal proper strong modules would otherwise be a transitive
block partition of equal block size dividing p. They must therefore be
singletons. Its automorphism group acts transitively on p vertices and has
p-adic valuation one, since it contains a p-cycle and embeds in S_p.

Compare

    U_p=Q[A,...,A],
    V_p=Q[B,A,...,A].                                    (F3)

Changing the first child preserves the labelled cyclic triples: a triple
inside the child retains its cycle status; one with two vertices inside
is transitive; one meeting three blocks depends only on Q. All children
are canonical maximal proper strong modules by primeness of Q. Formula F1 gives

    |Aut(U_p)|=3^p |Aut(Q)|,
    |Aut(V_p)|=3^p |Stab_Aut(Q)(0)|=|Aut(U_p)|/p.

This also exhibits the information loss as the coset space
`Aut(Q)/Stab(0)` of size p: a changed child marks one previously unmarked
quotient port. The example uses no unproved assertion about the full
automorphism group of a circulant tournament. QED.

The prime thresholds are therefore 12,20,28,44,... for p=3,5,7,11,... .
They concern a *change* in p-adic symmetry inside a fixed C3-fibre, not
the first occurrence of p as a tournament automorphism divisor. For example
C3 already has automorphism order three at order three.

## F3. Order twelve is the first failure, and its bad fibre is unique

At n<=11 the automorphism order is C3-invariant. At n=12, up to hypergraph
isomorphism there is exactly one realizable cyclic-triple hypergraph whose
tournament realizations have different automorphism orders. It is

    H=C3(C3[A,A,A]).                                      (F4)

**Proof.** F2 excludes every odd p below order12, hence excludes any order
change at n<=11. Alternatively, an odd-order quotient group fixes both
large-child ports when there are at most two. If variation occurs at12,
child automorphism orders are invariant by the smaller-order result; a
nontrivial change in Gamma requires at least three children of size at
least four. The total size forces exactly three children, all of size four,
and the prime quotient is C3. Its group is C3, so its colour stabilizer has
order three if all three child colours agree and one otherwise.

There are exactly four isomorphism classes of four-tournaments. Their C3
structures are: empty (the transitive class), one triple plus an isolated
vertex (A and B), and the two-triple structure (the unique strong class).
This tiny classification can be checked by all 64 labelled orientations;
the checker does so independently using all 24 permutations. Thus child
colours can change only in the A/B fibre, and for a stabilizer change all
three children must be of that hypergraph type. This proves uniqueness of
F4. The construction in F3 at p=3 attains it. QED.

Here is the full lift classification of F4; juxtaposition lists the three
four-vertex children around the directed triangle.

| Ordinary tournament class | Aut order | Orbit under Aut(H) | H(T) | tc(T)=H/Aut |
|---|---:|---:|---:|---:|
| AAA | 81 | 16 | 462915 | 5715 |
| AAB | 27 | 48 | 462915 | 17145 |
| ABB | 27 | 48 | 462915 | 17145 |
| BBB | 81 | 16 | 462915 | 5715 |

There are 128 labelled lifts: choose the two root orientations, independently
orient each of the three internal triangles, and independently choose A/B
in each child. This is also `2^4 (2!)^3` from the inherited fibre formula.
The root orientations are isomorphic after permuting blocks; for each
number of B children all placements are isomorphic. Different numbers of
B children are not isomorphic because canonical modules and their oriented
types are preserved. Hence the four listed classes are exhaustive.

The hypergraph has 67 triples: one inside each triangle and all 4^3 triples
meeting each block once. Cross-block pairs have codegree four; pairs inside
a block have codegree zero or one. Thus pair codegree recovers its three
blocks intrinsically. The isolated vertex inside each four-block is also
intrinsic; the three remaining vertices can be permuted arbitrarily.
Consequently

    Aut(H) = S3^3 semidirect S3,    |Aut(H)|=1296.          (F5)

Direct orbit-stabilizer gives the table's orbit sizes and the exact mass

    sum_[T in R(H)/Aut(H)] 1/|Aut(T)| = 128/1296 = 8/81.

This distinguishes class counting, labelled lift counting and symmetry
weighting. In particular C3 does not determine the fixed-Hamilton-path
metagraph multiplicity tc, despite determining H. The Hamilton counts in
the table are FINITE-EXACT subset dynamic-programming counts, independently
recomputed for all128 labelled lifts; the automorphism orders also have a
direct degree/adjacency backtracking check independent of modules.

## F4. The stronger path-count survivor

For k>=1 let P_k(T) count ordered k-tuples of nonempty directed paths whose
vertex sets partition V(T). Then the **entire ordered path-cover profile**
P_k(T), not only H(T)=P_1(T), is invariant under C3-fibre moves.

**Proof by explicit bijections.** In a substitution, break each directed
path into its maximal runs inside child modules. At a linear node each path
visits any child at most once and the runs occur in the child order. After
permuting that order, permute the runs within each path into the new order;
keep every run's internal orientation. The result is a path on the same
vertex set, and the inverse is the reverse reordering.

At a prime node, dualize Q and reverse the sequence of runs in every path,
while keeping the internal orientation of every run. Each old transition
i->j becomes the valid new transition j->i. Again this is an involutive
bijection preserving the number and order of the covering paths and each
path's vertex set. Repeated child visits cause no issue: adjacent runs have
distinct child labels before and after reversal.

For a move below the root, collect all runs through the affected module,
in the order of the outer covering paths and then their positions. They
form an ordered path cover of the module. Apply the corresponding bijection
and put the new paths back at those run positions. Every outside vertex is
uniform on the module, so all entering and exiting arcs remain valid.
This lifts the local bijection through every ancestor. The cited tree-move
description therefore proves the full C3 invariance. QED.

The old H-invariance conclusion survives; its proof is made explicit at
the composition level needed by the induction. The checker independently
computes the full profile for the four representatives, using subset
Hamiltonian counts and a minimum-vertex set-partition recurrence. Equality
of H does not constrain the changed quotient colour stabilizer in F1.

## F5. Two adjacent type corrections

1. A C3-hypergraph determines the **vertex sets** of strong components and
   their C3-types, not their ordinary tournament isomorphism types. The
   strong twelve-tournaments AAA and AAB already refute the latter statement:
   each has one strong component, their C3 structures agree, and their
   automorphism orders prove that their component types differ.
2. The split of an S_n-isomorphism orbit into two A_n-orbits is not the
   chirality bit in a prime labelled C3-fibre. For the regular five-tournament
   with steps {1,2} modulo5, the anti-automorphism i->-i is the even
   permutation `(1 4)(2 3)`. Thus T and T* belong to the same A5-orbit while
   the labelled C3-fibre is `{T,T*}`. Its S5-orbit still splits because all
   tournament automorphisms have odd order and are even permutations.

Both distinctions have literal controls in the receipt. Global converse
does preserve automorphism order; independently changing oriented child
types is the different operation responsible for F3.

## Reproduction and connection contract

```sh
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
```

Both receipts must agree byte-for-byte. The standard-library checker
exhausts order4 orientations, all128 order12 lifts, all1296 predicted
hypergraph automorphisms, four independent tournament-automorphism searches,
all128 Hamilton counts, four complete path-cover profiles, and sharp prime
construction controls at p=3,5,7,11,13. For these primes it checks quotient
primality by all subsets, its full automorphism and point-stabilizer orders
by backtracking, and literal equality of the two expanded C3 hypergraphs.

Source: a labelled tournament and its canonical decorated modular tree.
Target: C3 structure and its lift fibre. Map: forget child orders and prime
orientation choices. Preserved: strong-module sets, C3 quotient types and
all path-cover counts. Lost: oriented child isomorphism colours and therefore
their quotient stabilizer. Sidecar: the actual colour action Gamma, not only
the quotient group or child automorphism orders. Decisive test: 81 versus27;
number-theoretic generalization: sharp4p onset for every odd prime.

The classification is complete for symmetry-order variation through12 and
for the first p-adic onset at every prime. It is not a classification of all
prime tournaments, all metagraph moves, or every C3-fibre above12. The next
question is the action of Aut(H) on the whole lift set, retaining stabilizers
and move incidences rather than only its cardinality.
