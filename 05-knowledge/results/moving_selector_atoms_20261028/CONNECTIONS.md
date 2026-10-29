# The common mechanism: retain the action needed by the next operation

2026-10-28. This account distinguishes inherited mechanisms, the present
specializations, and OPEN extensions. It is not a claimed new general theory
of group actions, automata, or Bruck loops.

## Source recovery and what changed

The [previous singleton minimum](../minimal_blue_selectors_20261028/THEORY.md)
already separates minimal behavior from literal reconstruction. Its moving
selector question is answered here for the9-pair, while all-size faithfulness
of the enlarged Hamilton-presentation carrier remains open. C2's obstruction
concerns the more general ordered-run interface and cannot be substituted
for that unanswered question.

The [forest-complement packet](../strong_atom_continuation_vitali_20261026/README.md)
already proves agreement reconstruction and composition along quotient walks,
including the six-vertex15/23 hostile. C3 does not rediscover that composition
law: it extracts the induced component action, its exact order fiber and
the base criterion for replacing joint order by marginal probes.

### Concurrent result: a selector that writes C3 tags forces reconstruction

The incoming [joint-port packet](../joint_port_germs_20261028/README.md)
distinguishes a read-only cursor from arbitrary-position C3 substitution.
Its latter menu already forces a discrete behavioral quotient on Hamilton
presentations. Here is an all-size audit of the tag mechanism, credited to
that packet rather than presented as a separate discovery.

For an n-vertex presentation, choose distinct old positions i,j. Repeated
C3 substitution builds a strongly connected tournament of any odd order
at a selected vertex. Replace i by such a module P of order
p=n+1+(n mod2), and j by a module Q of order q=2n+1. The final order is
4n+(n mod2); every substitution preserves a Hamilton presentation by using
the internal identity Hamilton path. Only the allowed C3 atom is needed.

A strongly connected module cannot overlap another module without one
containing the other. Indeed, if M overlaps Q, a vertex of M outside Q
sees all of Q uniformly. The module condition on M then forces every
vertex of Q outside M to see M intersect Q in one common direction. This
is a directed cut in Q, contradicting strong connectivity.

There are only n-2+p < q vertices outside Q. Thus any other module of size
q cannot be disjoint from Q, cannot overlap it, and cannot strictly contain
or be contained in it. Q is the unique module of that size. Contract it;
there are then n-1 < p vertices outside P, so P is likewise the unique
module of size p. Both modules and their order are intrinsically decoded
from the unlabelled final tournament. Their contraction recovers
(T;i,j), hence the direction of the chosen old arc.

Any unequal identity-Hamilton presentations disagree on some old arc.
Applying the same two-port construction exposes that disagreement, so the
behavioral quotient for this richer language has exactly
2^((n-1)(n-2)/2) states. Sufficiency is literal presentation storage;
necessity is the distinguishing construction. This counts presentations,
not A000568 isomorphism classes, and gives no new enumeration complexity
bound. The context order is a sufficient bound, not an optimal one.

We replayed its checker in ordinary and optimized Python:118 tag contexts,
all ordered pairs and Hamilton presentations at orders2..4, plus four
selected presentations with two pairs each at order5;944 active checks.
This finite check supports the implementation; the preceding argument
supplies the unbounded quantifiers. The module argument also explains why
calling this operation a mere moving observation would lose information:
the tournament itself has become a persistent register.

The [frame theorem N1](../strong_attachment_tail_compiler_20261028/THEORY.md)
already proves LCA transport on free frames. The new requirement is operational:
the chosen anchors must kill the component stabilizer before individual
responses recover every ordering. A group base and a stabilizer chain are
standard notions; the [GAP reference](https://gap-system.github.io/gap/doc/ref/chap43_mj.html)
documents the underlying permutation-group machinery. We prove the scoped
if-and-only-if condition directly, without importing a group classification.

The [centralizer folded DP](../session90_even_graphs_odd_descent_20261029/oddness/CENTRALIZER_ORBIT_PATH_DP.md)
already reduces64 to24. C5 computes the coarsest equitable partition for its
terminal counting task, reducing further to16 and measuring precisely what
new named tests cost. This is an application of standard refinement, with
the whole finite quotient exported, not a runtime claim or an A000568 formula.

## Why tree compression and the sparse even graph are compatible but different

THM-4835
[`euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md`](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md)
already contains the tree end-edge map. Session90's
[tree codec](../session90_even_graphs_odd_descent_20261029/trijection/TREE_COMPRESSION.md)
makes the inverse and transport explicit. Fixing an oriented tree leaves
C(n-1,2) free coordinates, mapped bijectively to a tournament on its n-1
edges. The tree and its marked symmetry action remain part of the input.

The [enlarged sparse Royle-even codec](../session90_even_graphs_odd_descent_20261029/trijection/SPARSE_EVEN_CARRIERS.md)
instead stores the complete tournament in a larger uncoloured graph. Its
tail/head gadgets recover every arc and exclude extra automorphisms; distinct
rigid tree components adjust the graph's Euler characteristic to1 without
introducing component permutations. The resulting graph has e=v-1 and the
same full automorphism group as the tournament. Its group is odd, so its
Royle edge-reversal character is trivial.

The common preservation contract is exact marked isomorphism and actual
port transport. The operation correspondence must still be supplied:

* For the tree chart, C6 proves that complement off the frozen tree becomes
  converse on the edge tournament. Re-marking changes the chart.
* For the enlarged graph codec C with decoder D, a tournament operation F
  can be transported as C F D on the codec's image. This tautological
  conjugation preserves behavior but does not make F an ordinary graph
  operation or improve its computational cost.

Ordinary graph complement is not the transported native tournament complement.
A graph with N vertices and N-1 edges has complement with(N-1)(N-2)/2 edges,
which is not N-1 for N>=5. The enlarged codec therefore leaves its e=v-1
family under ordinary graph complement. Similarly, contracting all tail/head
gadgets forgets exactly the port choices and gives K_n for every source T.
These are concrete reasons to retain decoding and realization data, rather
than infer composition from equinumerosity or a shared cycle rank.

The tree action calculation separates another pair of notions. Its parity
sheet is defined relative to vertex labels and is exchanged by an odd
relabeling. The order-parity sheets in C3 are invariant under the particular
odd joint automorphism group. Neither is ordinary Eulerian degree parity,
and neither is the Royle character of an arbitrary graph.

Incoming Session91's [guarded-flip analysis](../session91_packing_golden_carriers_20261029/operators/GUARDED_FLIP_LEAKAGE.md)
provides a useful operator counterpart. For its free flip operator B and
feasibility projector P, PB^2P-(PBP)^2=PB(1-P)BP records excursions through
forbidden states. Its exact rank formula is proved for independent-set
constraints. We do not transfer that rank formula to Hamilton legality:
the constraint spaces differ. The common lesson, now visible in C6, is
that a free translation calculation does not determine guarded evolution.

## Nonassociativity records a projected phase, not a faulty composition law

The incoming [transport and Bruck-loop proof](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md)
uses the classical construction x circle y=sqrt(x y^2 x) on a uniquely
2-divisible twisted subgroup. The primary source is
[Vojtechovsky, Three lectures on automorphic loops, Proposition2.3](https://www.cs.du.edu/~petr/data/papers/three_lectures_on_automorphic_loops.pdf).
Its exact order27 occurrence on wreath9 is inherited, not a new classification
of loops or of all nine-vertex tournaments in this packet.

Actual permutation maps compose associatively. Pair-local square-root choices
can produce a nontrivial centralizer correction around a triangle. After
discarding that correction, the projected frame coordinates have a
nonassociative operation. A chosen common base frame instead gives flat
relative maps by g_i^(-1)g_j. This explains why both the Bruck description
and the LCA description are useful: they retain different parts of the same
specified transport problem.

The new even-parity holonomy check is a cheap hostile. Its two3-cycles have
no parity signal, although their action on named ports is nonidentity.
Thus retaining the single order bit from C2 cannot replace the actual
centralizer phase. Conversely retaining every phase is unnecessary for the
16-state folded observer, which deliberately ignores those port identities.

The literal Paley27 construction provides a different finite-field witness.
Frobenius fixes the prime-field anchors but rotates three non-prime-field
ports. Individually selected phases can match a transposed order one port
at a time; no common phase or full tournament automorphism matches the entire
ordering. Once a non-prime-field port is retained, the displayed C3 phase is
fixed. This is an exact map and a lost joint coordinate, not an identification
of Paley27 with the wreath tournament or its Bruck loop.

## The next tests

1. **OPEN: faithfulness of K_scan on Hamilton presentations.** The7- and
   padded9-pairs separate, and wreath9's39 presentations separate. Search
   sources for which both endpoint stabilizers are nontrivial; triviality
   on either complement sheet already implies reconstruction. A general
   ordered-run counterexample does not answer this narrower question.
2. **OPEN: minimum behavioral order information under supplied-run
   quotient-walk contexts.** The component-order orbit is faithfully
   sufficient, and C2 shows its two classes are observable in one fiber.
   Determine whether every distinct fiber class admits such an uncoloured
   assembly/complement context that distinguishes it, without freely
   inserting persistent selected tags, or find a genuine behavioral merger.
   For the richer arbitrary-position C3 language, the incoming reconstruction
   theorem already settles the presentation minimum.
3. **OPEN: selective frame refinement.** When an attachment requests a
   port, refine only by its stabilizer; retain common witnesses at each
   LCA. Compare this with full-base storage on the13-vertex hostile and
   the81/27 pair. The proof obligation is compatibility of joint witnesses,
   not just equality of individual orbit counts.
4. **OPEN: guarded tree action.** C6 completely describes the unguarded
   tree-complement translations. Determine which of those translations
   admit legal Hamilton realizations along a specified path of states.
   This returns the exact realization guard omitted by the easy linear span.

The reusable move is to minimize the observer after the operation language
is fixed, and separately identify the faithful interface needed to change
that language. The concrete evidence now includes the sixteen-state count,
the one-bit run-order obstruction, and the free-base criterion. This does
not transfer an unproved convergence, coloring, or Jacobian conclusion to
another problem merely by analogy. LRC(14) remains open.
