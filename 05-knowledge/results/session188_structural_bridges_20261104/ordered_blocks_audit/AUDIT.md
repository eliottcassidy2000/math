# Independent audit: physical cut ports and covered laminar families

Status: PASS. Independent normal and optimized runs completed successfully
with byte-identical receipts.

Target: ../covered_laminar/ORDERED_BLOCKS.md.

## Findings

1. The substitution profile lemma is valid for every number of ordered
   path slots. Compressing maximal runs retains both slot order and exterior
   names, and reading slots in order canonically numbers all module runs.
   The one-run Hamilton coefficient is exactly H of the contracted quotient,
   which is positive. This justifies the ambient lower bound; no invalid
   exact scalar multiplication of a general contrast is used.

2. Reversing every maximal module run in every ordered path slot is a
   bijection for a module converse. Its restriction to every fixed physical
   subset remains valid because the intersection of a module with that
   subset is a module in the induced tournament. A legal sequence therefore
   preserves all ordered-cover profiles hereditarily. This does not mean
   arbitrary simultaneous overlapping original modules preserve them.

3. The ordered-quotient injection is genuinely injective: each source path
   slot has at most one run of any block. Reordering these runs along a fixed
   Hamilton order of the occupied target quotient retains their internal
   directed orders; the physical block names recover the original slot.
   Thus it proves every q_h inequality, including zero and empty-slot
   boundary cases. No independence of endpoints is assumed.

4. A nonidentity permutation has a cut with a reversed cross-cut pair. The
   same physical y creates a directed triangle in the target extended
   quotient, whereas the source extended quotient is transitive. The
   x->y convention and the actual old-name incidences have the correct
   signs. The strict bound 2 times the product of block Hamilton counts
   follows from at least three target quotient Hamilton paths versus one
   source quotient path.

5. In a covered laminar family, maximal proper selected children of B are
   disjoint and cover B. Quotient SCCs are modules and group those children
   into an ordered join. Every selected descendant is inside one such group;
   B alone reverses their order. Within each group, all changes are a legal
   module-reversal sequence, hence preserve its complete q-profile. This
   supplies precisely the hypotheses of the physical cut theorem.

6. Maximality of B prevents an enclosing selected operation from changing
   its exterior incidences. Every other maximal selected set is disjoint.
   Lifted ports see every outside selected module uniformly, so those
   outside reversals preserve the literal query exactly, rather than merely
   modulo four. They cannot cancel the strict contrast.

## Terminology clarification

The hypothesis that B is "nonstrong" means that T[B] is **not strongly
connected**. It does not mean that B is not a strong module in the Gallai
modular-decomposition sense. These two notions must not be interchanged.
The parent was notified of this wording clarification. The mathematical
proof explicitly uses strong connectivity and is sound with that meaning.

## Independent controls

verify.py imports no producer code. It counts ordered path covers by
literal permutations and choices of cuts: every backwards adjacency forces
one cut, and each forward adjacency may or may not be cut. This is
independent of the producer's subset-cover convolution. Literal two-port
Hamilton counts use memoized endpoint recursion.

The checker tests all five nonidentity permutations of each three-block
choice from TT1, TT2, TT3, C3 with total old order at most six. It checks the
physical cut, quantitative strict bound, complete q domination through
seven extended vertices, and an actual strong C3 ambient lift in every
case. Each target block is also replaced by its converse, exercising the
profile-preserving replacement clause.

It separately enumerates every labeled tournament through five vertices,
every nonempty laminar family of its nontrivial original modules, and every
physical induced subset. Each intermediate operation is checked to be a
current module; complete hereditary q equality is checked by literal
word/cut enumeration. For every fully covered maximal selected module with
non-strongly-connected child quotient, it constructs the SCC block cut,
checks their target reversal and profile preservation, and evaluates the
actual physical query including all other selected reversals.

The completed independent controls comprise 130 ordered-block permutations
and 130 strong ambient lifts, 16,882 nonempty laminar families through five
names, 38,454 verified current-module reversals, 512,534 physical induced
subset cover comparisons, and 1,064 fully covered maximal-module cut tests.

The all-order theorem rests on the audited proof, not these finite checks.
The audit does not extend the conclusion to strongly connected maximal
modules or arbitrary legal inversion sequences.

Reproduction:

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
