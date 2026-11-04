# Full collisions require prime quotients at every highest changed node

## 1. Setting and dependencies

Let U be obtained from a named tournament T by edge-XOR of an arbitrary
family of modules of the original T. The family need not be laminar and
its original modules need not remain modules at intermediate steps.

We use two proved Session188 results:

1. ../module_decomposition/PROOF.md: the original strong-module tree remains
   a common recursive modular partition of T and U. At an originally prime
   quotient the output quotient is unchanged or conversed. At an originally
   transitive quotient it can change arbitrarily. Every physical induced
   subset's complete ordered-path-cover profile is coordinatewise
   nondecreasing. If H(T)=H(U), all these hereditary profiles are equal,
   and every originally transitive quotient remains transitive.

2. ../covered_laminar/ORDERED_BLOCKS.md: any nonidentity permutation of an
   ordered join of nonempty physical blocks is detected by a physical cut
   port, even when each block is replaced by an equal-q-profile block, in
   any fixed ambient quotient. Its strict lower bound uses positive
   substitution coefficients, not an assumed exact scalar multiplication.

Here a strong module means one that overlaps no other module; this is
distinct from a strongly connected tournament. We use the original
strong-module tree, not an asserted canonical strong tree of U.

## 2. Highest changed quotients

At each non-leaf original strong node M, its fixed original children form
a modular partition on both sides. Call M changed when the two quotients
on those same physical children differ. Call it highest changed when no
proper ancestor is changed. Let Hset be all highest changed nodes.

The nodes of Hset are pairwise disjoint: original strong nodes are nested
or disjoint, and an ancestor-descendant changed pair cannot both be highest.
There can be more than one highest changed node. If T!=U, Hset is nonempty,
because every old edge belongs to the quotient of the least original
strong node containing both its endpoints.

Every changed edge lies inside one of Hset's nodes. Indeed, its least
common ancestor is a changed quotient node and has a highest changed
ancestor. Thus outside the interiors of these disjoint nodes all old arcs
agree, including incidences between any such node and its exterior.

## 3. The reduction theorem

**Theorem.** Suppose H(T)=H(U). If any highest changed node M has an
originally transitive child quotient, then F_T != F_U. Consequently, if
T!=U but the full named two-port responses coincide, every highest changed
quotient is prime and is conversed.

**Proof.** Write M's original children in quotient order as
A_1 -> ... -> A_k. Scalar equality and dependency 1 imply that the target
child quotient is another transitive order. Because M is changed, this is
a nonidentity permutation. Crucially, for every child A_i,

q_h(T[A_i]) = q_h(U[A_i]) for every h.

This is a hereditary consequence of scalar equality, not an assumption
that a child's internal operation is a single reversal. In fact equality
holds for the complete q-profile on every fixed subset of A_i.

Temporarily form V by replacing only T[M] by U[M] and leaving all old names
outside M as in T. The exterior incidences are unchanged, so M is a module
with the same rooted ambient quotient R. Apply dependency 2. It supplies
the same physical ports on T and V, with both ports having M's exterior
incidences, for which

F_V(a,c)-F_T(a,c) >= 2 H(R) product_i H(T[A_i]) > 0.        (1)

For every other highest changed node N, its original and target complete
q-profiles are equal by hereditary scalar equality. Moreover N remains a
module in the actual port extension. It is an old module on both sides;
the old arcs between disjoint modules M and N are uniform, and both new
ports have M's incidences to N. Thus the ports also see N uniformly.

The substitution profile identity says that replacing a module by a
same-named equal-q-profile tournament preserves the actual Hamilton count
exactly. Replace each of these other highest nodes in V by its U interior,
one at a time. The nodes are disjoint and their boundaries do not change,
so every replacement remains valid. The resulting tournament is U, and
the same physical query has F_U(a,c)=F_V(a,c). Hence (1) separates T and U.
There is no cancellation by simultaneous changes at other highest nodes.

If a highest changed node is not linear, dependency 1 says its quotient is
prime. A changed prime quotient can only be its converse. This proves the
necessary condition for full collision. Finally full response equality
itself implies H(T)=H(U), since the physical query a=c=0 has x,y preceding
all old vertices and satisfies F_T(0,0)=H(T). QED.

## 4. Exact scope

Prime quotients of order three are included: the directed triangle is
prime, so the first unresolved quotient size is three, not five. There
is no omitted small-prime exception.

For example, take C3[TT2,C3,TT1] and reverse its first two children while
preserving the parent quotient. Those two children are original strong
modules and give two highest changed nodes, one linear and one prime. The
linear node's cut separates the full change even though the prime node
also changes. The verifier retains this explicit multiple-node control.

By contrast, in C3[TT2,TT2,TT2] with the cyclic quotient conversed, the
root is a highest changed prime node of size three. If the whole tournament
is conversed, its pair children also have changed linear quotients, but
they are descendants of the changed prime root and are not highest.
This reduction alone does not isolate those descendant changes. Some
such examples are already separated by other Session188 theorems; being
in the residual is not evidence of a collision.

The result is a necessary reduction, not a full-response collision or
injectivity theorem. It applies to simultaneous XOR of original modules.
It does not assert the same persistent original-tree action for an arbitrary
sequence in which later modules need not have been original modules.

## 5. Reproducibility

The exact checker enumerates every labeled tournament through five old
vertices and every original-module XOR family, retains the scalar-equal
cases, constructs the original strong tree directly from all modules,
and identifies every highest changed quotient. It checks child-profile
equality, the prime-converse restriction, and the literal cut query for
every highest changed linear node, including all other simultaneous changes.
It also retains six-vertex explicit multiple-node and prime-three controls.

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
