# An infinite family with growing anti-involution alignment memory

**Status:** the family formulas and nonflatness for every k≥2 are PROVED below.
`check_wreath_family.py` gives FINITE-EXACT controls through81 vertices.
The wreath-product symmetry mechanism is standard; no priority is claimed.
The useful consequence here is its exact anti-involution/centralizer accounting
and its explicit transport counterexamples at every scale.

## Family and intrinsic modules

Let W_0 be the singleton and W_k=C3[W_(k-1),W_(k-1),W_(k-1)]. It is the
lexicographic tournament on strings of k ternary digits: at the first differing
digit, orient according to0→1→2→0. Put N=3^k. It is regular of outdegree(N−1)/2.
Let t_k negate every digit modulo3. It is an anti-involution with unique fixed
vertex0...0.

Write A_k=Aut(W_k). For k≥1 the three top-level blocks have size s=3^(k−1).
They can be recovered directly from W_k: two vertices are in the same block
if and only if their out-neighborhood bit rows have symmetric difference
of size less than s. For two vertices in the same block this difference is
internal, even, and at most s−1. For vertices in different blocks it has size
exactly2s: the third block contributes s, and the two remaining blocks
contribute(s−1)/2 and(s+1)/2. Hence every automorphism preserves the three-block
partition, and its induced permutation is a rotation of C3. Conversely any
three child automorphisms and any macro rotation extend. Therefore

    A_k = A_(k-1)^3 semidirect C3,
    |A_k| = 3^a_k,       a_k=(3^k−1)/2.

This proof establishes the full group, not merely an exhibited subgroup.

## Exact centralizer and mark counts

Let C_k=C_(A_k)(t_k). A centralizing automorphism cannot have a nonzero macro
rotation, since t_k conjugates that rotation to its inverse. In its three
child actions, the action on block0 must centralize t_(k-1), the action on
block1 is arbitrary, and the action on block2 is forced by conjugation.
Thus

    |C_k| = |C_(k-1)| |A_(k-1)| = 3^c_k,
    c_k = (3^k−2k−1)/4.

The expression is an integer for every k≥0. All anti-involutions are conjugate,
so their number is

    i_k = |A_k|/|C_k| = 3^((3^k+2k−1)/4).

The useful factorizations are

    |A_k| = N |C_k|²,
    i_k = N |C_k|.

A_k is transitive on vertices. Each anti-involution has one fixed vertex, and
conjugation transports that vertex. Consequently **exactly |C_k|
anti-involutions fix any given named vertex**.

Examples:

| k | N | Aut order | centralizer order | anti-involutions | marks fixing each vertex |
|---:|---:|---:|---:|---:|---:|
|0|1|1|1|1|1|
|1|3|3|1|3|1|
|2|9|81|3|27|3|
|3|27|1594323|243|6561|243|

These are formulas for this specific iterated family, not universal formulas
for regular tournaments or arbitrary tournament modules.

## A procedural hostile example at every k≥2

Take the three anti-involutions s,t,u from the nine-vertex control in
`TRANSPORT_AND_HOLONOMY.md`. In W_k choose the nine-vertex submodule whose first
k−2 digits are all zero. It is fixed setwise by t_k. Replace the restriction
of t_k to that submodule by s, t, or u respectively, retaining t_k outside.
Each new map is an anti-involution: the cross-module arcs only depend on the
unchanged leading digits. The product of two such marks is identity outside
the module. Its odd square root is consequently identity outside too.
Therefore the triangular holonomy is exactly the order9 holonomy on that
module, extended by identity. It remains nontrivial for every k≥2.

This is a constructive infinite family, not extrapolation from the cases
checked through81 vertices.

## Information consequence, with the observation model explicit

Suppose a representation must specify the actual anti-involution, while its
fixed vertex is already known. In this family there are |C_k| possibilities,
so a fixed-length binary representation that distinguishes them all requires
at least

    ceiling(log_2 |C_k|)
      = ceiling(((N−2log_3 N−1)/4) log_2 3)

additional bits. This is a counting bound for that exact target. It does not
apply to a decision that ignores the mark, or prevent a short formula for the
family as a whole. Individual choices and a description of the family are
different information tasks.

Thus “designate a node for its orientation” is insufficient at increasing
scales: at order27 a fixed named node still permits243 anti-involution frames.
The three-vertex base happens to hide this issue because its corresponding
number is1. A global base mark plus exact relative maps resolves the transport
problem; simply erasing every base mark does not.

The centralizer is a second, separate issue: even after the anti-involution
has been fixed, |C_k| ordinary automorphisms preserve it. Named-port operations
may need their actual action, rather than just the centralizer order. This is
the same *type of information distinction* as prior clone-merger and LCA
accounts, but the concrete carrier here is the explicitly proved tournament
wreath action, not an asserted isomorphism with every other problem.

## Reproduction

    python3 check_wreath_family.py
    python3 -O check_wreath_family.py

Both produce identical reports with12336 active checks. They verify tournament
regularity, the intrinsic top-block distance criterion, the anti-involution,
the order-exponent recurrences, and the embedded nonflat witness for
k=0,...,4 (orders1,3,9,27,81). They do not enumerate the full groups at27 or81;
those group orders follow from the proved block-recovery argument.
