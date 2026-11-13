# Recover a complement-line spine from two unsigned wedge codes

This is an inverse form of the inherited
[path-complement trijection](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md),
using the [two-slot wedge observer](../session320_second_moments_and_native_contexts_20261113/WEDGE_RECONSTRUCTION.md).
It concerns literal lines, not arbitrary single-arc flips, and does not
replace the remaining global quotient-enumeration problem.

Give the vertices labels and consider the line graph L(K_n), whose
vertices are unordered pairs of labels. For adjacent pair-vertices e,f,
the wedge bit of a tournament T is 1 if its two edges point in opposite
directions at their shared vertex (one in, one out), and 0 otherwise.
This is exactly the positivity bit for simultaneous use in two optional
slots. For two aligned tournaments T,U, let D_e be 1 when their
directions disagree on e. Then
```
w_T(e,f) XOR w_U(e,f) = D_e XOR D_f.
```
The reference directions cancel. Integrate this difference along any
spanning tree of L(K_n); D is recovered up to replacing every bit
by its complement. For valid input this agrees on all other line-graph
edges. For abstract input the cycle sums give a necessary and sufficient
consistency test for the difference word.

Now try the two complementary supports for **agreement**. To be a
literal path-complement line one must have exactly n−1 edges, form a
spanning undirected path, and have w_T=1 at every adjacent wedge of
that path. The last test forces T to orient the whole path coherently.
Choose the polarity of U (U itself or its converse) that agrees on
that support. Then the two tournaments agree exactly on this directed
Hamilton path and disagree on every other edge, which is the native
complement operation. Conversely any actual aligned line passes these
tests. This recovers the unique marked spine from aligned unsigned data
for n>=5, since binom(n,2)>2(n−1). At order four the two path supports
can both work: take directed paths 0→1→2→3 and 2→0→3→1, which
partition all six pairs, and complement outside either path.

A useful cheap filter precedes the path connectivity test. If the
support is a line with start s and end t, then under the successful
alignment and polarity,
```
out_T(v)+out_U(v) = n−1 + [v=s] − [v=t].
```
Indeed off-spine arcs contribute 1 per incident pair, while retained
spine arcs contribute either 2 or 0. This is necessary, not sufficient;
the wedge-difference test gives the exact missing constraint.

To enumerate ordinary weighted darts between two isomorphism classes,
one can sum aligned incidence tests over all bijections between fixed
representatives, using the actual orientations rather than accepting
either polarity. The resulting count divided by
|Aut(T)| |Aut(U)| is the directed dart multiplicity, since every
Hamilton word giving the target yields |Aut(U)| alignments, and source
path markings form free Aut(T)-orbits. For converse-merged counts one
must aggregate the distinct ordinary classes in the corresponding
fibres, not count both copies of a self-converse class. This is an exact
algebraic constraint, not an efficient procedure for testing all
bijections at high order. The canonical Fano families provide a
different constraint by giving explicit regions with no internal lines.
