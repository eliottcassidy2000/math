# Functional graphs, graph linearity, and an inverse-limit coloring obstruction

2026-10-20. Elementary deductions connecting the two source blueprints to the session's LFSR and tournament results. These are statements about specified graphs and topologies, not a new resolution of a major coloring conjecture.

## 1. Two meanings of linearity meet, but do not coincide

For any function T:X→X, form a simple undirected graph with edges {x,T(x)} whenever x≠T(x). Collapse duplicate edges. Each finite subgraph is a pseudoforest: orient each retained edge using one of its generating arrows. The orientation has outdegree at most one, so every connected finite subgraph has at most as many edges as vertices and hence at most one cycle.

For any finite edge set E, its minimum number of monochromatic edges in a two-coloring is exactly its number of odd cycles. Remove one edge from each odd cycle and bipartition the remaining components; conversely each odd cycle requires a monochromatic edge. The cycles are vertex-disjoint and have length at least three, giving

    edgeCost(E) <= |V(E)|/3.

Thus every functional graph is 1/3-linear in the attached graph blueprint's sense, even if the function is nonlinear. A directed triangle makes the constant sharp. This is a simple special case with a direct 3-coloring, not a new application of the blueprint's difficult countable-coloring theorem. For a permutation whose finite orbit lengths are powers of two, its simple graph has no odd cycles and is bipartite. Fixed-point loops must be deleted: a graph with a retained loop cannot be properly colored.

## 2. A concrete Borel coloring of the tournament clock

Use the session's completed formal tournament space

    W = product_(n>=0) F_2[T_n],   U = I+S,

where S sends each tournament class A to the class s(A) obtained by adjoining a source. Source deletion proves s is injective. Give W the product topology, each finite degree space discrete. Fix an ordering of classes in each degree.

For w≠0, let r be its least occupied degree and A its first occupied class in that degree. Define c(w) to be the coefficient of s(A) in w; set c(0)=0. The least occupied degree and its coefficients do not change under U, because S raises degree. The coefficient at s(A) changes by the coefficient at A, namely one. Injectivity of source adjunction excludes other contributions. Therefore

    c(Uw) = 1-c(w),  w≠0.

This gives an explicit Borel proper two-coloring of the simple orbit graph. It is locally constant on W\{0}: the finite prefix through degree r+1 determines both A and c(w). Its Borel property at the one remaining point follows by adjoining {0} to the appropriate Borel color class. It uses actual coordinates, not a nonconstructive choice of an origin in every odometer orbit.

## 3. Continuity fails precisely at the accumulating fixed point

No proper coloring of this graph into any discrete palette can be continuous at zero. If it were, some basic zero-cylinder, fixing degrees through N to zero, would have constant color. Choose w to be a single tournament basis vector in degree N+1. Both w and Uw=w+Sw lie in that cylinder, and Sw≠0, so they are adjacent distinct vertices with the same color, a contradiction.

The zero vertex itself is isolated in the simple orbit graph because U0=0. The obstruction is topological accumulation of edges near zero, not an uncolorable self-loop. Every finite truncation is discretely two-colorable, and the inverse limit has the explicit Borel coloring above, yet continuous proper coloring of the full compact space fails. This distinguishes combinatorial colorability, Borel regularity, and continuous regularity without changing the underlying update rule.

## 4. A useful proof target suggested by the distinction

When importing a finite-state argument into an arithmetic or geometric inverse limit, specify which information must survive: finite defects, measurable choices, continuous coordinates, or algebraic module structure. The three losses documented in this session are different:

- cycle graphs can forget LFSR primary multiplicities;
- pair correlations can forget four-root extension costs;
- finite proper colorings need not assemble into a continuous limiting coloring.

A candidate cross-field reduction should identify an invariant that survives its chosen maps before using finite evidence to infer a global theorem. The examples here provide exact stress tests for such a reduction.

## Checks and review

The retained checker independently canonicalizes tournaments through degree four, builds source-adjunction columns, and checks all 512 formal states. It verifies the proposed coloring on all 496 directed nonfixed edges; ordinary and optimized Python outputs agree. The global topological conclusion is proved above, not extrapolated from this finite check. Independent mathematical review confirmed the arguments and emphasized that the tournament orbit graph itself is 0-linear, because it is bipartite. “Orbit graph” here always means immediate-update edges {w,Uw}, not a clique on each full orbit.
