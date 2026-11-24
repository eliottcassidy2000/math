# Exactly seven source/probe identifications are possible in the binary-Gram observer

Let A be the adjacency of a hypothetical Moore(57,2) graph on m=3250
vertices. Its standard binary feature matrix is X=A+I, with exact
Gram G=XXᵗ. Using one tournament vertex for each source and a
separate vertex for each probe, prescribe source→probe exactly where
the binary feature is 1. This note considers the same observer with
arbitrary binary features realizing the exact same Gram; sources
and probes must each be distinct within their own set, but a
source can also be a probe. Self-coordinate bits can be defined
separately (the construction uses zero there). The sharp number of
vertices used by the union is 2m−7, conditional on the Moore input.

## Exact construction

The reference graph has an induced path x_0,...,x_5. To extend
an induced path with current endpoint x_r for r≤4, exclude the
predecessor and any neighbor of x_r that is also a neighbor of
an earlier path vertex. There are at most 1+r excluded
choices, since distinct vertices have at most one common
neighbor. No other earlier path vertex can already be adjacent
to x_r, by the induction hypothesis. Degree 57 suffices.

The endpoints of the induced path are nonadjacent and have
one common neighbor, necessarily outside it. This closes a
seven-cycle, and any chord of this cycle would make a triangle
or four-cycle; hence it is induced. Label its vertices v_i
cyclically modulo seven. Identify the tournament source with
original label v_i with the probe that was to read column
v_(i+2) of X. Among the seven shared vertices, the bit required
in row i and probe j is then

    X_(v_i,v_(j+2)) = [j−i ∈ {−1,−2,−3} mod 7].

Thus the two requirements on a distinct unordered pair are
complements, exactly as required of tournament arrows, and
the required diagonal is zero. On any other unordered pair
of distinct vertices, at most one direction is already
prescribed by the source/probe read. Fill every unprescribed
arrow arbitrarily. This physically realizes the same full
binary matrix X, with seven identities between its source
and probe sets. In particular all Gram entries are unchanged.

## Optimality even with a different binary factorization

From the Moore identity A²+A=56I+J one has
G=(A+I)²=57I+A+J. Its diagonal is 58, its off-diagonal
entries are 2 on graph edges and 1 on graph nonedges.
It has rank m: A has eigenvalues 57,7,−8 and adding I
does not give an eigenvalue zero. Therefore any exact
binary feature factorization needs at least m probes.

Now suppose eight vertices were both sources and probes.
Their distinct-pair feature bits form a tournament on those
eight vertices. For each pair, its common-outneighbor count
inside this tournament is at most its full binary Gram entry.
It cannot exceed two; attaining two forces that source pair
to be adjacent in the reference graph. Accordingly the graph
of pairs attaining two must be triangle- and square-free,
since it is a subgraph of the corresponding Moore-induced
graph. Allowing nonzero separately defined self-coordinate
bits cannot weaken these inequalities.

The short accompanying checker rules out that eight-vertex
tournament exactly. Every tournament admits a Hamilton path,
so one may label its vertices in path order and exhaust the
2^21 remaining arcs. Of those anchored instances, 42,548
have all common-outneighbor counts at most two, and none of
their forced adjacency graphs is triangle-free (indeed
none is square-free either). No deduplication is used or
needed. Any larger overlap would give such an eight-subset
by restriction. Thus there are at most seven shared
vertices, and the probe rank bound gives ≥2m−7 in total.

This result concerns the exact binary Gram observer. It does
not rule out shorter encodings that only preserve a threshold,
use nonbinary features, or read adjacency with a different
decoding rule.

Replay:
```sh
c++ -O3 -std=c++17 moore-gram-overlap8.cpp -o /tmp/overlap8
/tmp/overlap8
```

## The probe lower bound survives exact distance-only observation

If the requirement is only to preserve the exact Hamming distances,
112 on edges and 114 on nonedges, the same probe lower bound m
still holds, even for arbitrary integer feature coordinates.
The relative source vectors f_i−f_m (i<m) have a Gram matrix
fixed by polarization from the distances alone. For the displayed
features its determinant is

    m · 64^1729 · 49^1520.

Indeed the difference coordinate basis e_i−e_m has ordinary
Gram determinant m and spans the zero-sum source subspace, on
which G has eigenvalues (7+1)^2=64 and (−8+1)^2=49 with the
stated multiplicities. Thus m−1 probes would force the matrix of
relative feature vectors to be square and integral, so its Gram
determinant would be a perfect square. The displayed determinant
is not a square, since m=3250 is not one and the other factors
are squares. Fewer probes already fail rank. The seven-vertex
overlap upper bound above concerns the exact Gram; arbitrary
distance-only factorizations need not have overlap inner
products at most two, so that particular argument does not
certify an optimal union size for them.
