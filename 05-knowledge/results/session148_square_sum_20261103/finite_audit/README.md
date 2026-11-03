# Independent square-sum growth and Hamilton audit

Status: FINITE-EXACT, with the general birth formula inherited from the
[Session 147 growth report](../../session147_square_sum_20261103/growth/REPORT.md).
Let G_n be the simple graph with vertices 1 through n, where distinct i,j
are adjacent exactly when i+j is a square. The birth of n is the set of new
edges {k^2-n,n} with n < k^2 < 2n. Thus its size is
floor(sqrt(2n-1))-floor(sqrt(n)). The upper bound is strict: a self-loop at
k^2=2n is not an edge. Counting by largest endpoint gives cumulative edges.

For n=1,...,15 the birth counts are 0,0,1,0,1,1,1,1,1,1,1,1,2,2,2;
the edge totals are 0,0,1,1,2,3,4,5,6,7,8,9,11,13,15. Independent
pair enumeration and square-sum enumeration agree with the birth formula
for all n=1,...,500. The checker also computes oriented Hamilton paths by
exact subset dynamic programming for n=1,...,25. Counts are ordered paths,
so reversal is counted separately except for the one-vertex graph.

At n=15 a witness is 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9. The graph
has 15 edges and this path uses 14; its one unused edge is {1,3}. Vertices
8 and 10 are 9-1 and 9+1, but their graph neighborhoods differ: N(8)={1}
and N(10)={6,15}. This observation alone asserts no symmetry.

Reproduce with `python3 -B verify.py` or `python3 -B -O verify.py`. The
verifier uses explicit failure gates, tests three independent edge counts,
checks each path edge and vertex, and rejects a repeated vertex, a broken
edge, a removed required edge and a false square-boundary birth formula.
`results.json` records the full finite output, including first 25 primes;
prime-index alignment is descriptive and not evidence of a prime theorem.
