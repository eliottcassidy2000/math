# Audit notes

The proof is uniform for undirected guards or reversal-closed directed
path banks at `n>=2`. Nonempty-bank scope is required for `Umax`; empty
banks have the separate constant polynomial. The empty subset is NOT
paired by reversal, but its contribution `2^N` meets the claimed bound
because `Umax>=1`. Leaf-half profile doubling applies only for `j>=1`.

Exact content uses integer translation in both directions; consequently
the gcd is unchanged by the binomial-moment basis shift. The presence of
`M0=2^N` rules out an odd factor in the content. Neither bound nor exact
content is inferred solely from maximum response.

The checker uses explicit exceptions, not assertions. It verifies all
74 small guards, 902 nonempty compatible subsets, literal tournament
histograms, profile-evenness, 112 leaf-half cases and both rejected scope
controls. Frozen n25/n26 source profiles give exact moment valuations,
primitive coefficients and rational probabilities. Normal and optimized
runs match. Isolated-copy portability includes only this packet and its
two authenticated source snapshots.

Independent reviewer accepted the mathematical proof and independently
confirmed n25/n26 contents. Its requested nonempty-subset quantifier in
the leaf paragraph has been applied. The wider independent audit remains
in the parent session's separate review packet.
