# Exact validation and scope

Run from the workspace root:

    python3 session36-runner/check.py > session36-runner/check-normal.out
    python3 -O session36-runner/check.py > session36-runner/check-optimized.out
    cmp session36-runner/check-normal.out session36-runner/check-optimized.out

Both completed successfully and outputs byte-match. The checker uses explicit exceptions, not Python assertions. It executes 92,927 always-active checks. There are no floating-point truth tests.

The exact controls comprise:

1. Ten convolution matrices, including a 196-cell base7/depth2 matrix. Rank is certified by a matching modular lower bound and explicit linearly independent rational null vectors. The predicted rank is M-M/q+1.
2. Six exhaustive binary-carrier universes totaling 152,064 carriers. Every nontrivial collision has a zero-mean q-periodic difference; both carriers have at least q components. Thresholds 2,3,4,7,4,8 are attained.
3. Sharp disjoint periodic hostiles at bases 2,3,4,7, including 343 components. Every retained overlap equals 1/(4b), and the next probe gives 1/(2b) versus 0. All phase breakpoints lie on the tested grid, so piecewise linearity upgrades the exact knot calculations to full phase responses.
4. An explicit rational finite-phase bank that cannot distinguish two disjoint one-component intervals of common length 3/686.
5. Exact wall partitions for tight AP rows {1,...,m}, m=1,...,20. They have only the predicted reduced-grid isolated points, with no safe positive-length cells. The m=13 carrier consists of six points.
6. Integer cyclotomic controls of mean plus first 1,2,3 complex moments on all 4,096 twelve-cell carriers. Every fiber containing a carrier with at most C components is a singleton.
7. 9,600 exact divisibility checks for prime7 valuation-bounded probe banks.
8. 65,024 binary-carrier pairs for the rational uniform-noise estimate in NOISE.md. The exact supremum of each response error is attained at tested phase knots. Empty/full carriers are checked with their separate mean bound.

These finite checks corroborate the symbolic proofs. They are synthetic observer models and tight-AP controls, not a census of the physical 165-row LRC frontier. No claim of a new native owner/word packet, universal component bound, physical arrival, or LRC closure is made.

Frozen checker SHA-256:

    c98c2c956ce06adf52fc255ef23e7f3fd2009a040bd9b49985addfff9d9a072a

Frozen receipt SHA-256:

    e7dfcb8cca9f7bfef66a41988d61c19a461d764d1cd6a9cc0bb9910da74000a2

`SOURCES.json` records exact input paths, byte counts, and SHA-256 hashes for thirteen read-only canon/frontier sources.
