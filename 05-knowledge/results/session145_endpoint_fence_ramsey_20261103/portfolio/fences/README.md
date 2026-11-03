# Fence geometry and cycle-sharing products

Read [PROOF.md](PROOF.md). This lane proves two scoped reductions:

1. A supported unit-fence configuration whose **split graph is a pseudoforest** has total field area at most n/4, with equality only for disjoint unit squares.
2. In the endpoint-only **product** problem, every cactus graph has an exact maximum given by balanced regular-polygon cycle partitions. At twelve edges that is 6+3 sqrt(3). The existing Friedman theta(5,5,2) example has product (89+26 sqrt(13))/16, strictly larger. Its exact coordinates, all unit lengths, and legal intersections are certified here. This recovers and explains a known construction, not a new record or a global optimum.

T-junctions prevent exporting the product cactus reduction to the fence problem. An exact five-fence cactus fixture has seven split edges and two non-unit triangular bases. [CORRECTIONS.md](CORRECTIONS.md) records the discarded exploratory extension. The no-through cactus corollary is safe.

## Verification

From this directory, Python standard library only:

    PYTHONDONTWRITEBYTECODE=1 python3 verify.py verification.json
    PYTHONDONTWRITEBYTECODE=1 python3 -O verify.py verification.optimized.json
    cmp verification.json verification.optimized.json

Both runs check 4,014 conditions and produce byte-identical receipts. Universes: all 1,099 labelled simple graphs of orders 1 through 5; all 1,885 integer partitions with parts at least 3 and total 3 through 30; all balanced candidates through n=50; and four exact planar geometry fixtures. Geometry uses rational arithmetic in Q(sqrt(13)), with a positive vertical scale factored out. Products are enclosed by rational intervals from Machin's identity and alternating sine/cosine series, then rounded outward. No floating point or optimization-removable asserts enter a certificate.

The geometric fixture engine covers connected or non-nested disconnected arrangements; it is not advertised as a general hole-aware face engine. The theorem itself includes nesting by a separate proof. Invalid crossing/overlap checks and the forbidden transfer of a T-junction fence to product mode are hostile controls.

Primary-source and repository pins are in [SOURCE_PINS.json](SOURCE_PINS.json). The source checkout was read-only and no publication was performed by this lane. Polygonal isoperimetry and graph facts are classical. The specific scoped synthesis and independently reconstructed coordinates are new to this recovered lane, with no historical novelty claim.
