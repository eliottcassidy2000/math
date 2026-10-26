# Decorated tournament composition: a complete ternary realization stratum

[CLASSIFICATION.md](CLASSIFICATION.md) gives the proofs and exact scope.

For n=4·3^d, one fixed labelled cyclic-triple hypergraph has ordinary tournament realizations with every automorphism order

    3^(3^d), 3^(3^d+1), ..., 3^((3^(d+1)-1)/2).

The exact number of types at each order is given by P_0(z)=2 and

    P_(d+1)(z) = (P_d(z)^3 + (3z-1) P_d(z^3))/3.

This is a full ordinary-isomorphism classification of the stated C3 realization strata, not merely two witnesses. The binary palette meets every ordinary orbit but does not contain every labelled realization. All path-cover counts agree; ordinary marked-Hamilton-path orbit counts vary by the inverse automorphism factor. The full metagraph still needs actual path-port data.

Canonical modular classification and the 12-vertex/prime-valuation results are credited to inherited and concurrent work in [PROVENANCE.json](PROVENANCE.json), with frozen source snapshots. The all-depth polynomial spectrum is the present continuation, without a global literature novelty claim.

## Replay

No external dependency or repository checkout is required:

```bash
python3 -B check_classification.py > /tmp/session68-classification.normal.json
python3 -B -O check_classification.py > /tmp/session68-classification.optimized.json
cmp /tmp/session68-classification.normal.json /tmp/session68-classification.optimized.json
cmp /tmp/session68-classification.normal.json receipt.normal.json
```

Both modes pass **29,077 exact checks**. The receipt records all finite universes: 1,099 labelled tournaments through order5; eight order12 palette members and all128 labelled realizations of their C3 hypergraph; all512 depth2 colors under81 literal permutations; explicit depth3 cyclic types; polynomial checks through depth5; and noncanonical-partition/prime-versus-simple-group hostiles. All-order statements rely on proofs, not finite extrapolation.
