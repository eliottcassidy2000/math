# Session 30: parity propagators and exact core compression

## Scope and outcome

This packet changes no source engine or repository files. It transforms the 31 audited weight-one and weight-two decorated-core kernels from the sibling `session30-decorated-engine` packet. It is an exact reformulation of that finite expansion, not a proof of its uniform asymptotic remainder.

Main exact counts:

- Four old channels: 2,364 nonzero symbolic monomials.
- Four parity channels: 553 nonzero symbolic monomials.
- After graph relabeling and commutation of parallel-edge factors: 263 contraction monomials.
- Full dense tensors: 31,396 slots; per-edge palette tensors: 12,970 slots.
- Every coefficient with an odd number of odd channels vanishes exactly.

The 553-term objects are literal transformed coefficient tensors. The 263-term objects are only contraction-equivalent tensors on skew-matrix inputs. These are different kinds of reduction and must not be confused.

## Exact rational identities

Let X=zM, with I-X^2 and I+X^2 invertible, and define

G=-X^2(I+X^2)^(-1), H=X^2(I-X^2)^(-1),
O=X(I-X^2)^(-1), P=X(I+X^2)^(-1).

The old propagators satisfy

R+=H+O, R-=H-O,
Q++=-H+P-O, Q--=-H-P+O, Q+-=-G.

Proof: substitute the definitions into R+=X(I-X)^(-1), R-=-X(I+X)^(-1), and Q_ab=-R_a R_b (I-X^2)(I+X^2)^(-1). Each difference is the zero scalar rational function. All factors are rational functions of the same matrix, hence commute, so the rational identities apply to matrices on their common domain. `check.py` verifies all five with exact SymPy cancellation.

Under M -> -M, G and H are even and O and P are odd. If M is skew-symmetric, G and H are symmetric and O and P are skew-symmetric; this remains transpose symmetry, not conjugate-transpose symmetry, for complex z. Hence odd-channel diagonals vanish.

## Exact tensor verification

`transform.py` applies the integer linear substitution independently at every edge slot, collecting coefficients with Python arbitrary-precision integers. It preserves each original automorphism factor/denominator. It certifies that no odd-O/P coefficient remains in any of the 31 transformed tensors. This is a finite exact certificate, not a numerical test or a general all-weight theorem.

An additional exact support property holds across these tables: all 37 loop-edge occurrences have palette {G,H} or a subset, and all 35 bridge-edge occurrences have palette {G,P} or a subset. At weight one every bridge uses G alone. The theta graph uses {G,H,O}, with no P. `check_structure.py` checks bridges by actual edge deletion, including multiplicity, so a member of a parallel bundle is not misclassified as a bridge. These support statements are certified for the present 31 tables only.

## Graph quotient proof

For a core with a fixed ordered edge list, rename the summed vertex indices by any automorphism of its uncolored multigraph. Reorder the resulting edge factors into the original edge order, allowing parallel-edge factors to commute. Reversing an edge contributes a sign precisely when its channel is O or P. Thus each colored monomial equals its canonical colored representative times that sign. `quotient.py` accumulates exact coefficients with these signs. It also checks whether a representative has both signs, in which case its contraction is identically zero; no such extra cancellation occurred here. The core's original denominator is retained without modification. This proves contraction equivalence of the 553- and 263-term forms for skew M.

## Numerical verification

- Four source fixtures, including transitive, regular, and almost-regular cases; real and complex z: 248 individual core comparisons for the quotient evaluator, maximum absolute error 2.67e-15.
- The palette-compressed dense evaluator independently makes the same 248 comparisons, maximum error 2.23e-15.
- Three deterministic weighted skew controls of orders 5,8,12, at real and complex z: another 186 core checks comparing sparse quotient, compressed dense, and global reversal against the source engine. Maximum error 4.45e-16.

These are floating-point implementation checks in addition to the exact algebraic certificate. They are not proofs of numerical conditioning near poles.

## Performance, including failed approaches

Sparse termwise contraction of all 263 quotient monomials was slower than the existing optimized four-basis engine: source/sparse median time ratios were 0.278, 0.443, 0.102 at n=10,20,40. Therefore a 4.27-fold symbolic reduction is explicitly not a 4.27-fold runtime reduction.

`parity_engine.py` instead compresses the palette axis separately at each edge and preserves the source engine's special K4 and doubled-square contractions. A first generic-greedy contraction plan was pathologically slow; the final implementation explicitly permits up to 2^24 intermediate elements during path search. This is an intermediate-size planning cap, not a complete process-memory guarantee.

Warm seven-repeat single-thread BLAS medians for weight two on deleted cyclic regular tournaments:

| n | Source seconds | Compressed parity seconds | source/parity |
|---|---:|---:|---:|
|10|0.00453|0.00406|1.11|
|20|0.01172|0.01004|1.17|
|40|0.07926|0.07489|1.06|

This is a modest observed benefit, sensitive to hardware, planner, matrix size, and caching. It is not a robust universal speed guarantee. Compilation/path-search time is excluded. The comparison preserves existing specialized contractions and includes propagator construction in both timed routes. `benchmark_hybrid.py` records the full measurements.

## Reproduction

From the common parent directory, with Python, NumPy and SymPy available:

    python3 session30-parity-basis/transform.py
    python3 session30-parity-basis/quotient.py
    OPENBLAS_NUM_THREADS=1 python3 session30-parity-basis/check.py
    OPENBLAS_NUM_THREADS=1 python3 session30-parity-basis/benchmark_hybrid.py
    OPENBLAS_NUM_THREADS=1 python3 session30-parity-basis/check_structure.py

The evaluator depends on the sibling audited engine and intentionally does not replace it. Exact JSON files preserve the complete coefficient and graph data. The sparse-quotient JSON is not suitable as a drop-in literal tensor identity for arbitrary non-skew inputs.

## Interpretation and next proof question

The even propagators and antisymmetric currents provide a clean algebraic separation of reversal parity. They do not establish a physical electrical-network equivalence. The most focused next structural question is whether the observed loop/bridge channel restrictions admit an all-weight local proof, and whether blockwise parity contraction can exploit them without the overhead seen in naive monomial summation. No claim about uniform z=1 error control or the Collatz bounds follows from this packet.
