# Independent score-gate and permutohedron review

**Verdict: PASS.** Final proof SHA-256: `583c1d2d112141c6e998e58369b077dde00952e6ab3c5c9458069734b98c9fc5`.

## Mathematical review

- The cut-slack identity is exact: a subset's Landau slack counts its outgoing cross arcs. Tight subsets are precisely terminal unions of strong components. The labeled score vector therefore determines the ordered labeled SCC partition and its minimal product face, but not internal arc orientations or Hamiltonian counts.
- The gate score map sends the entire closed core score permutohedron into the larger polytope's relative interior for a **nonempty** core and `0<t<=1`. Its exact minimum is unnormalized facet slack `t`, not Euclidean distance. The empty-core case was explicitly excluded after review.
- The full affine slice at gate scores is larger than the embedded core permutohedron. The relaxed vertex vector `(-1,1,2,...,n-2,n)` and all proper lower-sum bounds are correct. Strict containment requires `n>=2`; the singleton is a degenerate equality case.
- Directed-triangle reversal preserves scores and replaces the gate cycle count `H(T)` by `H(T-v)`. The forced segment and exact gate classification are valid. A mixed core vertex removes every gate. There is no universal strict count decrease: transitive cores are the stated equality boundary.
- The directed-four-cycle reversal also preserves every labeled score. Its forced entrance and exit yield the precise two-component endpoint-refined cover formula. The gate survives exactly when the original `u` is a sink and `v` a source. Strongness is explicit.
- The literal `HC=7` and `HC=21` witnesses have the stated matrices, scores and path/cycle counts. These are **cycle** values of ungated tournaments, not forbidden Hamiltonian-path values. They prove score-fiber loss of the gated class, not unrestricted spectrum completeness or an arithmetic coincidence.
- Historical attributions and corrected score/fiber/face distinctions are scoped carefully. No priority claim is made for Landau inequalities, triangle score moves, or earlier non-score-determination results.

## Producer replay

The final verifier was copied into an independent replay directory and run normally and under `python3 -O -B`. Both stdout files match frozen `RESULT.json` byte for byte:

- `verify.py`: `d774bc0e9dd853c2e11a1828b44ac5e298c847a301344d2896b182f6c9e4381c`
- `RESULT.json`: `451899d4228324b3815b2a0d02e8e18e200c1878ecccaee28ac09295e0d5a440`

The producer's scope includes all 1,099 labeled cores of orders 1..5, all 5,405 pointed triangle reversals, all 10,650 specified four-cycle switches, full score fibers through order five, the stated slice vertices, and direct independent permutation controls for both literal gap examples. The exact finite order-five ambiguity threshold is not inferred merely from an example.

Every entry of `SOURCES.json` was independently checked against its Git blob at the immutable reference commit.

## Independent implementations

None of these scripts imports producer functions:

1. `independent_score_gate.py` uses adjacency matrices, BFS, subset-DP counts and explicit gate search. It verifies all 1,099 cores and 5,405 triangle reversals, including singleton boundaries. Its minimum ungated strong order of five is explicitly restricted to orders at least three; the singleton convention is not hidden.
2. `independent_four_switch.py` reconstructs both supplied cycle reversals from their rows, verifies unchanged labeled scores, and checks the two literal witnesses by three separate cycle counts: subset DP, all rooted cyclic permutations, and endpoint-permutation/two-path convolution. It also checks Hamiltonian-path counts by subset DP and full permutations, and strongness and absence of gates directly. The resulting pairs are `HC 11->7, H 111->91` and `HC 37->21, H 389->293`.
3. `independent_gate_geometry.py` uses exact rational coordinates. It checks 612 core-vertex/time cases for `n=1..5` and `t=0,1/3,1/2,1`, and all 872 relaxed-slice vertices for `n=2..6`. Exact minimum Landau slacks agree with the proofs.

Each independent script was run normally and under optimized mode, with identical JSON outputs. The scripts and both result variants are approved for publication beside this review. All-order conclusions rely on the reviewed proofs; the censuses are reproducible controls.
