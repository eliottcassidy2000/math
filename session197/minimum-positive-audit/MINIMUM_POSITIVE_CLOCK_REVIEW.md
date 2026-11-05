# Independent review: globally least positive-unit ternary amplitude

Status: PASS. No correction or remaining blocker was found.

## Frozen packet

Reviewed source directory: session197-prime-atoms/minimum-positive-clock.

- PROOF.md: 819dfdfc8ea3ec095a82fbc530ddd0bc8ae90ad6fdb15aeb7ae3f364b07dca96
- verify.py: 2e22bc6890c8de8078a6b278b9bbc5542a4c855babf2a1c6d70d351ccab22172
- verify.cpp: 182ec1301f7c8ebbf4af764b4dd760e38f760470150966196a1d465275150640
- receipt-normal.json: d1b6732b677c756cd338fff26800eeb98bb6e05ffd3c65f71beb9ed7777feee6
- strong-bank.jsonl: 15fb4316f3555451e003a7666f75849a4a6d08cd83b11d32ff9027a537cb0b25

The theorem is exactly scoped to transitive cloning of one fixed marked root. The least positive amplitude in A*3^m+1 is72, even if the identity is required only eventually. This does not assert a minimum for unrelated graph operations or solve the full Hamilton-count spectrum.

## Proof audit

The positive root-skeleton formula forces the largest run index to be3. Hence A is the integral coefficient a3, the identity extends from a tail to every positive multiplicity, and the exact finite conditions are a1=3A+1, a2=3A, a3=A. Oddness of a1 makes A even. The two clone-count equations in the packet are algebraically equivalent to these conditions after the exact matching-rank-two reduction.

The matching/run-degree lemma was checked in its named dependency MINIMUM_TERNARY.md. Its lower-bound construction preserves the chosen matching edge as an ordered subsequence during insertion; cutting before the first outneighbor and after the last inneighbor gives valid outer fragments. The upper bound follows from vertex-disjoint crossing edges.

The SCC reduction works with the positive sign. External component counts divide both3A+1 and9A+1, whose gcd is1 because it divides2 and both are odd. A positive-amplitude source cannot have a singleton root SCC. Therefore only neutral singleton padding may lie outside its strong root core.

For A<72, evenness gives A<=70 and H<=211. The exact inequality211^3=9,393,931<9,765,625=5^10, combined with Busch's stated strong-tournament bound, forces order at most10. The official theorem and its uniform consequence were directly checked during the earlier negative-minimum review; the positive packet uses that same consequence correctly.

The maximal-proper-strong-subset deletion proof is valid, including the two-vertex complement case. Reversing this deletion chain gives complete mixed one-vertex extension generation. Induced-count monotonicity justifies the H cap and also the C++ early rejection on a completed subset count. In the forward recurrence, every proper subset is processed before the current mask, so the subset count used for that rejection is complete. No graph above the cap can later fall below it.

Exact degree-group canonicalization enumerates every isomorphism-compatible ordering. The finite bank therefore covers every possible smaller-amplitude strong core. The attaining nine-vertex source has the displayed rank-two vector(217,216,72), proving the all-order72 ray rather than just its first few values.

## Replays

The source verifier was run normally and with Python optimization. Both outputs reproduce the frozen receipt byte for byte. Its bank contains2,288 strong graphs and its independent Python stage checks2,642 relevant rank-two roots, with30 first-equation survivors and no full survivor.

The generated inventories are1,1,6,35,353,1074,644,174 at orders3 through10. As a separate low-order coverage control, the reviewer used the previously generated uncapped complete tournament bank through order7, recounted all530 representatives at orders3 through7, and filtered by strongness and H<=211. The resulting exact canonical sets match the new strong bank at every such order, including all353 order7 classes. This control starts from all tournaments, including weak ones, rather than from the capped strong-deletion construction.

## Independent root counter

A separately written C++ checker uses an exterior-gap polynomial recurrence, with no explicit root-revisit state. For each order of the nonroot vertices, a gap contributes either the direct arc or one root insertion; boundary gaps optionally contain the root. Subset accumulation of these literal gap factors yields the full coefficients through degree4.

This checker independently examines all2,642 relevant roots in the saved bank, checks their ordinary Hamilton coefficient, verifies the zero fourth coefficient, and tests the coefficient equations directly. It reproduces30 first-equation survivors and zero full hits. Optimized and unoptimized assertion-enabled builds give identical output.

Independent files under audits/minimum-positive are:

- gap_check.cpp: 22d2479f25231c6fc85168ab4ddd39efde30207fd5cd3bf56595c699741cbf37
- gap.json: 6dfa3f3ae4bea274c8f17ac7b5b2a93af9769b25f7e9460eaf574189d5bdfc0a
- order7.json: 677526178f6ec2890d7c21f08406533552e87c8c527ec5c8c78b7a4a0253d0d2

The proof supplies the all-order reduction. The complete finite computation then resolves exactly the bounded universe required by that reduction. No finite-to-infinite extrapolation is used.
