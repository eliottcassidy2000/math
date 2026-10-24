# Validation and scope

## Exact finite universe

The checker does not extrapolate theorem proofs from finite tests.

- Enumerate every conjunction of Boolean unary/two-literal clauses by finite intersection closure of their truth masks, then deduplicate solution families. For n=1,2,3,4, obtain respectively 4,16,166,4170 families, including the empty family.
- For every accepted ordered endpoint pair of every such family, compile switch implications, compute SCCs, and compare every literal endpoint-interval state with the reconstructed SCC ideals. Across n≤4, there are 144,470 accepted ordered endpoint pairs.
- For all directed graphs whose edges go from higher to lower labels on n≤4 vertices (1,2,8,64 presentations), all C2² charge assignments of total zero, and each declared integer-margin pattern, compare the full accepted interval DAG, its cover DAG, all neutral acyclic partitions, and SCC-merged neutral partitions. There are 12,561 charge/margin instances. Isomorphic posets can have multiple presentations; these are explicitly presentation counts, not unlabelled poset counts.
- Margin patterns are all-zero, alternating +1/−1 when n is even, and (1,0,...,0,−1) when n>1, deduplicated. These are not all possible integer-margin vectors. Costs are deterministic positive, nonuniform block slot and clock sums. The comparison retains simultaneous Pareto pairs, not separate minima.
- Check every accepted comparable ideal pair against its relative induced zero-defect ideal interval, and verify the one-SCC accessibility iff every SCC defect is zero boundary.
- Check chain hostiles for all n=2,...,64 and exhaust their literal Boolean cubes through n=12. Check paired lifts for logical n=2,...,10 and exhaust their full physical Boolean cubes through logical n=6. Larger ranges use the explicit implication-prefix characterization, not hidden 2^128 enumeration.
- Verify D(C2²)=3 by all 64 three-term sequences and a two-term zero-sumfree witness. The all-size obstruction is proved separately in THEOREM.md.

The checker uses `need` raising explicit exceptions, never removable Python assertions. Normal and optimized output must be byte-identical to the pinned receipt.

## Positive and hostile controls

Positive controls include the empty/full endpoint pair, full Boolean cubes, neutral original SCC blocks, arbitrary acyclic neutral packetizations, and exact cycle-merging repair. Hostiles include a ternary table not representable by the interval-ideal model, a neutral family not union-closed, a neutral subset that is not an ideal, and the four-slot directed packet-contraction cycle. These show which stronger shortcuts are actually false.

`portable_check.py` makes a real filesystem copy into a fresh temporary directory, changes the working directory away from the copied packet, verifies manifest hashes, reruns normal and optimized checks, and compares their bytes with the receipt. It separately changes a source file, changes an inherited input, deletes the theorem, and adds either an ordinary unmanifested file or a nested extra MANIFEST.json in disposable copies; each must be rejected. No original source is mutated by these controls.

## Scope and review

The arithmetic is finite exact integer/bitwise arithmetic with no floating point. All source and charge models here are abstract. No theorem is promoted from a finite run; no native source/target tournament, Hamiltonian path, electrical bank, or all-order receiver supply is inferred. A separate independent audit is requested from the parent session; this packet does not label its own checker as independent review.

The bounded startup read used clean repository commit d0187b4caa6d (full hash in SOURCEPINS.json), with its cached origin/main equal to HEAD. This lane did not fetch, publish, or modify that checkout. The optional agent inbox check could not run because that checkout has no `.machine-id`; mathematical work was unaffected. Parent session coordinates live-source freshness and publication.
