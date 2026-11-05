# Independent audit of the global minimum ternary clone amplitude

Result: PASS. The least positive amplitude for the specified single-root transitive-clone identity A*3^m-1 is 60, including when the identity is imposed only eventually. This audit does not extend the claim to other physical compilers.

Frozen inputs and SHA-256:

- `session197-prime-atoms/MINIMUM_TERNARY.md`: `45d183b1df1f030aff1bf755bd2b5b30ce441e8a5c13e01b99b6ef8613957c54`
- `session197-prime-atoms/verify_minimum_ternary.py`: `2ad41156b91eeb868170ee688cdd08f1ba2065d00948a812200ea990e1c9e27f`
- `session197-prime-atoms/verify_minimum_ternary.cpp`: `fd35f6b118fd2d1e651db1e2c9d3690ac4a89d7f5bfe8ba18bfb19c2ba08f000`
- `session197-prime-atoms/receipt-minimum-ternary.json`: `703249b318f4203d43233590cb4cf8e722c3c018ffa5d5f6aec7d5f8102d472a`

## Mathematical audit

I verified the run decomposition and the converse direction of the matching formula, including the delicate extension of one matched edge: insertion preserves its endpoint order, so the first-outneighbor/last-inneighbor cuts are ordered correctly; the discarded prefix and suffix have exactly the required in/out types. Distinct exponential terms then force matching number two and the stated exact coefficients, even from an eventual identity. The two clone equations are sufficient only after that matching bound, as the proof correctly requires.

The SCC reduction is valid. The component containing a nontrivial root remains strong after cloning. Every other component's Hamilton count divides both 3A-1 and 9A-1, whose gcd is one for even A. A nontrivial strong component has at least three paths, so all discarded components are singletons.

I directly read the specifically cited official Busch paper at https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf . Theorem 1 supplies the claimed safe inequality 5^(n-1)<=H^3 for strong tournaments; its piecewise bound is stronger. Thus the integer comparison excluding order at least eleven under H<=173 is sound.

The strong-deletion lemma is correct, including its two-outside-vertices case. Mixed one-vertex extensions of a strong graph are strong, and the lemma ensures that such extensions generate every strong graph. Induced Hamilton-count monotonicity makes the cap safe along an entire deletion chain. Degree-group minimization is exact isomorphism canonicalization, not a proxy invariant.

## Code and replay audit

I inspected both source implementations. The C++ cap checks the completed Hamilton count of an induced subset, so its early exit is justified by induced monotonicity. Its arrays cover every generated or cloned order, and the bounded intermediate sums are safely within integer range. The exhaustive mixed-mask loops, degree-class permutations, matching DP, root residue filter, and clone checks implement the proved search. Python uses a different matching algorithm and backward endpoint recurrence, checks every retained graph, and independently enumerates all labelled small strong graphs.

My full normal and Python-optimized runs reproduce the supplied receipt byte-for-byte, including 1,466 retained strong representatives, 1,624 candidate roots, 20 survivors of the first clone equation, and zero survivors of both equations. The deterministic bank hash is `03f1d7a24d36fc81a5b0e88e9c637668fe16d0582e9001200dba46af3874283e`.

As an additional independent coverage check, I filtered the separately generated complete uncapped bank of all 456 unlabelled order-seven tournaments. Its strong representatives with H<=173 are exactly the same 351 rows as the strong-deletion generator's order-seven bank. This checks a different construction of that entire inventory.

The attaining source and its all-order run proof agree with the separately audited three-port construction. All frozen hashes remained unchanged. No correction or remaining blocker was found.
