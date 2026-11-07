# Independent review: numeric14/15/27 and the H27 majority threshold

**Verdict: APPROVED in the stated scope.** The exact order14 ordinary and
converse-merged counts, degree histograms, minimum orders, existing
order14 H=15^6 witness and construction-rate bound are correct. The
additional fixed-H27 first-majority threshold is exactly23, with a
majority at every subsequent order. These are not global spectrum claims.

## Proof review

The complete low-count strong inventory is correctly reduced to order7 by
the cited Busch inequality and5^7>27^3. Completeness of the456-class input
atlas and Busch's theorem remain inherited/cited dependencies, not newly
proved by a limited sample. The exact strong inventory, ordered component
factorizations, singleton-gap counts and converse action are correctly
distinguished.

At even order14 neither H15 nor H27 admits a self-converse SCC word. Thus
the quotient degree is H/Aut on every node considered. Independent literal
word enumeration gives:

- H15:110 ordinary classes,55 merged nodes, degrees {3:5,5:41,15:9};
- H27:210 ordinary classes,105 merged nodes, degrees {1:28,9:65,27:12}.

In particular the28 leaves are only4/15 of the H27 stratum. The asymmetric
order7 atoms means Aut=1; one is self-converse as a strong atom, but it
cannot yield a self-converse even-order padded class. No histogram is
silently taken inside the induced fixed-H subgraph.

The minimum-order claims follow from complete smaller spectra plus explicit
witnesses. Dominance sums of the order14 atom give the stated7/3 upper
rate and mixed-prime construction; neither optimality nor a positive linear
lower bound is asserted.

For the majority addendum, the palindrome subtraction and merged-class
denominator are correct at both parity boundaries, including orders6,7,9.
The even/odd majority polynomials and their two-step differences are exact.
All signs through22 are negative, while O(23)=144 and E(24)=672. Their
positive differences thereafter prove the all-order tail rather than
extrapolating a finite plot. The22nd-order comparison is particularly tight:
280/561 is below one half, whereas336/660=28/55 is above it at23.

## Independent controls

`independent_numeric_audit.py` reads only the bundled atlas. It independently
examines1,346 induced matrices, counts H with a backwards subset recurrence
and canonicalizes the relevant strong cores by ALL permutations, without
degree pruning or producer imports. It recovers all13 atoms, their
automorphism orders and converse pairings. It then literally enumerates
the padded SCC words and quotients by reversal plus component converse,
recovering both complete order14 histograms. A separate backwards DP and
all-source reachability verify H=11,390,625 and strongness of the order14
witness.

From this audit directory, run
`python3 independent_numeric_audit.py --packet ..` and the same command
with `python3 -O`. Both outputs match `independent_numeric_audit.json`.

`independent_majority_audit.py` enumerates341,055 four-gap compositions
through order60, forms their actual reversal orbits, and checks the exact
majority signs. Its normal and optimized outputs both match
`independent_majority_audit.json`. The reviewed parity-polynomial argument,
not these finite tests, proves the unbounded tail.

Both producer checkers were separately replayed in isolated copies in
normal and optimized Python. All output pairs agree byte for byte and
match their respective frozen producer receipts.

Reviewed SHA256 values:

- PROOF.md:142f434af13cc656cc2062a9c5e3deb03b3a7a239f4219c2d9ed4c8341b2eb1d
- verify.py:9d9ba09da6faf5ca37a30f0e06c19480e9dca5721b0e6f3c1974b5b62ad89c6e
- RESULT.json:68b2cb905578662e64d7c134912c983195149e579fafa5ecf1149d11d0ad1c74
- MAJORITY_THRESHOLD.md:603273de810bb4773b91f4b09f90ac65721cdbd2c80c6700f6f18fc9041d27ad
- verify_majority.py:720b7c15d9d6821e6a0df03883a47755a5c01003fb38ccbdda90248e4656983c
- MAJORITY_RESULT.json:f76abe1d05ee811213fd1ee7bff99db0edb7eb57eb6fccce02dddf205db20cd1

## Boundaries

The threshold23 concerns the proportion of leaves inside the fixed-H27
paired-black stratum. It is not a threshold for realizing all large odd
Hamiltonian counts, and it says nothing about arbitrary graph connectivity.
The input atlas is an inherited complete order7 class source, not a fresh
all-labelled order14 census. Native labelled degree, multiplicity degree
in the paired-black quotient and distinct-neighbor degree remain different.
