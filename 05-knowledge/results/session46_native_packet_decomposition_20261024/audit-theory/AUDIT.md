# Review of the conditional packet-decomposition theory

Status: PASS within the explicit native-input hypotheses. This is an independent
proof review, supplementary finite controls and producer replay, not an assertion
that a native supply exists whenever a finite-group label can be written down.

## Proof review

1. The common-target current criterion is valid because the actual source rows
   and two target rows are disjoint, the target marked charges agree, and every
   varying residual constraint has been included. Equal group cardinality alone
   is never used to identify actual labels.
2. Minimal balanced subtrade extraction yields a conformal partition. Every
   partial sum remains between the endpoint bits. Under the explicit disjoint
   normalized-resource and joint-pin hypotheses this authenticates each native
   intermediate; arbitrary unresolved factor constraints would invalidate the
   simplified two-equation theorem.
3. Pairing a primitive balanced trade produces a minimal zero-sum sequence,
   giving support at most twice the Davenport constant. The converse is correctly
   rejected by a cross-pairing hostile. The p-group proof uses the actual group
   algebra augmentation ideal and coefficient of the identity, proving D=132
   for (C2)^4+C128. The upper bound 264 is not native sharpness or a constant
   across tournament orders.
4. The partial-sum buffer construction cancels both margin and charge at each
   event, using independently authenticated actual buffer slots. Primitive
   matching makes its intermediate charge positions distinct. It returns the
   buffer and so is a reversible-borrowing theorem, not a consumptive rule
   forbidding every later return to a released port.
5. Final independent review caught a normalization ambiguity before publication.
   The corrected full incidence version writes the actual perturbation as
   Delta L=2^b F and tests the DIVIDED F against ell_b, or equivalently tests
   actual Delta L against ell_0. Applying both doublings is explicitly forbidden.
   All native image certificates were already using the correct factor.
   The full incidence version retains all physical rows. Dickson's lemma gives
   finite conformally minimal lattice atoms, including modular constraints;
   subtracting one atom reduces the integer l1 norm. Endpoint coordinate bounds
   then imply feasibility of every partial sum. The proposed one-hot factor
   encoding is exact for integral assignments with all allowed tuples retained,
   not a claim about a fractional local relaxation or an efficient compiler.
6. The multi-target cyclic margin hostile correctly prevents a charge-only
   universal atom bound. It is explicitly abstract, not asserted native supply.
7. Integer lattice generation alone is correctly separated from bounded-box
   connectivity. This is consistent with the separately audited full28 native
   counterexample in the move-lattices packet.

## Independent arithmetic controls

check.py directly multiplies elements of the group algebra over F2 for C2, C4,
C2+C4 and C2+C2+C4. For every augmentation power it computes the span by binary
elimination and compares its dimension to the truncated-polynomial degree count.
The last nonzero power is one-dimensional and the next is zero in every case.
This independently exercises the nilpotency mechanism rather than merely
checking the printed Davenport formula.

For the actual abstract group Q, it exhausts all multiplicity choices in the
explicit proposed 132-term minimal-zero sequence. Only the empty subsequence
and the full length132 sequence sum to zero. This is not a native port claim.
Additional controls check fixed-pair versus cross-pair primitivity, a generated
integer lattice with no feasible binary generator, and every binary state of
multi-target cycles of lengths3 through9. Forty explicit checks pass in both
normal and optimized modes with identical output.

## Replay and provenance

The theory producer was independently rerun in normal and optimized Python,
with identical receipts matching its frozen output: 965633 explicit checks
plus the inherited original native checker with18321 checks. This replay is
labelled as replay; it does not replace the separate edge-level electrical
filtration or native audit.

Every manifest-listed theory input was hash checked. The referenced all-size
proofs are present in THEOREM.md, with classical antecedents credited in
REFERENCES.md. The review found no unsupported unconditional native supply,
complete-atlas rewiring, arbitrary-target 264 bound or hidden lift capacity.

From the published audit-theory directory:

    python3 -B check.py --packet ../decomposition-theory --out independent.normal.json
    python3 -O -B check.py --packet ../decomposition-theory --out independent.optimized.json
    python3 -B ../decomposition-theory/check.py --repo ../../../.. --out normal.json
    python3 -O -B ../decomposition-theory/check.py --repo ../../../.. --out optimized.json

Normal and optimized receipts must agree byte-for-byte in each pair.
