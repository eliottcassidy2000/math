# Core dependency status clarification

**PROVED, self-contained core.** `PARITY_INDUCTION.md` supplies an elementary
induction proving the classical tournament oddness theorem and pair-weight
integrality directly from the unhalved cut identity and same-order arc
reversal. Oddness also supplies path existence because counts are finite
nonnegative integers. Together with the existing block-insertion proof,
the core theorem needs no externally imported mathematical theorem.

This is an additive clarification to the cut-response packet published at
checkpoint `f5076a605f36777c190499e26c00f9271795acc3`, as identified by the
publishing parent. The original README, every program, and all verification
receipts remain untouched. The original external-dependency wording is
superseded by this note; no mathematical claim is retracted.

The optional triangle-packing corollary retains its cited odd-cycle
collection formula dependency. Classical attribution to Rédei is retained;
no novelty claim is made. All-order spectrum completeness remains open.
