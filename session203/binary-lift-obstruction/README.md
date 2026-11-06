# Bounded-offset binary-lift obstruction

Status: PROVED / independently audited in its stated certificate-format scope.

For any fixed offset bound c, finite seeds followed by physical O/Z/J/K1/K2 word lifts B_next=2B+epsilon, |epsilon|<=c, have only O(log M) distinct scalar outputs below M. They eventually reduce to finitely many Z rays. Therefore this natural bounded-shift binary-lifting method cannot certify an odd tail or F(8,137).

This is not a counterexample to F and does not rule out other radices, multi-input constructions, target-dependent rebuilding, or the finite-register presentation. The latter remains a possible route to occupancy.

PROOF.md contains the all-order proof. audit/REVIEW.md is the independent review, and audit/SOURCE_PROOF.md preserves its exact reviewed candidate. PROMOTION.json records the status-only promotion and import-path-only changes to the two independent check scripts; their originals are retained. SOURCE_SHA256SUMS is the original source freeze, not the release manifest. Its PROOF.md entry matches audit/SOURCE_PROOF.md after promotion. Dependency proofs retain their source-relative reference wording. ROOTED_ORDER_BOUND.md Section4 supplies the explicit weak-tail factorization; PROVENANCE.json pins all source copies.

Run from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/independent.py
    python3 audit/weak_factorization.py
    sha256sum -c SHA256SUMS

The first two outputs match NORMAL.json and OPTIMIZED.json. The independent scripts use assertions and should run without -O; outputs match audit/INDEPENDENT.json and audit/WEAK_FACTORIZATION.json. Their controls include 1,628,172 short word tests over 2,640 repeated-gate states, and all 772 positive-start labeled inputs through order five, including 201 weak inputs. These checks do not replace the all-order proof.

No navigation or global index is changed.
