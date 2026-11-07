# Low-run padding cannot supply all targets

[Proof](PROOF.md), [independent review](audit/REVIEW.md). The statement applies to any fixed finite bank of nonempty strong cores in the same two-sided padding construction. For every fixed 0<=rho<1/4, targets admitting padding with at most rho log2(N) one-runs form a set of size at most X^(H2(2rho)+o(1)) below X, hence zero density. A fixed run budget gives only a polynomial in log(X) number of targets.

This does not exclude finite-state languages allowing linearly many runs, does not prove that every individual target needs many runs, and does not settle high-complexity occupancy. The accompanying reversal-product identities retain their logarithmic correction and explicitly credit the inherited one-apex mechanism.

From this directory run normal Python3, without -O (the algebra checkers intentionally use assertions):

    python3 verify.py
    python3 audit/check_entropy.py
    sha256sum -c SHA256SUMS

The first output matches VERIFY.json, the second audit/ENTROPY_CHECK.json. Their checks concern identities and binary-run counts, not a new spectrum scan. The local core-state tables are byte-identical to the [audited Session216 table](../../session216_overnight_constructions_20261107/three-core/CORE_STATES.txt). The original proof and source-input hashes are preserved under audit/; the release proof changes only its audit-status label.
