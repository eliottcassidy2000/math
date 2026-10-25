# Outer hourglass K3: current exact result

**PROVED:** geometric Picard number 19 for every outer-deletion K3 and 58 for
the seven-root hourglass. See [EXACT_RANK_COROLLARY.md](EXACT_RANK_COROLLARY.md).
The seven deletion K3s form exactly two complex isomorphism classes,
distinguished by ranks 20 and 19. The elliptic fibration has Mordell–Weil rank one; the exhibited section is
nontorsion, but generator/primitiveness is not asserted.

[LOWER_BOUND.md](LOWER_BOUND.md) is the independently closed rank-at-least-19
proof snapshot. Its former two-value remainder is resolved by the independent
upper-bound packet, not by genericity. The exact-corollary file is current.

Run `python3 check_lower_bound.py` and `python3 -O check_lower_bound.py`.
`LOWER_BOUND_SHA256SUMS` covers the lower proof. `SHA256SUMS` freezes this
combined packet. `probe_counts.py` is an exploratory independent counter,
not a stand-alone geometric upper-bound proof.
