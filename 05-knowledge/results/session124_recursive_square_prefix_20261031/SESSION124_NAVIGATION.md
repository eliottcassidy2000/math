# Session124 checkpoint: explicit square-sum prefixes

**VERIFIED / independently audited.** The [exact theorem and minimal proof](README.md#theorem) show that, for every l≥1 and every N≥(107·49^l−1)/2, a square-sum Hamiltonian cycle on[1,N] contains the original unscaled P_l, a prefix of one explicit infinite square-sum ray. This is not arbitrary prescribed-path embedding.

The core uses only the [existing B-base hashes](BASE_DEPENDENCIES.json),98 new reusable paired transitions, and98 new terminal-cycle states. [Portable independent reproduction](INDEPENDENT_REPRODUCE.md) checks these dependencies without old transition tables or producer imports. Wider simultaneous-copy corollaries are not included in this minimal checkpoint.

**OPEN.** Tournament infinite completeness—the claim that every positive odd integer other than7 and21 occurs as a raw directed Hamiltonian-path count—remains OPEN. This checkpoint does not prove that claim, change its scope, or substitute an automorphism-normalized count. Session124's separate tournament work continues.

This is a package-local navigation note. Maintained global navigation is intentionally untouched pending reconciliation against fresh main.
