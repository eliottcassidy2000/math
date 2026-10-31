# Pinned Session126 READY observer proof excerpt

Source: ../primegame/fractran_ca/ACTUAL_RADIUS_ONE_CA_REPORT.md, section5.
This is an inherited observer theorem, not a bundled CA implementation.

## 5. Observer theorem and exact macro duration

**Theorem.** If the CA starts in the canonical READY encoding E_o(N) of any N≥1, it reaches another canonical READY encoding E_o(T(N)) after finitely many CA steps, where T is the fixed ordered PRIMEGAME integer update. No intermediate time before that endpoint is READY. Iterating these endpoints recovers the complete FRACTRAN orbit.

**Proof.** The READY transition reads the squarefree denominator guards exactly and selects their least enabled index. In a valid configuration there is one head, so the radius-one rule performs exactly the described finite-head tape move without producing ERROR. Every decremented coordinate was positive when selected. No coordinate occurs twice in the selected fraction's operation list, so it remains positive until its own decrement is executed.

A unary scan reaches the first 0 after exactly e right moves. Incrementing that boundary extends a unary block by one. For a decrement, the predecessor of the first 0 is its final1; clearing it shortens the unary block by one. All other tracks stay unchanged. RETURN finds the persistent origin, so every operation ends there. A zero-height increment stays at the origin; a height-one decrement clears the origin bit of that counter without moving left. These boundary cases, and RETURN's marker test, prevent the valid head from crossing left of the origin. Each operation preserves unary shape, the unique head, the marker and the spectator bits.

The final counter vector is exactly the initial vector plus the selected fraction's exponent difference. The last operation alone restores READY, yielding the canonical encoding of T(N). This proves the claim and supports induction over any number of macro updates. ∎

**Exact cost.** An operation on a counter with initial height e uses exactly **2e+1** CA steps, for both increment and permitted decrement. The READY choice itself costs one step. For selected fraction i, let D_i be the union of its numerator and denominator prime supports, and k_i=|D_i|. Distinctness of the touched coordinates gives

`τ_i(e) = 1 + k_i + 2∑_{j∈D_i} e_j`.

Since k_i≤5 and M=max_j e_j,

**τ_i(e) ≤ 10M+6 = O(M+1)**.

This is a per-FRACTRAN-update bound, independent of the spectator factor's bit length. It is not a bound on total prime-generation time. The timing formula and the observer relation are explicitly checked in the executable tests, not just stated in prose.

