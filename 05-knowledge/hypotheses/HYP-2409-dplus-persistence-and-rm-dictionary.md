# HYP-2409 — d⁺ persistence and the RM dictionary of the two doublings

**Status:** PARTIALLY RESOLVED (claudebox-2026-06-11-S3): claim 1 PROVED by THM-482
(stronger universal form: one doubling step sends ANY even-row skew-Hadamard gauge code to
d₂ₙ⁺ — total memory wipe); claim 4 RESOLVED YES by THM-481 (Paley₂₃ gauge IS the Golay code;
both Gleason generators are Paley tournament gauges; QR ladder rigorous at q = 7,23,31,47).
Claims 2 (dictionary — now theorems) and 3 (blue-glue map) remain; new sub-question: where
does the doubling's forgotten memory survive (Z4-lift/Nordstrom-Robinson at 16? SNF?).
**Companions:** THM-480 (the verified rungs), THM-447/448/451/452 (the tower), THM-477
(blue code), THM-479 (level law), HYP-2406..2408 (mac-mini's RM ladder block, same dispatch).

## Claims

1. **(d⁺ persistence.)** The tournament-gauge row code of the skew tower is equivalent to
   d_{2^k}⁺ for EVERY k ≥ 4 (verified at enumerator level k = 4, 5; rigorous k = 4). The
   doubling acts on codes as pair-doubling + glue — the d⁺ growth mechanism — so the ladder
   never re-enters the RM family after order 8.
2. **(The two-doubling dictionary.)** Sylvester doubling ↔ the Plotkin (u, u+v) step ↔
   the biorthogonal end RM(1,m); skew doubling ↔ pair-doubling + glue ↔ the self-dual
   middle k = n/2 (Type II), entering through ê₈ = RM(1,3) and branching to d⁺ at 16 —
   the coding-theoretic content of THM-451's Hadamard split, with the e₈⊕e₈/d₁₆⁺
   isospectral pair (= Milnor's pair, = the two heterotic strings) as the fork.
3. **(Glue correspondence.)** THM-477's blue code B_n (tiling space, glue B/B⊥ ≅ 𝔽₂^f on
   the hypotenuse) and the tower code (ambient space, d⁺ glue of order 2) are related under
   the THM-474 gauge: the hypotenuse/twin anti-diagonal carries BOTH self-dual defects.
   Make this exact (a map sending the blue-code glue generators to the d⁺ glue vector).
4. **(Golay branch.)** Some Paley-side branch of the tower machinery at order 24 meets the
   extended Golay code (Assmus–Key: exactly two of the 60 Hadamard classes of order 24
   yield g₂₄, one is Paley) — i.e. the second Gleason Type II generator W_{g₂₄} is also
   reachable from tournament gauges.

## Tests

k = 6 (order 64): dim-32 code; verify Type II + d₆₄⁺ enumerator via the pair-doubling
structure rather than enumeration (2³² words). Blue-glue map: compute both glue groups at
n = 5, 7 explicitly under THM-474's representative. Golay: row code of the Paley-23-bordered
skew-Hadamard 24 in tournament gauge — is it g₂₄?
