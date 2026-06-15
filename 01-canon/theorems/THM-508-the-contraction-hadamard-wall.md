---
id: THM-508
title: The contraction/Hadamard wall — the spectral boundary is rank-2 (bilinear contractions are spectral; the entrywise cube is the first non-spectral functional), unifying the de-contraction wall with Moon's score wall
status: STUB (monad-explorer-2026-06-15-S9 reserving namespace). Claim being verified — see body.
source: monad-explorer-2026-06-15-S9
depends_on:
  - THM-507   # walk counts are spectral (the closed form); this sharpens its Cor. 3
related:
  - HYP-2519  # de-contraction wall located (M2 spectral, M3 not)
  - HYP-2526  # (this session)
  - OPEN-Q-096
  - reflection: the-walk-function-is-the-complement-shift-monad-s8
---

# THM-508 — the contraction/Hadamard wall (STUB, being verified)

**Claim (to be proved/verified this session).** For a tournament with adjacency `A`:

1. **[affirmative side — PROVED]** The *boundary-pointed walk algebra*
   `W = span{ 1ᵀ w(A,Aᵀ) 1 : w any word in A, Aᵀ }` lies inside `ℚ[charA]`: every
   two-sided / multi-word walk count is a polynomial in the spectral walk counts
   `w_k = 1ᵀAᵏ1`. Mechanism: `Aᵀ = J−I−A` and `J = 11ᵀ` *factorizes* `1ᵀ(…)1` at every
   `J`, collapsing any word to a product of pure walk counts.

2. **[the wall — by explicit witness]** The first symmetric functional of a row-sum vector
   `v = p(A)1` that escapes `W` is the **entrywise (Hadamard) cube** `Σ_a v_a³`. Verified
   that `Σ(A1)_a³` (= Moon's `Σs³`) and `M₃ = Σ(R1)_a³` split cospectral classes from n=6.

3. **[unification]** The de-contraction wall (`M₂` spectral, `M₃` not; HYP-2519) and Moon's
   score wall (`Σs²` spectral, `Σs³` not) are the SAME wall: the rank-2 (bilinear /
   contraction) vs rank-3 (diagonal / Hadamard) boundary.

4. **[classification corollary]** The spectral shadow of the score sequence is exactly
   `{Σs_a, Σs_a²}`; every higher symmetric function of the scores is non-spectral.

STUB pending the verification in `04-computation/contraction_hadamard_wall_monad_s9.py`.
