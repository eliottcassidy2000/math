# THM-672 — The descent-window torsion-occupancy theorem: on k ∈ [15,28], blocking REQUIRES torsion classes; composite all-unit descents are NEVER blocked; wall primes are blocked iff every ±-class is hit

**Status:** CLAIMED (stub; proofs + machine verification landing this session)
**Source:** mac-mini-2026-07-09-S65 (cont. 2)
**Depends on:** THM-668 (pair-sum ruler theorem; certificates C1/C2). Targets the large-domain
statement "covering + sums > Q₀ ⟹ C1/C2 fires" (HYP-5730 obligation (2)).

## Claims (to be proved this session)

Fix the descent modulus `k` with `15 ≤ k ≤ 28` (so danger = `{0, ±1}` mod k), residue multiset
`R = {v_l mod k}` with no zero (else the descent is dead).

1. **Master ledger.** With occupied ±-classes `C` and per-class `g = gcd(class, k)`:
   blocked ⟹ `Σ_C w(g) ≥ k − 1`, where `w(1) = 2`, `w(g) = g − 1` for `g > 1`.
2. **Unit-pigeonhole theorem.** If `k` is COMPOSITE and all residues are units:
   `Σ ≤ 2·(φ(k)/2) = φ(k) < k − 1` ⟹ NEVER blocked (band solvable, C2 fires whenever `k`
   divides a pair sum).
3. **Wall-prime characterization (exact).** For prime `k ∈ {17, 19, 23}`:
   blocked ⟺ `0 ∈ R` or `R` hits EVERY ±-class mod `k`.
4. **Per-k torsion-occupancy necessary conditions** (from 1): e.g. blocked mod 16 ⟹ some
   `v ≡ 8 (16)`; mod 18 ⟹ some `v ≡ 9 (18)`; mod 21 ⟹ some `v ≡ ±7 (21)`; mod 22 ⟹ some
   `v ≡ 11 (22)`; mod 26 ⟹ some `v ≡ 13 (26)`; mod 27 ⟹ some `v ≡ ±9 (27)`; mod 25 ⟹ some
   `v ≡ ±5, ±10 (25)`; full mechanical list per k. The expensive classes are the SMALL-TORSION
   points of Z/k — blocking requires torsion occupancy.

Verification script: `lrc14_torsion_occupancy_macmini_S65cont2.py` (landing; cross-checks against
the exhaustive blocking census of S65 cont.).
