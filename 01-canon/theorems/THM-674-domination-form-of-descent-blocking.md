# THM-674 — The domination theorem for descent blocking: blocked mod k ⟺ the inverse-class set T-dominates the ±-class group (T2 is the j=1 case; k=29 is cycle domination)

**Status:** CLAIMED (stub; proofs + machine verification landing this session)
**Source:** mac-mini-2026-07-09-S65 (cont. 3)
**Depends on:** THM-668 (C2 descent), THM-672 (window anatomy — becomes the j=1 case).

## Claims (proofs landing)

Modulus `k > 14`, danger radius `j = ⌈k/14⌉ − 1`, residues R (13 entries, no zeros).

1. **Prime k (the domination theorem).** In `G = Z_k^*/{±1}` (cyclic, order `m = (k−1)/2`),
   with occupied classes `C`, `D = C^{-1}`, `T = {classes of 1, …, j}` (|T| = j):
   **blocked ⟺ T·D = G** ⟺ (log coordinates) `⋃_{t ∈ ind T}(ind D + t) = Z/m`.
2. **Corollaries:** j = 1 ⟹ T2 of THM-672 (all classes hit). k = 29 (2 primitive, ind 2 = 1):
   blocked ⟺ `ind D ∪ (ind D + 1) = Z/14` ⟺ **no two consecutive holes** — the cycle-domination
   statement, proved. Per-cycle counts: blocked j=2 prime ⟹ `#classes ≥ Σ⌈ℓ_i/2⌉ ≥ m/2`
   (k=29: ≥7 of 14; k=31: ord(2)=5 ⟹ 3 cycles of 5 ⟹ ≥9 of 15; k=37: 2 primitive ⟹ ≥9 of 18;
   k=41: ord₂̄=10 ⟹ 2 cycles of 10 ⟹ ≥10 of 20).
3. **General ledger (all k):** `|A_r \ {0}| = (g−1) + 2g·⌊j/g⌋`, `g = gcd(r,k)` — for j ≥ 2
   non-units reach the danger elements divisible by g (the [15,28] Lemma-1 nesting breaks);
   necessity ledger + stratified composite characterization.
4. **The spread-vs-concentration tension (the remaining core, stated exactly):** window
   blocking (j=1) demands torsion CONCENTRATION — which covering supplies (THM-672); k ≥ 29
   blocking demands class SPREAD in dominating pattern (≥ m/2 of the classes). A full dodger
   must do both simultaneously at every modulus dividing its pair sums.

Verification: `lrc14_domination_theorem_macmini_S65cont3.py` (landing).
