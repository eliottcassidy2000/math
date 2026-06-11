# THM-480 — The skew tower's tournament-gauge row code is the d⁺ ladder (ê₈, d₁₆⁺, d₃₂⁺)

**Status:** VERIFIED exact (orders 8, 16, 32; order-16 identification RIGOROUS, order-32 at
enumerator level) + mechanism proved at the generator level. claudebox-2026-06-11-S2.
**Provenance:** dispatch = RM duality (RM(r,m)⊥ = RM(m−r−1,m); "codes with k = n/2 are
self-dual"). **Companions:** THM-447/448 (the tower), THM-451 + hall_class_pin_kps1 (the
order-16 Hadamard split, tower-16 ≅ had.16.3 by exact witness), THM-477 (the blue code —
the OTHER self-dual-defect code of this dispatch, living in tiling space), HYP-2409.
**Literature:** Assmus–Key (codes of the 5 Hadamard classes of order 16: dims 5,6,7,8,8);
M. Hall, Combinatorial Theory §17.3; Armario–Frau arXiv:1311.2637 (self-dual codes from
skew-type matrices at the EW bound); errorcorrectionzoo (augmented Hadamard = RM(1,m)).

## Setup

Tower: H₈ = border(Paley₇), H_{2m} = [[H, H], [−Hᵀ, Hᵀ]] (THM-447 doubling at the
±1 level; skew-type is preserved: H + Hᵀ = 2I). The TOURNAMENT GAUGE is this skew
normalization M = I + S itself — no row/column renormalization. The binary code
C(H) = row span of (J − H)/2 over 𝔽₂ (1 ↔ entry −1).

## Statement

1. **(k = n/2 self-dual, every level.)** C(H₈) = ê₈ = the extended Hamming [8,4,4]
   = **RM(1,3)** — exactly the r = 1 self-dual Reed–Muller point RM(r, 2r+1).
   C(H₁₆) is a [16,8] doubly-even self-dual (Type II) code; among the only two such codes
   it is **d₁₆⁺** (indecomposable: the weight-4 support graph is connected), not e₈⊕e₈.
   C(H₃₂) is a [32,16] Type II code with the **d₃₂⁺** weight enumerator
   (A₄ = 120 = C(16,2); distribution 1, 120, 1820, 8008, 45638, …) — in particular NOT
   RM(2,5) (which has A₄ = 0, A₈ = 620).
2. **(Contrast: the Sylvester gauge.)** Sylvester rows span only RM(1,m) (with 𝟙: dim
   m+1, the biorthogonal/low end). The tournament gauge pumps the row space from dim m+1
   to dim n/2 = 2^{m−1} and lands it exactly on the self-dual middle — the dispatch's
   "k = n/2 ⟹ self-dual" realized by skewness.
3. **(Mechanism: doubling = pair-doubling + glue.)** In [[H, H], [−Hᵀ, Hᵀ]] the top block
   rows binarize to (b, b) — coordinate-PAIR-doubled words, the defining d-code mechanism —
   while the bottom rows supply glue (Hᵀ = 2I − H makes them complement-plus-point-mass
   words: the e_i twin corrections are the top-RM-degree point indicators). So the tower
   code grows the d⁺ way, and the ladder ê₈ → d₁₆⁺ → d₃₂⁺ (→ conjecturally d_{2^k}⁺, HYP-2409).
4. **(The split at 16, explained at code level.)** THM-451/kps1: the tower leaves the
   Sylvester Hadamard class exactly at order 16, with invariants unable to separate
   had.16.3 from had.16.4 (witness search needed 10⁷ nodes to exclude 16.4). Those two
   classes carry precisely the two Type II [16,8] codes — e₈⊕e₈ (had.16.4, decomposable)
   and d₁₆⁺ (had.16.3, indecomposable) — which share one weight enumerator
   (1 + 14y + y²)², y = x⁴: the coding avatar of the Milnor isospectral pair E₈⊕E₈ / D₁₆⁺
   (and of the two heterotic string theories). The tower picks the indecomposable D₁₆⁺ /
   SO(32) branch; decomposability of the weight-4 support graph separates in milliseconds
   what the monomial witness search separates in 10⁷ nodes.

## Verification (script skew_tower_dplus_code_cbx2.py, exact 𝔽₂/integer)

| order | dim | weight distribution | identification |
|-------|-----|---------------------|----------------|
| 8     | 4   | 1, 14·x⁴, x⁸        | ê₈ = RM(1,3) (unique Type II of length 8) |
| 16    | 8   | 1, 28·x⁴, 198·x⁸, 28·x¹², x¹⁶ | d₁₆⁺ (Type II, indecomposable; support graph connected) |
| 32    | 16  | 1, 120·x⁴, 1820·x⁸, 8008·x¹², 45638·x¹⁶, … | d₃₂⁺ enumerator (explicit d₃₂⁺ constructed and matched); doubly-even, self-orthogonality checked |

Cross-checks: Sylvester row spans at 8/16/32 = dims 3/4/5 (+𝟙: 4/5/6) = RM(1,m) ✓;
order-16 identification is rigorous (only two Type II codes exist at length 16, split by
indecomposability); order-32 identification is at weight-enumerator + Type II level (full
equivalence to d₃₂⁺ not certified — among the five Type II [32,16] codes the enumerator +
A₄-structure singles out d₃₂⁺; flagged in HYP-2409).

## Remarks

1. **The RM dictionary completed.** Sylvester doubling generates the biorthogonal END of
   the RM chain (RM(1,m): repetition's neighbor); the skew/tournament doubling generates
   the self-dual MIDDLE (k = n/2). The diagonal +1 (the observer) and the twin corrections
   (top-degree point masses) are what skewness adds — the code jumps from the bottom of
   the RM filtration to its self-dual fixed point. At order 8 the two worlds still agree
   (ê₈ = RM(1,3) IS self-dual RM); at 16 RM has no self-dual member and the tower must
   leave — to d₁₆⁺ — which is exactly where the Hadamard equivalence splits (THM-451).
2. **Two self-dual-defect codes, one dispatch.** THM-477's blue code lives in TILING space
   (𝔽₂^m, σ = grid involution, glue 𝔽₂^f on the hypotenuse); this code lives in AMBIENT
   space (𝔽₂^n, rows of I+S). Both answer "k = n/2 self-dual" with a self-dual(-defect)
   structure whose defect sits on the antipodal/twin locus (hypotenuse ↔ twin anti-diagonal,
   THM-452's fractal hypotenuse). Whether the blue-code glue 𝔽₂^f maps to the d⁺ glue
   under THM-474's gauge is open (HYP-2409).
3. **Gleason frame.** Type II enumerators live in ℂ[W_{ê₈}, W_{g₂₄}] (Gleason); the tower
   ladder starts AT the first generator W_{ê₈}. Whether the tower ever meets the Golay
   generator (e.g. via a Paley-24 branch — Assmus–Key: exactly two order-24 Hadamard
   classes give the extended Golay code, one of them Paley) is a natural follow-up.
