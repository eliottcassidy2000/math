# THM-450 — The Doubling Family Trichotomy: projector, projector, Hadamard mixer (a Clifford/Cayley-Dickson classification)

**Status:** PROVED (block algebra, with the forced-shape derivation) + VERIFIED exhaustively
(all 8+64 labeled tournaments n=3,4 × all 16 members; symmetry identities exhaustive).
Adversarially re-verified (`verify_E_family_kps1.out`).
**Source:** kind-pasteur-2026-06-09-S1 (branch E + verifier). Resolves HYP-2337.
**Related:** THM-447, OPEN-Q-044/045, THM-250..253 (rapidity), `lrc-cayley-dickson-tower-s387.md`.

## (1) The forced shape

Any 2×2-block doubling with blocks in span{M, I} that yields a valid tournament dominance
matrix for ALL tournaments T must have the form
```
M' = [[aM, bM + βI], [bM − βI, dM]],   a,b,d,β ∈ {±1}    (16 members)
```
(diagonal blocks: zero diagonal forces no I-term; skew-symmetry forces equal off-diagonal
M-signs and opposite I-signs; β ≠ 0 because twins need an arc). CAVEAT: blocks outside
span{M, I} (e.g. J) give further doublings such as the ordered sum T⊕T — the classification
is within the natural ansatz.

## (2) The Kronecker form and the square law

Write M' = P⊗M + Q⊗I with P = [[a,b],[b,d]], Q = [[0,β],[−β,0]]. Then
```
M'² = P²⊗M² + {P,Q}⊗M + Q²⊗I,   Q² = −I₂,   {P,Q} = β·(tr P)·J
```
**Chebyshev members:** M'² is block-diagonal for all T ⟺ tr P = 0 ⟺ {P,Q} = 0, and then
M'² = I₂⊗(2M² − I) — the Chebyshev T₂ law of THM-447.
**Skew-Hadamard preservation ⟺ tr P = 0 as well** (proof: S'S'^T = I − M'²; M² = (1−n)I gives
M'² = (1−2n)I iff tr P = 0; else the off-diagonal block ±2(b(1−n)I + βM) ≠ 0). The preserving
members are exactly those whose sign matrix P is, up to symmetry, the 2×2 Hadamard
[[1,1],[1,−1]] (det P = −2). **The 2×2 Hadamard matrix is the doubling's DNA.**

## (3) The three orbits

Under the symmetries (copy swap, global arc reversal, pre-composition with op), the 16 members
form exactly 3 orbits, with exact geometric meaning in the swap eigenbasis (e = (1,1),
f = (1,−1)):
```
ORBIT A (size 4):  T[K₂]    P = ee^T  = 2·projector on the swap-SYMMETRIC line   (det 0, rank 1)
ORBIT B (size 8):  D_skew   P = H₂    = the MIXER exchanging the two lines       (det −2)
ORBIT C (size 4):  SCblow   P = ff^T  = 2·projector on the swap-ANTISYMMETRIC line (det 0, rank 1)
```
P_A·P_C = 0, P_A + P_C = 2I. The family IS {projector, projector, mixer}.

## (4) Spectral pencils

spec(M') = ⋃_λ spec(λP + Q); char-poly product law det(xI − M') = ∏_λ [x² − (tr P)λx +
(det P·λ² + 1)] (verified all 16 members × 7 samples).
- **Orbit B (mixer):** lifts ±√(2λ²−1) per eigenvalue — Chebyshev (THM-447).
- **Orbits A, C (rank-1):** the two lifts of each eigenvalue multiply to exactly 1
  (det(λP + Q) = 1): a UNIMODULAR/boost pair — the rapidity-lattice theme of THM-251/252.
  Matrix relation: M'² − (tr P)(I₂⊗M)M' + I = 0 (exhaustive n=4).

## (5) Cayley-Dickson reading

(P/√2, Q) with {P,Q} = 0, (P/√2)² = +1, Q² = −1 is a real representation of Cl(1,1) ≅ M₂(R).
The CD twist jz = z̄j for purely-imaginary z (skew M) is LITERALLY {P,Q} = 0. D_skew = the CD
step with P = H₂ (Sylvester convention), Q = +J (counterclockwise j); β = −1 is the conjugate
orientation −j; a↔d is the sheet swap; negation is the opposite algebra — the 8 orbit-B members
are 2×2×2 copies of ONE Cayley-Dickson step. The doubling tower (THM-448) is thus a CD tower
on tournaments; cf. the repo's CD tower at n = 2^k + 1 — this one lives at 2^k − 1.

## (6) H-extremality (n=5, all 12 classes; doubles at n=10)

SCblow argmax in ALL 12 classes (T[K₂] ties only at regular T5: 15565); means
SCblow 15179.7 > D_skew 6906.3 > T[K₂] 3953.0. **The unique skew-Hadamard-preserving orbit is
NOT the H-maximizing one.** T[K₂]/SCblow orbits are fully op-equivariant (PROVABLE: inside
those orbits g3 = g1∘g2); the D_skew orbit is the ONLY non-op-equivariant doubling — its H
splits exactly on op-paired classes (189/333 at n=4; 1809/3561, 2817/5289 at n=5).

## (7) The consecutive-circulant coincidence (NEW LEAD → HYP-2344, renumbered from 2338)

All three orbits send C₃ to H = 45 = the n=6 maximum (scores (2,2,2,3,3,3)), but to TWO iso
classes: T[K₂](C₃) ≅ SCblow(C₃), while D(C₃) is distinct. More: T[K₂](U_n) ≅ SCblow(U_n) for
the consecutive circulants U₃ = C₃, U₅ = C₅({1,2}) = regular T₅, U₇ = C₇({1,2,3}) (iso verified
at orders 6, 10, 14; U₇: H = 24,855,901 both). But Paley T₇ SPLITS them: H(K₂) = 24,589,929 ≠
H(SCblow) = 24,453,597 — and so does the third regular-7 class (24,392,325 ≠ 24,485,541). The
criterion is consecutive-circulant structure, NOT regularity/DRT/SC (all refuted; evidence
caveat: n=7 is the only order with a discriminating competitor so far).

## Scripts

`block_doubling_classification_kps1.py`, `regular7_doubling_Htest_kps1.py`,
`consecutive_circulant_iso_kps1.py`, `verify_E_family_kps1.py` (+ .out).
