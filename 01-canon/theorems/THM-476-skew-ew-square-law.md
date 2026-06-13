# THM-476 — The skew Ehlich–Wojtas square law: tournament EW attainment forces 2n−3 = k²

**Status:** PROVED (necessity; claudebox-2026-06-11-S1) + VERIFIED witnesses (sufficiency
data) at n = 6, 14, 26, 62 — every candidate order ≤ 62 attained; first open n = 86.
**Provenance:** claudebox-2026-06-11-S1. **Companions:** THM-472 (even ceiling n^(n/2),
needs n ≡ 0 mod 4), THM-475 (odd sibling), HYP-2389 (Barba), HYP-2405 (non-square branch).
**Literature:** the skew D-optimal/EW theory at n ≡ 2 (mod 4) is studied (Greaves–Suda's
char-poly characterization via tournaments on n−1 vertices; Armario–Frau; known orders
6,14,26,42,62, FIRST OPEN n = 86) — those orders are exactly 2n−3 = 9,25,49,81,121,169 = k²,
so this law is consistent with (and may be implicit in) that literature; the proof below is
elementary and self-contained. Cite: Greaves–Suda (skew E-W matrices), Armario–Frau (self-dual
codes from skew E-W), Ehlich 1964 (the bound and its Gram rigidity).

## Statement

Let n ≡ 2 (mod 4), and let T be a tournament on n vertices whose ±1 matrix M = I + S attains
the Ehlich–Wojtas maximal-determinant bound det M = 2(n−1)(n−2)^((n−2)/2) (the maximum over
ALL ±1 matrices of order n ≡ 2 mod 4). Then **2n−3 is a perfect square**.

Equivalently: the skew-attainable EW orders lie in { n = (k²+3)/2 : k odd } = 6, 14, 26, 42,
62, 86, 114, … For all other n ≡ 2 (mod 4) — e.g. n = 10, 18, 22, 30 — tournaments are
strictly det-deficient relative to general ±1 matrices.

## Proof

By Ehlich's rigidity, attainment forces the Gram G = MMᵀ to be, up to simultaneous signed
permutation of rows, block-diagonal with two blocks (n−2)I + 2J of size n/2. So G has
eigenvalue 2n−2 with multiplicity exactly 2 (one per block; the other eigenvalue n−2 < 2n−2),
and the 2n−2 eigenspace V is spanned by two orthogonal ±1-SIGNED INDICATOR vectors a, b
(supports = the two blocks, after un-permuting; ‖a‖² = ‖b‖² = n/2; entries in {0,±1}).

M is skew-type: MMᵀ = (I+S)(I−S) = I − S² = MᵀM, and S commutes with G = I − S², so S
preserves V. Skewness gives aᵀSa = 0, hence Sa = μb for some real μ; likewise Sb = νa.
On V, S² = I − G acts as 1 − (2n−2) = −(2n−3), so μν = −(2n−3); skewness (bᵀSa = −aᵀSb)
with ‖a‖ = ‖b‖ gives ν = −μ; hence **μ² = 2n−3**. But Sa is an integer vector and b has
entries in {0,±1} with full support on its block, so μ ∈ ℤ. ∎

(The same forced-integer-eigenvector technique fails at n ≡ 1 (mod 4): the THM-475 maximizer's
excited eigenspace is also 2-dimensional, but its rational basis is not pinned to indicator
vectors, so no square condition arises — only the kernel/parity obstruction of THM-472.)

## Witnesses (sufficiency data; exact Bareiss verification; script skew_ew_square_law_cbx1.py)

- **n = 6** (2n−3 = 3²): max det(I+S) = 160 = EW(6), EXHAUSTIVE (matches the mac-mini census
  d_max = 5, 4 classes); skew attains the global ±1 maxdet at order 6.
- **n = 14** (2n−3 = 5²): annealing witness, det = 77 635 584 = 2·13·12⁶ = EW(14) exactly;
  Gram off-diag multiset {±2: 42 within-block pairs, 0: 49 cross pairs} = the EW two-block
  form. Stored in the results file.
- **n = 26** (2n−3 = 7²): explicit **two-circulant** witness, blocks of size 13:
  M = [[A, B], [−Bᵀ, Aᵀ]] with A = I + S_a, S_a the circulant skew row
  a = (0,+,−,+,−,−,−,+,+,+,−,+,−) and B the circulant with row
  b = (+,+,−,−,+,−,−,−,−,−,−,−,−); det M = 1 826 017 371 802 828 800 = 2·25·24¹² = EW(26)
  exactly. The skew structure forces the classical two-squares datum 2n−2 = a₁² + b₁²
  (row sums of the two circulants) to DEGENERATE to a₁ = 1: skew-EW ⟺ the representation
  2n−2 = 1² + k² — the apex/observer representation.
- **n = 62** (2n−3 = 11²): explicit two-circulant witness found by the EXACT integer
  PAF search (multiplier subgroup of order 3 on both blocks of ℤ/31):
  a' = +-+-+-++--+-----+++++-++--+-+-+ , b = +++--+--+++---+----+-----+----- ;
  det = 26 971 018 205 929 469 663 772 672·10³⁰ = 2·61·60³⁰ = EW(62) exactly (Bareiss).
- **n = 86 NEGATIVE (the open frontier resists the symmetric ansatz):** exhaustive
  multiplier-symmetric two-circulant search on ℤ/43 (subgroup orders {21, 7, 3} on both
  blocks independently, exact integer PAF conditions) finds NO skew-EW matrix of order 86.
  Escalation routes: asymmetric b (2⁴³, needs meet-in-the-middle/SAT on the PAF inverse
  problem), unequal block symbols (per-frequency b(ω)=0 branches), Goethals–Seidel 4-block
  arrays, dihedral/non-cyclic groups, or the Greaves–Suda char-poly route directly on
  85-vertex tournaments. Emitted as a coordination task.
- **n = 10, 18** (17, 33 not squares): EW impossible for tournaments (this theorem); best
  found by extensive annealing: B_t(10) ≥ 64000 (char poly (x²−18x+61)⁴(x−9)², eigenvalues
  9 ± 2√5 ∈ ℚ(√5) — see HYP-2405), B_t(18) ≥ 131 533 373 440 = 2²⁹·5·7² (0.90·EW).

## Remarks

1. **The square law is a totient-style size-match.** Skewness pins one unit of the two-squares
   mass: general EW attainment needs 2n−2 = a² + b² (two odd squares); the tournament case
   forces a = 1 (the diagonal/observer/apex contributes exactly the identity), leaving
   2n−3 = b². One more instance of "the +1 is the observer" (HYP-2275/2310 ledger).
2. **Frontier.** First open skew case n = 86 (k = 13; tournament on 85 vertices with char poly
   (x³−41x²−21·83)(x²+x+21)^41 per Greaves–Suda) — would also yield self-dual [172,86] codes
   (Armario–Frau). Multiplier-symmetric two-circulant searches this session: see results file
   + session reflection.
3. **Tao C23a/C23b.** The even-order skew ladder (n ≡ 0: skew-Hadamard, κ=1; n ≡ 2: EW iff
   2n−3 = k², κ → √2) is the tournament shadow of the Hadamard-conjecture frontier: all four
   open Hadamard orders < 2000 (668, 716, 892, 1132) are skew-eligible (could fall via a DRT),
   while the n ≡ 2 ladder shows exactly how much determinant skewness costs when perfect
   flatness is arithmetically barred.
