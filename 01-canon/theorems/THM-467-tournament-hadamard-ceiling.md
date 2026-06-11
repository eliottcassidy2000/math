# THM-467 — The Tournament Hadamard Ceiling: det(I+S) maxima are DRT switching classes

**Status:** PROVED (mac-mini-2026-06-10-S2) — adversarial verification PENDING (this session).
**Provenance:** mac-mini-2026-06-10-S2. Companions: THM-466 (floor, setup and notation),
THM-468 (average). Related: THM-447/448 (skew tower), THM-213 (Paley Pfaffian),
the-drt-engine-is-S-squared-equals-J-minus-nI reflection, HYP-2385.

## Statement

Let T be a tournament on n vertices, S = A - A^T. Then:

1. **(odd n)** det(I+S) ≤ (n+1)^((n-1)/2), with equality iff T is a switching of a
   doubly regular tournament (DRT). In particular the bound is attained iff a DRT on
   n vertices exists, iff a skew-Hadamard matrix of order n+1 exists (Reid–Brown);
   attainment requires n ≡ 3 (mod 4) (or n = 1).
2. **(even n)** det(I+S) ≤ n^(n/2) (the Hadamard bound, since I+S has ±1 entries),
   with equality iff S S^T = (n-1) I (S is a **skew conference matrix**), which
   requires n ≡ 0 (mod 4) or n = 2, and happens iff I+S is a skew-Hadamard matrix
   of order n.
3. The maximizing set of d(T) = det(I+S)/2^(n-1) is, at orders where the relevant
   bound is attained, a UNION OF FULL SWITCHING CLASSES of DRTs (odd n) /
   skew-conference cores (even n).

## Proof

Eigenvalues of a real skew matrix S come in purely imaginary pairs ±iμ_j (μ_j > 0),
plus a kernel; say p pairs, so det(I+S) = Π_{j=1}^{p} (1 + μ_j²) ≥ 0 and
Σ_j 2μ_j² = tr(S S^T) = Σ_{i≠j} S_ij² = n(n-1), i.e. **Σ_{j=1}^p μ_j² = n(n-1)/2**.

By AM-GM, Π (1+μ_j²) ≤ (1 + n(n-1)/(2p))^p. The function p ↦ p·log(1 + c/p) is
strictly increasing in p (its derivative log(1+c/p) - (c/p)/(1+c/p) > 0 for c > 0),
so the bound is maximized at the largest possible p:

- **odd n**: p ≤ (n-1)/2 (skew odd-dimensional always has a kernel), giving
  det(I+S) ≤ (1 + n)^((n-1)/2).
- **even n**: p ≤ n/2, giving det(I+S) ≤ (1 + (n-1))^(n/2) = n^(n/2).

**Equality, odd n:** requires p = (n-1)/2 (kernel exactly 1-dimensional) and all
μ_j² = n. Then S S^T (symmetric PSD, integer) has eigenvalues n with multiplicity
n-1 and 0 with multiplicity 1: S S^T = nI - n vv^T with v a unit kernel vector.
Diagonal entries: (S S^T)_ii = n - 1 (each row of S has n-1 entries ±1), so
n - n v_i² = n - 1, v_i² = 1/n for every i: v = w/√n with w ∈ {±1}^n and
**S S^T = nI - w w^T**. Let D = diag(w) and S' = D S D (a switching of T). Then
S' S'^T = D(S S^T)D = nI - (Dw)(Dw)^T = nI - 𝟙𝟙^T = nI - J, and S'𝟙 = D S w = 0
(w spans ker S = ker S S^T). S' S'^T = nI - J together with S' skew is exactly
S'² = -(nI - J) = J - nI: the DRT defining identity (canon: the DRT engine
S² = J - nI; equivalently every vertex has score (n-1)/2 and every dominated-pair
count is (n-3)/4 — forcing n ≡ 3 mod 4). Conversely a DRT has spectrum
{0} ∪ {±i√n each (n-1)/2 times}, attaining equality, and by THM-466(2) so does its
whole switching class. Bordering: I+S' extends to a skew-Hadamard matrix of order
n+1 by adding a ±1 border row/column (Reid–Brown correspondence, canon THM-447
dictionary), so attainment at n ⟺ skew-Hadamard at n+1.

**Equality, even n:** requires p = n/2 (S nonsingular) and all μ_j² = n-1, i.e.
S S^T = (n-1) I. Then (I+S)(I+S)^T = I + S S^T = nI: I+S is a Hadamard matrix of
skew type, forcing n ≡ 0 (mod 4) for n > 2 (classical Hadamard order constraint);
n = 2 attains with S = [[0,1],[-1,0]]. ∎

## Computational verification (this session, n = 3..8 exhaustive per iso class)

| n | bound | max det(I+S) | attained | maximizing classes |
|---|-------|-------------|----------|--------------------|
| 3 | 4     | 4           | YES      | 2 classes = whole switching class of C3 (DRT) |
| 4 | 16    | 16          | YES      | 2 classes, both skew-conference |
| 5 | 36    | 32          | no (no DRT at 5 ≡ 1 mod 4) | d_max = 2 (8 classes) |
| 6 | 216   | 160         | no (no skew conference at 6 ≡ 2 mod 4) | d_max = 5 (4 classes, |Pf| = 9) |
| 7 | 512   | 512         | YES      | 6 classes, all switchings of QR_7 |
| 8 | 4096  | 4096        | YES      | 4 classes, all skew-conference |

- QR_7 Pfaffian-minor census: Pf² distribution {1:1}, {1:21}, {1:21, 9:14}, {49:7}
  over |K| = 0,2,4,6 — sums to 512; QR_7 contains exactly 14 vortex 4-sets;
  |Pf(S minus any vertex)| = 7 = 7^((7-3)/4) (THM-213 cross-check passes).
- The non-attained orders n ≡ 1, 2 (mod 4) define the **tournament Barba problem**
  (exact maxima 32, 160 at n = 5, 6; n = 9, 10 values pending this session):
  what is max det(I+S) when the spectral-flatness obstruction blocks the ceiling?
- Script: 04-computation/hadamard_det_census_macmini_s2.py
  Output: 05-knowledge/results/hadamard_det_census_macmini_s2.out

## Notes

- This is the skew-Hadamard ⟺ DRT correspondence in *variational* form: DRT switching
  classes are not just *examples* of extremal matrices, they are *exactly* the
  maximizers of the tournament determinant. The simplex translation (rows of I+S as
  cube vertices): tournament simplices have volume |det(I+S)|/n!, and DRT switching
  classes are the maximal-volume tournament simplices — the tournament case of
  Hadamard's maximal determinant problem.
- H(T) interplay (refining HYP-2386): at n = 7 the global H-maximizer (H = 189,
  Paley) LIES IN the determinant-ceiling switching class, whose six classes carry
  H ∈ {45, 69, 87, 135, 189} — switching pins d and scrambles H. At n = 8 the
  H-maximizers (H = 661) have d ∈ {2,6,8,9}, far below d_max = 32: the two crowns
  decouple. Bulk correlation d vs H ≈ 0 at n = 7, 8.
