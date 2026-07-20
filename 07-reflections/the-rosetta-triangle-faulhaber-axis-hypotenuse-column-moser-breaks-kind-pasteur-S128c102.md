# The Rosetta triangle: a Faulhaber third axis, the hypotenuse column, and Moser breaks all the way down

**kind-pasteur-2026-07-20-S128c102** (HYP-8165) · owner: the triangle
{1},{2,1},{3,3,1},{4,6,5,1},{5,10,14,9,1},{6,15,30,37,17,1},{7,21,55,101,99,33,1},
"a third perspective on the triangular numbers... fibonacci... powers of 2...
mosers circle problem... 2n+1 and 2^x+1... n·2^x+1... 1,3/2,11/6,25/12,137/60...
and 1,5/2,29/3,109/12,1079/60... mine past work extensively."

## 1. Identification status (everything checked exactly; OEIS-new)

The triangle, its row sums (1,3,7,16,39,106,317), its NE-diagonal sums
(1,2,4,7,12,21,37,68), and the second series' numerators are ALL ABSENT from
OEIS (curl-verified). No constant/k/n/m-coefficient linear recurrence of support
≤ 3 exists. The structure that DOES hold, exactly:

- **Columns are Faulhaber power sums with rare +1 deviations.** With m = n−k+1
  (the term count), T(n,k) = Σ_{j=1}^{m} j^{k−1} for ALL entries except exactly
  three: T(6,4) = 36+1, T(7,4) = 100+1, T(7,5) = 98+1. So col 1 = n, col 2 =
  triangular, col 3 = square pyramidal — and the **penultimate column is the
  two-term power sum 1 + 2^{n−2}: the repo's HYPOTENUSE LAW / transitive-SC-
  neighbor H**, whose prime entries 2, 3, 5, 17 (→257) are THM-871's Fermat
  rungs. The first deviation appears in ROW 6 — the same row where Moser's
  circle count first breaks from 2^{n−1}.
- **The THIRD PERSPECTIVE named:** klein's figurate atlas (T1532) has the
  gonality axis (triangular→square→pentagonal) and the dimension axis
  (triangular→tetrahedral→simplex). This triangle threads the triangular
  numbers on the **EXPONENT (Faulhaber) axis**: Σj⁰ = n, Σj¹ = triangular,
  Σj², Σj³, … — the third classical family through T_n, here with the owner's
  deviation twist.

## 2. The diagonal law and its Moser break (the Fibonacci + 2^m answer)

NE-diagonal sums D_m = 1, 2, 4, 7, 12, 21, 37, 68, 130, 255…:

    D_m = F_{m+1} + r(m),   r = 0, 0, 1, 2, 4, 8, 16, 34, 75, 166, …

Pascal's diagonals give pure Fibonacci; this triangle gives **Fibonacci plus a
residual that is powers of 2 for m ≤ 7 and then BREAKS (34, 75, 166) — the
Moser phenomenon reproduced in the diagonals** (the law D = F_{m+1} + 2^{m−3}
holds m ≤ 7, fails at m = 8 by +2). Both owner hints (Fibonacci AND powers of
2) land in one identity, and the break is the Moser hint. The residual's eight
consecutive terms 1,2,4,8,16,34,75,166 match **A045648** (chiral n-ominoes in
(n−1)-space, one cell labeled) — computed from deviation-safe entries only
(D₁₀ needs no uncertain entry; r(10) = 166 ✓). If that law CONTINUES, it forces
**T(8,4) = 218 = 225 − 7** (a Moser-deficit −7!); the competing completions
give 225 (pure), 226 (step-dev), 227 (Fibonacci-dev). **One number from the
owner — row 8, column 4 — decides the triangle's universe.**

## 3. The two series

- Series 1 = the harmonic numbers H_n exactly (1, 3/2, 11/6, 25/12, 137/60).
- Series 2: the as-given reading (29/3) makes term₄ NEGATIVE (−7/12) — not a
  partial-sum series — so the 29/6 reading is forced. Then: **Σ_k T(n,k)/k
  matches EXACTLY for n = 1, 2, 3 (1, 5/2, 29/6)** and misses by 1/6 (n=4) and
  13/15 (n=5); no single integer modification of one triangle entry can
  reconcile either miss (proved: k·miss is never an integer for k ≤ n). Also
  ΣC(n,k)/k = Σ(2^k−1)/k matches 3/5 with different misses (6/12, 192/60). So
  series 2 = a harmonic-weighted transform of something STRICTLY RICHER than
  the displayed triangle — the second discriminating question for the owner.
- The Proth frame: series 1 = the x = 0 harmonic edge (Σ1/k), and the
  candidates for series 2 are harmonic transforms of the n·2^x+1 table's
  columns (Σ(2^k−1)/k = ΣC(n,k)/k is the x-summed edge) — "2n+1 and 2^x+1"
  interpolated by n·2^x+1 = the Proth numbers, whose two boundary families are
  the repo's H-spectrum odds (Ham-path counts = {odd}∖{7,21}) and the
  hypotenuse tower 1+2^d.

## 4. The repo crosswalk (the "mine past work" deliverable)

- Penultimate column = 1 + 2^{n−2} = H of the transitive tournament's big SC
  neighbor (CLAUDE.md principal line); its primes = THM-871's Fermat rungs.
- The two FORBIDDEN H-spectrum values {7, 21} both appear canonically: 7 =
  row-3 sum = 2³−1 and 21 = T(7,2) = C(7,2) — the triangle contains the
  H-spectrum's two holes at distinguished positions.
- klein T1532's atlas gains the exponent axis; T1533's "Moser deficit exactly 1
  at the 32nd region" rhymes with the +1 deviations; A045648's chirality onset
  = a NEW Moser-family sequence for the atlas (2-powers until dimension allows
  chirality — the break MECHANISM is named there: chiral pairs appear).
- Row sums 1,3,7,16,39,106,317 (OEIS-new, growth → ×3) and alternating sums
  1,1,1,2,1,0,7 are logged for future identification.

## 5. Honest ledger

- The triangle's generating rule beyond "power sums + three +1's" is NOT
  determined by the given data; three completions + the A045648-diagonal
  completion are registered with their row-8 discriminants (225/226/227/218).
- The A045648 match is 8 consecutive terms including irregular 34, 75, 166 —
  strong but correlational; A248890 shares terms through 166 and diverges at
  the next (370 vs 374): the residual r(11) = D₁₁ − 144 = 296 + T(8,4) − 144
  discriminates everything at once.
- Series-2 identification is OPEN (3/5 exact two ways, misses characterized,
  single-entry fixes excluded).
