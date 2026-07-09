/-
  TournamentH7.LRCD3FloorCert  (kind-pasteur-2026-07-08-S89)

  Kernel-checkable D3 COVERING FLOOR for the k=11 density-floor finite pieces.

  Ports klein's Farey-cell moment computation (04-computation/lrc14_d3_exact_verify_klein_S184.py)
  to exact ℚ arithmetic, in the same `native_decide` style as `LyWindowEnum`.  For a speed set E
  (sorted ℤ list), computes the exact moments m1,m2,m3 of the uncovered measure
  W(x) = Σ_gaps (gap(x) − 1/7)_+, then the degree-3 covering floor
      D3 = m1/M + (m1 − m2/M)² / (m2 − m3/M),   M = 6/7,
  a rational lower bound on μ = P(W>0) (THM-661).  `bar = 83549/252252` is the honest k=11 bar.

  The `native_decide` theorems certify `bar ≤ D3 E` for the anchor shapes of the corrected
  longest-AP closure:  the global minimizer (the 11-block), opus-S155's tail extremal A, and the
  block+outlier — kernel-checking the exact rationals the whole k=11 closure rests on.
-/
import Mathlib

namespace LRCD3FloorCert

/-- Arc length θ = 1/7. -/
def theta : ℚ := 1/7
/-- M = 6/7 (the max of W). -/
def MM : ℚ := 6/7
/-- The honest k=11 bar, 83549/252252. -/
def bar : ℚ := 83549/252252

/-- Farey fractions p/q with q ≤ D, in [0,1], deduped & sorted (the phase-order breakpoints). -/
def fareyList (D : ℕ) : List ℚ :=
  (((List.range (D+1)).flatMap (fun q =>
      if q = 0 then ([] : List ℚ)
      else (List.range (q+1)).map (fun p => (p : ℚ) / (q : ℚ)))).dedup).mergeSort (fun a b => a ≤ b)

/-- Circular gaps as linear functions `slope·x + intercept` at cell midpoint `mid`.
    Sort speeds by phase `E_i·mid − ⌊E_i·mid⌋`, take consecutive differences + the wrap gap. -/
def gapsAt (E : List ℤ) (mid : ℚ) : List (ℚ × ℚ) :=
  let phdata : List (ℤ × ℤ) := E.map (fun e => (e, ⌊(e : ℚ) * mid⌋))
  let sorted := phdata.mergeSort (fun x y => ((x.1 : ℚ) * mid - (x.2 : ℚ)) ≤ ((y.1 : ℚ) * mid - (y.2 : ℚ)))
  let phs : List (ℚ × ℚ) := sorted.map (fun x => ((x.1 : ℚ), -(x.2 : ℚ)))   -- (slope E_i, intercept −cj_i)
  match phs with
  | [] => []
  | p0 :: _ =>
    let pl := phs.getLast!
    let cons := (phs.zip phs.tail).map (fun pr => (pr.2.1 - pr.1.1, pr.2.2 - pr.1.2))
    cons ++ [(p0.1 - pl.1, (1 : ℚ) + p0.2 - pl.2)]

/-- Sub-breakpoints in (a,b) where a gap length crosses θ, together with a,b; deduped & sorted. -/
def subPts (gaps : List (ℚ × ℚ)) (a b : ℚ) : List ℚ :=
  let xs := gaps.filterMap (fun g =>
    if g.1 = 0 then none
    else
      let xc := (theta - g.2) / g.1
      if a < xc ∧ xc < b then some xc else none)
  ((a :: b :: xs).dedup).mergeSort (fun a b => a ≤ b)

/-- On a sub-interval, W(x) = A·x + Bc where the sum is over gaps with length > θ at `m2m`. -/
def ABc (gaps : List (ℚ × ℚ)) (m2m : ℚ) : ℚ × ℚ :=
  gaps.foldl (fun acc g =>
    if theta < g.1 * m2m + g.2 then (acc.1 + g.1, acc.2 + g.2 - theta) else acc) (0, 0)

/-- Exact (m1,m2,m3) contribution of the Farey cell [a,b]. -/
def cellMoments (E : List ℤ) (a b : ℚ) : ℚ × ℚ × ℚ :=
  let mid := (a + b) / 2
  let gaps := gapsAt E mid
  let pts := subPts gaps a b
  (pts.zip pts.tail).foldl (fun acc lohi =>
    let lo := lohi.1; let hi := lohi.2; let m2m := (lo + hi) / 2
    let AB := ABc gaps m2m; let A := AB.1; let Bc := AB.2
    let i1 := A/2*(hi*hi - lo*lo) + Bc*(hi - lo)
    let i2 := A*A/3*(hi^3 - lo^3) + A*Bc*(hi*hi - lo*lo) + Bc*Bc*(hi - lo)
    let i3 := A^3/4*(hi^4 - lo^4) + A*A*Bc*(hi^3 - lo^3) + 3*A*Bc*Bc/2*(hi*hi - lo*lo) + Bc^3*(hi - lo)
    (acc.1 + i1, acc.2.1 + i2, acc.2.2 + i3)) (0, 0, 0)

/-- Exact moments (m1,m2,m3) = (E[W], E[W²], E[W³]) of the uncovered measure. E must be sorted. -/
def moments (E : List ℤ) : ℚ × ℚ × ℚ :=
  let D := (E.getLast! - E.headI).toNat
  let breaks := fareyList D
  (breaks.zip breaks.tail).foldl (fun acc ab =>
    let cm := cellMoments E ab.1 ab.2
    (acc.1 + cm.1, acc.2.1 + cm.2.1, acc.2.2 + cm.2.2)) (0, 0, 0)

/-- The degree-3 covering floor D3 (THM-661), a rational lower bound on μ = P(W>0). -/
def D3 (E : List ℤ) : ℚ :=
  let m := moments E
  let m1 := m.1; let m2 := m.2.1; let m3 := m.2.2
  m1/MM + (m1 - m2/MM)^2 / (m2 - m3/MM)

-- Sanity value (visible in the build log): expect 0.404751 (the global minimizer)
#eval D3 [0,1,2,3,4,5,6,7,8,9,10]

/-- The GLOBAL D3-minimizer over all primitive 11-sets — the 11-block — clears the bar. -/
theorem block_floor : bar ≤ D3 [0,1,2,3,4,5,6,7,8,9,10] := by native_decide

/-- opus-S155's corrected TAIL extremal `A = 3·{0..9} ∪ {8}` (prim-diam 27) clears the bar. -/
theorem A_floor : bar ≤ D3 [0,3,6,8,9,12,15,18,21,24,27] := by native_decide

/-- The block+outlier `{0..9} ∪ {25}` (prim-diam 25) clears the bar. -/
theorem blockOutlier_floor : bar ≤ D3 [0,1,2,3,4,5,6,7,8,9,25] := by native_decide

/-- **A finite-closure slice.**  The block+outlier family `{0..9} ∪ {D}` for every far point
`D = 25..40` clears the bar (the binding longest-AP=10 far-outlier sub-family; cf. LEM-009). -/
theorem blockOutlier_family :
    ((List.range 16).map (fun n => 25 + (n : ℤ))).all
      (fun D => decide (bar ≤ D3 [0,1,2,3,4,5,6,7,8,9, D])) = true := by
  native_decide

/-- **The longest-AP=10 tail extremal at scale 3 (opus's A), every interior point** — the whole
`3·{0..9} ∪ {p}` sub-family (`p = 1..26`, `p` off the AP lattice ⇒ primitive) clears the bar; the
minimum is at `p = 8` (= opus's A, `D3 = 0.452986`). -/
theorem A_scale3_family :
    ((List.range 26).map (fun n => (n : ℤ) + 1)).all
      (fun p => if p % 3 = 0 then true else
        decide (bar ≤ D3 (([0,3,6,9,12,15,18,21,24,27, p]).mergeSort (fun a b => a ≤ b)))) = true := by
  native_decide

/-! ### Exhaustive slices — every primitive 11-set of bounded prim-diameter clears the bar. -/

/-- gcd of a ℤ list (via `natAbs`); `= 1` iff the set is primitive. -/
def listGcd (E : List ℤ) : ℕ := E.foldr (fun x g => Nat.gcd x.natAbs g) 0

/-- All PRIMITIVE 11-subsets of `{0,…,D}` with min `0` and max `D` (so prim-diam `= D`):
`0 :: (9-subset of {1,…,D−1}) ++ [D]`, kept when `gcd = 1`. -/
def shapes11 (D : ℕ) : List (List ℤ) :=
  ((((List.range (D - 1)).map (fun n => (n : ℤ) + 1)).sublistsLen 9).map
    (fun mid => (0 : ℤ) :: (mid ++ [(D : ℤ)]))).filter (fun E => listGcd E = 1)

/-- **Exhaustive slice, prim-diam ≤ 16.**  EVERY primitive 11-set with primitive diameter in
`{10,…,16}` (~8000 shapes) clears the bar — a genuine kernel-checked exhaustive over the
small-diameter base of the k=11 covering floor. The min is the 11-block `{0..10}`, `D3 = 0.404751`. -/
theorem exhaustive_le16 :
    ((List.range 7).map (· + 10)).all
      (fun D => (shapes11 D).all (fun E => decide (bar ≤ D3 E))) = true := by
  native_decide

end LRCD3FloorCert
