/-
COLLATZ EXPENSE RECORDS CERTIFICATE (opus-2026-10-06-S16).

The S9 note (05-knowledge/results/collatz_expense_diophantine_20261005.md) defines, for a first-descent segment
type (l, A) with 2^A > 3^l, the least admissible rate q(l, A) = ceil(l / (A - l log2 3)).  The worst type at
length l has A = A_min(l), the least A with 2^A > 3^l, and q_max(l) = q(l, A_min(l)).  S9 listed the records of
q_max for l <= 320 (exact power comparisons for the record lengths; q beyond l = 94 from 50-digit floats).

This file re-derives every q_max(l), l = 1..320, from ONE isolation of theta = log2 3, in the style of the
eleven-square packing formalization (github.com/Queuingtheorydotcom/11SquaresFormalized), which isolates its
degree-8 algebraic endpoint u in (9/25, 37/100) and then argues by exact arithmetic:

  (B) bracket   2^PLO < 3^QLO  and  3^QHI < 2^PHI,  i.e.  PLO/QLO < theta < PHI/QHI,
                with PLO/QLO = 16785921/10590737 and PHI/QHI = 301994/190537 consecutive convergents
                of theta (width 4.96e-13);
  (A) for each l:  2^(A-1) < 3^l < 2^A   (A = A_min(l));
  (Q) f(t) = l / (A - l t) is increasing on t < A/l, so f(PLO/QLO) < q-value < f(PHI/QHI); with
      dlo = A QLO - l PLO and dhi = A QHI - l PHI (both > 0),
      (q - 1) dlo <= l QLO  and  l QHI <= q dhi   give  q - 1 < l/(A - l theta) < q,  so q = q_max(l).

The glue "(B) and (Q) imply q = ceil(l/(A - l theta))" is the two-line monotonicity argument above (real
analysis, not formalized here).  Everything decidable is checked: Mathlib-free, Nat-only, native_decide (the
same trust model as LrcExtremalCert.lean and the packing formalization: Lean kernel + native compiler).
Mirror: 04-computation/experiments/collatz_expense_records_exact_20261006.py.
-/

def PLO : Nat := 16785921
def QLO : Nat := 10590737
def PHI : Nat := 301994
def QHI : Nat := 190537

/-- the root-isolation step: PLO/QLO < log2 3. -/
theorem bracket_lo : 2 ^ PLO < 3 ^ QLO := by native_decide
/-- the root-isolation step: log2 3 < PHI/QHI. -/
theorem bracket_hi : 3 ^ QHI < 2 ^ PHI := by native_decide

/-- least A with 2^A > 3^l (for l >= 1, 3^l is never a power of two). -/
def aMin (l : Nat) : Nat := Nat.log2 (3 ^ l) + 1

/-- ceil of the lower bracket value f(PLO/QLO) = l QLO / dlo. -/
def qCert (l : Nat) : Nat :=
  let A := aMin l
  let dlo := A * QLO - l * PLO
  (l * QLO + dlo - 1) / dlo

def certOK (l : Nat) : Bool :=
  let A := aMin l
  let q := qCert l
  let dlo := A * QLO - l * PLO
  decide (2 ^ (A - 1) < 3 ^ l) && decide (3 ^ l < 2 ^ A) &&
  decide (l * PLO < A * QLO) && decide (l * PHI < A * QHI) &&
  decide ((q - 1) * dlo ≤ l * QLO) && decide (l * QHI ≤ q * (A * QHI - l * PHI))

/-- every length 1..320 is certified by the single bracket. -/
theorem all_lengths_certified : ((List.range 320).map (· + 1)).all certOK = true := by native_decide

/-- running records of q_max over lengths 1..L, as (l, A_min(l), q_max(l)). -/
def recordsUpTo (L : Nat) : List (Nat × Nat × Nat) :=
  let step := fun (acc : Nat × List (Nat × Nat × Nat)) (l : Nat) =>
    let q := qCert l
    if q > acc.1 then (q, acc.2 ++ [(l, aMin l, q)]) else acc
  (((List.range L).map (· + 1)).foldl step (0, [])).2

/-- the S9 record table: lengths are the denominators of the upper best approximations of log2 3. -/
theorem records_320 : recordsUpTo 320 =
    [(1, 2, 3), (3, 5, 13), (5, 8, 67), (17, 27, 306), (29, 46, 804), (41, 65, 2480), (94, 149, 6951),
     (147, 233, 13984), (200, 317, 26668), (253, 401, 56382), (306, 485, 207489)] := by native_decide

/-- theta-free cross-check for the first four records: q = least m with 3^(m l) < 2^(m A - l)
    (monotone in m because 3^l < 2^A). -/
def directOK (l A q : Nat) : Bool :=
  decide (3 ^ (q * l) < 2 ^ (q * A - l)) && !decide (3 ^ ((q - 1) * l) < 2 ^ ((q - 1) * A - l))

theorem direct_small_records :
    [(1, 2, 3), (3, 5, 13), (5, 8, 67), (17, 27, 306), (29, 46, 804), (41, 65, 2480)].all
      (fun r => directOK r.1 r.2.1 r.2.2) = true := by native_decide
