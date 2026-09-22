/-
Lane lean_paste_audit_w6 (2026-09-22). Core Lean 4.30.0, NO imports, NO Mathlib, NO sorry, NO axiom.
Scratch file: NOT a build target of CollatzBlueprintAudit. Compile with
  lake env lean 04-computation/lean/standalone/collatz_mod6_20260922_w6_lean_paste_audit.lean
Vertices of Q_n carry the VALUES 1..n (the pasted Fin-n version uses value = index+1).
-/

-- square test by bounded search (values x+y <= 2n <= 64 here)
def isSq (m : Nat) : Bool := Nat.any (m + 1) (fun k _ => k * k == m)

-- (1) the pasted adjacency, on Fin-style indices: value = i+1, NO x != y guard
def adjPaste (i j : Nat) : Bool := isSq ((i + 1) + (j + 1))

-- index 1 has value 2 and 2 + 2 = 4 = 2^2: the pasted relation is reflexive at that vertex
theorem adjPaste_not_irreflexive : adjPaste 1 1 = true := by decide
-- so the pasted `loopless := sorry` is unprovable (a SimpleGraph needs Adj irreflexive)
theorem pasted_loopless_fails : ¬ (∀ i, i < 14 → adjPaste i i = false) := by decide

-- corrected adjacency on values 1..n: x != y and x + y a square
def adj (x y : Nat) : Bool := (x != y) && isSq (x + y)

theorem adj_irrefl_upto_32 : ∀ x, x < 33 → adj x x = false := by decide
theorem adj_symm_upto_32 : ∀ x, x < 33 → ∀ y, y < 33 → adj x y = adj y x := by decide

-- (2) reachability by iterating a Boolean closure. `S : Nat → Bool` marks values 1..n.
def stepR (n : Nat) (S : Nat → Bool) : Nat → Bool :=
  fun v => S v || Nat.any n (fun u _ => S (u + 1) && adj (u + 1) v)

def closure (n : Nat) : Nat → Nat → Nat → Bool
  | 0, s => fun v => v == s
  | k + 1, s => stepR n (closure n k s)

-- reach n s v : v reachable from s in Q_n (n-1 closure rounds suffice for n vertices)
def reach (n s v : Nat) : Bool := closure n (n - 1) s v

def connected (n : Nat) : Bool := Nat.all n (fun v _ => reach n 1 (v + 1))

-- number of components = number of values v in 1..n with no smaller value reaching them
def numComponents (n : Nat) : Nat :=
  Nat.fold n (fun v _ acc =>
    if Nat.any v (fun u _ => reach n (u + 1) (v + 1)) then acc else acc + 1) 0

-- session lead probe, machine-checked: Q_14 connected, Q_13 and Q_12 not
theorem square_sum_fourteen_connectivity : connected 14 = true := by decide
theorem square_sum_thirteen_disconnected : connected 13 = false := by decide
theorem square_sum_twelve_disconnected : connected 12 = false := by decide
-- the pasted comment "merger of the 3 historical components" is off: 3 at 12, 2 at 13, 1 at 14
theorem components_twelve : numComponents 12 = 3 := by decide
theorem components_thirteen : numComponents 13 = 2 := by decide
theorem components_fourteen : numComponents 14 = 1 := by decide
theorem components_fifteen : numComponents 15 = 1 := by decide
theorem components_sixteen : numComponents 16 = 1 := by decide
theorem components_seventeen : numComponents 17 = 1 := by decide

-- (3) the "Delta = 4 prime ladder": the statement is TRUE ...
theorem delta_four_prime_step_invariance_list :
    ([3, 7, 11, 17].all fun p => [7, 11, 15, 21].contains (p + 4)) = true := by decide
-- ... and `trivial` also closes the Prop form in core Lean (it tries `decide`):
theorem delta_four_prime_step_invariance_prop :
    ∀ p, p ∈ [3, 7, 11, 17] → p + 4 ∈ [7, 11, 15, 21] := by trivial
-- ... but the targets 15 and 21 are composite, and 3,7,11,17 is not an AP (17 - 11 = 6)
def isPrime (n : Nat) : Bool := 2 ≤ n && Nat.all n (fun d _ => d < 2 || n % d != 0)
theorem fifteen_not_prime : isPrime 15 = false := by decide
theorem twentyone_not_prime : isPrime 21 = false := by decide
theorem ladder_not_ap : 17 - 11 ≠ 7 - 3 := by decide

-- (4) `is_positive_sheet : True` carries no information: the pasted structure is inhabited
-- for every n with EVERY choice of the descent data, e.g. sheet-free trivial parameters.
structure PastedSignSpecificCertificate (n : Nat) where
  L : Nat
  K_L : Nat
  B_L : Nat
  is_positive_sheet : True
  descent_inequality : 3 ^ L * n + B_L < 2 ^ K_L * n

-- inhabited for every n >= 1 with L = 0, K = 1, B = 0 (no orbit is mentioned at all)
def pastedCertificateAlwaysInhabited (n : Nat) (h : 0 < n) : PastedSignSpecificCertificate n :=
  { L := 0, K_L := 1, B_L := 0, is_positive_sheet := trivial,
    descent_inequality := by simp; omega }
-- and the field `is_positive_sheet` is a proof of `True`, i.e. a subsingleton with one value:
theorem is_positive_sheet_unique (c d : PastedSignSpecificCertificate 7) :
    c.is_positive_sheet = d.is_positive_sheet := rfl

-- correct minimal shape (compression_and_lean_audit sec. 7, cited, restated here, not re-proved):
def iterate (f : Nat → Nat) : Nat → Nat → Nat
  | 0, x => x
  | k + 1, x => iterate f k (f x)
def collatz (n : Nat) : Nat := if n % 2 = 0 then n / 2 else 3 * n + 1
def DescentCertificate (n : Nat) : Prop :=
  ∃ t K L B : Nat, 0 < t ∧
    2 ^ K * iterate collatz t n = 3 ^ L * n + B ∧
    3 ^ L * n + B < 2 ^ K * n
-- witness n = 3, t = 6, K = 4, L = 2, B = 5 (orbit 3,10,5,16,8,4,2)
theorem descent_three : DescentCertificate 3 := ⟨6, 4, 2, 5, by decide, by decide, by decide⟩
-- the certificate is sheet-blind by construction: the same shape with the 3n-1 map
def collatzMinus (n : Nat) : Nat := if n % 2 = 0 then n / 2 else 3 * n - 1
def DescentCertificateMinus (n : Nat) : Prop :=
  ∃ t K L B : Nat, 0 < t ∧
    2 ^ K * iterate collatzMinus t n = 3 ^ L * n + B ∧
    3 ^ L * n + B < 2 ^ K * n
-- witness n = 3 on the minus sheet: 3 -> 8 -> 4 -> 2 (t = 3, K = 2, L = 1, B = -1 is NOT a Nat,
-- so take t = 3, y = 2: 2^K * 2 = 3^L * 3 + B with K = 3, L = 1, B = 7: 16 = 9 + 7 < 24)
theorem descent_three_minus : DescentCertificateMinus 3 := ⟨3, 3, 1, 7, by decide, by decide, by decide⟩

#print axioms square_sum_fourteen_connectivity
#print axioms components_thirteen
#print axioms delta_four_prime_step_invariance_prop
#print axioms pastedCertificateAlwaysInhabited
#print axioms descent_three_minus
