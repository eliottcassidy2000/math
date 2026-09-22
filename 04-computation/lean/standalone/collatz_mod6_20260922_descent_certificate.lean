/- A correct minimal Lean statement of the descent certificate (scratch; NOT part of
   the CollatzBlueprintAudit package, which already proves the equivalence
   `collatz_reachesOne_iff_globalTallyMargin` for the true orbit counters).
   Core Lean 4.30.0 only; no Mathlib. -/
def iterate {α : Type} (f : α → α) : Nat → α → α
  | 0, x => x
  | k + 1, x => iterate f k (f x)

def collatz (n : Nat) : Nat := if n % 2 = 0 then n / 2 else 3 * n + 1

/-- One start has a descent certificate: a positive time `t` at which the orbit value `y`
    satisfies the exact affine identity `2^K y = 3^L n + B` AND the margin `3^L n + B < 2^K n`.
    Both `K, L, B` are existentially bound but tied to the actual orbit through the identity. -/
def DescentCertificate (n : Nat) : Prop :=
  ∃ t K L B : Nat, 0 < t ∧
    2 ^ K * iterate collatz t n = 3 ^ L * n + B ∧
    3 ^ L * n + B < 2 ^ K * n

/-- The global (still OPEN) target: every `n > 1` has one. This is a definition, not a theorem. -/
def GlobalDescent : Prop := ∀ n, 1 < n → DescentCertificate n

/-- The certificate is nonvacuous: `n = 3` is certified at `t = 6` (orbit 3,10,5,16,8,4,2),
    with `K = 4` halvings, `L = 2` tripling steps and carry `B = 5`: `16*2 = 9*3 + 5 = 32 < 48`. -/
theorem three_certified : DescentCertificate 3 :=
  ⟨6, 4, 2, 5, by decide, by decide, by decide⟩

/-- `n = 7` is certified at `t = 11` (orbit 7,22,11,34,17,52,26,13,40,20,10,5 reaches 5 < 7):
    `K = 7` halvings, `L = 4` tripling steps, carry `B = 73`: `2^7*5 = 640 = 81*7 + 73 < 896`. -/
theorem seven_certified : DescentCertificate 7 :=
  ⟨11, 7, 4, 73, by decide, by decide, by decide⟩
