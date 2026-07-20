/- GMC2.MomentBasics — the circular-pair moment functional, kernel-pure start
   (boxeph-2026-07-20-S177). Scope per S173: formal (a,b)-monomial sums over Int,
   E(monomial a b) = if a = b then 2^a * a! else 0; first lemmas: unbalanced
   vanishing (definitional) and the S171 concrete span fact E2 = 4 for
   P = aZ^2 + Z + W (a decide-checked instance of the closed span).
   Build wiring into the lake project = next session (flagged, not hidden). -/

def fact : Nat → Nat
  | 0 => 1
  | n+1 => (n+1) * fact n

def wt (a : Nat) : Int := (2 ^ a * fact a : Nat)

/-- a monomial c * Z^a * W^b -/
structure Mon where
  c : Int
  a : Nat
  b : Nat

def Emon (m : Mon) : Int := if m.a = m.b then m.c * wt m.a else 0

def E (p : List Mon) : Int := (p.map Emon).foldl (·+·) 0

theorem Emon_unbalanced (m : Mon) (h : m.a ≠ m.b) : Emon m = 0 := by
  simp [Emon, h]

/-- product of two monomials -/
def mmul (x y : Mon) : Mon := ⟨x.c * y.c, x.a + y.a, x.b + y.b⟩

def pmul (p q : List Mon) : List Mon := p.flatMap (fun x => q.map (mmul x))

/-- the S171 span instance: P = Z^2 + Z + W (a = 1 slice): E[P^2] = 4. -/
def Pspan : List Mon := [⟨1,2,0⟩, ⟨1,1,0⟩, ⟨1,0,1⟩]

theorem E_Pspan_sq : E (pmul Pspan Pspan) = 4 := by decide

/-- charge of a monomial -/
def charge (m : Mon) : Int := (m.a : Int) - (m.b : Int)

theorem Emon_of_charge_ne (m : Mon) (h : charge m ≠ 0) : Emon m = 0 := by
  apply Emon_unbalanced
  intro hab
  apply h
  simp [charge, hab]

/-- S171 parity-fake death, formalized: 2P = Z^2 + 2Z + 2W - W^2 (the (1/2,-1/2)
    candidate scaled by 2) survives E2 = E3 = 0 but dies at E4:
    E[(2P)^4] = 16 * (-96) = -1536 != 0. Kernel-checked closure of the 4-charge span. -/
def Pfake : List Mon := [⟨1,2,0⟩, ⟨2,1,0⟩, ⟨2,0,1⟩, ⟨-1,0,2⟩]

theorem E_Pfake_sq_zero : E (pmul Pfake Pfake) = 0 := by decide

set_option maxRecDepth 4096 in
theorem E_Pfake_four_ne : E (pmul (pmul Pfake Pfake) (pmul Pfake Pfake)) = -1536 := by decide
