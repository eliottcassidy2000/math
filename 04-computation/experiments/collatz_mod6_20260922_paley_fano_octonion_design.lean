/- Lane paley_fano_octonion_design: core-Lean (no Mathlib) `decide` checks on the Paley tournament T_7.
   Run: lean 04-computation/experiments/collatz_mod6_20260922_paley_fano_octonion_design.lean -/

/-- Paley arc on Z/7: i -> j iff (j - i) mod 7 is a nonzero quadratic residue {1,2,4}. -/
def arc (i j : Nat) : Bool :=
  let d := (j + 7 - i % 7) % 7
  d == 1 || d == 2 || d == 4

def range7 : List Nat := [0, 1, 2, 3, 4, 5, 6]

/-- unordered triples a < b < c of Z/7 -/
def triples : List (Nat × Nat × Nat) :=
  range7.foldr (fun a acc => (range7.foldr (fun b acc2 => (range7.foldr (fun c acc3 =>
    if a < b && b < c then (a, b, c) :: acc3 else acc3) [] ) ++ acc2) []) ++ acc) []

def cyclic (t : Nat × Nat × Nat) : Bool :=
  let (a, b, c) := t
  (arc a b && arc b c && arc c a) || (arc a c && arc c b && arc b a)

def numArcs : Nat := (range7.foldr (fun i acc => (range7.foldr (fun j acc2 =>
  if arc i j then acc2 + 1 else acc2) 0) + acc) 0)

def numCyclic : Nat := (triples.filter cyclic).length
def numTriples : Nat := triples.length

/-- lines {s, s+1, s+3} carry the orientation s -> s+1 -> s+3 -> s -/
def dev013Oriented : Bool := range7.all fun s =>
  arc s ((s + 1) % 7) && arc ((s + 1) % 7) ((s + 3) % 7) && arc ((s + 3) % 7) s

/-- lines {s, s+1, s+5} carry the orientation s -> s+1 -> s+5 -> s -/
def dev015Oriented : Bool := range7.all fun s =>
  arc s ((s + 1) % 7) && arc ((s + 1) % 7) ((s + 5) % 7) && arc ((s + 5) % 7) s

/-- every ordered pair (i,j), i != j, is an arc in exactly one direction -/
def isTournament : Bool := range7.all fun i => range7.all fun j =>
  (i == j) || (arc i j != arc j i)

theorem paley7_is_tournament : isTournament = true := by decide
theorem paley7_arcs : numArcs = 21 := by decide
theorem paley7_triples : numTriples = 35 := by decide
theorem paley7_cyclic_triples : numCyclic = 14 := by decide
theorem paley7_dev013 : dev013Oriented = true := by decide
theorem paley7_dev015 : dev015Oriented = true := by decide

#print axioms paley7_cyclic_triples
#print axioms paley7_dev013
