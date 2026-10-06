/-
SIX-SEVEN CERTIFICATE (mac-mini-2026-10-06-sixseven).

Finite facts about the Paley tournament P_7 (i -> j iff j - i mod 7 is in QR_7 = {1, 2, 4}) used by
05-knowledge/results/sixes_and_sevens_20261006.md, checked by native_decide in the trust model of the
eleven-square packing formalization (Lean kernel + native compiler), Mathlib-free:

  (1) blocking: deleting any 5 of the 21 arcs leaves a Hamiltonian path, deleting the 6 arcs at a vertex
      leaves none, and exactly 63 six-arc sets block (the Hall bound N - 1 = 6 of a regular tournament);
  (2) THM-4524: in P_7 - 0 every one of the 15 arcs lies on an odd number of Hamiltonian paths;
  (3) the point-stripping ladder: |Aut(P_7)| = 21 and |Aut(P_7 - 0)| = 3 (the maps x -> x, 2x, 4x).

Python mirror: 04-computation/experiments/sixseven_20261006_structure.py (sections E and G).
-/

/-- arcs of P_7 as pairs (i, j), listed in a fixed order (21 of them). -/
def arcsP7 : List (Nat × Nat) :=
  (List.range 7).flatMap fun i => (List.range 7).filterMap fun j =>
    if i != j && (((j + 7 - i) % 7 == 1) || ((j + 7 - i) % 7 == 2) || ((j + 7 - i) % 7 == 4)) then some (i, j) else none

/-- out-neighbourhood bitmask of vertex v in the arc list A. -/
def outMask (A : List (Nat × Nat)) (v : Nat) : Nat :=
  A.foldl (fun m e => if e.1 == v then m ||| (1 <<< e.2) else m) 0

/-- number of Hamiltonian paths of the digraph A on the vertex set given by bitmask `all` (vertices < 7),
    by dynamic programming over subsets: cnt S v = number of paths covering exactly S and ending at v. -/
def hpCount (A : List (Nat × Nat)) (all : Nat) : Nat := Id.run do
  let outs := (List.range 7).map (outMask A)
  let mut cnt : Array Nat := Array.replicate (128 * 7) 0
  for v in List.range 7 do
    if (all >>> v) % 2 == 1 then
      cnt := cnt.set! ((1 <<< v) * 7 + v) 1
  for S in List.range 128 do
    for v in List.range 7 do
      let c := cnt[S * 7 + v]!
      if c != 0 then
        let ov := outs.getD v 0
        for w in List.range 7 do
          if (ov >>> w) % 2 == 1 && (S >>> w) % 2 == 0 && (all >>> w) % 2 == 1 then
            let T := S ||| (1 <<< w)
            cnt := cnt.set! (T * 7 + w) (cnt[T * 7 + w]! + c)
  let mut tot := 0
  for v in List.range 7 do
    tot := tot + cnt[all * 7 + v]!
  return tot

def hasHP (A : List (Nat × Nat)) : Bool := hpCount A 127 != 0

/-- all k-element sublists (as index lists) of [0, n). -/
def choose : Nat → Nat → List (List Nat)
  | _, 0 => [[]]
  | 0, _ + 1 => []
  | n + 1, k + 1 => (choose n k).map (fun s => s ++ [n]) ++ choose n (k + 1)

def deleteIdx (A : List (Nat × Nat)) (D : List Nat) : List (Nat × Nat) :=
  (List.range A.length).filterMap fun i => if D.contains i then none else A[i]?

theorem p7_has_21_arcs : arcsP7.length = 21 := by native_decide

/-- (1a) every 5-arc deletion leaves a Hamiltonian path (20349 sets). -/
theorem p7_survives_five : (choose 21 5).all (fun D => hasHP (deleteIdx arcsP7 D)) = true := by native_decide

/-- (1b) the star at vertex 0 (all 6 arcs at 0) blocks every Hamiltonian path. -/
theorem p7_star_blocks :
    hasHP (arcsP7.filter fun e => e.1 != 0 && e.2 != 0) = false := by native_decide

/-- (1c) exactly 63 of the 54264 six-arc sets block. -/
theorem p7_63_minimum_blocking_sets :
    ((choose 21 6).filter (fun D => !hasHP (deleteIdx arcsP7 D))).length = 63 := by native_decide

/-- P_7 - 0: vertices 1..6 (bitmask 126). -/
def arcsP7minus0 : List (Nat × Nat) := arcsP7.filter fun e => e.1 != 0 && e.2 != 0

/-- (2) THM-4524 at N = 6: every arc e of P_7 - 0 lies on an odd number H(T) - H(T - e) of Hamiltonian paths. -/
theorem p7minus0_all_arcs_odd :
    let H := hpCount arcsP7minus0 126
    arcsP7minus0.length = 15 &&
    arcsP7minus0.all (fun e => (H - hpCount (arcsP7minus0.filter (· != e)) 126) % 2 == 1) = true := by
  native_decide

/-- all permutations of a list. -/
def perms : List Nat → List (List Nat)
  | [] => [[]]
  | x :: xs => (perms xs).flatMap fun p => (List.range (p.length + 1)).map fun i => p.take i ++ [x] ++ p.drop i

/-- p (a list giving the images of the vertices `dom`, in order) preserves the arc relation of A. -/
def isAut (A : List (Nat × Nat)) (dom : List Nat) (p : List Nat) : Bool :=
  let img := fun v => p.getD (dom.idxOf v) 99
  A.all fun e => A.contains (img e.1, img e.2)

/-- (3) |Aut(P_7)| = 21 (the Borel subgroup x -> ax + b, a in QR_7, of PSL(2,7)). -/
theorem aut_p7 : ((perms (List.range 7)).filter (isAut arcsP7 (List.range 7))).length = 21 := by native_decide

/-- (3) |Aut(P_7 - 0)| = 3 ... -/
theorem aut_p7minus0_card :
    ((perms [1, 2, 3, 4, 5, 6]).filter (isAut arcsP7minus0 [1, 2, 3, 4, 5, 6])).length = 3 := by native_decide

/-- ... and they are the multiplications x -> x, 2x, 4x (the split torus, i.e. the Collatz map T on its
    trivial cycle {1, 4, 2} read through the code QR_7; image lists of the vertices 1..6). -/
theorem aut_p7minus0_are_multiplications :
    ((perms [1, 2, 3, 4, 5, 6]).filter (isAut arcsP7minus0 [1, 2, 3, 4, 5, 6])).all
      (fun p => p == [1, 2, 3, 4, 5, 6] || p == [2, 4, 6, 1, 3, 5] || p == [4, 1, 5, 2, 6, 3]) = true := by
  native_decide

-- axiom audit: each theorem should rest only on its own native_decide axiom (no sorryAx)
#print axioms p7_survives_five
#print axioms p7_star_blocks
#print axioms p7_63_minimum_blocking_sets
#print axioms p7minus0_all_arcs_odd
#print axioms aut_p7
#print axioms aut_p7minus0_card
#print axioms aut_p7minus0_are_multiplications
