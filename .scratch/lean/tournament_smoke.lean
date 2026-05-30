universe u

structure Tournament (V : Type u) where
  edge : V -> V -> Prop
  irrefl : forall v, Not (edge v v)
  total : forall u v, u ≠ v -> edge u v \/ edge v u
  antisymm : forall u v, edge u v -> Not (edge v u)

theorem tournament_no_loop {V : Type u} (T : Tournament V) (v : V) :
    Not (T.edge v v) := by
  exact T.irrefl v

theorem tournament_oriented_pair {V : Type u} (T : Tournament V) {u v : V}
    (h : u ≠ v) : T.edge u v \/ T.edge v u := by
  exact T.total u v h

theorem tournament_no_reverse {V : Type u} (T : Tournament V) {u v : V}
    (h : T.edge u v) : Not (T.edge v u) := by
  exact T.antisymm u v h

theorem lean_bootstrap_smoke : True := by
  trivial

inductive Two where
  | left
  | right
  deriving DecidableEq

open Two

def twoEdge : Two -> Two -> Prop
  | left, right => True
  | _, _ => False

def twoTournament : Tournament Two where
  edge := twoEdge
  irrefl := by
    intro v
    cases v <;> simp [twoEdge]
  total := by
    intro u v h
    cases u <;> cases v <;> simp [twoEdge] at *
  antisymm := by
    intro u v huv
    cases u <;> cases v <;> simp [twoEdge] at *

theorem two_left_beats_right : twoTournament.edge left right := by
  trivial

theorem two_right_does_not_beat_left : Not (twoTournament.edge right left) := by
  simp [twoTournament, twoEdge]
