/-
  TournamentH7.RootPackets — endpoint boundaries and closed root packets

  This module packages the type-A root telescoping theorem from `RootSigns`
  into two small reusable records:

    • `TypeA.RootWalk`: an open vertex walk with source, target, and middle;
    • `TypeA.RootPacket`: a closed vertex walk based at one vertex.

  The point is intentionally modest.  Open walks have boundary `e_source -
  e_target`; closed packets have zero total root.  Existing tournament
  `DirectedCycle`s can therefore be viewed as closed type-A root packets
  before any later support, disjointness, or Omega algebra is added.
-/

import TournamentH7.Cycles
import TournamentH7.RootSigns
import Mathlib.Data.List.OfFn

namespace TypeA

/-- A finite open walk packaged by its endpoints and intermediate vertices. -/
structure RootWalk (n : ℕ) where
  source : Fin n
  target : Fin n
  middle : List (Fin n)

namespace RootWalk

variable {n : ℕ}

/-- The vertex list associated to an open root walk. -/
def vertices (W : RootWalk n) : List (Fin n) :=
  W.source :: W.middle ++ [W.target]

/-- The total type-A root carried by an open root walk. -/
def rootTotal (W : RootWalk n) : RootSpace n :=
  walkRootSum W.vertices

/-- The root total of an open walk is exactly its endpoint boundary. -/
@[simp] theorem rootTotal_eq_boundary (W : RootWalk n) :
    W.rootTotal = root W.source W.target := by
  exact walkRootSum_append_single W.source W.target W.middle

/-- A closed open-walk record has zero total root. -/
theorem rootTotal_eq_zero_of_closed (W : RootWalk n) (h : W.source = W.target) :
    W.rootTotal = 0 := by
  rw [rootTotal_eq_boundary W, h]
  simp

end RootWalk

/-- A closed root packet: a base vertex plus the middle of a closed walk. -/
structure RootPacket (n : ℕ) where
  base : Fin n
  middle : List (Fin n)

namespace RootPacket

variable {n : ℕ}

/-- The closed vertex list attached to a root packet. -/
def vertices (P : RootPacket n) : List (Fin n) :=
  P.base :: P.middle ++ [P.base]

/-- The total type-A root carried by a closed packet. -/
def rootTotal (P : RootPacket n) : RootSpace n :=
  walkRootSum P.vertices

/-- Closed root packets carry zero total type-A root. -/
@[simp] theorem rootTotal_eq_zero (P : RootPacket n) :
    P.rootTotal = 0 := by
  exact walkRootSum_closed P.base P.middle

/-- A closed packet as an open walk with equal endpoints. -/
def toRootWalk (P : RootPacket n) : RootWalk n :=
  { source := P.base, target := P.base, middle := P.middle }

@[simp] theorem toRootWalk_rootTotal (P : RootPacket n) :
    P.toRootWalk.rootTotal = 0 := by
  simp [toRootWalk]

end RootPacket

end TypeA

namespace Tournament

open TypeA

variable {n k : ℕ} {T : Tournament n}

namespace DirectedCycle

/-- The distinguished base vertex of a directed cycle, using index `0`. -/
def rootPacketBase (C : DirectedCycle T k) : Fin n :=
  C.toFun ⟨0, by
    have hk : 3 ≤ k := C.three_le
    omega⟩

/-- The remaining vertices of a directed cycle after the distinguished base.

For a cycle of length `k`, this is the list of vertices at indices
`1, ..., k-1`.  Closing the list by returning to `0` gives the root packet
walk `0,1,...,k-1,0`. -/
def rootPacketMiddle (C : DirectedCycle T k) : List (Fin n) :=
  List.ofFn fun i : Fin (k - 1) =>
    C.toFun ⟨i.val + 1, by
      have hk : 3 ≤ k := C.three_le
      have hi : i.val < k - 1 := i.isLt
      omega⟩

@[simp] theorem rootPacketMiddle_length (C : DirectedCycle T k) :
    C.rootPacketMiddle.length = k - 1 := by
  simp [rootPacketMiddle]

/-- A directed cycle, viewed only as its closed type-A root packet. -/
def toRootPacket (C : DirectedCycle T k) : TypeA.RootPacket n :=
  { base := C.rootPacketBase, middle := C.rootPacketMiddle }

/-- The root packet associated to a directed cycle has zero total root. -/
@[simp] theorem toRootPacket_rootTotal (C : DirectedCycle T k) :
    C.toRootPacket.rootTotal = 0 := by
  exact TypeA.RootPacket.rootTotal_eq_zero C.toRootPacket

end DirectedCycle

end Tournament
