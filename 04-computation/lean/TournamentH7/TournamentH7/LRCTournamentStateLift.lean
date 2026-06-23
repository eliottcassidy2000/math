/-
  TournamentH7.LRCTournamentStateLift

  A small formal endpoint for the tournament route to LRC(14).

  HYP-2908 says the missing theorem should lift a primitive LRC14 bad atom to a
  tournament-conflict packet with activity-two value 7.  This file does not
  prove that lift.  It records the exact closure theorem: any such lift is
  impossible, by the already formalized tournament theorem `H_ne_seven`.
-/

import TournamentH7.H7

namespace LonelyRunner

/-- A proposed activity-two packet lift into the tournament category.

`packetValue` is kept as a separate field because the intended LRC application
first constructs a packet independence value, then proves it agrees with the
tournament Hamiltonian-path count `H`.
-/
structure TournamentStateLift where
  n : ℕ
  T : Tournament n
  packetValue : ℕ
  packetValue_eq_H : packetValue = Tournament.H T
  packetValue_eq_seven : packetValue = 7

/-- No activity-two packet lift with value `7` can land in the tournament
category. -/
theorem no_tournament_state_lift (L : TournamentStateLift) : False := by
  have hH : Tournament.H L.T = 7 := by
    calc
      Tournament.H L.T = L.packetValue := L.packetValue_eq_H.symm
      _ = 7 := L.packetValue_eq_seven
  exact Tournament.H_ne_seven L.T hH

/-- State-lift closure form for a named bad-atom predicate.

To close a concrete LRC14 bad-atom class by the tournament route, it is enough
to construct a `TournamentStateLift` from every bad atom. -/
theorem not_bad_of_tournament_state_lift {Bad : Prop}
    (hLift : Bad → TournamentStateLift) : ¬ Bad := by
  intro hBad
  exact no_tournament_state_lift (hLift hBad)

/-- Equivalent closure form when the lift is already phrased as a tournament
with `H = 7`. -/
theorem not_bad_of_H_eq_seven_lift {Bad : Prop}
    (hLift : Bad → ∃ n : ℕ, ∃ T : Tournament n, Tournament.H T = 7) :
    ¬ Bad := by
  intro hBad
  rcases hLift hBad with ⟨n, T, hT⟩
  exact Tournament.H_ne_seven T hT

end LonelyRunner
