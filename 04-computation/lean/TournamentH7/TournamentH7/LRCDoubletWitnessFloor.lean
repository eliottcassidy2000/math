/-
  TournamentH7.LRCDoubletWitnessFloor -- arithmetic import boundary for the
  genuine-wide doublet rho*/witness-floor scout.

  The exact measurement engine and finite scout live in
  `04-computation/lrc14_rhostar_gw_doublet_claudeopus_0622s4.py` with stored
  output `05-knowledge/results/lrc14_rhostar_gw_doublet_claudeopus_0622s4.out`.

  This module records only the final rational floor comparisons reported there:

    * audited doublet `rho*` floor: `2/147`;
    * audited doublet witness floor: `1066/2205`;
    * `2/147` is strictly above the old consecutive `rho*` floor `1/84`;
    * `1066/2205` is strictly above the THM-530 witness floor
      `14249/252252`.

  It is not a proof that the scout covers all configurations.  It is a
  Lean-facing arithmetic checksum for the imported finite evidence.
-/

namespace LonelyRunner
namespace DoubletWitnessFloor

/-- Reported global minimum of the audited genuine-wide doublet `rho*` scout. -/
def rhoFloorNum : Nat := 2
def rhoFloorDen : Nat := 147

/-- Reported global minimum of the audited genuine-wide doublet witness scout. -/
def witnessFloorNum : Nat := 1066
def witnessFloorDen : Nat := 2205

/-- The older consecutive `rho*` floor `1/84`, used as a comparison point. -/
def consecRhoFloorNum : Nat := 1
def consecRhoFloorDen : Nat := 84

/-- The THM-530 admissible witness floor `14249/252252`. -/
def thm530WitnessFloorNum : Nat := 14249
def thm530WitnessFloorDen : Nat := 252252

/-- The audited doublet `rho*` floor is positive. -/
theorem rhoFloor_positive : 0 < rhoFloorNum := by
  native_decide

/-- The audited doublet witness floor is positive. -/
theorem witnessFloor_positive : 0 < witnessFloorNum := by
  native_decide

/-- `2/147 > 1/84`: the audited doublet `rho*` floor is above the consecutive
floor reported in the THM-527 route. -/
theorem rhoFloor_above_consecRhoFloor :
    consecRhoFloorNum * rhoFloorDen < rhoFloorNum * consecRhoFloorDen := by
  native_decide

/-- `1066/2205 > 14249/252252`: the audited doublet witness floor is above the
THM-530 admissible witness floor. -/
theorem witnessFloor_above_thm530Floor :
    thm530WitnessFloorNum * witnessFloorDen <
      witnessFloorNum * thm530WitnessFloorDen := by
  native_decide

/-- The reported witness floor is also above the reported `rho*` floor. -/
theorem witnessFloor_above_rhoFloor :
    rhoFloorNum * witnessFloorDen < witnessFloorNum * rhoFloorDen := by
  native_decide

/-- Packaged checksum for the imported doublet floor scout. -/
theorem audited_doublet_floors_positive_and_separated :
    0 < rhoFloorNum ∧
      0 < witnessFloorNum ∧
      consecRhoFloorNum * rhoFloorDen < rhoFloorNum * consecRhoFloorDen ∧
      thm530WitnessFloorNum * witnessFloorDen <
        witnessFloorNum * thm530WitnessFloorDen ∧
      rhoFloorNum * witnessFloorDen < witnessFloorNum * rhoFloorDen := by
  exact ⟨rhoFloor_positive, witnessFloor_positive, rhoFloor_above_consecRhoFloor,
    witnessFloor_above_thm530Floor, witnessFloor_above_rhoFloor⟩

/-! ### Axiom audit -/

#print axioms rhoFloor_positive
#print axioms witnessFloor_positive
#print axioms rhoFloor_above_consecRhoFloor
#print axioms witnessFloor_above_thm530Floor
#print axioms witnessFloor_above_rhoFloor
#print axioms audited_doublet_floors_positive_and_separated

end DoubletWitnessFloor
end LonelyRunner
