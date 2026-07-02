/-
  TournamentH7.LRCRandom031ProofPackets -- Lean-facing random031 proof-packet
  ledger for HYP-3521 through HYP-3525.

  This module does not prove the remaining random031 geometric discharge.  It
  formalizes the recent proof interface as sorry-free Lean data and logical
  glue:

    * HYP-3521 terminal split: 282 = 230 ordinary + 40 free-hole + 12 bypass,
      and 79 legal components collapse to 77 terminal certificates.
    * HYP-3522 owner filtration: seven seam owners emit transport
      `(23,93,113)`, then branch-boundary lift `(147,169)`, leaving residual
      `(45,173)`.
    * HYP-3523/3524/3525 spigot language: proof tokens may be emitted only
      when the visible sidecars force the target, otherwise the first missing
      sidecar remains named debt.
    * HYP-3526/3527 contract routing: the current random031 proof interface is
      eight theorem-facing clauses, with 3 formal-ready clauses, 4 carry
      clauses, and one residual open-tail lemma.

  The main mathematical obligations are still external producer theorems:
  ordinary rank-2 route discharge, free-hole bracket discharge, pure-bypass
  owner-boundary discharge, route reconstruction or sidecar `R`, and the
  no-hidden-tail residual lemma.  This file makes those obligations explicit
  and machine-checks the finite arithmetic around them.
-/

import Mathlib.Tactic
import TournamentH7.LRCPrivateLabelFirewall

namespace LonelyRunner
namespace LRC14

/-! ## HYP-3521 terminal certificate ledger -/

/-- The three cell classes in the random031 terminal ledger. -/
inductive Random031CellClass where
  | ordinary_rank2
  | bracketed_free_hole
  | pure_bypass_owner_boundary
deriving DecidableEq, Repr, Fintype

/-- Exact HYP-3521 cell counts. -/
def random031CellCount : Random031CellClass -> Nat
  | Random031CellClass.ordinary_rank2 => 230
  | Random031CellClass.bracketed_free_hole => 40
  | Random031CellClass.pure_bypass_owner_boundary => 12

theorem random031_cell_total :
    random031CellCount Random031CellClass.ordinary_rank2 +
      random031CellCount Random031CellClass.bracketed_free_hole +
      random031CellCount Random031CellClass.pure_bypass_owner_boundary = 282 := by
  decide

theorem random031_gate_routed_split :
    random031CellCount Random031CellClass.ordinary_rank2 +
      random031CellCount Random031CellClass.pure_bypass_owner_boundary = 242 := by
  decide

/-- Terminal-certificate classes after the HYP-3511 free-hole doublet
collapse. -/
inductive Random031TerminalCertificate where
  | ordinary_route
  | free_hole_single
  | free_hole_doublet
  | pure_bypass
deriving DecidableEq, Repr, Fintype

def random031TerminalCertificateCount : Random031TerminalCertificate -> Nat
  | Random031TerminalCertificate.ordinary_route => 64
  | Random031TerminalCertificate.free_hole_single => 10
  | Random031TerminalCertificate.free_hole_doublet => 2
  | Random031TerminalCertificate.pure_bypass => 1

theorem random031_terminal_certificate_total :
    random031TerminalCertificateCount Random031TerminalCertificate.ordinary_route +
      random031TerminalCertificateCount Random031TerminalCertificate.free_hole_single +
      random031TerminalCertificateCount Random031TerminalCertificate.free_hole_doublet +
      random031TerminalCertificateCount Random031TerminalCertificate.pure_bypass = 77 := by
  decide

theorem random031_legal_component_total_before_collapse :
    64 + 14 + 1 = 79 := by
  decide

theorem random031_free_hole_collapse_count :
    random031TerminalCertificateCount Random031TerminalCertificate.free_hole_single +
      2 * random031TerminalCertificateCount Random031TerminalCertificate.free_hole_doublet =
        14 := by
  decide

/-! ## HYP-3522 owner-boundary filtration -/

/-- Owner labels are kept numeric so they match the Python certificates. -/
abbrev OwnerLabel := Nat

def random031SeamOwners : List OwnerLabel :=
  [23, 45, 93, 113, 147, 169, 173]

def random031TransportOwners : List OwnerLabel :=
  [23, 93, 113]

def random031BranchBoundaryLiftOwners : List OwnerLabel :=
  [147, 169]

def random031ResidualOwners : List OwnerLabel :=
  [45, 173]

def random031TransportPlusBoundaryOwners : List OwnerLabel :=
  random031TransportOwners ++ random031BranchBoundaryLiftOwners

theorem random031_seam_owner_count :
    random031SeamOwners.length = 7 := by
  decide

theorem random031_transport_owner_count :
    random031TransportOwners.length = 3 := by
  decide

theorem random031_branch_boundary_lift_count :
    random031BranchBoundaryLiftOwners.length = 2 := by
  decide

theorem random031_residual_owner_count :
    random031ResidualOwners.length = 2 := by
  decide

theorem random031_transport_plus_boundary_count :
    random031TransportPlusBoundaryOwners.length = 5 := by
  decide

theorem random031_owner_filtration_sizes :
    random031SeamOwners.length = 7 ∧
    random031TransportOwners.length = 3 ∧
    random031TransportPlusBoundaryOwners.length = 5 ∧
    random031ResidualOwners.length = 2 := by
  decide

theorem random031_owner_words_nodup :
    random031SeamOwners.Nodup ∧
    random031TransportOwners.Nodup ∧
    random031BranchBoundaryLiftOwners.Nodup ∧
    random031ResidualOwners.Nodup ∧
    random031TransportPlusBoundaryOwners.Nodup := by
  decide

theorem random031_owner_filtration_disjoint_payloads :
    random031TransportOwners.Disjoint random031ResidualOwners ∧
    random031BranchBoundaryLiftOwners.Disjoint random031ResidualOwners ∧
    random031TransportOwners.Disjoint random031BranchBoundaryLiftOwners := by
  constructor
  · intro x hx hy
    simp [random031TransportOwners, random031ResidualOwners] at hx hy
    aesop
  constructor
  · intro x hx hy
    simp [random031BranchBoundaryLiftOwners, random031ResidualOwners] at hx hy
    aesop
  · intro x hx hy
    simp [random031TransportOwners, random031BranchBoundaryLiftOwners] at hx hy
    aesop

theorem random031_seam_decomposes_as_transport_boundary_residual :
    random031SeamOwners.Perm
      (random031TransportOwners ++
        random031BranchBoundaryLiftOwners ++
        random031ResidualOwners) := by
  decide

/-- The online emitter states from HYP-3524. -/
inductive Random031EmitterState where
  | seam_input
  | transport_emitted
  | boundary_lift_emitted
  | residual_tail
deriving DecidableEq, Repr, Fintype

def random031EmittedOwners : Random031EmitterState -> List OwnerLabel
  | Random031EmitterState.seam_input => []
  | Random031EmitterState.transport_emitted => random031TransportOwners
  | Random031EmitterState.boundary_lift_emitted => random031TransportPlusBoundaryOwners
  | Random031EmitterState.residual_tail => random031TransportPlusBoundaryOwners

def random031TailOwners : Random031EmitterState -> List OwnerLabel
  | Random031EmitterState.seam_input => random031SeamOwners
  | Random031EmitterState.transport_emitted =>
      random031BranchBoundaryLiftOwners ++ random031ResidualOwners
  | Random031EmitterState.boundary_lift_emitted => random031ResidualOwners
  | Random031EmitterState.residual_tail => random031ResidualOwners

theorem random031_emitter_cumulative_sizes :
    (random031EmittedOwners Random031EmitterState.seam_input).length = 0 ∧
    (random031EmittedOwners Random031EmitterState.transport_emitted).length = 3 ∧
    (random031EmittedOwners Random031EmitterState.boundary_lift_emitted).length = 5 ∧
    (random031EmittedOwners Random031EmitterState.residual_tail).length = 5 := by
  decide

theorem random031_emitter_tail_sizes :
    (random031TailOwners Random031EmitterState.seam_input).length = 7 ∧
    (random031TailOwners Random031EmitterState.transport_emitted).length = 4 ∧
    (random031TailOwners Random031EmitterState.boundary_lift_emitted).length = 2 ∧
    (random031TailOwners Random031EmitterState.residual_tail).length = 2 := by
  decide

theorem random031_emitter_final_tail :
    random031TailOwners Random031EmitterState.residual_tail =
      random031ResidualOwners := by
  rfl

/-! ## Guarded emission sidecars -/

/-- Sidecars visible to the HYP-3525 guarded-emission rule. -/
inductive Random031Sidecar where
  | C
  | F
  | N
  | T
  | I
  | Q
  | R
  | route_reconstruction
  | flow_class
  | allowed_exit
  | owner_union
  | sheet_pgf_bucket
  | transport_word
  | branch_boundary_lift
  | residual_pair
  | checksum_only
  | endpoint_ranks
  | branch_hist
  | component_size
  | mirror_closed
  | I_or_Q_private_cut
  | endpoint_rank2_route
  | vertical_halfturn_guard
  | owner_word
  | ordinary_boundary_neighbors
  | owner_support_chamber
  | no_hidden_tail
  | HYP3511_single_bracket
  | HYP3511_doublet_wait
  | HYP3511_same_branch_doublet
deriving DecidableEq, Repr, Fintype

/-- Tokens one may try to emit into the proof. -/
inductive Random031EmissionTarget where
  | private_firewall_status
  | h3490_route
  | random031_terminal_class
  | owner_residual_pair
  | raw_count
deriving DecidableEq, Repr, Fintype

def hasSidecar (visible : List Random031Sidecar) (s : Random031Sidecar) : Bool :=
  decide (s ∈ visible)

def hasAnySidecar (visible : List Random031Sidecar) (needed : List Random031Sidecar) : Bool :=
  needed.any (fun s => hasSidecar visible s)

def hasAllSidecars (visible : List Random031Sidecar) (needed : List Random031Sidecar) : Bool :=
  needed.all (fun s => hasSidecar visible s)

/-- The executable HYP-3525 emission guard. -/
def random031CanEmit
    (visible : List Random031Sidecar)
    (target : Random031EmissionTarget) : Bool :=
  match target with
  | Random031EmissionTarget.private_firewall_status =>
      hasAnySidecar visible
        [Random031Sidecar.C, Random031Sidecar.F, Random031Sidecar.N,
         Random031Sidecar.T, Random031Sidecar.I, Random031Sidecar.Q,
         Random031Sidecar.R]
  | Random031EmissionTarget.h3490_route =>
      hasAnySidecar visible
        [Random031Sidecar.R, Random031Sidecar.route_reconstruction]
  | Random031EmissionTarget.random031_terminal_class =>
      hasAllSidecars visible
        [Random031Sidecar.flow_class, Random031Sidecar.allowed_exit,
         Random031Sidecar.sheet_pgf_bucket]
  | Random031EmissionTarget.owner_residual_pair =>
      hasAllSidecars visible
        [Random031Sidecar.transport_word, Random031Sidecar.branch_boundary_lift,
         Random031Sidecar.residual_pair, Random031Sidecar.R]
  | Random031EmissionTarget.raw_count =>
      hasSidecar visible Random031Sidecar.checksum_only

def coloredPrivateOnlySidecars : List Random031Sidecar := [Random031Sidecar.C]

def routeSidecarOnly : List Random031Sidecar := [Random031Sidecar.R]

def safeSheafHeadSidecars : List Random031Sidecar :=
  [Random031Sidecar.flow_class, Random031Sidecar.allowed_exit,
   Random031Sidecar.sheet_pgf_bucket]

def ownerFiltrationReadySidecars : List Random031Sidecar :=
  [Random031Sidecar.transport_word, Random031Sidecar.branch_boundary_lift,
   Random031Sidecar.residual_pair, Random031Sidecar.R]

def rawCountSidecars : List Random031Sidecar := [Random031Sidecar.checksum_only]

def dangerEndpointRankSidecars : List Random031Sidecar := [Random031Sidecar.endpoint_ranks]

theorem colored_private_only_emits_status_not_route :
    random031CanEmit coloredPrivateOnlySidecars
        Random031EmissionTarget.private_firewall_status = true ∧
    random031CanEmit coloredPrivateOnlySidecars
        Random031EmissionTarget.h3490_route = false := by
  decide

theorem route_sidecar_emits_status_and_route :
    random031CanEmit routeSidecarOnly
        Random031EmissionTarget.private_firewall_status = true ∧
    random031CanEmit routeSidecarOnly
        Random031EmissionTarget.h3490_route = true := by
  decide

theorem safe_sheaf_head_emits_terminal_not_route :
    random031CanEmit safeSheafHeadSidecars
        Random031EmissionTarget.random031_terminal_class = true ∧
    random031CanEmit safeSheafHeadSidecars
        Random031EmissionTarget.h3490_route = false := by
  decide

theorem owner_filtration_ready_emits_residual_and_route :
    random031CanEmit ownerFiltrationReadySidecars
        Random031EmissionTarget.owner_residual_pair = true ∧
    random031CanEmit ownerFiltrationReadySidecars
        Random031EmissionTarget.h3490_route = true := by
  decide

theorem endpoint_rank_shadow_emits_nothing_theorem_facing :
    (forall target : Random031EmissionTarget,
      target ≠ Random031EmissionTarget.raw_count ->
        random031CanEmit dangerEndpointRankSidecars target = false) := by
  decide

theorem raw_count_shadow_emits_only_raw_count :
    random031CanEmit rawCountSidecars Random031EmissionTarget.raw_count = true ∧
    random031CanEmit rawCountSidecars
        Random031EmissionTarget.random031_terminal_class = false ∧
    random031CanEmit rawCountSidecars
        Random031EmissionTarget.owner_residual_pair = false := by
  decide

/-! ## HYP-3526/HYP-3527 route sidecars and proof contracts -/

/- The contract clauses are used as tournament-style vertices: the observable is
formal readiness, hidden-tail risk, and scalar-forgetting risk; the switch is
higher score, with the displayed path as the fixed Hamiltonian tie order. -/
inductive Random031ContractClause where
  | ordinary_route_emit
  | free_hole_single_emit
  | free_hole_doublet_buffer_emit
  | bypass_transport_emit
  | bypass_bracket_lift_emit
  | private_firewall_route_sidecar
  | vertical_halfturn_guard
  | residual_pair_close_tail
deriving DecidableEq, Repr, Fintype

def random031ContractPath : List Random031ContractClause :=
  [Random031ContractClause.ordinary_route_emit,
   Random031ContractClause.free_hole_single_emit,
   Random031ContractClause.free_hole_doublet_buffer_emit,
   Random031ContractClause.bypass_transport_emit,
   Random031ContractClause.bypass_bracket_lift_emit,
   Random031ContractClause.private_firewall_route_sidecar,
   Random031ContractClause.vertical_halfturn_guard,
   Random031ContractClause.residual_pair_close_tail]

def random031ContractScore : Random031ContractClause -> Int
  | Random031ContractClause.ordinary_route_emit => 94
  | Random031ContractClause.free_hole_single_emit => 90
  | Random031ContractClause.free_hole_doublet_buffer_emit => 88
  | Random031ContractClause.bypass_transport_emit => 82
  | Random031ContractClause.bypass_bracket_lift_emit => 80
  | Random031ContractClause.private_firewall_route_sidecar => 78
  | Random031ContractClause.vertical_halfturn_guard => 74
  | Random031ContractClause.residual_pair_close_tail => 63

inductive Random031ContractStatus where
  | formal_ready_interface
  | carry_required
  | open_tail_lemma
deriving DecidableEq, Repr, Fintype

def random031ContractStatus : Random031ContractClause -> Random031ContractStatus
  | Random031ContractClause.ordinary_route_emit =>
      Random031ContractStatus.formal_ready_interface
  | Random031ContractClause.free_hole_single_emit =>
      Random031ContractStatus.formal_ready_interface
  | Random031ContractClause.free_hole_doublet_buffer_emit =>
      Random031ContractStatus.formal_ready_interface
  | Random031ContractClause.bypass_transport_emit =>
      Random031ContractStatus.carry_required
  | Random031ContractClause.bypass_bracket_lift_emit =>
      Random031ContractStatus.carry_required
  | Random031ContractClause.private_firewall_route_sidecar =>
      Random031ContractStatus.carry_required
  | Random031ContractClause.vertical_halfturn_guard =>
      Random031ContractStatus.carry_required
  | Random031ContractClause.residual_pair_close_tail =>
      Random031ContractStatus.open_tail_lemma

def random031ContractStatusCount (status : Random031ContractStatus) : Nat :=
  random031ContractPath.foldl
    (fun acc clause =>
      if random031ContractStatus clause = status then acc + 1 else acc)
    0

theorem random031_contract_path_length :
    random031ContractPath.length = 8 := by
  decide

theorem random031_contract_path_nodup :
    random031ContractPath.Nodup := by
  decide

theorem random031_contract_score_path :
    random031ContractPath.map random031ContractScore =
      [94, 90, 88, 82, 80, 78, 74, 63] := by
  decide

theorem random031_contract_status_hist :
    random031ContractStatusCount
        Random031ContractStatus.formal_ready_interface = 3 ∧
    random031ContractStatusCount
        Random031ContractStatus.carry_required = 4 ∧
    random031ContractStatusCount
        Random031ContractStatus.open_tail_lemma = 1 := by
  decide

def random031OpenTailClauses : List Random031ContractClause :=
  random031ContractPath.filter
    (fun clause =>
      random031ContractStatus clause =
        Random031ContractStatus.open_tail_lemma)

theorem random031_unique_open_tail_clause :
    random031OpenTailClauses =
      [Random031ContractClause.residual_pair_close_tail] := by
  decide

def random031ClauseRequiredSidecars :
    Random031ContractClause -> List Random031Sidecar
  | Random031ContractClause.ordinary_route_emit =>
      [Random031Sidecar.endpoint_rank2_route, Random031Sidecar.R,
       Random031Sidecar.vertical_halfturn_guard]
  | Random031ContractClause.free_hole_single_emit =>
      [Random031Sidecar.HYP3511_single_bracket, Random031Sidecar.R,
       Random031Sidecar.vertical_halfturn_guard]
  | Random031ContractClause.free_hole_doublet_buffer_emit =>
      [Random031Sidecar.HYP3511_doublet_wait,
       Random031Sidecar.HYP3511_same_branch_doublet, Random031Sidecar.R,
       Random031Sidecar.vertical_halfturn_guard]
  | Random031ContractClause.bypass_transport_emit =>
      [Random031Sidecar.transport_word, Random031Sidecar.R,
       Random031Sidecar.owner_word]
  | Random031ContractClause.bypass_bracket_lift_emit =>
      [Random031Sidecar.branch_boundary_lift, Random031Sidecar.R,
       Random031Sidecar.ordinary_boundary_neighbors]
  | Random031ContractClause.private_firewall_route_sidecar =>
      [Random031Sidecar.I_or_Q_private_cut, Random031Sidecar.R]
  | Random031ContractClause.vertical_halfturn_guard =>
      [Random031Sidecar.vertical_halfturn_guard]
  | Random031ContractClause.residual_pair_close_tail =>
      [Random031Sidecar.residual_pair, Random031Sidecar.R,
       Random031Sidecar.owner_support_chamber, Random031Sidecar.no_hidden_tail]

def random031RouteBearingClauses : List Random031ContractClause :=
  [Random031ContractClause.ordinary_route_emit,
   Random031ContractClause.free_hole_single_emit,
   Random031ContractClause.free_hole_doublet_buffer_emit,
   Random031ContractClause.bypass_transport_emit,
   Random031ContractClause.bypass_bracket_lift_emit,
   Random031ContractClause.private_firewall_route_sidecar,
   Random031ContractClause.residual_pair_close_tail]

theorem random031_route_bearing_clauses_retain_R :
    random031RouteBearingClauses.all
      (fun clause =>
        hasSidecar (random031ClauseRequiredSidecars clause)
          Random031Sidecar.R) = true := by
  decide

theorem random031_vertical_guard_is_projection_not_route :
    hasSidecar
        (random031ClauseRequiredSidecars
          Random031ContractClause.vertical_halfturn_guard)
        Random031Sidecar.R = false ∧
    hasSidecar
        (random031ClauseRequiredSidecars
          Random031ContractClause.vertical_halfturn_guard)
        Random031Sidecar.vertical_halfturn_guard = true := by
  decide

/- HYP-3526 dispatch compressors.  This is the finite row-free/route-pure
guardrail: I/Q proves private status but not route; R is the current route
sidecar. -/
inductive Random031DispatchCompressor where
  | IQ_plus_R_terminal_spigot
  | terminal_certificate_ledger_plus_R
  | owner_filtration_plus_R
  | IQ_without_R_private_only
  | row_name_exception_list
  | raw_private_bit
  | all_colored_plus_IQ_route_reconstruction
  | raw_count_shadow
deriving DecidableEq, Repr, Fintype

def random031DispatchPrivatePure : Random031DispatchCompressor -> Bool
  | Random031DispatchCompressor.raw_count_shadow => false
  | _ => true

def random031DispatchRoutePure : Random031DispatchCompressor -> Bool
  | Random031DispatchCompressor.IQ_plus_R_terminal_spigot => true
  | Random031DispatchCompressor.terminal_certificate_ledger_plus_R => true
  | Random031DispatchCompressor.owner_filtration_plus_R => true
  | Random031DispatchCompressor.IQ_without_R_private_only => false
  | Random031DispatchCompressor.row_name_exception_list => true
  | Random031DispatchCompressor.raw_private_bit => false
  | Random031DispatchCompressor.all_colored_plus_IQ_route_reconstruction => false
  | Random031DispatchCompressor.raw_count_shadow => false

def random031DispatchRetainsR : Random031DispatchCompressor -> Bool
  | Random031DispatchCompressor.IQ_plus_R_terminal_spigot => true
  | Random031DispatchCompressor.terminal_certificate_ledger_plus_R => true
  | Random031DispatchCompressor.owner_filtration_plus_R => true
  | _ => false

def random031DispatchLegalNow : Random031DispatchCompressor -> Bool
  | Random031DispatchCompressor.all_colored_plus_IQ_route_reconstruction => false
  | Random031DispatchCompressor.raw_count_shadow => false
  | _ => true

theorem random031_IQ_private_only_not_route :
    random031DispatchPrivatePure
      Random031DispatchCompressor.IQ_without_R_private_only = true ∧
    random031DispatchRoutePure
      Random031DispatchCompressor.IQ_without_R_private_only = false ∧
    random031DispatchRetainsR
      Random031DispatchCompressor.IQ_without_R_private_only = false := by
  decide

theorem random031_terminal_spigot_keeps_route_sidecar :
    random031DispatchLegalNow
      Random031DispatchCompressor.IQ_plus_R_terminal_spigot = true ∧
    random031DispatchPrivatePure
      Random031DispatchCompressor.IQ_plus_R_terminal_spigot = true ∧
    random031DispatchRoutePure
      Random031DispatchCompressor.IQ_plus_R_terminal_spigot = true ∧
    random031DispatchRetainsR
      Random031DispatchCompressor.IQ_plus_R_terminal_spigot = true := by
  decide

theorem random031_all_colored_plus_IQ_still_not_route :
    random031DispatchLegalNow
      Random031DispatchCompressor.all_colored_plus_IQ_route_reconstruction =
        false ∧
    random031DispatchRoutePure
      Random031DispatchCompressor.all_colored_plus_IQ_route_reconstruction =
        false := by
  decide

/-! ## HYP-3528 one-red-clause audit -/

inductive Random031AuditClause where
  | ordinary_route_emit
  | free_hole_single_emit
  | free_hole_doublet_buffer_emit
  | bypass_transport_emit
  | bypass_bracket_lift_emit
  | private_firewall_route_sidecar
  | vertical_halfturn_guard
  | residual_pair_close_tail
  | raw_count_shadow
deriving DecidableEq, Repr, Fintype

def random031AuditPath : List Random031AuditClause :=
  [Random031AuditClause.ordinary_route_emit,
   Random031AuditClause.free_hole_single_emit,
   Random031AuditClause.free_hole_doublet_buffer_emit,
   Random031AuditClause.bypass_transport_emit,
   Random031AuditClause.bypass_bracket_lift_emit,
   Random031AuditClause.private_firewall_route_sidecar,
   Random031AuditClause.vertical_halfturn_guard,
   Random031AuditClause.residual_pair_close_tail,
   Random031AuditClause.raw_count_shadow]

def random031AuditScore : Random031AuditClause -> Int
  | Random031AuditClause.ordinary_route_emit => 400
  | Random031AuditClause.free_hole_single_emit => 370
  | Random031AuditClause.free_hole_doublet_buffer_emit => 359
  | Random031AuditClause.bypass_transport_emit => 293
  | Random031AuditClause.bypass_bracket_lift_emit => 274
  | Random031AuditClause.private_firewall_route_sidecar => 269
  | Random031AuditClause.vertical_halfturn_guard => 249
  | Random031AuditClause.residual_pair_close_tail => -7
  | Random031AuditClause.raw_count_shadow => -217

inductive Random031AuditStatus where
  | finite_closed
  | guardrail_closed
  | open_proof_debt
  | checksum_only
deriving DecidableEq, Repr, Fintype

def random031AuditStatus : Random031AuditClause -> Random031AuditStatus
  | Random031AuditClause.ordinary_route_emit =>
      Random031AuditStatus.finite_closed
  | Random031AuditClause.free_hole_single_emit =>
      Random031AuditStatus.finite_closed
  | Random031AuditClause.free_hole_doublet_buffer_emit =>
      Random031AuditStatus.finite_closed
  | Random031AuditClause.bypass_transport_emit =>
      Random031AuditStatus.finite_closed
  | Random031AuditClause.bypass_bracket_lift_emit =>
      Random031AuditStatus.finite_closed
  | Random031AuditClause.private_firewall_route_sidecar =>
      Random031AuditStatus.guardrail_closed
  | Random031AuditClause.vertical_halfturn_guard =>
      Random031AuditStatus.guardrail_closed
  | Random031AuditClause.residual_pair_close_tail =>
      Random031AuditStatus.open_proof_debt
  | Random031AuditClause.raw_count_shadow =>
      Random031AuditStatus.checksum_only

def random031AuditStatusCount (status : Random031AuditStatus) : Nat :=
  random031AuditPath.foldl
    (fun acc clause =>
      if random031AuditStatus clause = status then acc + 1 else acc)
    0

theorem random031_audit_path_length :
    random031AuditPath.length = 9 := by
  decide

theorem random031_audit_score_path :
    random031AuditPath.map random031AuditScore =
      [400, 370, 359, 293, 274, 269, 249, -7, -217] := by
  decide

theorem random031_audit_status_hist :
    random031AuditStatusCount Random031AuditStatus.finite_closed = 5 ∧
    random031AuditStatusCount Random031AuditStatus.guardrail_closed = 2 ∧
    random031AuditStatusCount Random031AuditStatus.open_proof_debt = 1 ∧
    random031AuditStatusCount Random031AuditStatus.checksum_only = 1 := by
  decide

def random031AuditOpenClauses : List Random031AuditClause :=
  random031AuditPath.filter
    (fun clause =>
      random031AuditStatus clause = Random031AuditStatus.open_proof_debt)

theorem random031_audit_one_red_clause :
    random031AuditOpenClauses =
      [Random031AuditClause.residual_pair_close_tail] := by
  decide

theorem random031_stream_certificate_readout :
    79 = 77 + 2 ∧
    random031TerminalCertificateCount
        Random031TerminalCertificate.ordinary_route = 64 ∧
    random031TerminalCertificateCount
        Random031TerminalCertificate.free_hole_single = 10 ∧
    random031TerminalCertificateCount
        Random031TerminalCertificate.free_hole_doublet = 2 ∧
    random031TerminalCertificateCount
        Random031TerminalCertificate.pure_bypass = 1 := by
  decide

theorem random031_one_red_tail_geometry :
    random031TailOwners Random031EmitterState.transport_emitted =
      [147, 169] ++ random031ResidualOwners ∧
    random031TailOwners Random031EmitterState.boundary_lift_emitted =
      random031ResidualOwners ∧
    random031ResidualOwners = [45, 173] := by
  decide

/-! ## Abstract guarded-emission lemma -/

/-- A visible fiber has constant target if every legal hidden completion in the
fiber has the same target value. -/
def TargetConstantOnFiber
    {Hidden Target : Type}
    (fiber : List Hidden)
    (targetOf : Hidden -> Target) : Prop :=
  ∀ h1, h1 ∈ fiber -> ∀ h2, h2 ∈ fiber -> targetOf h1 = targetOf h2

/-- A proof token may be emitted when the fiber is nonempty and target-constant. -/
structure GuardedEmission
    {Hidden Target : Type}
    (fiber : List Hidden)
    (targetOf : Hidden -> Target)
    (printed : Target) where
  witness : Hidden
  witness_mem : witness ∈ fiber
  printed_eq : printed = targetOf witness
  constant : TargetConstantOnFiber fiber targetOf

theorem guardedEmission_correct
    {Hidden Target : Type}
    {fiber : List Hidden}
    {targetOf : Hidden -> Target}
    {printed : Target}
    (emit : GuardedEmission fiber targetOf printed)
    {h : Hidden}
    (hh : h ∈ fiber) :
    targetOf h = printed := by
  rw [emit.printed_eq]
  exact emit.constant h hh emit.witness emit.witness_mem

theorem guardedEmission_unique
    {Hidden Target : Type}
    {fiber : List Hidden}
    {targetOf : Hidden -> Target}
    {printed printed' : Target}
    (emit : GuardedEmission fiber targetOf printed)
    (emit' : GuardedEmission fiber targetOf printed') :
    printed = printed' := by
  rw [emit.printed_eq, emit'.printed_eq]
  exact emit.constant emit.witness emit.witness_mem
    emit'.witness emit'.witness_mem

/-! ## Conditional terminal packet interface -/

/-- Producer obligations needed to turn the random031 terminal ledger into an
actual LRC14 discharge for the concrete speed row `v`. -/
structure Random031TerminalProducers (v : Fin 13 -> Int) where
  ordinaryRoute : Random031TerminalCertificate.ordinary_route =
      Random031TerminalCertificate.ordinary_route
  freeHoleBracket : Random031TerminalCertificate.free_hole_single =
      Random031TerminalCertificate.free_hole_single
  pureBypassOwnerBoundary : Random031TerminalCertificate.pure_bypass =
      Random031TerminalCertificate.pure_bypass
  privateFirewallCompatible :
      (hyp3490RowDetail hyp3490HardPrivateOverlapRow).route =
        PrivateFirewallRoute.random031_private_firewall_hard_overlap
  verticalHalfturnGuard : True
  mreachFloor : (1 : Real) / 14 ≤ Mreach v

theorem random031_terminalProducers_to_mreach
    {v : Fin 13 -> Int}
    (packet : Random031TerminalProducers v) :
    (1 : Real) / 14 ≤ Mreach v :=
  packet.mreachFloor

theorem random031_terminalProducers_to_lonely
    {v : Fin 13 -> Int}
    (hv : ∀ i, v i ≠ 0)
    (packet : Random031TerminalProducers v) :
    ∃ t : Real, Lonely 14 v t :=
  lonely_of_Mreach_ge v hv packet.mreachFloor

/-! ## Axiom audit hooks -/

#print axioms random031_cell_total
#print axioms random031_terminal_certificate_total
#print axioms random031_seam_decomposes_as_transport_boundary_residual
#print axioms colored_private_only_emits_status_not_route
#print axioms owner_filtration_ready_emits_residual_and_route
#print axioms random031_contract_status_hist
#print axioms random031_unique_open_tail_clause
#print axioms random031_route_bearing_clauses_retain_R
#print axioms random031_IQ_private_only_not_route
#print axioms random031_terminal_spigot_keeps_route_sidecar
#print axioms random031_audit_status_hist
#print axioms random031_audit_one_red_clause
#print axioms random031_stream_certificate_readout
#print axioms random031_one_red_tail_geometry
#print axioms guardedEmission_correct
#print axioms guardedEmission_unique

end LRC14
end LonelyRunner
