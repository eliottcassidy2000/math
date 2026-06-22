# LRC14 Doublet Witness-Floor Formalization, Codex S81

This session converted the claude-opus S4 doublet floor scout into a small Lean
import boundary rather than trying to reprove the scout in Lean.  The new module
`TournamentH7.LRCDoubletWitnessFloor` records the exact reported minima:

- `rho*` floor `2/147`;
- witness floor `1066/2205`;
- comparison floor `1/84` from the consecutive rho* route;
- comparison witness floor `14249/252252` from THM-530.

Lean now proves the finite rational comparisons:

- `2/147 > 0`;
- `1066/2205 > 0`;
- `2/147 > 1/84`;
- `1066/2205 > 14249/252252`;
- `1066/2205 > 2/147`.

The important rigor boundary is explicit.  This is not a proof that the S4
Python scout covers every genuine-wide doublet configuration; it is a formal
checksum for the final rational values reported by that scout.  Coverage still
lives in `04-computation/lrc14_rhostar_gw_doublet_claudeopus_0622s4.py` and
`05-knowledge/results/lrc14_rhostar_gw_doublet_claudeopus_0622s4.out`.

## Formalization Status

The checksum is now wired into `LRCFourteenSkeleton.lean` as a finite DAG node
and into `TournamentH7.Verify` as an audit wrapper.  While running aggregate
Verify, a separate import gap surfaced: `Verify.lean` already referenced
`LRCQ6Contraction` wrappers but did not import `TournamentH7.LRCQ6Contraction`.
That import is now present, so the q6 wrappers and the new doublet wrapper both
resolve in the aggregate audit.

After pulling concurrent S24 work, the same pass also absorbed
`LRCGk8SingleFar.lean`.  I added a packaged theorem
`Gk8SingleFar.all_binding_checks` and wired it into the root import, the
skeleton, and `Verify`.  This makes the HYP-2829 far-count split Lean-visible as
an arithmetic checksum: bounded and single-far `L_yK8<10cap`, with single-far
below bounded, for k=8,9,10.

Verified targets:

- `lake build TournamentH7.LRCGk8SingleFar`;
- `lake build TournamentH7.LRCDoubletWitnessFloor`;
- `lake build TournamentH7.LRCFourteenSkeleton`;
- `lake build TournamentH7.Verify`;
- `lake build TournamentH7`.

The skeleton still has the intended direct `sorry` obligations:
`lonely_of_Mreach_ge`, `thm527_partA_density_pos_implies_reach`,
`thm527_partG_uniform_floor`, `doublet_Rtail_uniform_bound`, and
`gK8_concentration_extremality`.  The top-level `lrc14_from_thm527` depends
transitively on those open obligations.

## Tournament Analysis

The useful tournament vertices are not runners or arcs.  For this formalization
step, the vertices are proof-DAG/import-boundary objects:

- finite Python scout coverage;
- Lean rational checksum;
- skeleton DAG node;
- aggregate Verify wrapper;
- THM-527 witness-floor compactness bridge;
- LRC14 assembly.

The pairwise observable is predicate fidelity: which object preserves more of
the statement needed for LRC14 without importing unproved coverage.  The switch
orients an edge toward the object with fewer informal assumptions and a clearer
audit boundary.  The resulting path is transitive:

`finite scout coverage > Lean rational checksum > skeleton DAG node > aggregate Verify wrapper > compactness bridge > LRC14 assembly`.

This quotient preserves the distinction between arithmetic formalization and
configuration coverage.  It destroys row geometry, base/gap data, and the exact
enumeration logic, so it cannot replace the finite scout or close the analytic
floor by itself.

## Next Formal Targets

The most useful next Lean-facing step is not another rational checksum.  It is a
typed statement that separates coverage certificates from arithmetic readouts:
a structure carrying a finite-scout theorem, exact floor constants, and the
readout theorem that exposes those constants to the skeleton.  That would let
future exact-rational scouts plug into the same LRC14 DAG without blurring which
part is computational coverage and which part is arithmetic verification.
