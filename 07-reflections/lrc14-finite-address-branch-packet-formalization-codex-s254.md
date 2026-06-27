# LRC14 Finite-Address Branch-Packet Formalization

Source: codex-2026-06-27-S254.

This pass formalizes the current HYP-3083 frontier as a Lean interface in
`TournamentH7.LRCFiniteAddressBranchClosure`.  It does not prove LRC14, the
global packet normalizer, the covering-family bound, or the K33 state lift.
It does pin down what those missing statements must deliver if they are to fit
the cutting-edge proof spine.

The resulting Lean theorem shape is:

```text
CuttingEdgeBranchCoverage
  = every nonzero 13-speed row is discharged by
    q-witness / level-7 lift sieve / one-large-speed peeler
    or emits a low-apex top-balanced finite-address packet

CuttingEdgeBranchCoverage -> LRC14Statement
```

After rebasing over incoming S31af/S60, the early apex gate should be read in
its sharpened form: THM-573 closes rows with at least seven multiples of `7`,
subsuming the older THM-571 "at least seven multiples of 14" statement.  The
packet therefore retains both the covering datum "contains a multiple of 14"
and the new residual datum "`1..6` multiples of 7".

The packet itself joins the existing formal tracks:

- `LRCModularCuspLedger`: finite q-cusp principal part and q-Pochhammer tail.
- `LRCMedianCenterControl`: route-center packet and witness-floor payload.
- `LRCMomentDual`: covering-moment dual feasibility gives `p0 <= Ly`.
- `LRCTournamentStateLift`: K33/F7 exits must eventually build the state lift
  that THM-572 already closes.
- `LRCMreachConcrete`: the terminal floor feeds the concrete compactness
  theorem `Mreach >= 1/14 -> lonely time`.

## Small Investigations

### 1. The q-cusp sidecar can be made proof-bearing now

The new module packages the S247 `j` principal part into a `QCuspSidecar`.
The sidecar contains:

```text
principalPart
productTail
polarExitWord
hurwitzZeroPersistenceStatus
illegalInfinitePolarTailFlag = false
finitePrincipalPart
```

The checked lemma `qCuspSidecar_finite_principal_part` has no extra axioms.
This is exactly the modular-cusp lesson in formal form: a q-tail is allowed
only after the finite polar part is carried as data.  The Lean object is still
toy-sized, but the field list is the correct interface for a real HYP-2963 row.

### 2. The protected-branch condition is a field, not a slogan

The `ProtectedBranchCertificate` record retains raw and protected bridge
statuses plus the destroyed-coordinate sidecars that protect a quotient.  The
checked lemma `protectedBranch_no_raw_naked_bridge` is tautological in the
right way: any later graph audit must populate the field

```text
protectedBridgeStatus != raw_naked_bridge
```

before a contraction can enter the proof.  This turns the S250 audit rule into
a formal packet requirement.

### 3. The covering-moment bridge is already reusable

The new `CoveringMomentDualLedger` stores a dual polynomial `g` and exactly the
feasibility inequality needed by `LRCMomentDual.p0_le_Ly`.  The theorem
`coveringMomentDualLedger_p0_le_Ly` builds:

```text
slowμ(coverSet E).toReal <= Ly E g
```

from that ledger.  This means a future covering-family proof can focus on
producing exact feasible `g` rows and the remaining `Ly <= cap - delta`
certificate; the pointwise p0-to-moment reduction is already formalized.

### 4. The terminal floor is the narrow formal throat

`TerminalDischargeCertificate` carries:

```text
exit
witnessFloor
1/14 <= witnessFloor
witnessFloor <= Mreach v
```

The theorem `terminalDischarge_mreach` is just transitivity.  This looks
small, but it is the right formal throat: q-witness, gamma descent,
one-large-speed peeling, covering moment, C27 petal, AP/GW boundary, and
K33/THM-572 all have to end by producing this floor or a contradiction that
removes the row.

### 5. The hyperoperation grid is a scheduler into the packet

Incoming HYP-3087 fits cleanly if it is treated as an address scheduler, not as
another scalar proof.  A space-filling path through the `x+2` / `x*2` grid is
legal only when each visited operation cell retains:

```text
(p,q)
operation lane
danger deficit
endpoint owner
level-7 lift status
destroyed coordinate
terminal exit
```

Those are exactly the kind of fields the Lean packet expects.  The creative
synthesis is therefore not "prove by walking the grid"; it is "use the grid to
choose the next finite address, then require the packet to prove that no LRC
clock data was forgotten."

## What The Formalization Clarified

The remaining proof is not one missing scalar.  The formal interface says there
are exactly three unsolved producer theorems:

1. **Normalizer/coverage.** Every row after early gates must emit a
   `FiniteAddressBranchPacket`.
2. **Discharge.** Each terminal exit must produce the witness-floor fields,
   or, for K33/F7, the tournament state lift that contradicts THM-572.
3. **Family branch closure.** Parameterized tails must preserve the
   no-naked-bridge field, descend, or name new residual debt.

Incoming HYP-3085-gK8 adds a sharper producer target for item 2: the
covering-moment exit should likely be a low-order moment certificate, led by
the pairwise sector co-emptiness `S2` term and ultimately reduced to an exact
reflection-symmetric `3x3` Perron/eigenvalue bound rather than a full AP/GW
census.  Incoming HYP-3085-covering/K33 adds concrete O2/O3 shuttle rows:
nested-refinement rows should feed covering-moment discharge, while
cross-handoff rows should feed the THM-572 state-lift producer.  In packet
terms, both belong inside `CoveringMomentDualLedger`, protected-branch data,
and the terminal floor.

Incoming HYP-3087 adds the adjacent producer coordinate for item 1: the
hyperoperation grid should feed the packet normalizer by selecting
danger-weighted operation cells, but only with root `(p,q)`, level-7 status,
destroyed-coordinate label, and terminal exit intact.  This makes the grid a
front end for `FiniteAddressBranchPacket`, not a replacement for the packet.

Incoming S31ag/HYP-3088-HYP-3089 adds a compatible producer target for item 2:
the polynomial-method paper's Conjecture 7.1 at `k=13` can be read as a
uniform largest-lonely-arc floor.  In Lean packet terms, that is exactly a
candidate source of the `TerminalDischargeCertificate` inequalities
`1/14 <= witnessFloor <= Mreach`.  It does not remove the need for packet
coverage, but it clarifies that a largest-arc certificate and a moment-dual
certificate should land in the same formal throat.

The q-Pochhammer/Hurwitz idea is useful because it gives the same schema for
all three: infinite tails are legal only when generated by a finite address.
For modular functions this is the finite negative q-power principal part; for
Hurwitz/Markov/Pell it is continued-fraction/Markov/Pell seed data; for proof
graphs it is the sidecar that prevents a bridge from becoming naked.

## Assumption Challenge

Candidate vertices considered: runners, residues, multiples-of-14 blocks,
gaps, fixed circle sections, wall crossings, cover arcs, q-coefficients,
principal parts, Hurwitz seeds, endpoint owners, branch nodes, median centers,
moment-dual ledgers, K33 lifts, and proof obligations.

The Lean module chooses proof obligations and finite-address sidecars as
vertices.  This preserves the LRC predicate, exact scale, multiple-of-14
status, level-7 residual status, q-cusp finite polar debt, branch protection,
and terminal witness readout.  It destroys raw runner identity, raw
q-coefficients, and raw operation-grid traversal order only after they are
reconstructed, irrelevant, or named as residual debt.

Tournament path:

```text
global_packet_normalizer
> protected_branch_graph
> covering_moment_exit
> k33_state_lift_exit
> q_cusp_principal_part_guard
> hurwitz_seed_guard
> median_center_scheduler
> level7_lift_sieve_gate
> one_large_speed_gate
> q_witness_gate
> lean_sidecar_closure
> raw_scalar_shadow
```

Lean verifies this carrier set has `12` vertices.  The path is a work-order:
the first genuinely open theorem is still global packet coverage for the
primitive covering residual with `1..6` multiples of `7`, followed by the
HYP-3085 covering-moment/Perron certificate, the HYP-3085 covering/K33 shuttle,
and the K33 state-lift construction.

## Build Result

Target verified:

```bash
cd 04-computation/lean/TournamentH7
lake build TournamentH7.LRCFiniteAddressBranchClosure
```

The target builds.  Axiom output for the new frontier is as expected:
sidecar projection lemmas are axiom-free or foundational, while the final LRC14
conditional assembly depends only on the existing concrete `Mreach` bridge
foundation (`propext`, `Classical.choice`, `Quot.sound`), not on new analytic
assertions.
