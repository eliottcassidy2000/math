# LRC14 q-Cusp, Hurwitz, and Branch Closure Synthesis

Source: codex-2026-06-27-S253.

This pass expands HYP-3083 after reading the incoming q-Pochhammer,
Hurwitz/Markov/Pell, Lean q-cusp, median-hull, sixth-power, Robbins, and
branch-kernel work.  The useful synthesis is not a new analogy.  It is a
changed shape for the remaining LRC14 proof.

## The Shared Rule

The q-Pochhammer and modular-function cue says:

```text
full modular function + meromorphic at the cusp
  -> q-expansion has finitely many negative powers
  -> infinite positive tail is legal only after finite principal part is named
```

The Hurwitz cue has two relevant forms.

In the Diophantine Hurwitz-Markov-Pell lane, the visible extremal constant or
rare scalar is not the proof object.  The proof object is the finite address:
continued-fraction period, Markov depth, Pell unit, endpoint shell, and carry
residue.

In the complex-analytic Hurwitz lane, q-Pochhammer finite products are
nonvanishing on `|q|<1`, and zero/pole persistence prevents a limit argument
from creating hidden interior zeros unless the limiting object has the
corresponding divisor data.  So product tails also need a zero/pole ledger.

Robbins then gives the graph version:

```text
forgotten finite address
  -> naked dependency bridge
  -> illegal proof quotient
```

HYP-3082 verifies this on the bounded HYP-2963 bank: the raw scalar star has
five naked bridges; the protected graph has none after sidecars are attached.

Incoming S252b sharpens the shape of the source class.  After the q-witness
gate, THM-571 apex-majority gamma descent, and HYP-2906 one-large interval
peeler, the live residual is not arbitrary covering rows.  It is low-apex,
top-balanced covering rows with `1..6` multiples of `14` and no legal quotient
that forgets finite address.

## The Improved Argument

The LRC14 proof should now be attacked as a four-map theorem:

```text
low-apex top-balanced covering row
  -> finite-address HYP-2963 packet
  -> protected S250 branch graph
  -> terminal discharge or named residual debt
  -> formal witness readout M >= 1/14
```

This reconciles the older Lean-status document, THM-523, HYP-2963, S59, and
the new q-cusp/branch work.

THM-523 says the proof-critical rows are covering rows: if some small
denominator `q<=14` is missing, the direct `t=1/q` witness proves the bound.
S59 says the AP/GW census is the easy-side extremal, not the critical path.
S252b adds that THM-571 and HYP-2906 have already peeled off apex-majority and
one-large-speed cases, leaving the low-apex/top-balanced covering residual.
HYP-2963 says the bounded packet bank routes everything into q-witness, AP/GW,
C27 petal, K33 state lift, or covering moment, with no unknown bucket.  HYP-3082
says those routes form a protected graph once sidecars are attached.

So the remaining proof is not "find a scalar stronger than `M`, `q`, route, or
automatic word."  The remaining proof is:

```text
make the packet emitter global,
then discharge the covering and K33 exits,
while proving no family tail creates a new naked bridge.
```

## What Remains

1. **Global packet theorem.**  HYP-2963 is a bounded audit.  We need a theorem
   that every low-apex/top-balanced primitive covering residual emits the full
   finite-address packet: exact scale, q-threshold, endpoint owners, topology,
   C27/K33 labels, q-cusp principal part, arithmetic seed if present, and
   destroyed-coordinate exit.

2. **Covering-family proof.**  This is the real analytic bottleneck.  It
   should be phrased as a finite principal-part ledger: the polar debt is the
   set of killed small-denominator witnesses; the legal tail is p0/moment,
   gamma/apex-periodic, Haar/Ramanujan, or Node-3 equidistribution data; the
   exit is positive witness density feeding Part A.

3. **K33 state lift.**  THM-572 is formalized only after the lift exists.  The
   missing theorem is the finite-address K33 packet construction producing the
   tournament state with `H(T)=7`.

4. **Family no-bridge theorem.**  S250 is bounded.  A parameterized covering,
   K33, or q-cusp tail must preserve a reverse verification path, descend, or
   name new residual debt.  Otherwise it creates a Robbins bridge.

5. **Lean carrier record.**  The formal side should combine
   `LRCModularCuspLedger`, `LRCMedianCenterControl`, `LRCMomentDual`, and
   `LRCTournamentStateLift` through one finite-address branch-packet record.

## Concrete Next Pull

The best next computation is not another analogy scout.  It is to instantiate
one actual HYP-2963 covering-moment row as a q-cusp ledger:

```text
q_cusp_ledger_id
principal_part_order
polar_exit_word
partition_tail_certificate
log_derivative_divisor_channel
branch_id
bridge_status
terminal_covering_discharge
```

Then run the S250 bridge audit with that row's new sidecars visible.  If the
row still lands in the protected graph and the covering tail has a finite
principal part, the q-Pochhammer idea has become proof infrastructure.  If it
does not, the failed bridge names the exact coordinate the proof still lacks.

## Tournament Analysis

Vertices are proof obligations and finite-address carriers, not runners or
raw q-coefficients.

Pairwise observable: retained LRC predicate, finite cusp/arithmetic address,
destroyed-coordinate declaration, branch reverse-verification path, terminal
discharge strength, and formal payload readiness.

Switch: orient toward the carrier that retains more of those fields and loses
fewer coordinates.

Tie path:

```text
global_finite_address_packet_theorem
> protected_branch_graph_entry
> covering_family_p0_moment_discharge
> K33_THM572_state_lift_discharge
> q_cusp_principal_part_ledger
> Hurwitz_Markov_Pell_seed_address
> route_state_sidecar_closure
> Lean_formal_packet_record
> named_F7_residual_debt
> raw_scalar_or_tail_shadow
```

Assumption challenge: runners, gaps, fixed circle sections, section
boundaries, wall crossings, residues, cover arcs, Fourier modes, q-coefficients,
cusp principal parts, Hurwitz quadruples, branch kernels, endpoint owners, and
proof obligations were all considered.  Proof obligations and finite-address
sidecars are the useful vertices because they preserve exact scale,
open/boundary status, endpoint ownership, q-threshold, polar debt, bridge
status, and terminal exit.  Raw runner identity and raw q-coefficients are
destroyed unless another sidecar pays for them.
