# LRC14 proof-circuit past-work synthesis

This search pass changes the circuit-complexity prompt from analogy into a
proof-engine question:

```text
What is the fixed Boolean/proof circuit whose inputs are legal LRC sidecars,
and what exact input is missing when a proposed shortcut fails?
```

The older work is not peripheral.  It supplies gates.

This note extends the upstream HYP-3116 missing-input ledger by filling in the
past-work gate map and a first executable compiler sketch.

After rebase, HYP-3116 has a concrete base circuit: `12` inputs, `3` gates,
depth `3`, all inputs essential, and `0/10` tempting shortcuts closed.  Its
missing-input frequency is the practical target for this note:
`finite_address=10`, `observer_gluing=8`, `endpoint_owner=7`,
`uniformity=5`, `root_ear_sidecar=3`, `cocycle_exactness=2`, and
`state_lift=1`.  The past-work map below is a repair menu for those missing
coordinates.

## Map 1: Legacy Gates

HYP-2108 is the endpoint-cover circuit.  Its real content is simultaneous
resonance: a tight residual must place all component midpoints near compatible
integer centers at once.  That is a circuit predicate, not a bag of local
inequalities.

HYP-2112 is the output wire.  `Phi(C)` is the exact gap functional, and
`ker Phi` is the tight/worry set.  The compiler should aim to prove `Phi>0` or
name why a row remains in `ker Phi`.

HYP-2114/HYP-2115 are the fold gates.  Three-term folds are the load-bearing
objects; four-term energy is a translated shadow unless the hidden virtual sum
node is retained.  This looks like the same controlled-forgetting rule that
later appears as cocycles and observer sidecars.

HYP-2961/HYP-2963 supply the decision tree: every strict counterexample has to
fall through `qdiv`, open-front, unit-petal, K33, wide/decorrelation,
apex-multiple, covering, or source-kernel branches.  Read as circuit
complexity, those are input gates with a fixed order and named exits.

HYP-3016 supplies the branching-program warning: finite automaton states are
cheap but mix boundary and open rows unless a magnitude cocycle is retained.

HYP-3102 supplies the missing-input detector: every illegal quotient emits a
first obstruction class.

HYP-3107 supplies the proof bus: finite-address and observer-gluing packets are
the legal terminal outputs.

HYP-3109/HYP-3112 and HYP-3111/HYP-3115 supply newer sidecars: root curves,
ear payloads, relation walls, and the fixed input/depth/fan-in discipline.

## Map 2: Compiler Picture

The productive picture is:

```text
residual row
  -> labelled packet
  -> endpoint-cover P gate
  -> Phi exact gap wire
  -> fold/virtual-sum gate
  -> automaton magnitude-cocycle guard
  -> Lee-Yang root/ear sidecar
  -> relation-wall sidecar
  -> first-obstruction cocycle
  -> Lean finite-address or observer-gluing exit
```

This is not a demand that every proof literally execute in this sequence.  It
is a guardrail for shortcuts: a shortcut is legal only if it states which of
these wires it preserves, discharges, or intentionally bypasses with a stronger
certificate.

## Signals To Add

The next residual-row table should add:

```text
proof_circuit_input_basis_id
proof_circuit_missing_input_vector
uniformity_guard_status
endpoint_cover_P_gate
Phi_gap_output_wire
fold_gate_depth
hidden_virtual_sum_count
automaton_fiber_mixing_bit
magnitude_cocycle_height
first_obstruction_class
Lee_Yang_ear_payload_mean_level
root_motion_reconstruction_status
relation_wall_class
sidecar_fanin_profile
minimal_certificate_depth
gate_route_purity
terminal_exit_kind
```

The key derived object is `proof_circuit_missing_input_vector`.  A failed proof
attempt should not say "the scalar was inconclusive"; it should say which
required sidecar was absent.

## Wild Guesses Worth Testing

First, the same structure may be recurring under different names:

```text
3-term fold
hidden virtual sum
Lee-Yang one-runner ear payload
relation low-height wall
automaton magnitude cocycle
first obstruction class
```

These may all be "retained intermediate sum nodes" in different quotients.  If
so, the proof should try to define one common retained-node functor and show
that each existing sidecar is a projection of it.

Second, extremality may be less about the largest scalar and more about the
smallest legal proof circuit to a terminal exit.  The consecutive/AP-like rows
may be extremal because they compress into a shallow certificate only after
all hidden nodes are retained; nearby rows either open a positive `Phi` gap or
expose an obstruction earlier.

Third, the finite-bank classifier `apex7_error<=5` should be treated as a
probe for an input wire.  Its theorem-grade replacement might be:

```text
apex7_error small
AND root_motion_reconstructible
AND hidden_fold_nodes_retained
AND first_obstruction_class generated
=> finite-address or observer-gluing exit
```

That is a circuit statement.  The one-literal threshold is only the shadow.

## Tournament Check

The scout `lrc14_proof_circuit_past_work_scout_codex_20260627.py` uses proof
gates and sidecars as vertices.  It produces an acyclic fingerprint:

```text
score_hist={-3:1,18:1,19:2,20:1,21:1,22:1,24:1,25:1,28:1,31:1}
directed_3cycles=0
scc_sizes=eleven singleton components
hamiltonian_path_count=1
```

The priority path starts with `first_obstruction_cocycle_gate` and
`lean_frontier_obligation_bus`, then routes through the labelled packet
decision tree and `Phi`.  That ordering is a useful correction to the urge to
search for a small scalar circuit first.  The first thing to measure is what
information a proof map forgets.

Related: HYP-3117, HYP-3116, HYP-3115, HYP-3111, HYP-3112, HYP-3108, HYP-3107,
HYP-3102, HYP-3016, HYP-2963, HYP-2961, HYP-2114, HYP-2112, HYP-2108,
OPEN-Q-108.
