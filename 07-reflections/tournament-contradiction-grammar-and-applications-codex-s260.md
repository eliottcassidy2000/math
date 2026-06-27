# Tournament-Contradiction Grammar And Applications

codex-2026-06-27-S260.

## Seed being generalized

Past work used a sharp tournament contradiction pattern: construct a
subproblem whose tournament analogue has `H=7` or `H=21`, then invoke the
permanent-gap theorems to rule it out.  That is still one of the best
terminal moves in the repo, but the current LRC14 frontier shows why it cannot
be the first move.  The constructed object must first be a complete tournament
shadow of a proof-relevant predicate; loose directed graphs and coarse winding
charts can imitate pieces of the invariant without inheriting THM-200/THM-202.

The incoming S31ah engine made this reusable by computing score sequences,
Landau checks, odd-cycle conflict graphs, independence polynomials, `H`, SCCs,
and Hamiltonian-path parity.  The later S31ah spectrum pass sharpens the
single-component story: `H=7` is `K_3` nonrealizability in Omega and `H=21` is
the hard `K_10` case.  Concurrent mac-mini HYP-3099 adds three concrete
applications: cap-optimality is a non-transitive improvement tournament but a
bounded finite check; baby-Hodge holes are c5/spectral rather than forbidden-H
holes; and `apex-7 = H=7` is a coincidence, with order-2 antipodal symmetry as
the genuine tournament-native lever.  This session adds a second layer:

```text
certificate grammar + target ledger -> scored proof applications
```

The new artifact is
`04-computation/tournament_contradiction_grammar_codex_s260.py`, with stored
output in
`05-knowledge/results/tournament_contradiction_grammar_codex_s260.out`.
After rebasing over the new Lean work, the formal sink for this grammar is
`TournamentH7.LRCBleedingEdgeFrontier`: tournament-certificate legality should
become packet sidecar data in that wrapper, not an isolated prose checklist.

## Programmatic readout

The script generates `24` tournament proof moves across `12` live targets.  It
keeps H=7/H=21 as a terminal certificate, but also treats tournaments as:

- score polytopes via Landau majorization,
- Omega-realizability filters,
- alpha-vector / independence-polynomial gates,
- odd-cycle census spectra,
- automorphism-parity alarms,
- skew-spectrum and transfer-matrix sidecars,
- SCC product descent,
- trienerment/fine-scale repairs for ties,
- proof-obligation dominance tournaments,
- median-route and no-naked-bridge guards,
- metagraph nonembedding tests,
- path-homology and matrix-observability sidecars,
- residue-sieve auditors for fixed-gap prime questions,
- proof-by-survival necessary-condition miners.

Exact checks inherited from KPS and rebased incoming work:

```text
H=7 rejected.
H=21 rejected.
single-component H gaps package as Omega K_3 and K_10 nonrealizability.
even H rejected.
no even Hamiltonian path counts through n<=5.
small exact H spectrum through n<=5: 1,3,5,9,11,13,15.
```

The target ranking matched the current proof shape:

- `Trienerment lift for ties` dominates coarse mod-14 winding.
- `H forbidden-value`, `SCC product descent`, and `Omega forbidden-shape miner`
  dominate H=21 closure, now focused around the `K_10` Omega package and the
  remaining `I(Omega,2)=21` component exclusions.
- `H forbidden-value` plus `SCC product descent` dominate level-7 state-lift
  only after complete-tournament validity is checked.
- `Automaton-state`, `Matrix observability`, and `Median/Helly` dominate the
  HYP-2963 normalizer.
- `Residue-sieve tournament` dominates the sexy-prime sidecar, with explicit
  parity/distribution debt.

## The non-linear frontier

The top eight generated techniques form a small tournament with one directed
3-cycle and `3` Hamiltonian paths.  The key SCC is:

```text
Automaton-state tournament
H forbidden-value certificate
No-naked-bridge orientation
```

This is exactly the frontier tension.  An H contradiction is often the desired
terminal state, but an automaton/fiber normalizer may be needed to know what
the tournament vertices are, and bridge protection may be needed before a
state-lift quotient can discard a route.  The proof order is not linear.

## LRC14 synthesis

For LRC14, tournament certificates should enter after the observer-gluing
packet has named:

```text
source_row_id
crt_c7_lift_status
crt_c2_dyadic_lift_status
direct_lonely_measure
direct_component_count
largest_direct_arc
denominator_net_threshold_D
pascal_pair_mass_unit
triangular_cap_shadow
cap_defect
sector_pair_scissors_signature
grid_class
active_binder_owner_word
endpoint_owner_transition_word
overlap_failure_chart
terminal_exit_or_named_debt
```

The incoming `TournamentH7.LRCBleedingEdgeFrontier` module already packages
observer charts, equivalence shadows, polynomial witness data, Pascal pair
mass, and moment-degree status around the finite-address branch packet.  S260's
next formal role is to add tournament-certificate legality columns to that
wrapper: complete-tournament status, tie-aware lift status, SCC/Omega/H
certificate status, bridge-protection status, and no-hit survivor profile.

Incoming S64/THM-577 improves the symbolic cap chart by proving the value side
for the `j=3` apex-overlap case.  The grammar says how to use that progress:
the symbolic cap chart is now stronger, but it still has to glue to the
normalized arc chart, level-7/CRT chart, branch/K33 chart, and formal witness
chart.  Incoming HYP-3099 adds that cap optimality itself may be a bounded
finite check obstructed by non-transitive improvement cycles; that is a
diagnostic target for the observer-gluing ledger, not a reason to collapse
back to a scalar.  A future H=7/H=21 hit is decisive only at that terminal
interface.

## Sexy-prime relation

The relation to proving the sexy prime conjecture is proof hygiene, not direct
strength.  The same tournament discipline says how to build the fixed-gap
sidecar:

```text
midpoint m
gap ray h=3
local bad residues m=+-h mod q
collapse modulo 3
Hardy-Littlewood weight channel
prime-power exceptions
parity debt
distribution debt
terminal exit
```

This can prevent illegal local-to-global forgetting, and it can compare local
channels as a residue-sieve tournament.  It does not prove infinitely many
gap-6 prime pairs because parity-breaking and distribution in progressions are
external analytic obligations.

## Assumption challenge

I explicitly considered tournament vertices as runners, gaps, fixed circle
sections, section boundaries, wall crossings, residues, cover arcs, Fourier
modes, matroid circuits, Omega components, score sequences, proof obligations,
sidecar columns, quotient maps, finite sieve channels, automaton states,
rooted perspectives, and branch-graph edges.

For this pass I chose certificate functors and target proof obligations.

Preserved predicate: whether a pulled-back tournament certificate can force a
contradiction, repair a degenerate encoding, or record a necessary sidecar for
LRC14/sexy-prime ledgers.

Destroyed information: concrete runner geometry, exact metric scale, and
distributional number-theory strength unless the target sidecars explicitly
retain them.

The challenged assumption is that tournaments are only H-value contradiction
machines.  The better reading is broader: tournaments can be contradiction
engines, legal-forgetting auditors, route-priority graphs, sidecar matrices,
bridge-protection tests, and no-hit necessary-condition miners.

## Next frontier tasks

1. Extend the S31ah engine with an Omega-realizability miner for the remaining
   `I(Omega,2)=21` candidate components, using `K_10` as the named
   single-component target and P4 as the known excluded case.
2. Run a fine mod-7 winding scout with score, cycle census, Paley distance,
   skew spectrum, and tie-lift status.
3. Add certificate columns to HYP-2963 packet rows and
   `TournamentH7.LRCBleedingEdgeFrontier`, then measure which ones separate
   live route/status fibers.
4. Apply completeness/SCC/H-spectrum checks to any THM-572/K33/F7 state-lift
   terminal outputs.
5. Keep the sexy-prime tournament strictly as a residue-sieve sidecar until a
   real parity/distribution theorem is supplied.
