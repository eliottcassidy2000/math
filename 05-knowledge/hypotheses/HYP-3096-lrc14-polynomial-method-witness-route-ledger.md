---
id: HYP-3096
title: LRC14 polynomial-method witness-route ledger
status: SYNTHESIS / proof-obligation ledger; not a proof
source: codex-2026-06-27-S255
tangent: T1176
technique: LTI-237
tournament_technique: LTT-135
related:
  - HYP-3088
  - HYP-3089
  - HYP-3084
  - HYP-3085
  - HYP-3083
  - HYP-3003
  - HYP-3004
  - HYP-2866
  - HYP-2827
  - THM-573
  - THM-565
  - THM-530
  - OPEN-Q-108
  - arXiv:2604.23906
---

# HYP-3096: LRC14 Polynomial-Method Witness-Route Ledger

## Claim

The S31ag paper bridge should be treated as a proof-obligation ledger, not as
one more analogy.  For 13 speeds / 14 runners, the polynomial-method paper's
composite `k+1` obstruction and the project's LRC14 residual are the same
object:

```text
k+1 = 14 = 2 * 7
paper fallback: prime-factor lifts at c=2 and c=7
project route:  descent 14 -> 7 -> 2
```

The project already has the `c=7` side as THM-573, the level-7 lift sieve.  The
live residual is therefore the `c=2`/analytic side: primitive covering rows
with at most `6` multiples of `7`, after all q-witness and level-7 exits have
fired.

What remains is to replace the paper's prohibitive enumeration of
`I(k,p,1)` by a uniform witness-interval theorem:

```text
For every non-tight primitive covering 13-tuple S in the THM-573 residual,
prove

  mu(L(S)) >= m0 > 0
  number_of_components(L(S)) <= A0

where L(S) = {t : ||s_i t|| >= 1/14 for all i}.

Then max_arc(L(S)) >= m0 / A0, so every denominator d >= ceil(A0/m0)
has a grid witness a/d in L(S).
```

This is the project-strength version of Conjecture 7.1(13).  It is stronger
than the scalar statement `M(S) >= 1/14`, and it implies both the LR property
and the paper's large-denominator witness condition.  The reverse implication
should not be used without a compactness/finite-packet normalizer: LRC itself
only gives a real witness, while Conjecture 7.1 asks for a uniform rational
net after tight rows are removed.

## Ledger Form

Every surviving packet should carry these fields before it is allowed to
claim a terminal proof exit:

```text
source_tuple_or_family
primitive_gcd_status
covering_status
fourteen_free_q_witness_status
count_7_divisible
crt_c7_lift_status
crt_c2_dyadic_lift_status
level7_residual_leq6_status
exact_M_or_lower_bound
tight_boundary_status
non_tight_status
I_discrete_grid_status
continuous_bad_phase_measure
lonely_measure_floor
direct_1_14_component_bound
largest_lonely_arc_floor
denominator_net_threshold_D
hyperoperation_root_pq
sum_lane_p_plus_q
product_lane_p_times_q
power_stress_word
finite_principal_part_or_finite_bad_denominator_budget
preserved_lrc_predicate
destroyed_coordinate
required_sidecar_or_debt
terminal_exit
```

The new fields connect four active threads.

1. **Polynomial method / CRT:** `crt_c7_lift_status` and
   `crt_c2_dyadic_lift_status` record the prime-factor fallback that replaces
   the failed field interpolation over `Z/14`.

2. **Witness route:** `lonely_measure_floor`, `direct_1_14_component_bound`,
   and `largest_lonely_arc_floor` are the route from THM-530/THM-565 to
   Conjecture 7.1(13).

3. **Hyperoperation grid:** `hyperoperation_root_pq`, `p+q`, `p*q`, and the
   power-stress word keep the operation-address chart from T1169/LTI-233
   attached.  The grid is a scheduler only; the proof clock is still the
   danger deficit and terminal exit.

4. **q-cusp finite-principal-part discipline:** a legal rational-grid proof
   may leave only finitely many small-denominator exceptions, just as a
   full-level modular function meromorphic at the cusp has only finitely many
   negative q-powers.  Infinite bad-denominator tails are illegal unless they
   are generated, bounded, or routed to named residual debt.

## Relation To `I(k,p,1)`

The paper's `I(k,p,1)` is a finite-field/grid count of tuples with no witness
on the denominator-`p` grid.  In the LRC14 witness language, for a fixed tuple
this asks whether

```text
L(S) intersects (1/p) Z / Z.
```

A largest-arc lower bound bypasses the enumeration.  If
`max_arc(L(S)) >= ell0`, then every denominator
`p >= ceil(1/ell0)` has a witness, so the corresponding improper-grid count
is empty after that threshold for the packet family.  The remaining
finite-principal-part debt is the small range `p < ceil(1/ell0)`, which can be
handled by exact packet normalizers or named computation.  Thus the analytic
substitute for `I(13,p,1)` is not just positive total measure; it is positive
measure plus a uniform component bound for the direct `1/14` lonely set.

## S258 Starter Ledger: Direct Arcs Plus Pair Scissors

`04-computation/lrc14_observer_gluing_ledger_codex_s258.py` builds the first
small observer-gluing ledger over representative rows and writes
`05-knowledge/results/lrc14_observer_gluing_ledger_codex_s258.out`.  It
alternates this HYP's witness-route fields with HYP-3097's pair/Pascal scissors
fields on the same rows.

The exact direct `1/14` lonely-set scan gives three immediate checks.

```text
q-witness row {1,...,13}:       boundary-only direct grid exit, components=0
synthetic H7=7 row:             THM-573 terminal exit, largest arc=3/343
H7=6 boundary residual row:     live residual, components=42, largest=19/1372
apex family V=13:               live residual, components=24, largest=3/637
apex family V=200:              live residual, components=102, largest=3/9800
divisor-loaded B=8 row:         live residual, components=860, largest=1/82320
```

The back-and-forth conclusion is useful: direct arcs are computable packet
fields, but they are not yet a scalar theorem.  The divisor-loaded rows keep
positive direct measure while shattering the largest direct arc, exactly the
THM-575 warning.  Therefore the HYP-3096 target should be sharpened from
"prove a direct component bound for all residual rows" to:

```text
either prove a direct component/floor theorem on a finite bounded-apex packet
or pass through the normalized slow/ruler chart before using any denominator
net statement globally.
```

The same ledger records `H7_pair_shadow`, `even_pair_shadow`, mod-7 and mod-14
residue-count scissors signatures, and Farey `p+q,p*q mod 91` lanes.  Those
sidecars are not decorative: live rows with the same terminal status split
into `5` mod-7 scissors signatures in the sample, so the next packet-bank
ledger must keep these fields before asking a moment or branch chart to glue.

## S259 Lean Bridge

`TournamentH7.LRCObserverGluingLedger` lifts this ledger into a checked proof
interface.  `DenominatorNetNumerics` stores the exact largest-arc fraction and
reciprocal grid threshold, and the S258 samples are now named obligations:

```text
s258H7Eq6BoundaryObligation:      19/1372, D=73
s258DivisorLoadedB8Obligation:    1/82320, D=82321
```

The module does not turn these rows into proof exits.  It separates an
`ObserverGluingObligation` from an `ObserverGluingCertificate`; only the
certificate carries a `TerminalDischargeCertificate` and hence implies
`Mreach >= 1/14`.  The new theorem
`lrc14_from_observer_gluing_coverage` says the remaining global theorem can be
phrased as coverage by early gates or such certificates.

This sharpens the witness-route target.  The direct denominator-net field is
legal for bounded-apex rows only after a uniform component/floor theorem is
proved.  Divisor-loaded rows should instead route through normalized
slow/ruler coordinates, moment/Perron discharge, branch/K33 discharge, or
named finite-denominator debt.

## Proof Route Target

The current sharpened theorem target is:

```text
primitive 13-speed row S
  -> 14-free q-witness exit, or covering residual
  -> if count_7_divisible(S) >= 7, discharge by THM-573
  -> residual core count_7_divisible(S) <= 6
  -> if tight, route to AP/GW/dilation boundary debt
  -> if non-tight:
       prove mu(L(S)) >= m0
       prove components(L(S)) <= A0
       infer largest_arc(L(S)) >= m0/A0
       infer denominator-net witness for all d >= ceil(A0/m0)
       infer Conjecture 7.1(13) for the residual
       infer LRC14 for the residual
```

The open step is the direct component bound for `L(S)` at threshold `1/14`.
THM-565 already supplies a scale-separated arc-count bound for the
`maxgap > 1/7` witness object after peeling.  This HYP promotes the needed
bridge as a named proof obligation:

```text
direct_lonely_component_bound:
  THM-565 Framing A component control
  -> component control for L(S) itself,
     or a loss-controlled reduction from L(S) to the maxgap witness object.
```

## Tournament Analysis

Vertices are proof obligations and retained sidecars, not runners:

```text
largest_arc_denominator_net
direct_lonely_component_bound
lonely_measure_floor
crt_c7_level7_lift_THM573
crt_c2_dyadic_lift
continuous_I_substitute
finite_principal_part_bad_denominator_budget
hyperoperation_grid_address
polynomial_prime_field_packet
raw_I_table_enumeration
```

Pairwise observable:

```text
preserves_LR_predicate
reduces_CRT_debt
retains_exact_denominator_clock
retains_endpoint_or_component_topology
replaces_or_bounds_I_k_p_1
keeps_finite_bad_denominator_budget
names_destroyed_coordinate
provides_terminal_exit
```

Gauge: orient `A -> B` when `A` moves closer to a uniform witness-time theorem
while retaining more of the LRC predicate and destroying less packet data.
Tie Hamiltonian path:

```text
largest_arc_denominator_net
> direct_lonely_component_bound
> lonely_measure_floor
> crt_c7_level7_lift_THM573
> crt_c2_dyadic_lift
> continuous_I_substitute
> finite_principal_part_bad_denominator_budget
> hyperoperation_grid_address
> polynomial_prime_field_packet
> raw_I_table_enumeration
```

Fingerprint for this synthetic routing tournament: transitive, score
histogram `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}`, no directed cycles,
singleton SCCs, one Hamiltonian path, and zero edge flips against the tie path.

## Assumption Challenge

Candidate vertices considered: runners, residues modulo `14`, CRT factors,
polynomial interpolation points, grid witness denominators, finite q-series
principal parts, lonely intervals, component boundaries, hyperoperation cells,
proof obligations, and terminal exits.

The chosen vertices are proof obligations because the preserved predicate is
not "same residue" or "same operation-grid address"; it is existence of a
lonely witness at threshold `1/14`, strengthened to a denominator-net witness
for all large `d`.  The quotient destroys the individual runner geometry,
raw tuple enumeration, and exact finite-field table rows, but retains the CRT
lift status, lonely-set topology, denominator clock, and named terminal exit.

The challenged assumption is that Conjecture 7.1(13) is automatically
equivalent to the scalar LRC14 statement.  It is the right proof-strength
target for the project, but the equivalence needs the finite-union/compact
packet normalizer and the uniform largest-arc lower bound.  Until those are
proved, the safe implication is:

```text
uniform largest-arc route
  -> Conjecture 7.1(13)
  -> LRC14.
```

## Remaining Work

- Split the direct `1/14` component-bound obligation into bounded-apex direct
  packets versus large-apex normalized slow/ruler packets, as exposed by the
  S258 ledger.
- Make the `c=2` dyadic lift as explicit as THM-573 is for `c=7`.
- Turn the continuous substitute for `I(13,p,1)` into a theorem with a finite
  small-denominator budget.
- Attach the S258 ledger fields to HYP-2963 packet rows and outside-bank
  normalizer attempts: direct component data, CRT status, pair-scissors
  signatures, Farey lanes, and named gluing debt.
- Formalize the route as a finite-ruler statement: largest arc
  `>= ell0` implies witnesses in `(1/d)Z` for every `d >= ceil(1/ell0)`.
