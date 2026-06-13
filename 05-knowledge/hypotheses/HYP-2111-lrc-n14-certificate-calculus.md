---
id: HYP-2111
status: FORMALIZATION + extension program; no new proof of n=14
source: codex-2026-06-03-S580
related:
  - THM-398
  - HYP-2108
  - HYP-2107
  - HYP-2106
  - HYP-2105
  - HYP-2104
  - HYP-2103
  - HYP-2102
  - HYP-2101
  - HYP-2100
  - HYP-2099
  - HYP-2095
---

# HYP-2111: the n=14 LRC frontier is a certificate calculus with one common positivity residual

## Claim

The cutting-edge n=14 LRC work now has a useful formal shape:

```text
SpeedRow -> CertificateSection or NamedResidual
```

where the active certificate sections are:

```text
witness_1_over_n,
cheap_pair,
positive_measure,
ledger_failure,
residual.
```

The gates compose as a cascade product of conditional clearances.  A row exits
when some gate supplies a closed witness, an unblocked small-pair witness, or a
positive-measure certificate.  A row remains only if it lies in a named residual
that survives all restrictions.

HYP-2104 now sharpens the Cprime side: the Vitali handoff is exactly `n|v`,
and the Bprime/Vitali covering criterion proves the long-interval branch while
leaving only short-interval arc alignment.  The new synthesis is that the two
most important open residuals seem to be the same kind of object:

```text
HYP-2103/HYP-2104 all-short Vitali alignment residual:
  one AP of thin danger arcs covers the fixed safe set G(S')

HYP-2101 failed sheaf gluing residual:
  all cheap-pair local sections are shielded or anchored

Common target:
  endpoint-owner cover circuit -> private pivot / peel -> positive measure.
```

So the next extension should try to prove endpoint-cover-circuit positivity,
not merely enumerate more 64-class rows.

## Formal Gates

S580 records the following gates:

| Gate | Status | Output | Residual |
|---|---|---|---|
| `G0_no_multiple_n_clock` | proved | closed witness `t=1/n` | none |
| `G1_Cprime_reduction` | proved reduction | LRC follows from Cprime | Cprime |
| `G2_dominance_or_long_interval_dodge` | proved partial Cprime | positive measure | all-short Vitali alignment |
| `G3_all_short_Cprime_residual` | open | positive measure for small multiples | AP thin-arc alignment cover |
| `G4_paired_or_anchored_cheap_pair` | reduction + verified split | cheap-pair pinch witness | block-all positivity |
| `G5_fixed_boundary_owner_functor` | bounded scaffold | 64 fixed fibres with owner labels | realization independence |
| `G6_unit_lift_cheap_sieve` | bounded zero residual | cheap pair before exchange | simultaneous multi-unit no-cheap rows |
| `G7_apex_lift_certificate_sheaf` | bounded zero restriction residual | global section or positive measure | apex monodromy |
| `G8_endpoint_cover_circuit_positivity` | open bridge | positive measure from failed gluing | ownerless cover circuits |
| `G9_tie_wall_limit_functor` | open formalization | labelled tie-wall object | perturbation-direction ambiguity |
| `G10_certificate_calculus_closure` | proposed | all rows exit or residual SCC impossible | terminal residual acyclicity |

This makes THM-398 the global spine and HYP-2104 the measure/construction
handoff.  If a row has no multiple of `n`, `G0` certifies it by construction.
If it has a multiple of `n`, all remaining work is Cprime on the Vitali-covering
side.  The dominance/long-interval dodge proves the easy Cprime side; only
all-short arc-alignment rows remain.  HYP-2095/HYP-2101 add a second exit: if
the row has an unblocked small pair, it is done by a pinch certificate.

## S580 Tournament Analysis

The vertices are proof obligations/gates, not runners and not naked classes.

Pair observable:

```text
(open_rank, residual_weight, dependency_count, -evidence_strength, tie_order)
```

Switch:

```text
harder unresolved proof obligation beats easier certified gate.
```

Fingerprint from `04-computation/lrc_n14_certificate_calculus_s580.py`:

```text
score_hist: {0:1, 1:1, ..., 10:1}
directed_3_cycles: 0
sccs: 11 singleton SCCs
hamiltonian_path_count: 1
hardness top:
  G10_certificate_calculus_closure
  G8_endpoint_cover_circuit_positivity
  G3_all_short_Cprime_residual
  G5_fixed_boundary_owner_functor
  G9_tie_wall_limit_functor
```

This is still a transitive proof ledger.  Cycles are not expected until
certificate transports inside fixed fibres are made explicit.  If nontrivial
cycles appear there, they should be interpreted as apex monodromy or
incompatible endpoint-owner transports, not as Hamiltonian-path complexity in
the all-0 staircase warning.

## Extension Queue

1. **All-short Vitali alignment to endpoint cover.**  Enumerate small-multiple
   all-short rows with `v=14` and record which endpoint owners cover every
   component of `G(S')`.  If every AP-cover attempt becomes a THM-397 owner
   circuit that peels, this proves the n=14 Cprime residual left after HYP-2104.
2. **Fixed-fibre sheaf table.**  Attach local section data to AP, `V*`,
   transversal flips, and minimal `gcd 3/gcd 9` lifts: cheap pair, reduced sum,
   mod-7 owner, mod-27 shell, shields, anchors, and D/U/N private pivots.
3. **Residual acyclicity.**  Build the directed graph of S579 certificate
   transports.  Union-only `{ledger,ledger}->cheap` rows should be two-step
   composites, not anomalies.
4. **Component discrepancy.**  Use Erdos-Turan, Vitali density, or
   three-distance tools on `{k/(nw)}` against the interval components of
   `G(S')`.
5. **Level-3 cascade product.**  Treat the proof ledger itself as the repeated
   multiplication closure of runner clearances: gate factors are witness,
   cheap, positive-measure, ledger-failure, and residual.

## Assumption Challenge

Possible tournament vertices considered:

```text
runners, gaps, fixed sections, section boundaries, wall events, residues,
cover arcs, Fourier modes, certificate germs, and proof obligations.
```

The chosen quotient uses proof obligations.  It preserves the predicate needed
by this session:

```text
every row gets a witness, a cheap section, a positive-measure section,
or a named residual.
```

It destroys exact speed values, full round-class identity, and most
Hamiltonian-path complexity.  That is intentional.  The 64 fixed classes are a
routing surface, not the proof object until owner labels and sheaf restrictions
are attached.

## Honest Status

This is not a proof of n=14.  It is a sharper formalization of the current
frontier after THM-398/HYP-2103.  The useful extension is concrete: prove that
the terminal residuals are endpoint-cover circuits with private owners/pivots,
or prove the all-short Cprime residual analytically by discrepancy.

**See:** `04-computation/lrc_n14_certificate_calculus_s580.py`
(+ `05-knowledge/results/lrc_n14_certificate_calculus_s580.out`),
`07-reflections/lrc-n14-certificate-calculus-formalization-s580.md`,
THM-398, HYP-2104, HYP-2103, HYP-2102, HYP-2101, HYP-2100, HYP-2099, HYP-2095.
