# LRC14 Observer-Gluing Lean Frontier

*codex-2026-06-27-S259. Formal continuation of S258, HYP-3095, HYP-3096,
and HYP-3097.*

This pass adds the Lean module
`04-computation/lean/TournamentH7/TournamentH7/LRCObserverGluingLedger.lean`.
It does not prove LRC14.  It formalizes the bleeding-edge proof interface that
S258 exposed:

```text
observer charts + overlap discharge + direct-arc/pair-scissors fields
  -> terminal discharge certificate
  -> Mreach >= 1/14
```

The new conditional theorem is:

```text
ObserverGluingCoverage -> LRC14Statement
```

where `ObserverGluingCoverage` says every nonzero 13-speed row is either
discharged by an early gate or emits an `ObserverGluingCertificate`.  A
certificate is deliberately stronger than a ledger row: it must contain an
existing `TerminalDischargeCertificate`, so the module cannot smuggle in the
missing analytic proof.

## Task A: formal edge

The module introduces:

- `ObserverChart`: arithmetic lift, direct arc, normalized arc, pair scissors,
  cap/Pascal, moment/Perron, branch/K33, fine-scale tournament, formal witness,
  raw scalar shadow, and coarse mod-14 winding.
- `GluingDischarge`: reconstructs, dual-annihilates, fiber-constant, descends,
  terminal exit, or named debt.
- `ChartOverlapCertificate`: an audited overlap with a named destroyed
  coordinate.
- `DenominatorNetNumerics`: exact largest-arc numerator/denominator and the
  reciprocal grid threshold.
- `ObserverGluingObligation`: the nonterminal S258 row shape with direct
  component counts, pair shadows, mod-7/mod-14 scissors signatures, Farey
  lanes, coarse-winding status, overlap certificates, and remaining debt.
- `ObserverGluingCertificate`: an obligation plus terminal discharge.

Two S258 rows are now named as Lean obligation examples:

```text
H7=6 boundary residual:    largest arc 19/1372, threshold D=73
divisor-loaded B=8 row:    largest arc 1/82320, threshold D=82321
```

This is useful because it preserves the exact obstruction.  The first row is a
bounded-looking direct-arc candidate; the second is the THM-575 warning that a
raw direct largest-arc scalar shatters under divisor loading.

## Task B: synthesis and refinement

The formalization clarifies the current proof frontier:

```text
prove ObserverGluingCoverage
```

not merely "find a scalar lower bound."  The producer theorem should split
the live residual into at least four certificate families:

1. bounded-apex direct packets, where a direct component/floor theorem can
   justify the denominator-net field;
2. large-apex normalized slow/ruler packets, where THM-575 says raw time is
   the wrong coordinate;
3. moment/Perron packets, likely using the gK8/S2/reflection side of the
   cap bound;
4. branch/K33 packets, where active binders and endpoint-owner words must
   build the state lift or name residual debt.

The coarse mod-14 winding warning is also formal now:

```text
coarseWinding_degenerate_not_terminal
```

An antipodal-tie coarse winding chart is not a terminal proof carrier.  It may
remain as a shadow, but the proof must attach fine-scale mod-p data, packet
scissors, or an independent terminal discharge before using it.

## Relation to the sexy prime thread

The useful transfer is structural, not theorem-level.  Both LRC14 and fixed
gap prime pairs need ledgers that distinguish a count shadow from a proof
carrier.  For LRC14 the shadow is coarse H, pair/Pascal mass, or a direct arc
scalar; for sexy primes it is local residue survival or a Hardy-Littlewood
singular-series factor.  The shared discipline is:

```text
local shadow
  -> retained sidecars
  -> named destroyed coordinate
  -> terminal analytic/sieve input or explicit debt
```

The LRC module therefore suggests a `sexy_prime_pair_ledger` shape with
residue sidecars, parity/decomposition debt, distribution input, and terminal
exit fields.  It does not transfer an LRC proof to primes, but it does prevent
the same overclaim: local equinumerosity is not equidecomposition and is not a
lower-bound sieve.

## Tournament Analysis

This session again avoids using runners as tournament vertices.  I considered
the natural alternatives:

```text
runners, gaps, fixed circle sections, section boundaries, wall crossings,
residues, cover arcs, Fourier modes, matroid circuits, proof obligations,
observer charts.
```

The chosen vertices are observer charts because that quotient preserves the
predicate actually needed by the proof: the existence of a terminal route to
`Mreach >= 1/14`.  It destroys raw runner identity and some exact scale, so
the record forces the destroyed coordinate and discharge mode to be explicit.

Pairwise observable: chart `A` beats chart `B` when `A` preserves the LRC
predicate needed by `B` and retains a coordinate that `B` would otherwise
forget.  Tie path:

```text
formal_witness
> fine_scale_tournament
> branch_k33
> moment_perron
> normalized_arc
> direct_arc
> pair_scissors
> cap_pascal
> arithmetic_lift
> coarse_mod14_winding
> raw_scalar_shadow
```

The challenged assumption is the point of the module: coarse tournament H,
pair/Pascal mass, and direct largest arcs are not proof vertices until their
observer overlaps are certified or routed to named debt.

## Next pull

Populate `ObserverGluingObligation` from HYP-2963 packet rows and outside-bank
normalizer attempts.  For each row, require one of:

```text
bounded direct certificate
normalized slow/ruler certificate
moment/Perron terminal discharge
branch/K33 terminal discharge
explicit named residual debt
```

Then prove producer lemmas from those populated obligations into
`ObserverGluingCertificate`, keeping `ObserverGluingCoverage` as the global
target.
