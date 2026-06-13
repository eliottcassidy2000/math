---
id: HYP-2165
status: SUPPORTED by S609 bounded bridge computation; full lift theorem open
source: user-2026-06-03; codex-2026-06-03-S609
related:
  - HYP-2164
  - HYP-2163
  - HYP-2162
  - HYP-2161
  - HYP-2160
  - HYP-2101
  - HYP-2100
  - HYP-2099
  - HYP-2098
  - HYP-2097
  - HYP-2096
  - HYP-2095
  - THM-401
  - THM-406
---

# HYP-2165: the n=14 `Res_27` bridge separates parity scaffolds from owner certificates

## Claim

The useful next n=14 proof lemma is a bridge between two quotients:

```text
64 self-converse round classes
  -> C=27 unit/nonunit/owner certificate layer.
```

The first quotient carries HYP-2160's solved-face parity truth.  The second
quotient carries the HYP-2095/HYP-2099 certificate predicates: unblocked small
pair, positive measure, D/U/N obligations, and unit-spine slack owners.

S609 supports the following bounded bridge form:

```text
In the canonical C=27 unit-spine slack fibre through slack <=81,
every full D/U/N cover is certified by either
  (a) an unblocked small pair, or
  (b) positive strict safe measure.

The only measure-zero rows are AP and V*, and both use (1,13) at 1/14.
```

## S609 Evidence

`04-computation/lrc_n14_res27_fixed_bridge_s609.py` combines the S578 fixed
round scaffold with an extended `C=27` owner scan.

### Fixed-Class Parity Table

For the `64` self-converse round classes on `m=13` moving runners:

```text
all Hamiltonian-path counts odd: True
distinct Hamiltonian-path counts: 64/64
H min=1, H max=3711175
score-span histogram: {0:1, 2:20, 4:21, 6:12, 8:7, 10:2, 12:1}
SCC-count histogram: {1:63, 13:1}
```

This verifies HYP-2160's odd-parity transfer on the actual n=14 fixed table,
not only on small-model round tournaments.  It also warns against treating the
scalar Hamiltonian-path count as a loneliness certificate: the regularity gauge
and HP-scarcity gauge disagree on `2005/2016` pair edges.  Parity is the robust
solved-face invariant; scalar `H` is a class-only shadow.

### `C=27` Unit-Spine Slack Scan Through `81`

S578 scanned canonical slack through `42`.  S609 extends the canonical slack
fibre itself through `81`:

```text
slack rows:                         17550
full D/U/N quotient covers:          9506
clock-ledger failures:               8044
full covers with unblocked pair:   9504/9506
block-all positive-measure rows:      2/2
measure-zero full covers:              2
measure-zero + cheap pair:             2/2
open residual rows:                    0
```

The route histogram is:

```text
cheap pair + positive measure: 9502
floor via cheap pair:             2
block-all but positive-measure:   2
clock-ledger failure:          8044
```

The two floor rows are exactly:

```text
AP: slack=(3,6,9,12)
V*: slack=(3,6,9,24)
```

Both use the cheap pair `(1,13)` at `1/14`.  The two no-cheap block-all controls
are unchanged from S578 and are positive-measure:

```text
(3,9,12,42)
(3,12,27,42)
```

Thus the canonical bounded owner scan improves from S578's `531` full covers
through `42` to S609's `9506` full covers through `81`, with the same zero
residual.

## Interpretation

The 64 self-converse classes are a parity scaffold, not a proof certificate by
themselves.  They preserve round/converse/fixed-boundary information and force
Redei parity, but they forget the speed-owner labels that decide whether a
small pair is blocked.

The `C=27` owner layer restores the missing labels.  In the canonical fibre
through one full `3C=81` lift range, it has no live residual:

```text
measure-zero => cheap pair,
block-all    => positive measure.
```

Incoming HYP-2162/THM-407 supplies the adjacent shell-orbit reduction: after
the twisted-involution normalization, the `13` raw `C=27` shells fold to three
gcd strata `{1,3,9}`.  S609 is the owner-certificate complement to that result:
its shell and gcd signatures test whether the normalized strata still carry
the cheap-pair / positive-measure discharge.

Incoming HYP-2163 supplies the stronger efficient proof-pipeline claim up to
`n=14`, and incoming HYP-2164 supplies a sharper `Res_27` pinch-certificate
thread.  S609 should therefore be read as a bridge/compatibility artifact: it
explains where the fixed-class Redei parity table and the `C=27` owner
certificate ledger sit relative to those routes.

So the bridge theorem should not try to prove loneliness directly from the
unlabelled 64 class table.  It should prove that any n=14 fixed-boundary
realization either:

1. lifts to the `C=27` owner layer where the S609 dichotomy applies; or
2. fails to lift because some owner obstruction already exposes an unblocked
   small pair or positive-measure escape.

## Tournament Analysis

S609 uses two Tournament Analysis quotients.

Fixed-class vertices:

```text
vertices: 64 self-converse round classes
observable: score span, anti symmetry, exact Hamiltonian-path count
switch: regularity burden, with class-id tie path
fingerprint: transitive ranker; 0 directed 3-cycles; 2005/2016 edge flips
             between regularity and HP-scarcity gauges
```

Slack-row vertices:

```text
vertices: canonical C=27 full-cover slack rows
observable: route, positive measure, private D/U/N count, max speed
switch: residual > no-cheap positive > floor > cheap positive
fingerprint: transitive proof ledger; 0 directed 3-cycles; no open residual
```

This is intentional.  The current proof route is a transitive certificate
ledger, not a cyclic residual.  If a strict SCC appears, it should appear only
after arbitrary unit representatives, endpoint-owner fibres, or fixed-class
realization data are reattached.

## Assumption Challenge

Candidate vertices considered:

```text
runners,
gaps,
fixed circle sections,
section boundaries,
wall-crossing events,
residues,
cover arcs,
Fourier modes,
matroid circuits,
proof obligations,
self-converse round classes,
C=27 slack rows.
```

S609 chooses the last two because they preserve complementary predicates:

```text
fixed classes preserve solved-face parity;
slack rows preserve owner/certificate discharge.
```

The challenged assumption is that one quotient must do everything.  The
evidence says the n=14 proof should use a coimage bridge: parity lives in the
fixed-class quotient, but the loneliness certificate lives in the owner quotient.

## Honest Status

This is not a proof of LRC n=14.  It is a stronger bounded bridge target.

Open tasks:

1. Prove the canonical `C=27` slack dichotomy without the `<=81` bound.
2. Prove that every tight-boundary fixed-class realization can be functorially
   lifted to the owner layer, or that failed lifting exposes an immediate
   certificate.
3. Reattach arbitrary simultaneous unit representatives beyond the one-lift
   S578/S579 audits.

## See

`04-computation/lrc_n14_res27_fixed_bridge_s609.py`,
`05-knowledge/results/lrc_n14_res27_fixed_bridge_s609.out`,
`07-reflections/lrc-n14-res27-fixed-bridge-s609.md`,
`01-canon/theorems/THM-407-twisted-involution-shell-reduction-of-the-LRC-additive-residual.md`,
`07-reflections/lrc-twisted-involutions-laminar-flow-shells-s599.md`,
`04-computation/lrc_twisted_involution_flow_shells_s599g.py`,
`05-knowledge/hypotheses/HYP-2161-coimage-yoneda-2nm1-resonance-cancellation.md`,
`07-reflections/lrc-redei-solved-face-transfer-menu-s599.md`,
`05-knowledge/hypotheses/HYP-2099-lrc-n14-fixed-round-certificate-scaffold.md`,
`05-knowledge/hypotheses/HYP-2100-lrc-n14-unit-spine-exchange-cheap-pair-sieve.md`.
