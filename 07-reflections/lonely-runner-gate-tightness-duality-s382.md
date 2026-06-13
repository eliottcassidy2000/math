# Lonely Runner Gate Tightness Duality

codex-2026-05-31 S382

The sentence "at `n=14`, tightness requires no multiple of `14`, but a
counterexample requires one" is almost exactly the right shape.  The useful
version is:

```text
no 14-gate  => unit skeleton survives as safe boundary mass
14-gate     => unit skeleton is closed, but endpoint debt moves elsewhere
```

This is not just a computational vibe.  For the six unit points

```text
a/14,  a in {1,3,5,9,11,13},
```

a speed `v` forbids `a/14` exactly when

```text
||v a / 14|| < 1/14.
```

Since `a` is a unit mod `14`, this strict inequality is possible exactly when

```text
v == 0 mod 14.
```

So any set with no `14`-multiple leaves all six unit points safe.  If such a
set has forbidden length `1`, it is automatically boundary-only, not an open
cover.  This is the clean counterexample obstruction.

## What The Probe Adds

`lonely_runner_14_gate_tightness_s382.py` checked the local tight basin around
the initial segment.

The initial segment is boundary-only:

```text
(1,2,3,4,5,6,7,8,9,10,11,12,13)
```

It has no `14`-multiple and all six unit witnesses remain unprotected.

Replacing one speed by a `14`-multiple immediately protects the unit skeleton,
but in every checked case it creates positive gaps.  The pure gate replacement
ledger

```text
remove r, add 14q, q <= 12
```

gave `156/156` positive-gap outcomes.

The scan did find a second local boundary-only example:

```text
(1,2,3,4,5,6,7,8,9,10,11,13,24).
```

This is the nuance.  Tightness is not "the initial segment only"; it is
"unit-skeleton boundary behavior", and that behavior still avoids `14`-gates.

## The Search Split

The right proof/search split now seems very crisp:

```text
Branch A: no 14-gate.
  Unit points survive.
  This branch cannot contain an open-cover counterexample.

Branch B: at least one 14-gate.
  Unit points are protected.
  The search must account for new positive gaps / descendant endpoint debt.
```

This explains why speed-first disproof attempts keep feeling like they are
doing two incompatible jobs.  Tightness wants to preserve the Dirichlet unit
skeleton.  Counterexamples must kill it.  Killing it is cheap, but the debt
reappears at `98`, `182`, `196`, and related quotient layers.

The next useful theorem would be an endpoint-debt lower bound: after a
`14`-gate protects the six unit points, some descendant endpoint layer must
remain exposed unless a full protection cycle appears.
