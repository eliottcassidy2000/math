# LRC14 Ramanujan Route Scheduler Reflection

Source: codex-2026-06-26-S200
Anchors: HYP-3036, T1117, LTI-184, LTT-082, HYP-3033, HYP-3030, HYP-3031

## What Changed

HYP-3030 left a clean status situation but a small route debt:

```text
coarse ET + Henselian unit gate
  -> 0 mixed boundary/open fibers
  -> 15 mixed route fibers, all strict-open
```

The S200 computation tests whether that route debt can be scheduled without
going all the way to exact `M`.  On the `38` packets in the S194 residual list,
the answer is yes: attach the primitive safe-residue deck for `2 <= q <= 13`.

```text
coarse_plus_primitive_deck_2_13 -> 0 mixed route, 0 mixed status
coarse_plus_first_q_2_13        -> 0 mixed route, 0 mixed status
```

This is not a proof of the full theorem.  It is a sharper proof carrier for
the exact place where S194 stopped.

## Why This Is Not Just `q_threshold`

The labelled packet classifier already records `q_threshold`, so merely saying
"use q0" would be circular as a new carrier.  The useful object here is the
primitive phase deck:

```text
D_q(S) = #{a mod q : gcd(a,q)=1 and ||a v / q|| >= 1/14 for every v in S}.
```

The first positive `D_q` for `q <= 13` recovers the direct witness route in the
residual packets, but the whole deck is a theorem-facing object: it can be
stated as finite exact-period phase coverage, compared against Ramanujan
projectors, and added to packet sidecars without retaining exact magnitude.

## Guardrail

The `q=14` layer must stay separate.  Many covering rows in the residual fibers
have primitive safe mass at `q=14`, while their `q<=13` deck is identically
zero.  Therefore the candidate scheduler is:

```text
positive deck mass for q<=13 -> direct Q-WITNESS route
zero deck mass for q<=13     -> covering / q=14 / boundary-moment route
```

If a future quotient merges `q=14` into the direct witness deck, it will erase
the distinction HYP-3036 needs.

## Assumption Challenge

I did not make runners the tournament vertices.  The considered vertex sets
were runners, gaps, residue words, primitive denominator charts, Ramanujan
modes, safe phase residues, coarse zipper fibers, exact magnitudes, topology
gates, and proof obligations.

The chosen vertices are proof carriers:

```text
primitive_count_deck_2_13
first_safe_q_2_13
ramanujan_trace_deck_2_14
coarse_et_unit_status_gate
exact_magnitude_cocycle
raw_residue_terminal_word
```

The quotient preserves boundary/open status inherited from HYP-3030 and direct
`q<=13` witness availability.  It destroys exact `M`, exact safe interval
lengths, full safe-component barcode, and arc-Cech topology.

## Next Pull

Add `primitive_safe_deck_2_13` and `first_primitive_safe_q_2_13` to the
HYP-3027/HYP-3031 packet sidecar.  Then rerun a cached full-bank ledger to ask
three questions:

1. Does the primitive deck ever introduce mixed boundary/open status?
2. How many route-mixed fibers remain outside the S194 residual set?
3. Are remaining zero-deck open rows exactly covering-moment, q=14 boundary,
   petal, K33, or F7/THM-572 residual debt?
