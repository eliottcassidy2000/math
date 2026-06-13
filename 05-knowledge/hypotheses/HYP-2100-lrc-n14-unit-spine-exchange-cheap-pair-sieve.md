---
id: HYP-2100
status: SUPPORTED by bounded S578 exchange audit; full multi-unit/fixed-class proof open
source: codex-2026-06-03-S578
related:
  - HYP-2101
  - HYP-2099
  - HYP-2098
  - HYP-2097
  - HYP-2096
  - HYP-2095
  - HYP-2094
  - HYP-2090
  - HYP-2088
  - THM-397
---

# HYP-2100: unit-spine exchange first passes through the cheap-pair sieve

The HYP-2096 unit-spine exchange lemma should be split into two stages.

```text
U0. cheap-pair sieve:
    any large unit-shell representative that creates a full D/U/N lifted row
    already exposes a HYP-2095 unblocked small-pair witness, unless it lies in a
    still-smaller no-cheap residual.

U1. true exchange:
    only after U0 fails do we need to lower unit representatives to the
    canonical spine while preserving D/U/N private pivots.
```

S578 gives bounded evidence that the U1 residual is empty in the tested local
unit-lift regimes.

## Evidence

`04-computation/lrc_n14_unit_spine_exchange_s578.py` audits:

- all canonical slack quadruples through `42` with exactly one unit-shell
  representative replaced by any same-shell representative through `81`;
- selected slack layers AP, `V*`, first open-gap, and two non-cover controls
  with all one- and two-unit lifts through `81`, plus several all-shell extreme
  representative patterns.

The main result:

```text
one-unit-lift rows:              45045
one-unit-lift full D/U/N covers: 13169
cheap unblocked-pair witnesses:  13169
no-cheap residual:               0
```

For the named local/extreme stress:

```text
AP slack:             240/240 full covers cheap
V* slack:             336/336 full covers cheap
first open-gap slack: 240/240 full covers cheap
zero-slack control:     0 full covers
double-3 control:       0 full covers
```

The cheap witnesses are not random. In the one-lift scan, the most common
witness denominator is `14`:

```text
D=14: 7943 witnesses
D=26: 1206
D=20: 1067
D=40:  873
```

So the lifted unit representative tends to leave the old `n`-straddle or a
nearby small-pair pinch unblocked. This is exactly the HYP-2095 paired-or-
anchored route.

## Interpretation

The original exchange lemma asked:

```text
Can large unit-shell representatives be lowered to the canonical spine?
```

S578 suggests the sharper local question is:

```text
Why is a large unit representative ever needed if a full lifted row already
has an unblocked small-pair witness?
```

In the tested regime, every full lifted row is already done by HYP-2095. There
is no no-cheap row left for the exchange lemma to lower.

The two non-cover controls are also informative. Their canonical failures are
`D9`/`D12`-style composite obligations. Unit representatives are coprime to
`27`, so they cannot repair the `gcd 3` D-obligations that actually make those
canonical slack layers fail. This explains why the controls had zero full
covers.

## Next Lemma

The local proof target becomes:

```text
Unit-lift cheapness lemma.

Fix n=14 and C=27. Start from one representative in each unit shell and four
nonunit/zero-residue slack runners. If changing one unit representative within
its shell yields a full D/U/N quotient cover, then either there is an unblocked
small-reduced-sum pair, or lowering that representative loses no D/U/N
obligation that could matter in the strict sub-edge residual.
```

The bounded audit supports the stronger first branch for one-unit lifts:

```text
full D/U/N cover => unblocked small pair.
```

After that is proved, the hard exchange lemma only needs to handle simultaneous
multi-unit interactions that destroy all cheap pairs. That is exactly the kind
of row the 64 fixed-round certificate table from HYP-2097 should isolate.

S579/HYP-2101 refines this one step further: simultaneous lift objects should
be treated as an apex-lift certificate sheaf.  The denominator-14 apex chart is
large but not complete, side-chart cheap witnesses are real, and ledger-failure
restrictions must be legal local sections.

## Tournament Analysis

Vertices are exchange audit patterns:

```text
one-unit-lift all-slack,
AP slack,
V* slack,
first open-gap slack,
zero-slack control,
double-3 control.
```

Pair observable:

```text
(bad lowering count, below-floor count, no-cheap residual count,
 full-cover count, max representative)
```

Switch:

```text
harder = more bad lowerings, then more below-floor rows, then more no-cheap
residual, then more full covers.
```

Fingerprint:

```text
directed_3_cycles=0
sccs=6 singleton SCCs
hardness order starts with one_unit_lift_all_slack.
```

This remains a transitive proof ledger. The interesting cyclic object is still
expected inside the HYP-2097 64-row fixed-round certificate table, after the
cheap-pair sieve and unit-spine labels are attached.

## Assumption Challenge

Possible tournament vertices considered:

```text
runners, unit shells, exchange moves, D/U/N obligations, small pairs, endpoint
owners, slack rows, and fixed self-converse round classes.
```

Chosen quotient:

```text
exchange audit patterns.
```

Predicate preserved:

```text
whether large unit-shell representatives are needed after HYP-2095 removes
rows with unblocked small-pair witnesses.
```

Information destroyed:

```text
the full 64 fixed-class identity and arbitrary simultaneous choices of all nine
unit representatives.
```

Challenged assumption:

```text
that the exchange lemma must begin by lowering large unit representatives.
```

The audit says to sieve for unblocked small pairs first. In the bounded local
regime, that sieve empties the exchange problem.

## Files

- `04-computation/lrc_n14_unit_spine_exchange_s578.py`
- `05-knowledge/results/lrc_n14_unit_spine_exchange_s578.out`
- `07-reflections/lrc-n14-unit-spine-exchange-cheap-pair-sieve-s578.md`
