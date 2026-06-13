---
source: codex-2026-06-03-S578
status: exchange-lemma sharpening + bounded audit
tags:
  - lonely-runner
  - n14
  - unit-spine
  - exchange-lemma
  - cheap-pair
  - paired-anchored
  - private-pivots
---

# n=14 Unit-Spine Exchange: Cheap Pair First

I went into this pass expecting to find a real lowering problem: a large
unit-shell representative might cover some D/U/N obligation that the canonical
spine does not, so lowering it would need a careful private-pivot exchange.

The audit found something cleaner in the tested regime.

Once HYP-2095 is allowed to fire first, the residual disappears.

## What Was Tested

The new script `04-computation/lrc_n14_unit_spine_exchange_s578.py` tests two
bounded regimes.

First, it does the local exchange audit exhaustively:

```text
canonical slack quadruples through 42
one unit-shell representative changed within its shell through 81
```

That is `45045` one-unit-lift rows. Among them, `13169` are full D/U/N quotient
covers.

Every one of those `13169` full covers has a HYP-2095 unblocked small-pair
witness. So:

```text
full cover = 13169
cheap pair = 13169
no-cheap residual = 0
```

Second, it checks selected named slack layers with all one- and two-unit lifts
through `81`, plus all-reflected/all-shifted extreme unit patterns:

```text
AP slack:             240 full covers, all cheap
V* slack:             336 full covers, all cheap
first open-gap slack: 240 full covers, all cheap
zero-slack control:     0 full covers
double-3 control:       0 full covers
```

So the candidate exchange obstruction is not showing up. The large unit
representatives are not quietly preserving a hard no-cheap residual; they are
opening a small-pair pinch.

## The Useful Shape

The exchange lemma should not start with "lower the representative." It should
start with:

```text
Does this lifted unit row have an unblocked small pair?
```

In the tested local regime, yes, always.

The witness histogram explains why this is not accidental. In the one-lift
scan, the cheap witness denominators are led by:

```text
14, 26, 20, 40, 28, 22, 49, 16.
```

The `D=14` witnesses dominate. That is the old straddling-pair route in a new
guise: the unit lift perturbs the unit spine, but it usually leaves the
`n`-clock/small-pair route open.

## What The Controls Say

The two non-cover controls matter because they show the other side.

The canonical zero-slack and double-3 controls fail D-obligations like `D9` and
`D12`. Unit representatives are coprime to `27`, so they cannot repair the
gcd-3 part of those D-obligations. In the local/extreme unit stress they never
became full covers.

That points toward a proof sentence:

```text
Unit-shell changes can perturb unit and n-clock ownership, but they cannot
repair the composite D9/D12 debt. If they make the row full, they tend to
create an unblocked small-pair witness first.
```

The first half is arithmetic. The second half is now a precise local lemma to
prove.

## How This Advances HYP-2096

HYP-2096 said:

```text
normalize the unit-shell spine or prove an exchange that lowers it while
preserving private U pivots.
```

S578 says the front of that lemma should be:

```text
large unit representative
  -> cheap-pair sieve
  -> only then exchange/lower
```

In the bounded local regime, the last arrow is never reached.

That is genuine progress. It means the exchange lemma may be less about
constructing a delicate lowering homotopy and more about proving that local
unit-shell motion exposes the paired-or-anchored small-pair certificate from
HYP-2095.

## Remaining Gap

This does not prove the full unit-spine exchange lemma.

The audit does not enumerate all simultaneous choices of all nine unit
representatives, and it does not attach the HYP-2097 64 fixed-round class
identity. The dangerous row, if one exists, should be a simultaneous multi-unit
configuration that kills all cheap small pairs while staying in a fixed
self-converse round fibre.

That is a much smaller target than yesterday's vague "large representatives"
problem.

## Next Move

Build the next certificate table with these columns:

```text
fixed self-converse round class
unit shell owner pattern
cheap unblocked pair, if any
if no cheap pair: lost D/U/N obligations under lowering
composite D9/D12 debt
endpoint owner cells
exact witness or descent certificate
```

The operational slogan is:

```text
cheap pair first, exchange second.
```
