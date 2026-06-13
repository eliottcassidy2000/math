# LRC Small-Denominator Sieve and Slack S388

**Session:** `codex-2026-05-31-S388`
**Script:** `04-computation/lrc_small_denominator_sieve_s388.py`
**Stored output:** `05-knowledge/results/lrc_small_denominator_sieve_s388.out`
**Theorem:** THM-366
**Hypotheses:** HYP-1855, HYP-1856

The useful formalization this session is simple but high leverage:

```text
At level n, a primitive point a/m with m <= n can only be strictly
protected by speeds divisible by m.
```

This is THM-366.  It generalizes THM-360 from the top unit layer `a/n` to
every small denominator `m <= n`.  If no speed is divisible by `m`, then every
unit `a/m` is safe.  For `m<n` it is an open lonely witness; for `m=n` it is a
boundary lonely witness.

So every full-open-cover counterexample must be sieve-complete:

```text
for every m in {2,...,n}, some speed is divisible by m.
```

## What This Explains

The initial segment `{1,...,n-1}` is not mysterious in this language.  It hits
every denominator below `n` and misses exactly `n`.  That is why its surviving
witnesses are exactly the unit boundary points `a/n`.

At `n=14`, the initial segment leaves:

```text
missing=(14,), boundary t=1/14, unprotected=6, coreE=0
```

Dropping `13` and adding `14` pays the top gate but reopens denominator `13`:

```text
missing=(13,), open t=1/13, gap/th=0.045455
```

This is the clean whack-a-mole mechanism behind the recent gate/tightness
duality notes.

## The New Distinction

The important new distinction is:

```text
missing a small denominator gives an immediate rational witness;
satisfying the sieve only moves the problem into endpoint debt.
```

For example, replacing the dropped `13` by `lcm(14,13)=182` is sieve-complete,
but it is not more dangerous:

```text
n=14 drop 13 add 182:
  sieve-complete
  gap/th=0.065934
  unprotected=48
  coreE=0
```

The same pattern appears in the focused largest-repair rows:

```text
n=15 drop 14 add 210: gap/th=0.061905, unprotected=24, coreE=0
n=18 drop 17 add 306: gap/th=0.052288, unprotected=64, coreE=0
```

The transfer tables for `n=14,15,18` reinforce the point: the `lcm(n,d)` rows
remove the explicit missing denominator, but all remain positive-gap.

This is HYP-1855: sieve completion exports endpoint debt.

## How This Connects To Labelled Cycles

THM-365 says a counterexample needs a directed endpoint-protection cycle.
S386 says this is the circular arithmetic analogue of tournament good-cut
protection.  THM-366 now says that any such cycle has to pay every small
denominator gate.

That suggests the next formal object: cycle slack.  If

```text
e = (n*m + eps)/(n*u)
```

is protected by speed `p`, the strict integer condition is

```text
|p*(n*m+eps)-a*n*u| < u.
```

The arrow has positive slack

```text
u - |p*(n*m+eps)-a*n*u|.
```

HYP-1856 is that this slack cannot circulate around a genuine integer-realized
endpoint core without either paying a small-denominator gate or exporting debt
to a higher endpoint layer.  In tournament language, this would be the LRC
version of condensation order: the support graph may look cyclic, but the
labelled slack layer should force a descent or a private endpoint leaf.

## Updated Search Grammar

The proof/disproof branch is now sharper:

```text
not sieve-complete:
  THM-366 gives a rational witness immediately.

sieve-complete:
  inspect endpoint debt, private pivots, critical-radius surplus, and
  labelled cycle slack.
```

This also changes how to rank counterexample attempts.  A set satisfying the
small-denominator sieve is not automatically close to a counterexample.  It
has merely passed the first filter.  The real question is whether the endpoint
incidence matrix has a nonempty core and whether the labelled slack can close
around that core.

Current evidence says no sampled family gets that far.
