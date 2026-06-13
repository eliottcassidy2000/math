---
source: codex-2026-06-03-S571
status: synthesis + bounded audit
tags:
  - lonely-runner
  - antipodal-witnesses
  - summand-graph
  - multiplication
  - addition
  - odd-even
  - unit-action
  - nonunit-hole
  - spectral-gap
---

# The `2n-1` Antipodal Witnesses Are the Summand Graph Acted on by Units

The S553 witness family

```text
t_k = k / (2n - 1)
```

is not a separate trick from the summand graph.  It is the summand graph at the
odd node

```text
C = 2n - 1
```

with the multiplicative unit group acting on its incoming additive pairs.

This gives a clean dictionary:

```text
summand node C                 <->  modulus C = 2n - 1
incoming pair {a, C-a}          <->  antipodal residue pair
missed incoming pair            <->  possible 2/C witness
unit a mod C                    <->  invertible clock k = a^{-1}
nonunit a                       <->  composite-modulus blind spot
perfect transversal             <->  one choice from every summand pair
AP                              <->  all-lower transversal
```

The companion audit is:

```text
04-computation/lrc_antipodal_summand_units_s571.py
05-knowledge/results/lrc_antipodal_summand_units_s571.out
```

A follow-on bounded audit checks the second-gap separator:

```text
04-computation/lrc_second_gap_transversal_audit_s572.py
05-knowledge/results/lrc_second_gap_transversal_audit_s572.out
```

In primitive boxes through `k=6`, every row below `2/(2n-1)` is already
`n`-clock tight with `M(S)=1/n`, and every such row is a perfect antipodal
transversal.  The bounded flip-set menu is AP plus the known `{2}` sporadics,
so raw sumset minimality is not the exact separator; floor-tight transversal
structure is sharper.

## The Core Mechanism

For `C=2n-1`, a runner with residue `s` at time `k/C` has collar below `2/C`
only if

```text
k s == 0, +1, or -1 mod C.
```

If `a` is a unit and a speed set misses the antipodal pair `{a,C-a}`, choose

```text
k = a^{-1} mod C.
```

Then multiplication by `k` moves that missed additive pair to `{+1,-1}`, the
only nonzero bad residues.  So every runner avoids the bad pair and

```text
M(S) >= 2/C = 2/(2n-1).
```

Read this as:

```text
addition creates the pair shell {a,C-a};
multiplication by a unit moves the shell to the observer.
```

So the S553 lemma is a two-operation statement.  Addition supplies the missing
summand pair; multiplication supplies the inverse clock that makes the missing
pair visible.

## Why This Looks Like the Summand Graph

In the ordinary summand graph, node `C` has incoming pairs

```text
{1,C-1}, {2,C-2}, ..., {floor((C-1)/2), ceil((C+1)/2)}.
```

For odd `C`, these are exactly the antipodal pairs of the cyclic group
`Z/CZ`.  There is no midpoint.  The pair count is

```text
(C-1)/2 = n-1.
```

That is why the LRC speed count also fits perfectly: `n-1` runners can hit
every one of the `n-1` additive pair shells exactly once.  The residual after
S553 is therefore a perfect transversal: one speed in each summand pair.

The AP is not merely "consecutive."  It is the all-lower choice:

```text
{1,2,...,n-1}
```

from the `n-1` pairs `{a,C-a}`.  A flip-set is literally a choice to replace a
lower summand by its upper complement.

## Multiplication Versus Addition

Addition is dense and pair-complete:

```text
a + b = C
```

gives every antidiagonal shell.  Its one-input shadow is just the transitive
order `a<C`.  That is why addition creates many candidate pinch denominators
but, by itself, does not know which time exposes a witness.

Multiplication is sparse and selective:

```text
k * a == 1 mod C
```

exists only when `a` is a unit.  It turns an additive shell into a clock.  This
is exactly the divisibility/sieve side of the repo's add/multiply stack:

```text
addition:       candidate shells / pinch denominators / sumset room
multiplication: unit inverses / divisibility tests / CRT visibility
```

So a missing summand pair is useful only if multiplication can invert it.  This
is the precise sense in which multiplication tests the witness that addition
creates.

## Odd Versus Even

Odd `C` is clean:

```text
a -> -a
```

has no nonzero fixed point.  Every nonzero residue belongs to a two-element
antipodal pair.  The summand graph and the cyclic antipodal quotient are the
same object.

Even `C` has a fixed midpoint:

```text
C/2 == -C/2 mod C.
```

But the distinct-summand graph excludes the midpoint pair

```text
C/2 + C/2 = C.
```

This is the same shape as the LRC apex/half-turn obstruction: the even node has
a self-paired shell that is not a genuine distinct pinch pair.  Multiplication
by `2` creates the parity quotient, but it also creates the degenerate midpoint
that must be handled outside the odd antipodal story.

In short:

```text
odd  = free antipodal pairing, all additive shells are genuine
even = midpoint/apex fixed point, one additive shell degenerates
```

This is why `C=2n-1` is so natural for the spectral gap: it deliberately moves
to an odd completion where the midpoint is absent.

## Composite Odd Moduli: The Nonunit Hole

If `C` is prime, every antipodal pair is unit-visible.  S553's link sees every
missed pair.  The audit shows:

```text
C=11: 5 summand pairs, 5 unit-visible pairs, one unit orbit
C=13: 6 summand pairs, 6 unit-visible pairs, one unit orbit
```

If `C` is composite, the additive shells split by gcd strata.  Multiplication
by units preserves gcd, so it cannot move a nonunit missed pair to `{+1,-1}`.
For `C=27`:

```text
unit-visible pairs: 9
gcd-3 pairs:        {3,24}, {6,21}, {12,15}
gcd-9 pair:         {9,18}
```

The n=14 sporadic tight set

```text
V* = {1,2,3,4,5,6,7,8,9,10,11,13,24}
```

has exactly the composite-modulus fingerprint:

```text
missed_nonunit = {12,15}
doubled        = {3,24}
```

So `V*` is not mysterious.  It is an additive transversal defect that lives
inside the gcd-3 unit-action orbit.  The unit witnesses cannot see it because
`12` is not invertible modulo `27`.

Likewise the n=8 sporadics at `C=15` miss only nonunit pairs:

```text
{6,9} or {3,12}.
```

This explains the exact limitation of the antipodal-transversal reduction:

```text
prime C:     missed pair -> unit inverse -> witness
composite C: missed unit pair -> witness;
             missed nonunit pair -> residual gcd-stratum problem
```

## What This Says About `n=14`

For `n=14`, the spectral-gap modulus is

```text
C = 27 = 3^3.
```

The hard residue space is not all `13` antipodal pairs equally.  It is a
three-layer unit-action stack:

```text
gcd 1: 9 visible pairs
gcd 3: 3 nonunit pairs
gcd 9: 1 nonunit pair
```

The AP is the all-lower perfect transversal.  `V*` is a one-defect
non-transversal that shifts mass inside the gcd-3 layer.  This suggests that
the real reduced proof target is:

```text
prove the gap on perfect transversals;
prove nonunit gcd-stratum defects are boundary-only or clear by a second clock.
```

The second clock should probably be a lift that restores invertibility, for
example working modulo `p*C` or mixing the `n`-clock with the `2n-1` clock.  The
repo already hints at this: all known tight sets, including nonunit-hole
sporadics, are still certified at the `n`-clock `j/n`.

The S572 audit adds a useful warning: sumset minimality is sufficient for AP
tightness but not necessary for all floor-tight rows.  Some `{2}` flip-set
sporadics have positive sumset excess and still sit at the floor.  So the
transversal/flip-set coordinate is closer to the equality stratum than raw
Freiman excess.

## Summand Graph Similarities Worth Keeping

1. **Pair count parity.**  Summand node `C` has `(C-1)/2` pairs when `C` is odd,
   and `C/2-1` pairs when `C` is even because the midpoint pair is excluded.
   This is the same odd/even split as free antipodal pairing versus apex.

2. **AP as lower transversal.**  In both the summand graph and LRC residues,
   `{1,...,n-1}` is the canonical lower-half choice at odd node `2n-1`.

3. **Flips as upper summand choices.**  The S553 flip-set is not an arbitrary
   bitstring.  It is the record of which summand pairs choose the upper
   complement.

4. **Addition gives abundance, multiplication gives visibility.**  The sumset
   creates many shells; the unit group decides which shells can be clocked.

5. **Composite moduli create hidden layers.**  Nonunit summand pairs are present
   additively but invisible to multiplicative inversion.  This is exactly where
   sporadic tight families live.

## Tournament Analysis / Assumption Challenge

The vertices should be antipodal summand-pair shells, not runners:

```text
P_a = {a, C-a}.
```

The pairwise relation should not be a runner comparison.  A clean quotient is:

```text
P_a ~ P_b if gcd(a,C)=gcd(b,C)
```

with unit-action orbit fingerprints:

```text
orbit size, gcd stratum, missed/doubled occupancy, lower/upper flip
```

If a tournament is needed, use a loss-aware one:

- **Vertices:** antipodal pair shells `P_a`.
- **Pairwise observable:** which shell is more visible to multiplicative clocks,
  ordered by `gcd(a,C)` and orbit size.
- **Switch/gauge:** orient `P_a -> P_b` when `P_a` has smaller gcd, then larger
  unit orbit, then lower representative.
- **Tie Hamiltonian path:** increasing gcd, then increasing lower representative.
- **Fingerprints:** score histogram by gcd stratum, SCCs before/after quotient,
  edge flips after marking missed/doubled occupancy.

This tournament will be mostly transitive; that is the point.  Its failures or
marked edge flips are precisely the nonunit-hole data.  A runner-vertex
tournament would destroy the predicate we care about: whether a missed additive
pair is multiplicatively visible.

## Short Version

The `2n-1` witness family is the summand graph plus unit multiplication:

```text
addition builds antipodal pair shells;
multiplication by units turns visible missed shells into witnesses;
oddness removes the midpoint;
evenness creates the apex;
compositeness creates nonunit holes.
```

That is the clean addition/multiplication and odd/even understanding the repo
has been circling.
