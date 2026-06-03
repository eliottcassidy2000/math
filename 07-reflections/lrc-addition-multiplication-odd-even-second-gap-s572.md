---
source: codex-2026-06-03-S572
status: synthesis + exact bounded audit
tags:
  - lonely-runner
  - second-gap
  - antipodal-transversal
  - summand-graph
  - addition
  - multiplication
  - odd-even
  - n-clock
  - spectral-gap
---

# Addition, Multiplication, Odd, Even: What Each One Actually Does in the `2n-1` Gap Picture

**S573 correction.** The floor-only conclusion below `2/(2n-1)` was true only
for the original S572 small boxes.  The expanded S573 clock-blocker audit finds
open-gap lifted rows, including `(1,5,6,11,16,17)` with `M=5/33` at `n=7` and
two `n=8` rows with `M=3/23`.  The corrected theorem target is the
three-clock blocker ledger in HYP-2088, not a global claim that every
below-edge row is floor-tight.

The repo now has enough pieces to say something sharper than "AP is hard."

For the spectral-gap route around

```text
1/n  <  ?  <  2/(2n-1),
```

the four words in the user's prompt have distinct jobs:

```text
addition        = builds the odd summand shells {a, 2n-1-a}
multiplication  = decides which shells are visible by inverse clocks
odd             = removes the midpoint and gives genuine antipodal pairing
even            = reappears only at the floor, as midpoint/apex degeneracy
```

That division of labor is clearer than the older "AP versus non-AP" phrasing,
and S572 gives a bounded exact audit in its favor.

## The clean object is the odd node `C = 2n-1`

S571 already says the S553 witness family

```text
t = k / (2n-1)
```

is the summand graph at the odd node `C=2n-1`, acted on by multiplicative
units.  That means:

```text
addition:       shells P_a = {a, C-a}
multiplication: inverse clocks k = a^{-1} when gcd(a,C)=1
oddness:        no midpoint shell
```

So the off-floor part of the spectral-gap story is not fundamentally about
the physical circle first.  It is about whether the speed residues cover the
odd summand shells, and whether the missed shell is unit-visible.

## What S572 adds

`04-computation/lrc_second_gap_transversal_audit_s572.py` scanned the same
bounded primitive boxes used in the S570 witness-or-core audit:

```text
k=3, max_speed=20
k=4, max_speed=16
k=5, max_speed=13
k=6, max_speed=11
```

and asked only about the rows with

```text
M(S) < 2/(2n-1).
```

The result was rigid:

1. every such bounded row is already `n`-clock tight, with `M(S)=1/n`;
2. every such bounded row is a perfect antipodal transversal mod `2n-1`;
3. the bounded flip-set menu is exactly:

```text
AP                 -> empty flip-set
sporadic floor rows -> flip-set {2}
```

In particular, the true separator below the second floor does not look like
"sumset-minimality" alone.  One of the bounded sporadics,

```text
(1,3,4,5,9),
```

still lies below `2/(2n-1)` while having positive sumset excess.  So AP-style
small sumset is sufficient for floor-tightness, but it is not the right
characterization of the whole sub-gap branch.

The right bounded characterization is:

```text
sub-2/(2n-1)  =>  already floor-tight  =>  perfect odd-shell transversal.
```

## Why this answers addition versus multiplication

Addition and multiplication are not rival explanations.  They control different
cuts of the same object.

### Addition

Addition supplies the shell geometry:

```text
{1,C-1}, {2,C-2}, ..., {n-1,n}.
```

This is the summand graph side.  It explains:

- why there are exactly `n-1` antipodal pairs at modulus `2n-1`,
- why AP is the all-lower transversal,
- why flip-sets are the natural coordinates on the residual branch.

### Multiplication

Multiplication supplies visibility:

```text
k*a == 1 mod C.
```

This is what turns a missed shell into a witness clock.  If the shell is a unit
shell, multiplication can bring it to `{+1,-1}` and expose the `2/C` witness.
If the shell is nonunit, multiplication cannot see it inside `(Z/CZ)^*`, and
that is the composite-modulus hole.

So:

```text
addition says what shell is missing;
multiplication says whether the missing shell is witness-visible.
```

## Why odd and even are not symmetric here

The `2n-1` modulus is odd on purpose.

For odd `C`, every nonzero residue belongs to a genuine antipodal pair.  There
is no midpoint.  So the summand graph and the antipodal shell picture line up
perfectly.

For even moduli, the midpoint shell degenerates:

```text
C/2 == -C/2 mod C.
```

That is the old apex/half-turn defect again.  Distinct summands exclude the
midpoint pair, and LRC floor witnesses at `t=j/n` live exactly where this even
midpoint/apex issue comes back.

So the picture is:

```text
odd shell world (2n-1):   off-floor margin story
even clock world (n):     floor-tight equality story
```

That is the conceptual reason the proof wants both clocks:

```text
2n-1 clock for the gap,
n-clock for the floor.
```

## The improved proof route

The old slogan was:

```text
AP is the wall.
```

The better slogan is:

```text
perfect floor-tight transversals are the wall.
```

That is a real upgrade, because it explains the known sporadics instead of
treating them as noise.

The route now looks like:

1. Use S571:

```text
missed unit shell -> 2/(2n-1) witness.
```

2. Isolate the residual:

```text
perfect transversals
plus composite-modulus nonunit-hole defects.
```

3. Use S572's bounded evidence and the old S553 data:

```text
in the original small boxes, the only rows below 2/(2n-1)
were already n-clock-tight floor rows.
```

4. Finish by proving that any non-floor residual clears either by:

```text
a unit-shell witness,
a second clock on the nonunit hole,
or an endpoint-core contradiction.
```

## Short version

If I compress the whole session to one line, it is this:

```text
In the original S572 boxes, below the second floor was not "small sumset";
it was "already floor-tight perfect transversal."
```

And that explains the roles:

```text
addition       -> shell coverage
multiplication -> shell visibility
odd            -> honest antipodal pairing
even           -> midpoint/apex floor defect
```

That is the cleanest conceptual split I have for the `2n-1` route so far.
