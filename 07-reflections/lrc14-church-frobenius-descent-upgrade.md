# LRC14 Church-Frobenius descent upgrade

Source: codex-2026-06-13.  External source: Benjamin Church,
arXiv:2508.14876.  Companions: HYP-2445, HYP-2446, HYP-2443,
HYP-2444, HYP-2463, HYP-2464, HYP-2465, HYP-2469.

## What Changed Since HYP-2445

HYP-2445 imported Church's product-quotient paper as a scalar/support warning:
Shioda supersingularity is not the obstruction; diagonal forms on every
asymmetric partial Frobenius twist are the obstruction channel.

The repo now has a sharper LRC14 analogue.  HYP-2465 shows that Q27 blocking is
not a residue-name problem.  It is a set-cover problem over safe twist
obligations.  In the carry window `1..1092`, any primitive row retaining at
least nine speeds of `7*{1,...,12}` has a Q27 witness.

That means the proof target can be stated much less fuzzily:

```text
any no-Q27 row
=> below-nine-core, outside-window, or named exception.
```

Church gives the proof shape.  HYP-2465 gives the LRC-specific blade.

## The Real Analogy

The bad analogy is:

```text
1092 appears in both places, therefore something mystical is happening.
```

The good analogy is:

```text
scalar quotient passes,
retained side channel survives twists/fibers,
generic bad object descends,
only finite named exceptions remain.
```

In Church:

```text
Q = Shioda supersingularity.
S = diagonal symmetric forms on every asymmetric partial Frobenius twist.
D = rational/elliptic curves descend through partial Frobenius by a projection-degree drop,
    unless they lie in finite exceptional types.
```

In LRC14:

```text
Q = "plain q<=27 blocked" or raw no-witness status.
S = Q27 obligations plus resource coordinates:
    13-clock, deleted-core address, shell-27 class, divisor fiber, owner/Bprime support.
D = a row either enters a certified finite atlas, opens an exception, or loses resource rank.
```

This is why the plain shell is noisy but Q27 is stable.  The plain shell is a
scalar shadow; Q27 remembers enough of the support channel to make the finite
certificates bite.

## Where The Current Proof Stands

The certified part is now substantial:

```text
one-stranger family:       936 rows, all Q27/Bprime covered.
hard replacement hull:     77,520 rows, all Q27 covered.
two-stranger residuals:    877 plain-shell residuals, all Q27 covered.
near-core set-cover:       299 deletion cases through |D|<=3, all infeasible.
```

The two live holes are clean:

```text
below-nine-core:
  what happens when four or more 7-core speeds are deleted?

outside-window:
  what happens to speeds beyond 13*84?
```

That is much better than "LRC14 is hard."  It is a finite descent interface.

## The Below-Nine-Core Thought

Deleting core speeds increases Q27 obligations.  HYP-2465 shows that through
three deletions, `e+1` primitive additions cannot cover them all.  At four
deletions, a naive cardinality MILP may become less informative because five
additions can distribute across more resource types.

So the next computation should not only ask:

```text
can five additions cover all obligations?
```

It should ask:

```text
can five additions cover all obligations while also paying the 13-clock,
divisor-fiber, shell-27, owner, and low-clock invoices?
```

This is exactly the Church upgrade: not one scalar degree, but a retained
multi-channel descent rank.

## The Outside-Window Thought

The `1092` cutoff is not arbitrary in the repo: it is `13*84`, the one-stranger
carry window, and also the Hurwitz order `|PSL2(F_13)|`.  But the proof should
not lean on coincidence.  The real outside-window lemma should be operational:

```text
if v > 1092 appears in a would-be blocker,
then v either opens Bprime(any), dominates/transports an existing core speed,
or reduces through a divisor/carry fiber to a bounded representative.
```

That would play the same role as spreading out and partial Frobenius twists in
Church's proof: a global object is pushed back into the finite arena unless it
is already exceptional.

## What To Try Next

1. Typed `|D|=4` set-cover: split candidate additions by resource role instead
   of allowing an undifferentiated budget of five.
2. Bprime normalizer for `v>1092`: search for domination, owner-private, or
   carry-fiber reductions.
3. Support-transport curvature: compare obligations under `q -> 2q`,
   `q -> 7q`, and `27 -> 9 -> 3`.
4. Exception list as a formal object: AP, Vstar, nonprimitive 2AP, `q=91`,
   `q=161`, owner-private/Bprime, and low-clock exits.

The punchline: HYP-2469 does not claim the proof is done.  It says the proof now
has the right silhouette: finite certified blocks, two named open portals, and
a descent measure to attack them.
