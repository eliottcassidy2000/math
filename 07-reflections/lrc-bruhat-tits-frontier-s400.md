---
source: codex-2026-05-31-S400
status: exploratory formalization
tags:
  - lonely-runner
  - bruhat-tits
  - p-adic
  - endpoint-debt
  - product-law
---

# Bruhat-Tits Frontiers for LRC Endpoint Debt

The Bruhat-Tits angle finally gives a geometric place for "endpoint debt" to
live.

For a prime `p`, the boundary of the `PGL_2(Q_p)` tree is `P^1(Q_p)`.  A
rational endpoint `x=a/b` is a boundary point.  If `p` divides the denominator
more times than it divides the base denominator `n`, then the endpoint has
moved into the infinity cusp of the tree.

That extra depth is

```text
h_p(x;n) = max(0, v_p(b) - v_p(n)).
```

The raw endpoint count is a horosphere population.  The normalized
Bruhat-Tits mass is

```text
sum p^(-h_p).
```

Now the S398 product law has a shape:

```text
ArchGap * raw_debt = (ArchGap * p^h) * BT_frontier_mass.
```

On a clean branch, the visible real gap shrinks because the endpoint frontier
is translated deeper into the p-adic tree.

## Why n=16 Is Clean

The `n=16` dyadic ladders are exactly the single-tree lab.

```text
d=2:   h=1, raw=34,  BT=17,   gap*2^h=2/33
d=4:   h=2, raw=68,  BT=17,   gap*2^h=2/33
d=8:   h=3, raw=140, BT=35/2, gap*2^h=2/33
d=16:  h=4, raw=280, BT=35/2, gap*2^h=2/33
```

So the `35/34` jump is not a failure of the translation picture.  Translation
is perfect: `gap*2^h` stays fixed.  What changes is frontier mass.  The
frontier grows from `17` to `35/2`.

The residue profile says the same thing in a more arithmetic voice.  At `d=8`
the infinity-cusp unit residues modulo `16` split into counts

```text
1,7,9,15       -> 14 each
3,5,11,13      -> 21 each
```

At `d=16` the same pattern simply doubles.  That looks like a stable branch
measure after the phase tax has been paid.

## Why n=14 and n=18 Are Product Buildings

The mixed rows are not as pure.

For `n=14`, the seven-ladder sees `p=7` depth `1`, but the dyadic profile is
already split:

```text
p=2: depths 0 and 2
p=7: depth 1
```

For `n=18`, the `9`-ladder has

```text
p=2: depths 0 and 4
p=3: depths 2 and 3
```

So the correct object is not one tree.  It is a finite subset of a product of
Bruhat-Tits trees.  That is exactly what "adelic product law" should mean here.

## Proof Shape

The Bruhat-Tits proof would not try to chase every real interval directly.
It would prove that endpoint repairs preserve positive frontier mass.

The three outcomes become:

```text
real gap positive        -> ordinary lonely interval;
boundary endpoint debt   -> endpoint witness;
BT frontier mass positive -> debt cannot vanish at infinite p-adic depth.
```

For `n=16`, I would try to prove a dyadic theorem first:

```text
Any all-repaired 16-gate branch has nonzero dyadic frontier mass
after quotienting the base n=16 lattice levels.
```

Then the endpoint-cycle theorem THM-365 becomes a path statement in the tree:
a labelled protection cycle must translate around the cusp without losing
frontier mass.  If it tries to close with zero mass, it needs a cancellation
that should show up as the residue phase tax.

This feels like the right language for the user's sentence.  The real gap is
an Archimedean length.  Endpoint debt is p-adic distance into a cusp.  A
counterexample asks for both to be zero; the tree says the repair moved the
witness, not erased it.
