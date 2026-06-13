# LRC three-layer stack audit after S514

S513 proposed that a bad LRC row has to fail several gauges at once:
additive scarcity, endpoint debt, product-sum trouble, A000568 projection
trouble, and nonpeeling pressure.  S514 tests that conjunction directly on
selected n14/n16/n18 hard rows.

The stack does produce real tournament shape.  The full row tournament has:

```text
H=87, c3=10, SCCs=(8,1,1,1).
```

So the three-layer object is not just a sorted list.  It has cycles and a
large core, led by:

```text
n14 double-gate,
n16 gate,
n18 double-gate,
n18 gate.
```

The layer split is the useful part.  Operation-weighted danger and denominator
labels behave very differently:

```text
operation_weighted: H=1, c3=0
denominator_static: H=1, c3=0
operation vs denominator edge flips: 0.53
```

Each layer alone is a ledger, but they disagree enough to create shape once
runner dynamics enter.  The full stack is close to operation-weighted danger
(`0.07` edge flips) and farther from static denominator labels (`0.49` edge
flips), which says the current LRC pair-cell danger is doing most of the
motion.  The denominator stack supplies coordinates, not the moving proof by
itself.

The conjunction audit is the sharper result.  Using median flags for additive
scarcity, endpoint debt, product-sum collisions, A000568 odd survival, and
pressure, the sampled n14 and n16 rows hit four out of five flags, while n18
hits three out of five.  The missing flag is always the same:

```text
pressure(nontrivial SCC) = absent.
```

This is good news, in the very researchy sense of "the failure has a clean
shape."  The add/multiply/A000568 stack is not refuted.  It is telling us that
coarse row-time pressure SCCs are too easy to peel.  The pressure coordinate
has to be the labelled owner-compatible endpoint core from THM-380, not a
plain snapshot SCC.

So the next proof object is:

```text
runner dynamic gauges
+ pair-cell operation danger
+ denominator add/multiply/A000568 labels
+ owner-compatible endpoint pressure core.
```

A true counterexample-shaped row must keep all four layers bad at once.  The
current hard rows do not: their arithmetic and danger texture can look severe,
but the pressure layer still peels.
