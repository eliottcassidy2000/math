# LRC14 Two-Large Lift Opposition

**Source:** codex-2026-06-19-S25, HYP-2634 / T882.

The useful new fact is that the opposite bounded signs are not mysterious
floating-point drift.  They come from exact low-height integer relations that
survive in one lift and not the other.

The seed family is

```text
S_a = (1, a, 8, a+7, 15, 22),  a=2,...,6.
```

The finite HYP-2632 packet sees `a=2` and `a=4` as the same QR row with weight
`-25U`.  The reciprocal lift does not.  At `H=16`, `a=2` is positive while
`a=4` is negative.

The mechanism is a low-height defect sieve.  Every `S_a` has the universal
positive height-2 relation

```text
(1,1,-1,-1,-2,2),   defect 0.
```

But the `a=4` lift alone has larger negative exact relations:

```text
(-1,3,-1,-1,2,-1),   defect 2a-8,
(-4,4,-1,3,-1,-1),   defect 7a-28.
```

Those vanish only at `a=4`.  Shell totals show the effect immediately: the
`a=4` lift nearly cancels by height `4`, then the `h=8..12` block pushes it
negative; the `a=2` lift keeps a positive low-height reservoir.

Moving the repeated pair up one residue period removes the sign opposition in
the small consecutive-ladder scan.  So the next proof object is not just the
finite residue class or the Legendre character.  It is the integer lift offset
together with its low-height relation-defect polynomials.

For HYP-2633, this suggests a concrete summation-by-parts statistic: delete or
isolate low-height defect-zero motifs first, then prove equidistribution only
for the residual lift discrepancy.  The finite packet is the coefficient layer;
the defect sieve is the lift layer.

LRC(14) remains open.  The progress is a structural explanation of the first
opposite-sign bounded pair and a sharper target for the residue-lift lemma.
