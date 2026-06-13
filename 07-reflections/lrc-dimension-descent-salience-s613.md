# LRC Dimension Descent Salience (S613)

The reduction question sharpened immediately once I stopped treating `n` as the
size parameter.  In this proof family the visible size is `n-1` runners, but
the certificate size is largely controlled by the residue modulus

```text
C = 2n - 1.
```

That is why `n=12` can be a cleaner model than `n=11`: `23` is prime and `21`
is not.

## What Changed

The same quotient audit was run for `3 <= n <= 14`: unit-shell rows modulo
`C`, D/N gates, and pair-sum pinch classification.  It recovers the n=14
least-positive quotient exactly:

```text
unit_rows=340928, D/N survivors=27733, strict=27730, floor=3, below=0.
```

The descent from `14` to `11` is not monotone in ordinary size:

```text
14 -> 13: 27=3^3 to 25=5^2, survivor burden drops 27733 to 4883.
13 -> 12: 25=5^2 to 23 prime, survivor burden drops 4883 to 420.
12 -> 11: 23 prime to 21=3*7, survivor burden grows 420 to 3415.
```

That third line is the whole point.  Removing one runner can make the quotient
larger if it changes `2n-1` from prime to composite.

## Salient Structure

The proof should be decomposed by `C` type:

```text
prime C       -> pure unit-shell quotient
p^e C         -> one nonunit channel with depth e
squarefree C  -> several independent nonunit channels
lifted rows   -> carry fiber over the coimage
```

This reads the n=14 work as a depth-3 prime-power problem, not as a generic
13-runner problem.  The n=13 case is the same shape with depth 2.  The n=12
case is the prime reset.  The n=11 case is the squarefree-composite warning.

The carry-fiber theorem from S611 remains genuinely n=14: there is no honest
map `Res_27 -> Res_25 -> Res_23`.  The moduli are not nested.  What descends is
the obligation ledger seen by Yoneda-style probes: D clocks, N clocks, pinches,
owner caps, and floor atoms.

## Collapse Family Connection

The output also answers the user's p0-collapse concern in a satisfying way.
The sporadics are not a side quest:

```text
n=5 floor: (1,3,4,7)
n=6 floor: (1,3,4,5,9)
n=8 floors: AP[6->12] and (1,4,5,6,7,11,13)
n=14 floor: V*=AP[12->24]
```

They are all quotient floor atoms.  The phrase "2n-1 resonances = the
cancellation" now feels literal: the shell coimage decides which resonance
families survive down to the all-order `p_0=0` wall.  Additive-chain structure
is necessary because it supplies relation-lattice cancellation; it is not
sufficient because the cancellation has to land in the right `C` shell strata.

The incoming HYP-2169 mod-3 law fits inside this picture rather than replacing
it: the doubled atom `AP[(n-2)->2(n-2)]` appears exactly when `3|(2n-1)`.  That
explains the doubled sporadics at `n=8` and `n=14`, while the off-law chains at
`n=5`, `n=6`, and the second `n=8` floor atom remain genuine collapse-family
residuals.

## Proof Route

The most useful next reduction is not "prove n=13, then n=12, then n=11" as a
linear descent.  It is:

```text
1. Prove the prime-reset lemma for prime C.
2. Prove the squarefree composite product lemma for independent nonunit strata.
3. Prove the p^e tower lemma for nested nonunit strata.
4. Attach S611 carry-fiber conservativity for C=27.
```

That would put n=14 in the smallest genuinely difficult box: one prime channel,
depth 3, with a carry fiber.

## Assumption Challenge

I considered runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
proof obligations, shell orbits, and n-level summaries as tournament vertices.
For this session, n-level summaries were the correct vertices because the user
question was about proof-size descent.  This choice preserves the LRC predicate
after the D/N/pinch quotient but destroys the actual integer lift data; that is
why S611 remains separate.

## Takeaway

The slogan I would keep:

```text
dimension descent = residue-depth descent + collapse-family floor control.
```

Smaller n is only easier when `2n-1` becomes arithmetically simpler.
