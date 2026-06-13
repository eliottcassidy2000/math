# Lonely Runner Scalar Excision Formal Session S371

This session converted the scalar-ramp obstruction from a computational warning
into a theorem and a sharper finite target for the fourteen-runner case.

## Formal Step

THM-364 proves the scalar-ramp cell blocking identity.  If

```text
v_i = m i mod n,
```

then on any open micro-staircase cell with pattern

```text
b_i = floor(n*{i alpha}),
```

the shifted residues

```text
s v_i + b_i mod n
```

are exactly the initial-segment floor vector at the shifted time
`alpha+s*m/n`.  If all coordinates avoided `0` and `n-1`, THM-358 would force
the shifted time to be a unit endpoint `a/n`; but then `alpha` itself is an
`i=1` breakpoint, contradicting that the cell is open.

So scalar ramps are not exceptional bad examples for a micro-staircase lemma.
They are the Dirichlet equality case in residue coordinates, and any generic
witness statement must excise them first.

## Computation

`04-computation/lonely_runner_scalar_excision_s371.py` rebuilds the S364
`n=14` full cell system:

```text
patterns = 812
candidates = 14*812 = 11368
```

It reconstructs one representative open interval for each floor pattern and
checks the scalar midpoint identity for all shifts and all scalar multipliers.
The result is zero midpoint failures.

The same script scans all one- and two-coordinate deformations of scalar ramps.
No non-scalar full blocker appears in those neighborhoods.  The best one-defect
vectors all miss exactly `56` candidates.  The two-defect neighborhood is worse:
its best vectors miss `112`.

## The 56-Cell Target

The S364 best `n=14` non-scalar vector

```text
(8,2,10,4,12,13,0,8,2,10,4,12,6)
```

is not a generic near-blocker.  It is the scalar ramp `m=8`

```text
(8,2,10,4,12,6,0,8,2,10,4,12,6)
```

with one defect:

```text
coordinate 6: 6 -> 13.
```

Its `56` missed cells have a clean structure:

```text
shifts: odd s only, 8 cells for each s in {1,3,5,7,9,11,13}
widths: 14 cells each at 1/728, 1/882, 1/1176, 1/1386
minimum residue margin: always 1
base unique blocker: coordinate 6 for all 56
```

Thus the 56 cells are exactly the scalar cells uniquely protected by the
coordinate that was changed.  This is a better proof target than "search found
no blocker": prove that any non-scalar deformation of a scalar ramp exposes at
least one uniquely protected cell.

## Next Moves

1. Turn the radius-1/radius-2 scalar-neighborhood evidence into a general
   fragility lemma for scalar-ramp deformations.
2. Separate unit scalar ramps from nonunit ramps in the prime-grid lifted
   variables; nonunit ramps should be quotient/descent cases under THM-360.
3. Extend the exact search beyond scalar neighborhoods.  A global backtracking
   proof should choose exposed cells as constraints rather than search raw
   residue vectors.
4. Compare the four missed-cell widths with the public verifier's
   `I(13,p,1)` bad-set construction: those denominators may suggest the first
   useful pruning threshold for the initial cover search.

