# LRC14 Pi/Unital Flower Carriers

The prompt's pi facts are useful, but only after splitting the quotient.

`22/7` explains the flower observation because `1/pi` is close to `7/22`.
Turning by `1/pi` radians and then reducing modulo one radian gives a clean
model: step `+7` on `Z/22`, with `gcd(7,22)=1`, so all 22 families are visited.
The exact drift after 22 petals is only about `0.002817` radians from 7
radians.  But this is not a full-circle closure.  On the circle the rotation
number is `1/(2*pi^2)`, whose small-denominator behavior is different.

`cuberoot(31)` is the other side: numerically better than `22/7`, but less
directly a flower quotient.  Its useful reading is `31=27+4`.  For the q=3
geometric unital, the parameters are `2-(28,4,1)`, so `27=C=2*14-1` and `4`
is the block size.  That makes `cuberoot(31)` a mnemonic for a cubic shell plus
one block, not a substitute for exact LRC arithmetic.

The word "unital" is a trap unless typed.  Geometric unital means a pair-unique
block design.  Algebraic unital means identity/unit preservation.  Unit groups
are yet another thing.  All three touch LRC14 differently: the geometric q=3
unital suggests a finite pair-coverage carrier, the algebraic sense says
quotient maps must preserve the floor identity, and the unit-group sense
explains the `+7` step on `Z/22` and unit-visible C27 holes.

The carrier tournament puts these below exact `M`/Farey branch, C27 transfer,
bigraded relation signature, and Kpq/K33 incidence.  That feels right.  The
new ideas should help type the p>=3 owner branch, not reorder the proof.

Next useful attempt: after a low-gap atom is already labelled by exact
`M(S)`, C27 transfer, and K33 owner data, try to package the three-owner branch
as a q=3 unital-style pair-unique incidence packet, then see whether that
packet can feed the HYP-2908 forbidden-H7 state lift.
