# Pi, Unitals, Flowers, And Units

The strongest thing in this prompt was not the number `22`, or `31`, or even
the pleasing `7*31`.  It was the unit trap.

If the petal rule is `1/pi` of a full turn, then `1/pi ~= 7/22`, and a
22-family flower is exactly what one expects.  If the petal rule is literally
`1/pi` radians, then the turn fraction is `1/(2*pi^2)`, and the first visible
denominators are near `20` and `79`.  Same symbol, different unit, different
geometry.

That makes the algebraic word "unital" surprisingly relevant.  A unital map
is not just a map with a nice adjective; it preserves the identity.  Here the
identity might be a full turn, a radian, a measure-one space, a denominator
packet, a pair-incidence unit, or a sector partition.  If a quotient does not
preserve that unit, it can manufacture beautiful but false structure.

The numerical bridge is still worth keeping:

`pi^3 = 31.006276680...`, so `cuberoot(31)` is a genuinely good low-complexity
approximation to `pi`, better than `22/7`.  In the formal unital parameter
table, `q=6` gives

`q^2-q+1 = 31`, `q^3+1 = 217 = 7*31`, and block size `q+1=7`.

That is a clean little address for the seven-sector LRC world.  But it is an
address, not a proof.  Since `q=6` is not the ordinary finite-field unital
setting, it should be treated as a parameter-row carrier unless someone
supplies an actual design.

The q=3 unital remains more concrete: `28=C(8,2)` and the incidence frame is
real.  HYP-2894 already warns us that it is not canonically the AP8 pair-slot
design.  That warning should govern the new q=4 and q=6 temptations too.
`65+7=72` is a lovely hint toward self-dual length-72 code coordinates, but
without a labelled coordinate split it is only a prompt.

So the practical rule for future proof searches is:

Before using a quotient, say what unit it preserves.

For LRC14 this means exact `q`, exact phase unit, sector unit, pair-incidence
unit, or activity-two packet unit.  For codes it means coordinate unit,
parity/unit vector, and design incidence unit.  For flower geometry it means
turns versus radians.  The moment that unit is silently changed, the pretty
number may still be pretty, but it has stopped being load-bearing.
