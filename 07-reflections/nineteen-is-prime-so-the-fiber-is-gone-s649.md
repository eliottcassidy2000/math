# Nineteen is prime, so the fiber is gone

*S649 reflection. On being asked for the full Lonely Runner at nineteen instead of fourteen, and finding
that changing one number changed which half of the project the problem belongs to.*

The owner said: work the full conjecture, but nineteen, not fourteen. It sounds like a small edit — bump
the number by five — and it is the opposite of small, because fourteen is two times seven and nineteen is
prime, and almost everything I had built for fourteen lived on the fact that it factored. Fourteen I could
take apart. The Chinese Remainder Theorem splits the clock into a two-part and a seven-part, and the seven
is where the difficulty hides, and the doubling map on the seven is a little order-three rotation, two
neat three-cycles, the cube root wearing a small hat. I fibered the fourteen-runner problem over a
seven-runner base and the structure cooperated. None of that survives the move to nineteen. Nineteen has
no divisors to fiber over. And the doubling map, which on seven was a tame three-cycle, on nineteen is a
single eighteen-cycle, because two is a primitive root mod nineteen — it visits every nonzero residue
before returning, maximal mixing, no sub-structure to grab. The tool I sharpened for fourteen is exactly
the tool that does nothing at nineteen. So the honest first thing to say is that nineteen is harder, not
easier; the prime case is where the divisor methods go to die.

What I could still do is the corner that does not care about any of this. The consecutive configuration,
speeds one through eighteen, is lonely at one-nineteenth for the same three-line reason the fourteen case
was lonely at one-fourteenth: every runner sits at a fraction of nineteen, those fractions live between
one-nineteenth and eighteen-nineteenths, and anything in that band is at least one-nineteenth from the
nearest integer. I formalized it, and this time I did it for every n at once rather than for a single
number, so the consecutive family is now closed in the library — fourteen, nineteen, and all the rest are
one lemma. The machine checks that the canonical nineteen-runner config is lonely. It is a real proof of a
real nineteen-runner statement, with the same asterisk as before: it is the tight, friendliest corner, and
the full conjecture for all eighteen-speed sets is open and stays open.

The part that made the session feel like more than a relabeling was asking what nineteen is *for*, if not
fibering. And nineteen turns out to be one of the special numbers — just on the other side of the arc.
Seven was special on the divisor side: composite-friendly, the seven-clock, the doubling cube root.
Nineteen is special on the field side. It is a Heegner number; the imaginary quadratic field of square
root minus nineteen has class number one; it is four times five minus one, the rotation field for an
Eisenstein norm of five, the conjectural next rung of the chromatic-number tower after the Moser spindle's
minus eleven — square root minus three to square root minus eleven to square root minus nineteen, the
chromatic number of the plane climbing from three to four to five. And nineteen is the radius-two centered
hexagon, one plus six plus twelve, the Eisenstein lattice's second shell, with two-n-minus-one equal to
thirty-seven, the radius-three hexagon. So when the owner moved from fourteen to nineteen they did not just
pick a harder instance of the runner problem; they moved the problem from the two-adic, divisor, fiber half
of the arc to the cube-root, CM, Heegner half. Fourteen's leverage was its seven. Nineteen's leverage is
its square root of minus nineteen. The two frontier cases the owner has now handed me sit on the two seams
the whole project runs between, and the right attack on the full nineteen-runner problem is not a fiber
that is not there but the cyclotomic depth at the prime, with the Heegner field doing for nineteen what
the spindle field did for the fourth color. I proved the corner and mapped the country; the interior is
still the conjecture, and it is now wearing the cube root instead of the divisor.
