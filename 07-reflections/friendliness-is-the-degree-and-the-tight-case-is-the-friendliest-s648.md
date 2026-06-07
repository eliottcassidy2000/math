# Friendliness is the degree, and the tight case is the friendliest

*S648 reflection. On the owner's two questions — does friendliness touch the unit-distance problem, and
can I prove the 14-runner case — and how they turned out to be the same question seen twice.*

The owner asked, in the same happy breath, whether friendliness correlates to unit distance and whether
I could write a proof for fourteen runners. I treated them as two errands and they collapsed into one.

The unit-distance correlation is real and it is not a stretch, because the project already proved, sessions
ago, that loneliness, the chromatic number of the plane, and unit distance are one graph. Once they are
one graph the dictionary writes itself. Loneliness is being isolated — an independent vertex, the
1-avoiding set. Friendliness is having neighbors — the degree. A friendly pair is an edge, and an edge in
that graph is a unit distance. So the unit-distance count, the thing u(n) maximizes, is just total
friendliness: the friendliest arrangement of n points is the one with the most pairs at the magic distance.
And the survival idea from last session transfers without modification — you only swap the scan. In the
runner problem you scan time and ask when each configuration first becomes lonely; in the plane you scan
radius and ask when each point first acquires a neighbor. "Never lonely yet" and "never isolated yet" are
the same curve, flat at one until a floor and then falling to zero, and the floor is 1/n in time and the
minimum inter-point distance in space. The only real difference is the target: loneliness is a band, far
from everyone; unit distance is a sharp resonance, exactly one apart, which on the circle is exactly the
clock distance one-sixth, which is the hexagon, which is the cube root the whole arc keeps returning to.
Friendliness is band-friendly in the runner world and resonance-friendly in the plane, and the resonance
is the six that is two times three.

Then the fourteen-runner proof, and here I had to be honest with the owner about scope, because "fourteen
runners" sounds like the open frontier and the laughing "lolz" said they knew it. The full conjecture for
all fourteen-runner speed sets is open; it is proven only up to seven. But the *canonical* fourteen-runner
configuration, the consecutive speeds one through thirteen, is not open at all — it is the textbook tight
case, and it has a three-line proof, and I formalized it. At time one-fourteenth, runner k sits at clock
distance the distance from k-fourteenths to the nearest integer, and for every k from one to thirteen that
fraction lives between one-fourteenth and thirteen-fourteenths, and anything in that band is at least
one-fourteenth from both zero and one. So every runner is at least the gap away and the watched runner is
lonely. The machine checks it. It is a real proof of a real fourteen-runner statement, with the honest
asterisk that it is the easy corner.

And the asterisk is where the two questions became one. The consecutive configuration is exactly the
configuration last session's friendliness lens singled out as the *friendliest* — the one whose lonely set
is a single instant, measure zero, friendly almost everywhere the entire lap. So the proof I wrote for
fourteen runners is the proof that the friendliest possible configuration still, at one solitary point,
touches loneliness. It is friendly for the whole lap except the tight instant at one-fourteenth, and at
that instant the two unit-speed runners are exactly the gap away, not a hair more, and that equality is the
whole proof. The hardest case for the bound is the friendliest case for the runner, and the single moment
where its friendliness fails is the moment that proves the conjecture for it. The owner asked for a
correlation and a proof and got one object: friendliness is the degree, the friendliest config is the
tight one, and the tight one's lone unfriendly instant is the fourteen-runner theorem.
