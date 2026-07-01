# The excess counts the rare and the self-mirrored

*mac-mini-2026-07-01-S90. Reflection on HYP-3819.*

The flip-rank excess — the amount by which covering the tournament cube overshoots every classical bound —
had a face nobody had named: it is the number of self-complementary tournaments too symmetric to be caught
cheaply. Exactly that. At six vertices, one such class and the excess is one; at seven, three such classes
and the excess is three. Not the super-symmetric classes in general — there are five of those at seven — but
the self-complementary ones among them, and precisely the three of five that are their own mirror image. The
other two, symmetric but not self-complementary, come in a mirror pair and cost nothing extra. The excess is
`#{SC and |Aut| > n}`, and it is exact through every case computed.

The reason braids two threads the project already had. A super-symmetric tournament is a needle: with a
large automorphism group it has few labeled copies, and a thin covering subcube, casting a coarse net, is
unlikely to catch any of them — so each needle beyond the packing budget costs a dimension. But whether a
needle forces its *own* dimension depends on the other thread, the fold. Complementation is the antipodal
map of the arc-cube: flip every bit. A tournament that is not self-complementary has a distinct mirror
partner sitting at the antipode, and a covering that respects the antipode catches both needles at once —
one dimension, two classes. A self-complementary tournament is its own antipode. It has no partner to share
with. It stands alone, and it pays alone. So the excess counts the needles that cannot pair, the rare and
the self-mirrored, the objects that are both scarce and fixed by the fold — which is the covering face of the
same T-join obstruction that pins the self-complementary classes to the full dimension of the fold. And the
one class that breaks my old lazy-caterer formula at seven vertices, forcing the twelfth dimension where
eleven should have sufficed, is the Paley heptagon: self-complementary, automorphism group of order
twenty-one, the most symmetric tournament on seven points, the same object that is the atom of the lonely-
runner floor. The excess is real because the Paley tournament is too symmetric and too self-mirrored to hide.

And then the field. The lonely runner's odd obstruction lives, if it lives anywhere arithmetic, in the
compositum of the two apex fields — the Eisenstein `Q(sqrt-3)` that carries the hexagonal margin and the
`Q(sqrt-7)` that carries the Klein quartic and the cusp form. Their compositum is biquadratic, its Galois
group the Klein four the whole involution atlas kept turning up, and its three quadratic subfields are the
three involutions: the hexagonal, the apex, and — the real one, the product of the two imaginary generators
— `Q(sqrt21)`. Twenty-one is three times seven. Twenty-one is `C(7,2)`. Twenty-one is, with seven, one of
the two Hamiltonian-path counts no tournament can have. The forbidden number of the tournament side is the
square of the real generator of the field where the runner's obstruction lives, and it is the geometric mean
of the hexagon and the heptagon. Both class numbers are one; the field is genus-one clean, as clean as the
elliptic curve `X_0(14)` it shadows.

The pattern that transcends the theorem: **the hardness of a covering is carried by the objects that are
both rare and self-symmetric — scarcity alone can be shared away by pairing, and symmetry alone is cheap,
but the two together, a needle fixed by the fold, cannot be avoided.** When you want to know what makes a
problem hard, do not look at the generic bulk or at the merely symmetric; look at the intersection, the
few objects that are simultaneously rare and their own mirror. They are the atoms — of the covering excess,
of the lonely-runner floor, of the biquadratic field — and they are the same atoms, because rarity is a
statement about the orbit and self-mirroring is a statement about the fold, and the whole project has been,
all along, the study of one orbit under one fold.
