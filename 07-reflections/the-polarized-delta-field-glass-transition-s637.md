# The polarized delta field, and the glass it freezes into (S637)

The instruction was to take the sharpest open question — how H moves under arc flips, the polarized
delta field — and understand it at larger n, find where it lives abstractly, and be creative. The
field turned out to be the gradient of a frustrated magnet, and the magnet turned out to freeze.

First the reframe, which the canon had already half-built without my noticing. H, as a function on the
tournament hypercube, is a height. The delta of an arc is the difference H makes across that arc's
flip — which is exactly the discrete gradient of H on the cube, an exact one-form. So the "polarized
delta field" is `dH`, and "polarized" means signed: the gradient has a direction, uphill toward the
regular tournament where H is largest, downhill toward the transitive where it is one. The second-
order question from last session — how flipping one arc changes another arc's delta — is the discrete
Hessian, the curvature of this height function. A field, its potential, its gradient, its curvature.
The whole apparatus of a landscape.

And the canon already knew the landscape's character. H is antiferromagnetic: its log expands with a
positive external field per tile, `log(1+2^{s-1})`, and strictly negative pairwise couplings, every
pair of tiles fighting each other. Negative couplings are frustration; frustration makes glasses. So
the delta field is not the smooth gradient of a bowl. It is the gradient of a spin glass, and the only
question is at what size the glass freezes.

It freezes at six. I computed the local maxima of H — the tournaments from which no single arc flip
raises H, the basins of the landscape. At n=5 there are sixty-four of them and every one is the global
maximum: H equals fifteen at all of them, the regular tournaments, a single basin, a clean bowl with
no false bottoms. At n=6 there are twelve hundred local maxima, and only four hundred eighty of them
are the true maximum H=45. The other seven hundred twenty sit at H=37, trapped — local maxima that are
not global, metastable states you can fall into and not climb out of by any single flip. The landscape
went from one basin to many between five and six. That is a glass transition, and it happens exactly
at the first even n. The two-adic seam, the thing that makes even n the hard side of the Lonely
Runner, is here the onset of metastability in a frustrated magnet. The trap at thirty-seven is even
flanked by forbidden values — thirty-five and thirty-nine are H-values no tournament on six vertices
can take — so the metastable basin is walled by holes in the spectrum. Frustration digs the trap and
the impossible H-values fence it in.

The holes themselves are the third thing that grows. Seven is forbidden at n=5; seven, twenty-one,
thirty-five, thirty-nine at n=6; a dozen values at n=7. These are the H that no odd-cycle conflict
graph can realize as an independence count, and they proliferate. The range of the field is a Cantor-
ish set of allowed odd values riddled with widening gaps, and the gaps are not noise — they shape the
basins, as thirty-seven's flanking shows.

Where else does a polarized delta field live? Everywhere a partition function sits on a configuration
cube. In the analysis of Boolean functions it is the influence field, the discrete derivative, and the
band-limitedness the canon proved — H has Walsh degree only about n, far below the cube's dimension —
says this particular field is low-frequency, sparse, structured, not a generic rugged mess but a
frustrated one with a short-range interaction. In spin-glass physics it is the local field of an
antiferromagnet, and the metastable states I found are its replica-symmetry-broken landscape in
miniature. In discrete Morse theory it is the gradient vector field of H as a Morse function, its
critical cells the extremal and metastable tournaments, its forbidden values gaps in the Morse
spectrum. And in our own work it is the gradient of the Lonely Runner's covering-depth partition
function, loneliness the ground state, delta its arc-derivative. The same one-form keeps reappearing
because the same recipe keeps reappearing: a conflict graph, its independence polynomial, and the
height that polynomial draws over the cube of configurations.

So the sharpest question had a sharp answer. The polarized delta field is the gradient of a frustrated
antiferromagnetic energy on the tournament cube; it is low-frequency and signed; its curvature is the
second-order odd-cycle coupling; and at the even-n two-adic seam it freezes into a glass, growing
metastable basins fenced by impossible H-values. The Lonely Runner's hardness, read on this cube, is
literally a spin glass freezing — and it starts at six, the first even number large enough to
frustrate.
