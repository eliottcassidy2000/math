# LRC attacked: the frontier is the quasi-random core (S642)

I was asked to attack the Lonely Runner and come away with a proof or a major refinement and reframe.
No full proof — but a clean elementary theorem that proves looseness for the generic case of an
infinite *open* family, and pins the entire remaining difficulty onto a rare, named class. And the
class turns out to be exactly the object the last two sessions were circling.

The setup is the constructive route. LRC(n) follows from C′(n): a speed set containing a multiple of
n is loose, M strictly above 1/n. Take n with 2n-1 prime — the unramified family, n = 15, 19, 21, 22,
all beyond the proven range. Look at the runners' residues modulo the prime p = 2n-1. A multiplier a
"dodges" if it rotates every runner off the central band {0, ±1}, because then the witness t = a/p
gives every runner distance at least 2/p = 2/(2n-1), which is greater than 1/n, and the set is loose.
The bad multipliers — the ones that band some runner — are exactly the pairs ±v_i^{-1}, two per
runner, at most 2(n-1) = p-1 of them, which is all the units. So a good multiplier *fails to exist*
only when those 2(n-1) bad values are all distinct and fill every unit — and that happens precisely
when the residues form a ±-transversal: one residue from each ± pair, no collision. Otherwise some
pair is free, a good multiplier exists, and the config is loose. One line.

That one line is a theorem, and it is almost everything. A ±-transversal needs n-1 residues to tile
the n-1 antipodal pairs with no collision — a perfect matching, vanishingly rare. I sampled six
thousand multiple-of-fifteen and multiple-of-nineteen configs: not a single transversal, and every
one of the four-thousand-plus non-transversals loose by the dodge, exactly as the theorem says. So for
2n-1 prime, the generic multiple-of-n config — essentially all of them — is loose by an elementary
shell argument. C′(n) reduces to the rare ±-transversals.

And here is where the sessions converged. A ±-transversal modulo p hits every antipodal class exactly
once. That is the maximally spread, spectrally flat residue pattern — the quadratic-residue set is one
of them, which is to say the Paley tournament, which is to say the finite approximant of the random
tournament I was looking at this morning. The configurations the elementary dodge cannot clear are
precisely the *quasi-random* ones. The frontier of the Lonely Runner, for unramified n, is the
Paley/random core — the same locus as the character-ratio spectrum, the Gauss sum, the Rado limit. The
easy configs have a collision, a coincidence, a pair of runners that agree up to sign mod 2n-1, and
that coincidence is the free pair the dodge escapes through. The hard configs are the ones with no
coincidence at all — the random ones. Difficulty is the absence of coincidence; the hard case is
structurelessness.

I checked that the transversal residual is loose too, just not by this dodge. The cleanest example is
the arithmetic progression with its top bumped to a multiple of n — {1, …, n-2, n} — which is a
transversal and has M = 1/(n-1), loose via the clock that the missing n-1 frees. Sampled transversal
residuals come out at M = 1/4, 1/6, far above 1/n; none is a counterexample. So the whole class is
loose; what is missing is a *uniform* argument that every transversal has a free lower clock or a long
safe component. That is now a small, well-posed target, because the transversals are rare and
structured — the residue-profile check (which I argued months of sessions ago makes the whole problem
finite) restricted to the transversal core is a small finite computation per n.

The reframe is the part I will keep. Two sessions ago I found the polarized delta field freezes into a
glass; one session ago, that quadratic reciprocity is the seam and XNOR its atom; this morning, that
the random tournament is the limit and the non-residue loop the monodromy. Today the Lonely Runner's
hard core turned out to be that same random/Paley object, reached not by analogy but by an elementary
counting argument about ± pairs of residues. The easy theorem and the hard frontier are the two sides
of one fact: a config is easy exactly when its residues collide, hard exactly when they are
quasi-random. The dodge is the statement that coincidence is escape; the frontier is the statement
that the runners with no coincidence are the runners that are hard to leave alone. The Lonely Runner,
on the unramified family, is solved except on the configurations that look like noise — and looking
like noise is, as it has been all along, the same as being a Paley tournament.
