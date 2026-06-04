# The grid disproof and LRC: rotations are the whole story (S624)

The prompt sent me to the Erdős unit-distance problem at n=22, told me to find its deep relation to
the Lonely Runner and to tournaments, and to read the recent disproof of the grid conjecture *in our
context*. I want to record why these are not three tasks but one, and what the disproof is telling us
about where LRC is hard.

Start with the number that made it click. A compact 22-point patch of the triangular lattice has 49
unit distances. The optimum for 22 points is 60 (maybe 61 — that single edge is the open frontier).
Eleven of the sixty edges are simply not available to a lattice. Where do they come from? Rotations.
The triangular lattice is the Eisenstein ring `ℤ[ω]`, the CM field `ℚ(√−3)`, and the rotations that
overlay it with itself and manufacture new unit distances are the modulus-one elements `α/ᾱ` of that
field — angles like `arccos(13/14)`, bounded-height algebraic numbers on the unit circle. The optimal
unit-distance graphs are lattice plus a few of these rotations. The lattice gives you the bulk; the
rotations give you the excess that beats the trivial bound.

That is the entire content of the disproof, scaled up. Erdős conjectured the grid is the best you can
do — `n^{1+o(1)}`. The disproof says no: take a high-degree CM field, harvest its dense supply of
bounded-height modulus-one numbers (pigeonhole over ideal classes), keep the supply growing with an
infinite class field tower whose ramification stays bounded (Golod–Shafarevich, bounded root
discriminant), seed it with a split prime, and project a window of the resulting lattice to the plane.
"Differ by a modulus-one element" becomes "distance one." The grid loses for one reason: too few
rotations. The CM tower wins for one reason: many.

Now read that sentence again with the LRC dictionary in hand, because it is our sentence. The
rotations on the unit circle — the modulus-one algebraic numbers — are the **multiplier orbit**, the
**perspective key** the user has been circling for a year. Roots of unity, the literal lattice
rotations, are the witness orbit `(ℤ/m)*`; they are rare, which is exactly why a single shell is never
enough at the LRC frontier. The disproof's move — give up on roots of unity, use the *dense*
bounded-height modulus-one numbers in a CM field — is the move from one shell to the whole tower. And
the tower is the object I kept hitting at n=14: `2n−1 = 27 = 3³`, the 3-adic shell tower, the first
even case where doubling stops being shell-transitive (THM-407). That failure of transitivity is
ramification. The reason n=14 is the LRC frontier and the reason the disproof needs a *tower* of CM
fields rather than a single field are the same reason: an odd prime power ramifies, and a ramified
tower has structure a single level cannot see.

So the three objects line up exactly. The grid conjecture — "the rigid grid is optimal" — is the
unit-distance twin of "the AP is the unique LRC extremal." Both are false, and false the same way:
last year's sporadic tight LRC configs (the AP is *not* alone) are the small-scale version of this
year's grid disproof (the grid is *not* optimal). In each case the rigid object is beaten by twisting
it with the rotation group, and in each case the rotation group is the CM/multiplier symmetry. The
tournament side is the same tower once more — Cayley–Dickson `ℝ→ℂ→ℍ→𝕆`, the doubling that the canon
has been climbing — with `ℂ` the complex-multiplication step that the unit circle lives on.

What I take from it, as a direction rather than a theorem: the twisted-shell dodge's "free unit
±-pair" (THM-412) is a finite, ramified window lemma. The disproof projects a CM lattice through a
geometric window to count planar resonances; the shell dodge projects the clock through a multiplier
to count loneliness. The honest next step is to make that analogy a statement at `ℚ(ζ_27)` and see
whether the LRC extremal `M = 2/(2n−1)` (THM-415) is precisely the config whose residues form a
CM-rotation orbit — the LRC echo of "ℤ[ω] is the densest small unit-distance lattice." If it is, then
the constructive route to LRC(14) and the grid disproof are running the same engine, and the
bounded-root-discriminant control that makes the disproof uniform in `n` is the thing to port into a
uniform shell-dodge bound. That would be an actual attack, not an analogy.

The smaller size did help, exactly as asked: n=22 is small enough to *see* the eleven missing edges,
and seeing that the missing edges are rotations is what exposes the whole bridge.
