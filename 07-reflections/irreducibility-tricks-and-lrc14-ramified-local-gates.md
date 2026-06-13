# Irreducibility Tricks And LRC14 Ramified Local Gates

The most useful transfer from polynomial irreducibility is procedural.

For polynomials, "irreducible over `Z`" becomes manageable because we almost
never attack it at once. We first strip content. Then we look for a residue
prime where every nontrivial factor split dies. If that residue picture is
ramified or collapsed, we keep valuation data and use Eisenstein/Newton slopes.
If the row still survives, we use evaluation values, factor-capture, Cohn/Perron
dominance, or recombination constraints.

The LRC14 proof route now has the same shape.

Q27 is the residue face. HYP-2470 shows that Q27 alone is almost sufficient for
rows retaining eight core speeds, but it has two finite exceptions. The new
diagnostic says those exceptions are not random: both have 12 of 13 speeds in
the 7-ideal and exactly one primitive non-7 speed, and that primitive escape is
13-clock. That is the LRC version of a ramified Eisenstein packet: the residue
test does not empty the survivor set, but the valuation/carry face points to
the missing shells `31` and `33`.

So the next proof should not merely scan farther. It should prove a ramified
portal lemma:

```text
12 speeds in 7Z + one primitive 13-clock escape
=> q=31/33/41 witness or Bprime opening.
```

That would make HYP-2470 feel like a theorem rather than a census.

The concurrent Q31-fiber and dilated-band covering work slots neatly into the
same picture. A shell witness is an uncovered unit after covering `(Z/q)^*` by
the dilated danger bands `v^{-1}B_q`. That is the complement of the set-cover
language used here. HYP-2480 says "the added speeds cannot cover all the right
obligations once the ramified layer is kept"; the band-cover lens says "the
units must leak at the next band/fiber layer." The one-stranger evaders leaking
at `40/41` and the four-deletion exceptions opening at `31/33` are now two
faces of the same instruction: keep the first band-2/fiber layer as valuation
data.

The assumption challenge also sharpened. The tournament vertices are not
runners. For this session the proof-bearing vertices were local proof tricks
and denominator obligations. That preserves what matters: whether a hidden
lift survives. It destroys continuous timing geometry, so it must be paired
with Bprime/positive-measure checks, but it is the right quotient for this
piece of the proof.

The high-leverage path now looks like:

```text
set-cover dual certificates
-> ramified 7/13 portal lemma
-> below-eight-core survivor ledger
-> outside-window Cohn/Perron normalizer
-> Church-style descent theorem.
```

This is a real narrowing of the LRC14 proof search. The analogy with
irreducible polynomials is not decorative; it tells us which surviving local
state to keep when the scalar gate fails.
