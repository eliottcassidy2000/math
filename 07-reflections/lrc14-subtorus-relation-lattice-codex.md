# LRC14 Subtorus Relation Lattice

The useful reframing this session is that spread is not the object.  It is a
shadow cast by the real object: the affine relation lattice of the orbit
`x -> E*x` on the torus.

The smooth minorant identity makes this hard to unsee:

`int G(E*x) dx = (5/7)^k + sum_{n in Lambda_aff(E), n != 0} prod_i psi_hat(n_i)`.

That is the whole large-spread question in one line.  If the relation lattice
has many short primitive vectors, the signed correction can pull the measure
down.  If the relation lattice has only high-product vectors, the product
kernel should not have enough mass to erase the positive main term.

The exact bank made one subtle point very clear.  In one dimension, every
triple is affinely dependent, so merely counting support-3 relations is the
wrong invariant.  The carrier is coefficient product.  A triple `(a,b,c)` has
primitive relation `(b-c, c-a, a-b)` up to gcd, and the Fourier coefficient
sees those numbers.  Dense runs are dangerous because they generate many
small-product triples; raw spread is only indirectly related.

That explains the proof-carrier tournament.  `triple_decay` and `small_triples`
beat additive quadruples, runs, and spread.  Spread was the sink.  This does
not mean spread was useless; it found the door.  But the door opens onto
relation height, not metric diameter.

The live proof shape now feels cleaner:

1. Work with `G`, not the discontinuous max-gap indicator.
2. Split by relation height in `Lambda_aff(E)`.
3. Low height gives finite relation patterns.
4. High height gives a product-kernel tail bound.
5. Then pass the resulting continuous reservoir to the 14-colored CRT
   placement machinery.

This also tells us what not to overclaim.  Triple-only approximations already
failed in prior notes; support-4 and higher signed terms are real.  The
relation-height split should keep the full lattice tail, using low-support
statistics only as diagnostics and finite-pattern generators.

Assumption challenged: I did not use runners or arcs as vertices.  I considered
runner dominance, gap sections, universal centers, residues, Fourier modes,
short relations, subtori, and proof obligations.  For this lane, the quotient
that preserves the predicate is the relation lattice itself; it loses `G_P`
alignment and color-14 placement, so it must be paired with HYP-2593 to become
an LRC14 proof.
