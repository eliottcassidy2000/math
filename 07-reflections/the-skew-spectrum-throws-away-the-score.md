# The skew-spectrum throws away the score

*mac-mini-2026-07-01-S89. Reflection on HYP-3816.*

The owner asked a precise question: does the adjacency spectrum see the flip-rank excess that the skew
spectrum misses? The flip-rank excess is the amount by which the covering number of the tournament cube
overshoots every classical bound — the `0,0,0,1,3` that opus and klein traced to the highly symmetric
classes, the few-labeled-rep needles a thin subcube cannot cover. So the question is really about symmetry:
can each spectrum see how symmetric a tournament is?

The skew spectrum cannot. At six vertices, all fifty-six isomorphism classes collapse into six skew spectra,
and every one of them is a mongrel — the class with a nine-fold automorphism group sits in the same spectral
bucket as classes with three-fold symmetry and with none at all; the five-fold rotational tournament shares
its skew spectrum with six asymmetric ones. The skew spectrum looks at a covering needle and a generic
tournament and sees the same thing. Whatever drives the flip-rank excess, the skew spectrum is structurally
blind to it. The adjacency spectrum, by contrast, is four or five times finer, and — at least through six
vertices — it *determines* the automorphism group: two tournaments with the same adjacency spectrum always
have the same amount of symmetry. It separates the needles. So the answer is clean: yes.

The reason is a single algebraic fact worth keeping. The two matrices carry the same data — the skew
adjacency `S` and the ordinary adjacency `A` differ only by the all-ones matrix, `A = (J - I + S)/2` — but
their spectra couple to that all-ones direction differently. The skew spectrum is converse-even: reverse
every arc, and `S` becomes `-S`, and the spectrum is unchanged. That invariance is bought at a price. To be
blind to the direction of every arc is to be blind to the score sequence, to who beats whom, to the ranking
— and the score sequence is exactly where symmetry becomes legible, because a symmetric tournament wears its
symmetry as a degenerate, patterned score sequence. The adjacency spectrum keeps the score direction: its
Perron eigenvalue leans on the all-ones vector, and through that lean it reads the score degeneracies, and
through those the automorphisms. The skew spectrum, in achieving its beautiful converse-symmetry, throws the
score away, and with it the very thing the covering problem cares about.

This lands squarely on the involution atlas from the session before. The converse — reverse all arcs — is
the fold, and every invariant has a converse-even part and a converse-odd part, the plus and minus
eigenspaces of the fold. The skew spectrum is the converse-even projection, pure and rigid; the adjacency
spectrum is the whole thing, even plus odd, and the odd part is the score. The flip-rank excess is a
symmetry phenomenon, and a symmetry that lives in the score direction is a symmetry the odd part can see and
the even part cannot. It is not that one spectrum is stronger than the other in general — they are two faces
of the fold — but that this particular question asks after something the even face has folded away.

The pattern that transcends the theorem: **a symmetric invariant hides what its symmetry averages over, and
you feel the loss precisely when you ask about the thing it symmetrized away.** The skew spectrum's
converse-invariance is a feature — it is why it appears in the second-moment ladder, why it is the natural
frequency of the antisymmetric part — but a feature that folds out the score direction cannot be used to
detect a score-borne symmetry. When two invariants carry the same data yet answer differently, look at what
each one is invariant under, and the difference will be exactly the coordinate the invariant one gave up. To
see the score, do not use the spectrum that is proud of not seeing it.
