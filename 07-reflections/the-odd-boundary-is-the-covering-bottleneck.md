# The odd boundary is the covering bottleneck

*klein-2026-07-01-S77. A reflection on HYP-3810 — asking whether a parity obstruction explains a covering
hardness, and finding the two describe the same objects.*

Two threads of the last two weeks had stayed separate. One was combinatorial hardness: the flip-rank, the
minimum-dimensional subcube that covers all isomorphism classes, and its excess over the information floor —
the fact that covering is harder than packing, that at `n=6` you need a seventh bit where six should
suffice. The other was parity structure: the merged metagraph's blue lines forming a T-join whose boundary
is exactly the self-complementary nodes, each carrying an odd count. The owner asked whether the second
obstructs the first — whether the T-join boundary parity is *why* the low-dimensional covers fail on the SC
classes. Testing it, the two threads turned out to be one.

The setup is clean once you see it. The grid-symmetric tilings — the blue ones — are the fixed points of a
linear involution on the tile coordinates, so they are not just a set but a *linear subspace* `W`. Inside
`W`, the self-complementary classes appear as clusters, and every cluster has odd size; the
non-self-complementary classes do not appear in `W` at all. So `W` is partitioned into odd-sized SC pieces —
the odd boundary of the T-join, made concrete as an odd-cluster partition of a vector space. That is the
parity object.

Then I split the covering problem by class type and measured the excess separately. The result was
unambiguous: the excess lives on the SC side. At `n=5`, covering the eight SC classes alone is exactly as
hard as covering all twelve classes, and the SC subset overshoots its own information floor by one while the
non-self-complementary classes overshoot by zero. At `n=6`, the SC subset overshoots by two, the largest of
any subset. The hardness of low-dimensional covering is concentrated precisely on the odd boundary of the
T-join. The parity structure and the covering hardness are not two facts about the metagraph; they are one
fact wearing two descriptions, and both point at the self-complementary classes.

This also reconciles with the symmetry story from the flip-rank thread — that the classes forcing the excess
are the high-automorphism ones, the Paley heptagon and its kin, because high symmetry means few labeled
representatives and therefore few chances to be caught in a thin subcube. Those high-symmetry classes are
self-complementary; the self-complementary classes are the odd T-join boundary; the odd boundary is the
covering bottleneck. Three vocabularies — parity, symmetry, covering — and one set of objects. When several
independent invariants all single out the same subset, that subset is the real content, and the honest thing
is to say so plainly: the SC classes are where the metagraph is rigid, where the parity is odd, where the
symmetry is large, and where the cover must spend its extra dimensions.

I want to be careful about what is proved and what is observed, because the temptation is to say the parity
*forces* the excess. It does not, quite. The parity forces the handshake — that the number of SC classes is
even — and it forces the odd-cluster structure of `W`. But the *size* of the covering excess is not read off
from the parity alone; it is the symmetry-folding, the `|Aut|` arithmetic, and parity only tells you where
to look, not how far. So the precise claim is the modest and correct one: the T-join boundary parity is the
structural *marker* of the covering bottleneck, co-located with it exactly, and a genuine parity *lower
bound* on the SC flip-rank — a Fourier or character argument on the odd-cluster partition of the linear blue
subspace — is the next thing to try. The lesson is that a parity you can prove and a hardness you can only
measure can still be shown to live on the same objects, and that co-location is itself a result: it tells
you the hardness has an arithmetic cause, and it tells you which arithmetic to interrogate next.
