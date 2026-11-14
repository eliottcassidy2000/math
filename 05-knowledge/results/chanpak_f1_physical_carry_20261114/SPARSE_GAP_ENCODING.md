# Sparse arc encoding for positive multi-contact gaps

This strengthens the quantitative obstruction in the
[fixed-gap ray theorem](FIXED_GAP_RAYS.md). For an apex with
at least r≥2 outcontacts, an ordinary target H≤X and positive
rooted gap d=R₊−H≤D, put p=binom(r,2). The number of possible
ordinary values H is at most

    2^((D/(4p)) log₂ D + O(D)) · log(2X),

with an absolute constant in the O(D). In particular, for any
fixed c<4p, almost every odd H≤X admits no such witness with
d≤c log₂ X/log₂ log₂ X. Taking r=2 includes the whole
multi-contact class and gives c<4. This does not address one
contact or gaps that grow faster than this budget.

## Pivot order: backward arcs cost directed triangles

For a tournament T, choose a pivot vertex v and order the
in-neighbors recursively before v and the out-neighbors
recursively after v. An arc crossing backward between these
two groups forms a directed triangle with v, and all backward
arcs inside either group are handled at a later recursive step.
No triangle is charged twice. Thus this literal order has at
most t(T) backward arcs, where t(T) counts directed triangles.
By the [odd-cycle packing expansion](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14),
t(T)≤(H(T)−1)/2: each single triangle contributes two in
addition to the empty packing.

## Spend one common budget, averaging over every contact pair

Use a normalized finite template G → C → B from the fixed-gap
theorem, after extracting the free transitive suffix. Let m and q
be the path counts of G and B, treating an empty block as one.
Let k≥r be the actual contact count and P=binom(k,2).
The physical defect equations in that theorem give

    m q E_C ≤ d.

For an empty suffix q=1 and the inequality simply drops the
nonnegative full-end term. For a nonempty suffix the expression
is equal to d. Select a pair of core contacts {u,v} minimizing
its contribution E_uv, so E_uv≤E_C/P. Its AND and OR contractions
have path counts summing to at most E_uv; take the smaller,
A≤E_uv/2. With W=C−{u,v}, H(W)≤A, because all odd-cycle
packings of W survive in that contraction with the mark fixed.
Therefore

    m q A ≤ d/(2P),   and
    t(G)+t(W)+t(B) ≤ (m+q+A−3)/2 ≤ (m q A−1)/2.

For the last step use x+y+z≤xyz+2 for positive integers, which
follows by applying x+y≤xy+1 twice. The three pivot orders
together use at most d/(4P) backward arcs, hence at most
D/(4p). Choosing only one contraction or one arbitrary contact
pair would miss this saving.

## Count codes, then retain the exact integer tests

The fixed-gap proof bounds all normalized templates through D
by 3D+1 vertices, before the transitive tail. Encode three
ordered blocks G, W, B, the insertion of two named contacts
in the core, the contact mask, and incident orientations of the
two named contacts. These choices require O(D) bits. Across
blocks the arcs are prescribed. Within the three ordered blocks,
the remaining orientations are encoded by the locations of at
most floor(D/(4p)) backward arcs among O(D²) candidates.

With K=floor(D/(4p)), the logarithm of the binomial sum for K≥1
is at most (D/(4p))·(log₂ D+log₂ p+O(1))+O(log D) by the
standard binomial bound; the extra log₂ p term contributes
O(D) uniformly. For K=0 only the empty choice occurs.
Each valid template supplies at most O(log(2X)) outputs for D≤X,
because the only remaining freedom is a dyadic exponent; when
D>X the bound on the number of targets is trivial. Substituting
D=floor(c log₂ X/log₂ log₂ X) gives X^(c/(4p)+o(1)).

The same code is a search improvement: enumerate sparse arc
sets in place of all core orientations, then use the literal
defect equation and the exact power-of-two test. The encoding
itself is an overcount, not a new sufficiency criterion.
Run **node replay_sparse_gap_encoding.js** for small literal
pivot and pair-budget controls.

The pivot-order approach and a bound with coefficient 1/2
before D log₂ D came from a background pass. Minimizing over
both contractions and then over contact pairs gives the
coefficient 1/(4p) here.
