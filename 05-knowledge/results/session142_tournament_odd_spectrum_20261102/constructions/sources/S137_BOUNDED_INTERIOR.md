# Fixed height and interior budget have an exact dyadic normal form

**Status: PROVED structural normal form and inverse; FINITE-EXACT illustrative
base bank.** This strengthens the zero-interior face result in `THEORY.md`.
The unrestricted odd tournament Hamiltonian-count spectrum remains **OPEN**.
No literature-priority claim is made.

Use the actual source and coordinates `(h,F,q,g)` and `(C,L,g,t)` from
`THEORY.md`. A graph/source normalization below may reorder singleton SCCs;
it is NOT an induced deletion-only operation and is not claimed to preserve
old labels or arbitrary later re-marking of old vertices. It preserves the
entire four-state response and has an explicit actual graph recipe.

## 1. A selected-pair bound on internal singleton padding

The inherited S124 canonical weights satisfy

    L=sum_({u,v} subset S) b_uv.

The S136 positive path-cut kernel is

    R_uv=sum_P H(A-V(P)),

where P ranges over simple directed paths from u to v with a positive even
number of vertices. Its cut response equals the canonical weighted-cut
response. Comparing the coefficient of the product of the two cut-indicator
variables gives the exact identity

    b_uv=R_uv+R_vu.                                      (1)

This uses `session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md`,
Sections 1-2, and the S124 cut/interior identities pinned in `THEORY.md`.
It does not identify the whole directed R network with the orientation-gauged
canonical network; their circulation can differ.

Suppose selected vertices u,v lie in different SCCs, u before v, and r
singleton SCCs lie strictly between their SCCs. For every even-cardinality
subset of those r vertices, u followed by that subset in dominance order
followed by v is an even-cardinality directed path. If r>=1 there are exactly
2^(r-1) such paths. Every complementary H is at least one, so

    L>=b_uv>=2^(r-1),
    r<=1+floor(log_2 L).                                 (2)

The intermediate singleton vertices may themselves be selected; (1) places
no cut restriction on the interior of its paths. For r=0 no logarithmic
inequality is needed. If L=0 there cannot be two selected vertices at all,
by positivity of b_uv.

## 2. Unselected suffix blocks commute with a zero singleton

For any nonempty tournament B seen by x with empty cut, write
H=H(B) and E=the number of ordered two-path covers of B, allowing empties.
The inherited ordinal-block action specializes to

    (h,F,q,g) -> (Hh, HF+(E-H)q+Hg, Eq, 0).               (3)

For the singleton this is S0, with H=1,E=2. Direct multiplication shows that
(3) commutes with S0. Both orders give

    (Hh, HF+(2E-H)q+Hg, 2Eq, 0).                         (4)

This is a complete state identity on actual sources, not only scalar H.
It is implemented by reordering actual homogeneous SCC blocks after the last
component meeting S. No rearrangement of a selected block is justified by
(4).

## 3. Finite bases and the exact normal form

Put

    B(h,0)=h+1,
    B(h,L)=h+3+floor(log_2 L), L>=1.                      (5)

**Theorem.** Every actual source of height h and interior coordinate L has
the same four-state vector as either:

1. an actual source with parent order at most B(h,L); or
2. a source obtained from such an actual base by a positive number k of S0
   appends.

The base has the same h,L. In the second case, if its native state is
`(h,F0,q0,g0)`, then

    F=2^k q0+h-2L,
    q=2^k q0, g=0,                       k>=1.           (6)

The empty-cut case has F=h,q=g=L=0 and is constant; its small representative
is included separately. For all other cases q0>0.

**Proof.** Write A as its ordered SCC sum and first suppose S is nonempty.
Let B_first and B_last be the first and last SCC meeting S. Delete every
unselected singleton SCC before B_first: all those components have cut empty,
and their height-one prefix action is neutral before the first selected
component. Retain all nonsingleton SCCs.

Every singleton SCC after B_last is unselected. By (4), move all of them past
the other unselected suffix blocks into one final S0 tail. Remove that tail
for the base. This reordering and removal preserves h,L, and restoring the
actual S0 tail restores the entire source state.

If L=0, S is a singleton. Apart from the nonsingleton core, the base retains
at most the selected vertex when its SCC is singleton. If L>=1 and B_first
and B_last differ, all singleton components strictly between them are
bounded by (2), and at most two singleton endpoint components remain. If the
first and last selected component agree, it is nonsingleton (L>=1 requires
two selected vertices), so no singleton components remain in that interval.
Thus the base has at most its nonsingleton-core order plus1 when L=0, and
plus3+floor(log2 L) when L>=1.

The elementary inherited SCC bound from `THEORY.md` gives core order at most
h: each nontrivial strong component has H at least its order, and the sum of
its positive counts (each at least3) is at most their product h. This proves
(5). If S is empty, delete singleton SCCs and use a one-vertex parent if no
core remains; F=h is preserved and the same bound holds.

Finally the inherited S0 recurrence gives, for k>=1,

    F=2^k q0+F0+g0-q0.

The identity defining L gives F0+g0-q0=h-2L. This proves (6). Every normal
form uses actual graph deletion/reordering and a literal S0 tail. Conversely
each actual base and each such tail is a valid source with the same h,L. QED.

A padding sidecar can retain the deleted prefix singleton positions and the
original positions of the moved suffix singletons. It reconstructs the old
SCC word and hence the original labelled graph when labels are retained.
It is not needed for four-state right-local responses, but must not be silently
dropped if later operations refer to a particular old vertex.

## 4. Exact target inverse and growth consequence

For fixed h,L, enumerate all parent tournaments of order at most B(h,L),
all cuts, and retain exactly those with the requested h,L. This finite actual
bank, together with (6), decides the scalar source target N at all orders:

- accept the count F0 of any base;
- or accept when, for some base with q0>0,

      (N+2L-h)/q0 = 2^k for an integer k>=1.

A successful test yields the explicit base graph followed by k actual S0
vertices. The divisibility, positivity and positive-tail condition are all
necessary. No arbitrary numerical vector is admitted as a base. The algorithm
is not claimed efficient in h,L or their bit lengths; the point is an exact
finite bank and a dyadic rather than unbounded graph search in N.

For any finite set of pairs (h,L), the union of all its actual source H-values
therefore has O(log X) distinct elements up to X. This excludes every exact
spectrum construction confined to bounded parent height AND bounded interior
coordinate, without assumptions on operation count, finite phases or variable
radices. To evade this obstruction, h or L must be unbounded. This is a
necessary resource requirement, not a sufficiency theorem.

This result is distinct from S134's bounded marked-word run obstruction and
S124's fixed singleton-source fibres. It covers ALL actual parent tournaments
at the prescribed h,L, not one chosen word grammar, while allowing arbitrarily
many neutral prefix and unselected suffix vertices.

## 5. An explicit first positive-interior fibre

At h=1,L=1 the canonical bases are exactly

    (h,F0,q0,g0)=(1,1,3,1), (1,3,5,1).

They are the transitive parents with cut words 11 and101, respectively.
Indeed L=1 forces exactly two selected vertices, and their transitive pair
weight is1 precisely when they are adjacent or have one intervening vertex.
After deleting the neutral prefix and removing the unselected suffix these
are the two possible words. Hence the exact all-order scalar source language is

    {1,3} union {3*2^k-1,5*2^k-1 : k>=1}.                (7)

For example the actual H11 source lies at h1,L1, not on the impossible h1,L0
state with the same C. This is a source-profile classification, not a hole
claim for the full spectrum. Parent height1 here means the parent is transitive;
the full source may or may not be strong.

## 6. Verification boundary

`bounded_interior.py` checks the entire labelled parent universe of orders1..5,
every nonempty cut, the exact SCC normalization, the bound (5), the complete
restored native state, and direct adjacency counts after suffix reencoding.
It also checks the suffix commutation on all labelled blocks through order4,
the complete h1,L1 base universe through the proved bound4, and every cut on
transitive parents through order10 against (7). Normal, optimized and relocated
receipts must agree. These are finite implementation controls; the universal
statements are the path bound, actual SCC normalization and exact tail law.
