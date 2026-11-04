# Arbitrary-prefix capacity: changing the strong suffix is necessary

2026-11-03. **PROVED**, with a separate finite exact audit. This is an
obstruction to one attachment operation, not to the Hamiltonian-count spectrum.

The inherited [canonical zero-ready identity](../triangle_return_170807_20261103/ZERO_READY.md)
is

    H(Z(A;B,w)) = h H(A) + F P_A(2),
    h=H(B), F=E_w(B,2), A nonempty,

where B is strong, every A vertex dominates every B vertex, and the new root
r is beaten only by w. The first slot of E_w is nonempty and ends at w;
the second may be empty. Both slots of P_A may be empty individually.
All vertices, attachment relations, and the selected port are part of the
input. The scalar pair (h,F) evaluates this particular operation only.

## Sharp elementary capacity bounds

For every tournament A on n>=1 vertices,

    P_A(2) >= 2^n,
    P_A(2) >= 2 H(A) + 2^n - 2.                      (P1)

The first equality holds exactly when A is transitive. If A is not
transitive, the sharper bound is

    P_A(2) >= 3*2^(n-1),                             (P2)

with equality exactly for

    A = TT_a -> C3 -> TT_b,  a,b>=0, a+b=n-3.

Proof. Write H(empty)=1. Splitting vertices between two ordered slots gives
the exact identity

    P_A(2) = sum_{S subset V(A)} H(A[S]) H(A[V(A)\S]).

Every tournament, including every induced subtournament, has a Hamiltonian
path. For completeness, insert a new vertex in an existing directed path:
put it first if it beats the first vertex, otherwise after the last initial
segment of vertices that beat it; the next vertex, if present, is beaten by
it. Hence every summand is at least one. The empty and full subsets each
contribute H(A), proving both inequalities in (P1). A transitive tournament
has exactly one path on every subset. Conversely, a nontransitive tournament
contains a directed triangle, whose induced Hamiltonian count is three, so
some summand is strictly larger than one.

Fix a directed triangle D in a nontransitive A. There are 2*2^(n-3)
two-colour assignments in which all three vertices of D lie in one slot.
In each such assignment that slot has at least three Hamiltonian paths:
this follows by inserting its other vertices into each of the three paths
of D. The insertion procedure preserves their relative order, so the three
outputs are distinct. These assignments therefore contribute at least two
above the baseline one. This proves (P2).

If there is a second directed triangle E, colour E monochromatically and
some vertex of D\E with the other colour. Such an assignment makes D
nonmonochromatic but contributes more than the baseline. Consequently
equality in (P2) permits exactly one directed triangle. Every vertex outside
that triangle must uniformly dominate it or be dominated by it: a mixed
incidence would create another directed triangle using two of its vertices.
The first group dominates the second, since an edge in the other direction
would create another triangle using a vertex of D. Both groups are
transitive. Thus the indicated ordinal-sum form is necessary. Conversely,
ordered two-slot covers multiply under ordinal sum: concatenate the paths
in each slot in component order. Therefore this form has
P_A(2)=2^a*12*2^b=3*2^(n-1), proving sufficiency.

The relative-order-preserving insertion observation also proves the
monotonicity H(A[S])<=H(A) used in the argument. It does not assert induced
extension of a prescribed path with fixed endpoints.

## Complete exclusion at the actual 170807 state

The actual marked strong suffix in the inherited named graph has

    h=2195, F=42153, A=(50->51), r=39, w=31,
    H=2195+4*42153=170807.

Keep that exact B and w, but replace A by an arbitrary nonempty tournament,
of any size. If |A|>=4, (P1) gives

    H(Z(A;B,w)) >= 2195 + 16*42153 = 676643,

already above both 512423=3H+2 and 512425=3H+4. At smaller orders the
complete list of possible (H(A),P_A(2)) is

|A|H(A)|P_A(2)|H(Z(A;B,w))|
|---|---:|---:|---:|
|TT_1|1|2|86501|
|TT_2|1|4|170807|
|TT_3|1|8|339419|
|C3|3|12|512421|

None gives either required return. In particular the apparently natural
triangle replacement gives exactly 3H and misses the two increments. This
rules out *all* prefix replacements with this suffix, including nontransitive
prefixes; it is not merely the earlier obstruction to resynchronizing two
transitive clocks.

## General finite sieve and its boundary

For any fixed marked strong suffix (h,F), every candidate target N requires

    2^|A| <= (N-h)/F.                                (P3)

This is a necessary integer capacity condition, not a sufficient realization
test. For a nontransitive prefix, substitute 3*2^(|A|-1) for 2^|A|.
Thus arbitrary-size prefix replacement is a finite exact problem for each
fixed N. The equality objects in (P1),(P2) give positive controls and the
first nonlinear boundary. Changing B, its port, or its external incidence
escapes this sieve and is precisely the live direction here.

## Audit and provenance

`python 05-knowledge/results/zero_ready_strong_switch_20261103/check_prefix.py`
checks every labelled tournament of orders 1..5 (1099 graphs), including
transitive and single-triangle equality controls, against direct subset
Hamiltonian dynamic programming. It checks (P1), (P2), both equality
characterizations, and the four actual target values. This finite audit is
separate from the all-order proof above. It uses no spectrum exclusion as a
filter, and performs no write to inherited packets.

The ordinal-sum path-slot operation is inherited; the new use is its sharp
capacity/equality boundary at this named zero-ready source. No external
priority claim is made, and no conclusion about all other odd counts follows.
