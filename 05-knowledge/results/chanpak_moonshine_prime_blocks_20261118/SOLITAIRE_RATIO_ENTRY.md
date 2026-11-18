# Physical ratio changes at fixed height and a two-orientation inverse

**Status:** verified literal ratio changes, an exact multiplicative
target compiler, and a reversible marked dynamics on the 11:13 lane.
No assertion that every prescribed count lies in the compiled domains.

## What the three links suggest

[Eriksson–Jonsson–Sjöstrand](https://combinatorialpress.com/article/ojac/vol15/192.pdf)
replace the bottom-layer choice in Bulgarian solitaire by a prescribed set
of Young-diagram layers. This separates what stays fixed (total cards) from
what the rule reads (layers, not the card total). [Pham](https://arxiv.org/abs/2208.14496)
uses necklaces for the recurrent orbits and a reverse forest to study
transient levels and limiting rational generating functions. In ordinary
solitaire the inverse makes this extra information explicit: for a
partition of length p choose a part k as the newborn pile. Remove it,
add one to the p−1 remaining piles, and add k−(p−1) singleton piles.
This is a legal predecessor exactly when k≥p−1. Distinct chosen sizes
give distinct predecessors because the recovered pile count is k.
For example (3,2,1) has both itself and (4,2) as predecessors.
Keeping the scalar total of six does not tell us which inverse to take.

The linked [unrefereed collection](https://github.com/1337isnot1337/solved-conjectures)
has three useful warnings, not dependencies of our theorems. Its
[L-solitaire note](https://github.com/1337isnot1337/solved-conjectures/blob/main/l-solitaire-limit-shapes/PROOF.md)
constructs a two-cycle at fixed horizontal scale using a pick block
larger than half the deck; it distinguishes that boundary case from
diverging scaling. Its [necklace note](https://github.com/1337isnot1337/solved-conjectures/blob/main/pham-denominator-degree/PROOF.md)
reports a primitive period-ten example with reduced denominator degree
five, not a refutation of the separate common-denominator conjecture.
Its [chromatic-symmetric note](https://github.com/1337isnot1337/solved-conjectures/blob/main/li-coefficient-sign/PROOF.md)
uses an explicit non-spider tree and connected partition counts to
challenge a sign prediction: legal internal cuts, not the outside
shape label alone, decide the coefficient. We use these as inspiration,
not as claims of independent priority review.

The [Berggren tree](https://kconrad.math.uconn.edu/blurbs/linmultialg/descentPythag.pdf)
has the additional property that each non-root primitive triple
has a unique legal smaller parent. Our archived
[reflection cospan](../../../07-reflections/berggren-fibonacci-keller-dihedral-reflection-cospan-codex-20260816.md)
already separated a common reflection presentation from equality
of the underlying dynamics. Bulgarian solitaire does not literally
generalize the tournaments or primitive triples. The shared scheme
is a reversible normalization with an inverse-legality test. Here
we ask for a move at fixed H that changes a hidden ratio, followed
by an inverse on that richer state which has a descent rank.

## A physical ratio change at height 1063

Take the two [literal tape witnesses](ORDERED_PHASE_CARRY.md) with
rooted profiles (H,e,t,d)=(1063,1,1847,1) and (1063,1,1967,1).
Reverse **all the arcs in each** and keep the same roots.
Path reversal preserves H, but rooted weighted starts become
the old weighted ends. In increasing start order the two new
profiles are (1063,1,1105,1) and (1063,1,1177,1), directly
recounted on the two order-12 tournaments. Thus the physical
count-preserving move changes 1847:1967 to **1105:1177**,
not just its reciprocal or a common dyadic scaling.

The [finite reflected-word test](REFLECTED_INDEPENDENT_RETURNS.md)
accepts the new ratio. Use W=2072723099975 and its bit
reverse R=1945929503375, of width 41. They obey
1105W=1177R and have common debit 1834165701692.
One paired return therefore gives 2288524859771745
and keeps 1105:1177. The transfer proves the count, not
a full subset recount at order 53.

This move is available at an **absolute, target-selected**
class of heights. Given N=1063s and a concrete tournament
S with h(S)=s, place the same S before each order-12
witness by an ordinal join. Then H,e,t,d scale by s for
the original root. The old weighted ends instead scale by
the ordered two-cover count B(S): the path ending at the
root can take either an optional prefix in S or leave that
prefix to its complementary path. When the whole graph
is reversed, the new ordinary end count is zero for
nonempty S (the marked root was not a global start).
Both new start reads scale by the same B(S), so the new
ratio is still 1105:1177. Thus exact division by 1063 and
a witness for the quotient select a real move, not an
estimate of the correct height. The replay checks an S=C3
host at height 3189. This does not select other targets
and does not require s to be prime.

There are also two verified substitutions at different targets:

- At the prime **863**, switch from the old d=3 ratio
  11:13 to the d=1 ratio **385:509** by replacing the
  witnesses with a one-apex tournament over TT11 and
  its converse. The apex incontact mask in forward
  tape order is 1154 (the other converse orientation
  uses 1526). The profiles are (863,1,1155,1) and
  (863,1,1527,1). They admit the reflected addresses
  260099,196735 and return exactly at 300220007.
  This is a change of physical source at the *same*
  prescribed prime, not a local edit of the old
  order-11 graphs.
- At **29**, keep the old 13 branch (29,0,78,3)
  and replace the 11 branch by the tournament with
  masks (110,124,120,48,97,64,8), root 4, followed
  by an ordinal sink. The profile becomes (29,0,114,3),
  so the new ratio is **13:19**. Its reflected words
  57,39 return to 4373. The replay directly counts
  both resulting tournaments. This changes only one
  branch in the pair.

The same ordinal-prefix construction extends those exact
entries to heights 863s and 29s when s is physically
supplied. It is not an argument for arbitrary prime heights.

## A native two-orientation return, with its inverse

The old 11:13 pair at height 29 has weighted starts
(11·6,13·6), weighted ends (13·3,11·3), a common deletion
two-cover count of 96, d=3, and ordinary start counts both
3. Its first reflected return corresponds to the marked
state
```
(H,u,v,B)=(863,93,99,1536),
```
where starts have scales (11u,13u), ends (13v,11v),
and B is the common deletion two-cover count. Both
ordinary start and end reads are 3 at this point.

Under one retained-root letter b the weighted end read
changes by b B: if the new vertex has a back arc to the
root, a selected path ending at the root can use that
vertex only just before it. Its other vertices and its
complement form an arbitrary ordered two-cover of the
deletion. With an ordinary sink that contribution is absent.
The deletion gets an ordinal sink either way, so B doubles.
Also the ordinary rooted start count is unchanged for
order at least two: a path starting at the root must put
the added vertex last (it cannot use a back arc into the
already visited root). This appending preserves the exact
start count and sets the ordinary end to d when the last
letter is 1.

A reflected return with j preliminary zero letters
therefore has the full transfer
```
H' = H + 143·2^j u −24−3[j=0],
u' = 2^(4+j)u−3,
v' = v+B·2^j,
B' = 2^(4+j)B.
```
Converse both and reorder the two branches to exchange
u,v at no change to H or B. The two scales remain odd
since B is even. Thus at **every** returned marked state
one can choose either orientation before the next return.
The ratio 11:13 stays fixed here; the converse changes the
absolute growth scale independently. At 863 the choices 93
and 99 give absolute outputs 14135 and 14993, both checked
by literal order-15 subset counts.

For a fixed post-return orientation, the smaller predecessor
is chosen by j=v2(u'+3)−4. Recover u=(u'+3)/2^(4+j),
B=B'/2^(4+j), and v=v'−B·2^j, then subtract the exact
increment from H. Test integrality, positivity and entry
into the literal seed. One may try both predecessor
orientations, since the converse step has order two.
Along each backward macro the deletion two-cover count
strictly decreases, so this is a terminating marked inverse
forest. The reversible fixed-height moves alone do contain
two-cycles; the descent only applies after factoring them
out or allowing one choice per macro.

This example makes the denominator warning precise. If
one always uses the same orientation and j=0, the observed
scalar H after n returns is the sum of a constant, a term
linear in n and a multiple of 16^n; its reduced generating
denominator is (1−16z)(1−z)^2. The hidden v and B are
not needed to predict that scalar trace, yet a legal
converse makes v control the next absolute height. A small
observable denominator is not a certificate that the
underlying dynamics has no additional states.

**Remaining proof task:** identify a height-preserving
physical switch whose target domain reaches outside these
multiplicative classes, or a decreasing inverse that always
finds a legal entry. A legality test for a supplied marked
state is not a proof that an arbitrary absolute H has one.

Replay: `python3 ratio_entries.py` with the earlier modules
in the same packet. Small stated counts use literal subset
DP. Longer reflected continuations use the proved transfer.
