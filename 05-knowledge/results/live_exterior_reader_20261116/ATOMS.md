# One square, two retained states

Fix arcs 0→1, 0→3, 1→2 and 2→3, contributing minimum outdegrees
(2,1,1,0). Give arcs 02 and 13 their forward direction at flags 00,
and reverse them with flags x and y respectively. The two free edges are
a perfect matching; the four fixed edges form its complementary cycle.

| xy | outgoing rows | class | directed triangles |
|--|--|--|--|
| 00 | (14,12,8,0) | transitive | 0 |
| 01 | (14,4,8,2) | source over triangle | 1 |
| 10 | (10,12,9,0) | triangle over sink | 1 |
| 11 | (10,4,9,2) | strong | 2 |

Let P(t) count ordered t-slot path covers with empty slots, for integer
t≥1. The positive odd-cycle packing formula, applied to all slot
assignments, gives
```
P(t) = t^4 + 2(x+y)t².
```
Empty cycle packing has four singleton choices of slot. A directed
triangle has coefficient 2 and must be in one slot, leaving its
singleton free. No two triangles are disjoint. In particular
H=1+2(x+y), P(2)=16+8(x+y).

## Recursive invisibility, and how to break it

Ordinal sum multiplies P(t): in each slot concatenate the path
from each ordinal factor, with unique boundary decomposition.
Replacing each of k ordinal atoms independently by 01 or 10 yields
2^k **nonisomorphic** tournaments, with exactly the same P(t) for
every t, namely (t^4+2t²)^k. To recover the word from the isomorphism
class, list the strongly connected component sizes in ordinal order,
then group by each successive interval of four vertices: each atom is
(1,3) or (3,1). The partial sums 4j are component boundaries and
specify the same grouping intrinsically.

This invariance is useful when an exterior must be kept live
without changing what an unmarked reader measures. But two-slot
compression is insufficient for a selectively marked reader. At
root 1 the two chiral atoms have
(h,g,R,S)=(3,1,1,2) and (3,1,3,3), respectively. Thus
(Z,I)=(2,1) and (5,−1) in the reader coordinates Z=R+S−g,
I=h−2S+2g. The root-deleted tournaments are transitive; in
their respective transitive orders the root's out-neighbor words
are 010 and 101. A common synchronizing 0 followed by address
b and final 0 gives count difference
```
−2 + 6(2D−rev(b)).
```
Thus this selective exposure does not produce a unit change.
The example pinpoints which data the reciprocal construction has
to preserve, rather than postulating that common path-cover
moments suffice.

## Price and the old Fibonacci-box work

In H. Lee Price, [The Pythagorean Tree: A New Species, §2–3](https://arxiv.org/pdf/0809.4324),
the four box entries in window order are (p−q,q,p,p+q),
and the two half-angle tangents x,y satisfy xy+x+y=1.
The three perfect matchings on the entries are exactly the
three product pairings in the old
[fibonacci-box note](../../../07-reflections/fibonacci-boxes-farey-and-the-k4-matching-triple-opus-20260728.md).
The tournament atom selects one of those three matching types
to be free; no assertion that the count polynomial equals
the Pythagorean form follows.

The HAT identity has a more direct reader translation. For
p=b/D and q=rev(b)/D, the dominant count is
D(Z+gp)(2−q). With u=gp/Z and v=1−q this is ZD(1+u)(1+v).
At target 2ZD, Price’s HAT partner v=(1−u)/(1+u)
is exactly the exterior compensation curve. At any target, replace 2
by N/(ZD). The affine-coordinate reciprocal map says which
real exterior partner is needed; it neither proves that the two
addresses come from one word nor absorbs the remaining integral
defect.

## Chair and Monster scope

The old recursive-interfaces note uses the Chair44 mechanism
that unique coarsening divides an alleged integer translation
period by 2 at each level. The useful ingredient is not simply
that parents exist; it is an intertwining relation with a strictly
shrinking integer obstruction. A reader hierarchy needs either
that kind of finite descent for its exact defect, or a constructive
termination lemma. An infinite tree of nonempty coarse exterior
fibres does not prove exact placement. See
[recursive interfaces](../recursive_interfaces_20261016b.md).

The earlier [Monster/Baby Monster investigation](../monster_levels_square_sum_7_21_20261022.md)
distinguishes the exact level 1|2|6 incidence filtration from
the valuation game 4→2→5→3→7→13 (where 4 must be read as 13 to
close the word). Its concrete parallel is refiltering the same carrier
by changing the observer; the valuation word provides no count identity.

## Compression under bounded endpoint observation

The same bank survives all ordered-slot cover experiments prescribing
start or end vertices in at most r of its k atoms, as long as the marked
atoms' states are held fixed. Condition on the colored paths inside the
marked atoms, including satisfaction of the within-atom endpoint
constraints. In an unmarked atom each slot is forbidden if its prescribed
start lies in a later atom or its prescribed end in an earlier atom.
For an unmarked atom with s allowed slots, there are P(s) completions,
including P(0)=0. Multiply over unmarked atoms and sum over the
conditional choices. This is an exact endpoint-profile collision with
2^(k−r) nonisomorphic alternatives. It applies to these marked-slot
observations, not to selective contacts that reorient arcs inside
the unmarked atoms or let a slot revisit a block.

## Uniqueness and locality of the compressed ambiguity

Up to relabelling, the displayed square is the only two-free-arc
face containing all four classes. Label its transitive corner in
transitive order. Flipping any of its three adjacent arcs remains
transitive, so the two free arcs must be among 02,13,03. Flipping
03 by itself is strong. Flipping 03 together with either 02 or 13
gives score sequence (1,1,2,2), also strong at order four. This
rules out the two choices involving 03.

