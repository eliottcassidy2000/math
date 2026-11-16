# A rooted reopen after the square reset

For the specific but infinite d=3 branch of [the radix-four carry](../physical_radix_carry_20261116/RADIX4.md), let
v≥2, width r=v+3 and a=3·2^v. Then m=2^r+a=11·2^v. The reset
lands at (4m^2+1=121·2^{2v+2}+1). A different physical reader gives
exactly (4m^2+3) at every height of this branch.

Take the eight-vertex tournament A with outgoing row masks
```
(6,12,40,33,15,83,159,63)
```
and mark vertex 0. Let h=H(A), g=H(A−0), and let r₀ and s₀ denote
the sums over even-order and odd-order paths starting at 0 of
H(the induced complement), with H(empty)=1. An independent subset-DP
gives ((h,g,r₀,s₀)=(107,37,69,89)).

Append a new vertex z beaten by every old vertex except vertex 0,
which z beats. Keep 0 marked. The exact rooted extension formula
from [the earlier inverse notes](../inverse_fences_landing_20261116/INVERSE.md) is
```
(h,g,r₀,s₀) → (h+2r₀,g,r₀+s₀−g,r₀+s₀).
```
A direct explanation: the new odd-cycle pack contribution is
precisely the nonempty option of an odd cycle with z and an old
even-order rooted path; this adds 2r₀. In the rooted sums an old
path can append z unless it is the singleton at 0, which has
complement count g; appending changes the parity of the path.
A new rooted path cannot contain z except as its last vertex.
Without z on the path, z is a neutral sink in the complement since
0 has been removed. Deleting 0 in the whole graph likewise leaves
z as a neutral sink.

After its first use the two rooted sums differ by g, so thereafter
the even sum doubles at each step. Hence after n≥1 appends the
new full count is (h−2s₀+2g+(r₀+s₀−g)2^n=3+121·2^n).
Choose n=2v+2. The literal output is strong, since A is strong
and each newcomer reaches the root while each nonroot can reach
the newcomer. This gives a fourth consecutive count after the
three-state canonical carry on this infinite subfamily. The new
reader has one fewer vertex than the tape reset witness. Thus this
argument is a source replacement, **not** a vertex-retaining +2
operation on the original tape, nor a statement for all r,v.
At the first height v=2 the counts are 7741,7743,7745,7747.

## A vertex-retaining conversion with only three fixed old-arc changes

The original reset can also reopen without discarding even one old
vertex, by adding two vertices. Collapse its dyadic padding into the
literal unmarked transitive spacer: the reset at width r=v+3 and
masks a=3·2^v,b=5·2^v is the same as widths 4 and 3, masks a=3,b=5,
and spacer length n=2v+2. Call its nine remaining boundary vertices
```
L0,L1,L2,L3,R0,R1,R2,sentinel,head.
```
The five vertices L0–L3 and head beat every spacer vertex. Each
spacer vertex beats the other four. Following the reopen, interpret
the old boundary vertices as the roles
```
6,7,4,3,0,8,9,10,5
```
respectively, where roles 0–7 are the A seed and roles 8→9→10
are three neutral terminal sinks beaten by the entire other graph.
Add the two new boundary vertices that fill A roles 1 and 2.
Both newcomers beat every old spacer vertex. Within the entire
new boundary, use the orientations of A→TT3.

This keeps **every** old arc incident to the spacer, including its
internal arcs. Among the nine old boundary vertices, the role
assignment reverses exactly three arcs:
```
{R0,head}, {R1,head}, {sentinel,head}.
```
Thus the physical reset tournament is converted to the new one
by three fixed old reversals plus two additions, no matter how long
the spacer grows. The reopened tournament consists of the strong
rooted reopen core followed by the three neutral terminal sinks,
so their presence multiplies its count by one. This version is
not strong and does not show that a second physical reopen can
reuse the same edit motif. The role comparison is independently
checked in the replay.

## Next +2 without adding any vertices

There is also a strong eight-vertex seed B with outgoing rows
`(2,36,65,23,135,93,155,47)`, root 1, and profile
((h,g,r₀,s₀)=(247,31,121,152)). The same extension law gives
(5+242·2^j) after j≥1 appends. Use j=n−1 and pad with four
terminal sinks to match the order of the preceding A construction's
n appends and three sinks.

To do this without changing **any** arc incident to the first n−1
spacer vertices, reinterpret the old A roles 0 through 7 as the new
roles
```
1,4,2,0,3,6,7,5
```
in B. The old root remains the new root. Move the last old spacer
vertex to become the first of the four terminal sinks. This entails
one flip between that last vertex and the old root, plus four flips
inside the old A roles: `{0,1},{0,7},{1,6},{3,5}`. No vertices are
added or removed. Every other old edge survives. This yields
(121·2^n+5) at every branch height, or a total of five consecutive
values across the original radix-four carry and the two localized
reopens. The smallest example is 7741,7743,7745,7747,7749.

A finite bank of rooted readers of one odd slope q has a general
version of this locality. If their raw slopes are q·2^e with e
drawn from a finite nonnegative set, choose N>max(e) so each gets
N−e appends. Separate off a common suffix of N−max(e) appended
vertices, absorbing the others into its fixed boundary. Every
common vertex is beaten by all boundary core vertices except the
root, which it beats. Prior appended vertices lie on the incoming
side. Neutral prefix and suffix padding make the incoming and
outgoing boundary sizes respectively match across readers. Matching
these roles literally preserves all arcs incident to the common
spacer, so changing readers requires only a constant-size boundary
rewire at each step. This compiles a *finite* interval of known
intercepts, and makes no claim that the bank supplies every intercept
or repeats the low-digit refill at the next high address.
