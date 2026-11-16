# Reopening the square reset with rooted readers

Use v≥2, and put s=2v+2, m=11·2^v and N=4m²=121·2^s.
This is the width family r=v+3, a=3·2^v from the
[previous carry](../physical_radix_carry_20261116/RADIX4.md).
The old construction gave N−3→N−1→N+1. The reset at N+1
now genuinely reopens to N+3. In addition there is an exact
finite library of **38 consecutive odd counts, N−5 to N+69,
at every v≥2**. At v=2 the whole band is 7739,7741,…,7813.

## What changes in the reader

For a rooted seed with root x, record h=H(T), g=H(T−x) and
the weighted sums R,S over paths starting at x with even and odd
numbers of vertices respectively, giving each path P the weight
H(T−V(P)). Include the singleton x in S. Append a new vertex z
which beats x and loses to every other old vertex.
The exact update is

    (h,g,R,S) → (h+2R, g, R+S−g, R+S).

A new odd cycle through z uses exactly an even-order old root
path. For the two new root sums, a path through z ends there and
is any old root path except the singleton, with its parity toggled.
Paths not through z keep their complement weights (z is a sink
relative to the complement, since it does not contain x).
This proves the update by the odd-cycle-pack count. In particular
after ell≥1 appends the count is B+A·2^ell where

    A = R+S−g,          B = h−2S+2g.

The N+3 reader is a literal eight-vertex seed with outgoing rows
(6,12,40,33,15,83,159,63), root 0, and (h,g,R,S)=(107,37,69,89).
Thus A=121 and B=3. The old transitive seam stops at the high
square, whereas this reader carries the same growing coefficient
with a new constant term. It is not claimed that this reader
can be created by an induced extension.

## Exact local conversion (no changes touching the long spacer)

The reset has a marked head plus left tape L0…L3, long transitive
unmarked spacer X of length s, right tape R0…R2 and sentinel t.
All tape arcs point forward; the old marked sites are L1,L3,R0,R1
and t. In the new graph assign the nine old boundary roles
(L0,L1,L2,L3,R0,R1,R2,t,head) to roles
(6,7,4,3,0,8,9,10,5) in the rooted eight-vertex reader plus three
neutral sinks (roles 8,9,10). Add the reader vertices 1,2. Nonroot
reader vertices beat all X, X beats root 0, and X beats all three
sinks. Comparing the old and new prescriptions leaves every arc
touching X unchanged, and reverses just three old boundary arcs:
head vs R0, R1, and t. This exact conversion retains all old
vertices and adds two. It gives 7745→7747 at the first height.
[An even smaller next-step rewrite](LOCAL_MOVES.md) reaches
N+5 with five reversals and no additions.

## A whole bank on one preserved long spacer

readers.json lists (B,A,root,rows). Each row was independently
verified; for each B=−5,−3,…,69, A=121·2^j for 0≤j≤4. Giving
that seed s−j appends produces N+B. Thus this is an exact
constant-width band at unbounded scale, not an extrapolation of
a finite count census. Discovery scanned each order-8 atlas parent
with a new root; the resulting witnesses do not rely on atlas
completeness. For discovery, old rooted weighted sums at i are
R_i,S_i. If the new root beats old positions in Q, its r,s sums
are Σ_(i∈Q) S_i and g+Σ_(i∈Q) R_i. For h, each possible
root-cycle contributes the even old path count weighted by its
complement; this gives the exact weighted directed cut. The code
root_atlas.cpp replays that search, while verify.py directly
checks the displayed seeds instead.

These are physical changes on a common vertex set, not just
existence of unrelated graphs. Set aside the first old spacer
vertex as an extra neutral source, and repurpose R1 as the last
active append vertex. The new orientation for N+3 uses the same
two newcomers as above and nine old arc reversals instead of
three. Retain the first s−4 of the remaining old spacer vertices
as the fixed long core of the tail. The remaining three central
vertices plus R1 supply four following seats. There are now:

- nine reader seats, one fixed as root and eight preceding the
  fixed tail; any spare reader seat is a neutral source;
- s−4 common appended vertices, each beating the root and
  losing to all the other reader seats;
- four following seats, of which the first 4−j become more
  appended vertices and the rest are neutral sinks;
- the two remaining terminal singleton sinks.

Every reader has precisely s−j active appends under this scheme.
All incidences touching the common s−4 vertices agree, as do
incidences touching the two terminal sinks. Only the 13 mutable
reader/following seats need rewiring, so any next +2 step in the
band uses at most 78 old arc changes and no added vertices,
uniformly in v. The smaller special edits above are useful
as independent controls, not a minimum over all allowed rewires.

This reopens many high cycles after the reset. The band is only
a fixed number of counts around the ray N=121·4^(v+1).
It does not select an arbitrary high address or establish
unboundedly many consecutive carries within a single band.

The first reader also has a short two-apex description. Use the
increasing spine 202310 followed by s copies of 2, with symbol
2[a beats site]+[b beats site] and mutual arc a→b. In order
(spine sites 0…5, apex a, apex b), the seed rows are
(190,252,184,48,96,192,141,24). Relabel them to roles
(6,7,4,1,2,3,5,0) of the same H107 seed; the marked root is b.
Every following 2 site is exactly a rooted append. Thus the
two-apex tape is the same reader in simpler coordinates, not an
extra independent source of offsets.
