# Two genuinely different H=9 source mechanisms

The [H=3 right-fiber theorem](
../chanpak_h3_right_fiber_20261114/PROOF.md) separates the number
of possible offsets from the width of a consecutive window. At
H=9, merely using two triangles does not fix that problem. A
single five-vertex strong component does: a physical right fiber
gives four consecutive odds and an all-height positive family.

A right fiber fixes the source quantities (h,f,P,S): its ordinary
path count, two-cover count, total weighted starting-path count at
the apex outcontacts, and odd-minus-even part of that total.
Writing K for the apex path count, appending the same marked
right tournament rescales all K differences by its ordinary path
count and adds a common constant. This is the
[ordinal cut-state formula](
../chanpak_selective_target_20261114/CUT_STATE.md).

## Which H=9 sources are even possible?

Path counts multiply over strongly connected components. By
[Busch's strong-order minimum](https://doi.org/10.37236/1141),
a strong component with count nine has at most five vertices.
An exhaustive finite check through order five (in the replay)
gives precisely the two strong classes with outgoing masks
A=(14,28,24,16,1), B=(10,28,25,16,1); both have f=84. Thus
every H=9 base has either two triangle components, or exactly
one of A,B, with arbitrary singleton padding. The resulting
two-cover counts have disjoint odd parts: 9 and 21 respectively.

## Two triangles: a signed four-bit interaction, but no triple

Compress each marked triangle as its contact count r=0,1,2,3 in
two bits. The rest are the actual singleton contact bits. Fix
the resulting word W from P/9 and the width from f=9·2^w.
Choose two disjoint adjacent pairs, at i<j, to expand. Compare
the scalar K against the reference 9g(W), where g is the
transitive-spine apex count. If δ_i,δ_j denote the H=3 single-
replacement offsets calculated on W, the exact new offset is
[
    K-9g(W)=3(delta_i+delta_j)
     +[j=i+2],d(r_i)e(r_j),
]
with d=(0,−2,2,0), e=(0,−1,1,0).

For a direct derivation use a virtual two-bit block in each
triangle's place with its scalar transfer coefficients scaled
by 3. Its h, f, and starting total P then match the triangle.
Relative to it a triangle changes the start bias by d(r) and
the ending bias by e(r). The start bias resets after each
ordinal block, so its change can affect only the next block's
count update. If both replacements are adjacent that next
block's ending bias is changed by e as well; the product of
the two changes gives the sole cross term. All independent
one-triangle offsets multiply by the other block's h=3. This
also keeps the correct terminal bias: S=9 times the last bit
unless the last block is a triangle of contact count r, in
which case S=3r.

Each δ is even. Thus any offset not divisible by six must
come from a pair of contiguous mixed triangles. For strict
prefix x, p=x mod 2, and q=1 iff a following bit exists and
is losing, the four options are as follows.

| Four bits | Offset | Residue modulo 6 |
|---|---:|---:|
|0101|45x+14−3p−6q|+2|
|1010|−45x−40+3p+6q|+2|
|0110|−27x−14−3p+6q|−2|
|1001|27x+16+3p−6q|−2|

The two nonzero residue classes are **never within four of
one another** in a single word. For positive offsets compare
0101 with 1001. If 0101 comes first with prefix x, the other
cannot begin until at least three positions later, so its
prefix y is at least 8x+2. Its offset is at least 27y+10,
whereas the first is at most 45x+14. If 1001 comes first,
the next 0101 is at least two positions later with prefix
y≥4x+2; now the latter offset is at least 45y+5 and the
former at most 27x+19. Both separations exceed four. For
negative offsets compare the absolute magnitudes for 1010
and 0110. When 1010 comes first the earliest legal overlap
starts three places later with prefix y≥8x+5; compare
27y+8 with 45x+40. When 0110 comes first the earliest overlap
is two places later with y≥4x+1; compare 45y+31 with
27x+17. Again the separation exceeds four. Opposite signs
are separated even at prefix zero.

Every three consecutive odds in this chart would use both
nonzero residue classes modulo six, an impossibility.
Furthermore an adjacent pair must touch offset class zero,
and K's reference and the common right-reader shift are
multiples of three (h,P,S all divisible by three). If the
same right reader is nontransitive, it rescales every nonzero
gap by at least three instead. Thus the two-triangle mechanism
does not produce a prime/7p unit pair with both counts
coprime to three in a fixed right fiber. These are restrictions
on a particular chart, not on the H=9 source count in general.

## A single strong H=9 component supplies four consecutive odds

Use A from above. Form TT_k → A → TT_m. In the next table
the left and right columns are fixed-width binary addresses
of apex outcontacts on those transitive slabs. All four rows
have exactly (h,f,P,S)=(9,672,245,9), but different K.

| Left slab | Central contact mask | Right slab | K |
|---|---:|---|---:|
| 0 | 27 | 01 | 89 |
| 01 | 18 | 1 | 91 |
| empty | 9 | 101 | 93 |
| empty | 20 | 101 | 95 |

These are four literal choices, not an interpolation in the
right profile. The K values and the shared coordinates are
verified by subset-DP on the actual graphs. The four-way
choice is retained under every common transitive right
reader. For instance, taking a losing suffix of length m≥1
gives
[
        245·2^m +{-147,-145,-143,-141}.
]
The first strong window is 343,345,347,349. For arbitrary
transitive right reader reverse-complement address y the
increment is 9(g(y)−1)+245y+9(y mod 2). At y=45 the
window is 11249,11251,11253,11255; the first two are
7·1607 and a prime, respectively, with 1607 prime. The
replay checks primality by trial division and all four
actual tournament counts. After removing the first two rows’ neutral leading zero, connectivity follows since
the apex can reach the initial SCC and the final singleton
beats it for each of these suffixes.

There is a direct target inverse on the displayed dyadic
window: for each of four offsets c, require (N−c)/245
to be an integral power of two with exponent at least one.
It selects the actual row and suffix. This is sparse and
not a universal selector.

### Exact inversion of a prescribed strong H=9 right profile

For either of the two five-vertex cores R, store its 32 marked
profiles with (P_R,S_R,Q_R,T_R,K_R), where Q_R,T_R are the total
and signed weighted **ending** sums on the losing side (distinct
from the two-cover f_R). If f=84·2^b, the total number of padding
singletons is b. For a proposed right width m (0≤m≤b) and left
width k=b−m, let x,y be the two actual contact addresses. Then
[
  P=(84x+P_R)2^m+9y,qquad 0\le x<2^k,quad 0\le y<2^m.
]
Because 9 is odd one recovers the unique y as 9⁻¹P modulo 2^m.
Put u=(P−9y)/2^m. The central starting total lies in [0,75]
(the full weighted start total is f_R−h_R=75). Hence its
only candidate is u modulo 84; the quotient gives x. This
decodes the two pads and restricts directly to the atlas
entries with the right central P_R. The terminal bias is S_R
for m=0 and 9 times the final bit for m>0. Reject failed guards
and compute K for the survivors with the ordinal formula;
compare to the absolute target. No enumeration over pad
masks is needed. The inverse takes b+1 passes through a fixed
64-row profile atlas.

Even varying f at fixed (h,P,S) adds only redundant leading
losing padding to this decoder beyond O(log(P+1)) positions.
If u>0 then 2^m≤P, and leading zeros in x do not alter the
count or profile. If u=0 then P_R=x=0, and arbitrary leading
zeros in the right address add no effect since the source
prior to its first contact has no outgoing apex arc. The
m=0 boundary versus positive m must still be distinguished.

The strong core thus repairs the no-triple obstruction, but
the process still lacks a way to choose the right profile for
an arbitrary target. No global width-four upper bound for the
strong H=9 branch is asserted.

Run `python replay_h9.py`.
