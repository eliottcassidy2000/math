# The h=3 right fiber: as many offsets as positions, but never three consecutive odds

This follows the [minimum right-only switch](
../chanpak_f1_physical_carry_20261114/RIGHT_ONLY_SWITCH.md).
The marked puncture B has h=H(B), f=P_B(2), and P,S the total
and odd-minus-even weighted sums of paths starting at an
apex outcontact. A path is weighted by the Hamilton count of
its complement. The apex output is K=H(B+v). Fix the four
right coordinates (h,f,P,S).

**Theorem.** For h=3 and f=3·2^w (necessarily w≥2) there are
at most w−1 distinct K values at any fixed right coordinates;
this bound is attained at every w. Nevertheless at every
scale the longest consecutive odd interval inside one such
fiber has length at most two. The latter is sharp for w≥3.
An adjacent pair must include a distinguished scalar K0
defined below, and a prescribed K can be decoded by a
constant number of binary-prefix tests. The conclusions
about cardinality and no triple persist after any one
fixed ordinal right continuation (even one with a nontransitive
puncture). They do not compare different right profiles
or different continuations selected at the target.

## One triangle and a binary address with one pair expanded

The only nonsingleton SCC of an H=3 tournament is a
directed triangle: H multiplies over the ordinal SCC word,
and Busch's strong-order lower bound gives H≥5 from order
four on. All other vertices are singleton SCCs. The triangle
has two-cover count 12, total start contribution 3r if r of
its vertices are contacted, and signed start contribution r.
These follow from cyclic symmetry and the signed two-cover
identity, or directly by listing its directed paths.

Replace the marked triangle by the two bits of r, and retain
each singleton's actual contact bit. The resulting word W
has width w, since f=12·2^(w−2). Its binary address is P/3.
Indeed when a block R is appended, the new start total is
f_R P+h P_R; before the triangle h=1, after it h=3. A
singleton doubles and adds its bit (times h); the triangle
quadruples and adds r in the normalized address. Leading
zeros are known from f. Thus recovering a candidate puncture
is exactly choosing one adjacent pair of W to expand to a
triangle, with its r contacts unique up to triangle rotations.
If a singleton follows the triangle then S=3 times the final
bit. If the triangle is last then S=r, so the last mixed pair
belongs to a different S fiber from all earlier placements.
A final homogeneous pair has precisely the common S. This
gives the upper cardinality bound without an order census.

## The exact offset at a chosen pair

Let K0=3g(W), where g(W) is the path count of the apex over
the transitive word W. This is a reference scalar, not
always an available puncture. Let x be the strict prefix
address just before the pair and p=x mod 2. Let q=1 if the
following bit exists and is losing (0), otherwise q=0.
The offset from K0 is
[
 Delta_{00}=Delta_{11}=0,qquad
 Delta_{01}=2+3x-p-2q,qquad
 Delta_{10}=-4-3x+p+2q.
]

For a derivation, use the [ordinal cut transfer](
../chanpak_selective_target_20261114/CUT_STATE.md).
For a two-vertex transitive block the nonempty end
total and signed bias, into the uncontacted vertices, are
(3,1), (1,1), (2,0), (0,0) in the four binary cases.
Its augmented K is respectively 1,1,3,1. For the triangle
with r contacts the corresponding numbers are (3(3−r),3−r),
and its augmented K is 3,5,5,3. The actual and reference
scalars after the two bits have offsets respectively
0, 2+3x−p, −4−3x+p, 0. Their start totals and h agree, but
the triangle's signed start bias differs from the scaled
reference by 0,−2,+2,0. Only the next singleton can read
that difference, and it does so exactly when losing. This
proves the formula including the terminal boundary.

An alternating word starting at 1, with its last bit
duplicated (truncate to w−1 before duplicating), puts all
pairs in the common S fiber. The last is homogeneous and
the others are mixed. In this word all 01 prefixes are
positive and all 10 offsets are negative. Distinct
same-type pairs have strictly increasing absolute offsets
by the estimates below. Thus this word attains w−1 distinct
values. At width seven, W=1010100 gives
187, 209, 215, 219, 221, 233.

## All unit adjacencies must touch K0

The 10 offsets are at most −2. For consecutive occurrences
of the same 10 pair type in the word, if its strict prefix
changes from x to y then y≥4x+2. Hence the increment in
absolute offset is at least 3(y−x)−3≥9x+3. This is even,
so no two negative offsets are adjacent at distance two.

The 01 offsets are nonnegative. Here y≥4x+1 for two
successive occurrences of that type. If x≥1 the increase
is at least 9x; if x=0 and y≥2 it is at least 3, hence
at least four by parity. The only exceptional case is
x=0,y=1: the occurrences must be immediately consecutive
as pairs, since appending any bit after the first 01 makes
the prefix at least two. So the first occurrence has q=1
and its offset is zero. Thus no two strictly positive
offsets are adjacent either.

Nor can both −2 and +2 occur in one word. The −2 case
forces x=0,q=1 for a 10: the first 1 begins 100.
The +2 case either has x=0,q=0 at an 01 (the first 1
is followed by another 1) or x=1,q=1 at an 01 (the first
1 begins 1010). These conditions conflict. Therefore any
adjacent pair includes K0, and a three-odd run is impossible.
The width-three word 100 attains two, and any number of
leading zeros can be inserted without losing this pair
(recompute W at fixed widths). An appended fixed right
reader rescales differences by its positive odd h and
adds a common baseline by the cut transfer. It cannot
create a triple of consecutive odds.

### The width coordinate can be dropped

At fixed (h,P,S), even if f varies without bound, the same
no-triple bound holds. Indeed the fixed binary address P/3
determines W up to leading zeros. A leading zero does not
affect g, old pair offsets, or the bias guard. Once there
are two leading zeros, increasing their number adds only
more 00 placements with the same offset and bias. Thus it
suffices to use the binary expansion of P/3 with two leading
zeros (truncated to 00 for address zero). This also shows
that unbounded width alone does not make a fiber infinite
when its actual right-readable coordinates (h,P,S) are fixed.

In fact the raw no-triple proof allows all positions of a fixed
address even before the S guard. A common transitive suffix
preserves this stronger assertion at fixed (h,P), since the
candidate expansions are then positions in one longer word.

### Every adjacent pair on a right fiber includes a multiple of three

If a fiber has more than one output, S is either 0 or 3:
a terminal mixed triangle gives S=1 or 2 and is unique.
Both h and P are multiples of 3. Hence the additive shift
from any common right reader is divisible by 3, as is its
scaled image of K0. Every consecutive-odd pair contains that
image. In particular the minimum h=3 family cannot match
the earlier h=9 prime/7p unit-pair mechanism with both counts
coprime to three, even though it has unbounded offset
multiplicity as the right address varies.

## Exact target-dependent inverse inside the fiber

Compute W from P/3 and f, and compute K0 in a single scan.
Subtract the known right-reader contribution if one was
used, then divide by its h if integral, and subtract K0.
For a nonzero remaining offset, its sign determines
01 or 10. Try the four options for p and q; each gives the
only possible integer x from the displayed affine equation.
For x>0 its position is forced by its bit length and the
number of leading zeros of the fixed-width W. At x=0
there is at most one occurrence of the given mixed pair:
it must use the first 1. Verify the pair, its prefix and
q, and the terminal S check. For zero offset consider all
homogeneous pairs (a scan finds one, if any) and the at
most one leading 01 that can have offset zero. By replacing
the selected pair with a literal marked triangle, the
inverse constructs the cut without testing any unrelated
tournament or relying on residue matching.

A specified right profile at h=3 is now completely decoded.
The result rules out a tempting extension of the earlier
unit switch: adding more insertion positions creates
unbounded offset multiplicity but never enough local
width for even three consecutive odds. Selecting the
right profile from a completely unrestricted target,
or achieving wider windows by higher h or different
reader classes, remains a separate problem.

Use Arthur H. Busch's [strong-order bound](https://doi.org/10.37236/1141).
