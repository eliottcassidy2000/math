# A variable-coefficient reader with polynomially many exact windows

The previous rooted readers reopened the carry at fixed odd
leading coefficients. This note changes the coefficients as well
as proving an all-scale absolute coverage bound. The result is
**at least a constant times X^(2/3) disjoint three-count windows
below X**, and an explicit selector that brackets every large
target within O(sqrt X). Neither is full odd-count coverage.

## Three physical readers with a synchronized binary address

Use outgoing rows and root 3 in each of the following seeds.
Record h=H(T), g=H(T−root) and weighted even/odd root-start path
sums R,S as in the
[previous reader](../reader_reopen_collatz_20261116/BANK.md).

| Rows | h | g | R | S |
|:--|--:|--:|--:|--:|
| (2,4,25,19,35,15) | 23 | 9 | 24 | 26 |
| (6,4,24,19,35,15) | 25 | 9 | 24 | 26 |
| (10,4,17,22,35,15) | 29 | 9 | 23 | 27 |

Let Z=R+S−g and I=h−2S+2g. Before synchronization
Z=41 and the three I values are −11,−9,−7.
There are two actual appends:
- 0: the new vertex only beats the root (near-sink).
- 1: the new vertex loses to every old vertex (ordinary sink).

For 0 the state maps to (h+2R,g,R+S−g,R+S) as before.
For 1 h,g are unchanged and both R,S become R+S. Indeed a root path
ending at the new sink extends any old root path, switching parity;
paths not using it preserve their complement counts. Thus after
one common 0 append, Z=82, g=9, and the distinct I values persist
with the same R,S. After any common binary word the differences
of actual counts remain exactly 2 and 4. For each append of
bit e, the next value satisfies 2R'=Z'+g e. Therefore `Z→2Z+g e`,
and I decreases by Z+g p on a 1, is unchanged on a 0, where
p is the preceding bit. The initial common 0 supplies p=0.

Let b be the ordinary forward address of a subsequent word of n
bits and D=2^n. Then Z=82D+9b. In particular the slope visits a
**complete arithmetic progression of D values**, not a fixed
list of odd parts. Requiring the final bit to be 1 already
leaves D/2 distinct odd slopes. Append one final 0 to read
the actual count. All three outputs are strong: the strong
core contains the root, every ordinary sink reaches it
through the final near-sink, and all vertices of the old
core other than the root beat that last vertex.

Both arbitrary addresses and the last append are identical across
the three outputs. To change the reader at a fixed address,
only old seed arcs need change; the root labels are already
aligned, the seed orders are the same, and all appended arcs
are untouched. With the displayed labels, the first low step flips
just one seed arc and the second flips three. Thus every interval below has literal bounded
rewires joining its three witnesses.

## Exact count formula for any synchronized state

Starting with any current Z,I and previous bit zero, apply an
n-bit address b followed by one final 0 (not included in b).
Let R_n(b) be the n-bit reversal, and J(b) the sum of
2^(d−1)+[d=1] over two set bits at distance d in b. Then

    H(b) = I + Z(2D − R_n(b)) + g(2b − J(b)).

Indeed each 1 pays the current Z plus g times the preceding
bit out of I. The leading initial Z pays the reverse address.
Each earlier 1 contributes its place weight to the leading
coefficient at each later 1, plus the adjacency surcharge.
Together they contribute gJ(b). At the final 0 the count is
I_new+2Z_new. The same formula holds when no final 0 has yet
been appended if one just wants the updated I from these
subtractions.

We use 0≤J(b)≤b: for each set upper bit i the sum of
interactions with lower bits is at most
(2^i−1)+1=2^i, with no contribution for i=0.

## Counting distinct absolute anchors rather than programs

Use a prefix of m bits, freely chosen, then a common 0, then
a suffix of 2m bits and the final 0. Put L=2^m and D=L².
The intermediate common 0 resets the previous-bit flag and
doubles the prefix slope, so the L distinct prefix slopes are
Z_p=164L+18p, with 0≤p<L. Their baseline I_p can depend on p;
this will not affect any comparison within the same prefix.

At fixed prefix, if two different suffixes b,c collide, the
exact formula implies

    Z_p (R_(2m)(b)−R_(2m)(c))
       = 9(2b−J(b) − (2c−J(c))).

The right side has magnitude at most 18D, hence their reversed
addresses differ by at most (18/164)L. Reversal is a
permutation of all D addresses, so there are at most
D floor((18/164)L) unordered eligible pairs. For each pair
there is at most ONE prefix where it collides: the coefficient
of Z_p is nonzero, and I_p cancels. Averaging over the L prefixes
gives one with at most (18/164)D colliding pairs. Therefore
that prefix alone has at least (1−18/164)D distinct starts,
as losing k−1 starts from a fiber of size k costs at least k−1
pairs. With each start retain its two other reader counts.
A greedy choice of disjoint triples loses at most a factor
of three, since starting counts are odd.

All these graphs have count O(LD), uniformly in the bits:
the prefix baseline changes by at most 91L+11 and the
displayed formula gives the bound directly. Since D=L²,
there are Ω(X^(2/3)) disjoint exact triples below X. To get
the bound at all sufficiently large X, choose the largest
power of two L with the uniform bound cL³≤X. This is a lower
bound on **distinct absolute counts**, not presentation
multiplicity or residue coverage. Finding the good prefix
by finite collision enumeration is a valid decoder, not yet an
efficient method for choosing a prescribed count.

## A target-first bracketing algorithm in absolute distance

For a single address (no extra prefix) the same construction
has an elementary approximate product law. Let x=b/D and
y=R_n(b)/D. Then

    H(b)/D = (82+9x)(2−y) + O(n/D),

uniformly. To verify it expand bR_n(b)/D over marked pairs:
the increasing pairs give the geometric part of J; the
opposite-direction weights sum to at most n/2 and diagonal
weights to at most n/2, while the adjacency correction is at
most n. Thus J=bR_n(b)/D+O(n) with no assumption of independent
addresses at full precision.

The first floor(n/2) bits of x and y come from disjoint ends
of the address word. Given any point (x,y) in the unit square,
take those many initial digits for each coordinate, placing
them at their respective ends (use the all-ones grid limit at 1).
This deterministically chooses the physical graph whose
normalized count differs from (82+9x)(2−y) by O(D^(-1/2)).
The continuous expression covers [82,182]. For a desired
large N choose D dyadic so that N/D is in [88,176).
Given z in that range, if z≤164 set x=0,y=2−z/82; otherwise
set y=0,x=(z/2−82)/9. There is a fixed margin to the edges
of the allowed range. Applying this recipe to N−C sqrt D
and N+C sqrt D with a sufficiently large absolute C yields
two literal witnesses on either side of N, both within
O(sqrt N) of it. Each starts an exact three-count repair window.

This bank no longer sits on a finite set of sparse rays and
there is a constructive approximate position for every target.
The unresolved task is to spend the O(sqrt N) remaining error
by a physical signed repair, or to extend the reader window
until adjacent positioned anchors overlap. A constant-width
window at a nearby anchor is not exact coverage.
