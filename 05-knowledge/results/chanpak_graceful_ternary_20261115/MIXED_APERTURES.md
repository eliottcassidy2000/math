# A mixed-scale ternary tree with an exact target-first list decoder

There is a way to mix *some* square-sum apertures without the Pell
commutation collision. It uses the two host labels, not just
their height. This gives a denser exact window family than the
single-aperture ternary example, while retaining a rigorous sparse
support boundary.

For an absorber triple (a,b,c), a²+b²=c²+1, write S=x+y. Its move
is (x,y)↦(a²S−x,b²S−y). Thus S'=c²S and the normalized fraction
r=x/S lies after the move in the open interval

    ((a²−1)/c², a²/c²).

Take these three letters:

    (4,7,8), (7,4,8), (7,11,13).

Their intervals (15/64,16/64), (48/64,49/64), (48/169,49/169)
are pairwise disjoint. All have min(r,1−r) in (1/5,1/3).
Given the ordered host pair, its ratio determines the only possible
last letter, even when the height admits both apertures. Divide
S by the resulting c², then invert x to a²(S/c²)−x. Each step must
remain integral and positive, which is the exact square-sum
admissibility test; the height strictly decreases. Accept only at
the prescribed seed (1,3). Thus the host graph is an actual free
ternary tree. Disjoint ratio intervals, rather than the abstract
Berggren coefficients alone, are what make this lift invertible.

## Decode an absolute tournament count without enumerating the tree

Send each pair to the same literal three-center bank as in
[SQUARE_TREE.md](SQUARE_TREE.md). At each fixed sum, the ordering
theorem from that note and the [discriminant inverse](DISCRIMINANT.md)
decode a candidate absolute count uniquely to an ordered pair.
At positive depth, coarse bounds give 35 S²≤H_j≤700 S² (the wider
ratio interval used in the earlier proof also contains these
three new ones). The possible sums are 4·8^(2a)·13^(2b).
For a prescribed N, loop over O(log N) possible a and use the
constant multiplicative interval between sqrt(N/700) and
sqrt(N/35) to find the O(1) possible b per a. At each sum use
the discriminant inverse, then the ratio descent and exact seed
test. Every accepted witness is physical and exact. At most one
witness occurs for each sum; do not assume scalar injectivity at
different sums, nor test just a single preferred height. The
depth-zero seed is separate.

## Counting the gain and its limit

Let σ be the positive solution of

    2·8^(−4σ) + 13^(−4σ) = 1.

There are Θ(X^σ) physical windows up to height X **with host
multiplicity**. Indeed the height of each positive-depth window
is within fixed positive factors of its weight W, the product of
c^4 along its word. Stop the weighted tree at the first W>Y.
Assign each child probability c^(−4σ), with the two c=8 letters
distinct. The stopping leaves partition the infinite words and
have total mass 1. Each leaf weight lies in (Y,13^4Y], so there
are Θ(Y^σ) leaves. The full finite ternary tree before these
leaves has a fixed linear relation between leaves and internal
nodes. This counts all words with W≤Y and proves the assertion.
The sum decoder bounds the multiplicity of any scalar by
O(log X), so the number of distinct counts is between
Ω(X^σ/log X) and O(X^σ). In particular σ>log(3)/(4log(18)),
since both apertures are strictly smaller than 18.
The support is still sparse (indeed the broader square-sum
bound applies). The constructive improvement is permitting
variable scale without losing an exact inverse, not a proof that
all absolute targets can be assigned to one of these scales.
