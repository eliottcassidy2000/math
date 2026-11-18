# Strong blue probes that preserve the exact count as decodable information

For a tournament X write h=H(X), and let q=P_2(X) count ordered pairs of
directed paths partitioning X, with empty paths allowed. For any
nonempty X and any self-converse B (including empty), build the ordinal
sum a -> X -> B -> X^op -> b with singleton ends a,b. Reverse just a->b.
The result W(X,B) is strong: all other vertices lie on a directed
cycle a->v->b->a. Reflect a and b, exchange the labeled X and X^op,
and use an anti-involution of B. This is an anti-involution of W.
When B has a given blue path, a blue frame of W is given by a, any
Hamiltonian path of X, the blue frame of B, the reflected reversed
first path, b. Empty and small B used here have such a frame, and the
existence for general B also follows by pairing paths and Redei.

All Hamiltonian paths avoiding b->a must consume the blocks in
ordinal order. Their number is h² H(B). Every path using b->a has a
pre-path ending at b and a post-path starting at a. Restrict them to
each interior block, obtaining its ordered two-slot covering. The
choices are independent: every slot concatenates the selected pieces
in ordinal order, and a and b can never be reused. This is a bijection.
Therefore

    H(W(X,B)) = h² H(B) + q² P_2(B).

No uniform-converse assumption on X is needed. With B empty and
singleton respectively, the two counts are Y0=h²+q² and Y1=h²+2q²,
so h=sqrt(2Y0-Y1) and q=sqrt(Y1-Y0). Thus two Hamiltonian count probes
of size at most 2|X|+3, **both strong and blue**, recover the exact
arbitrary tournament count. This is a computational reduction, not
a spectrum inclusion. If desired, use B=TT3 and C3 (all three in a
directed Hamilton frame). Their profiles are (1,8) and (3,12):
the 12 includes 6 nonempty split coverings and 6 two-slot coverings
with one empty side. The counts Y,Z now come from **equal orders**
2|X|+5, and the probes differ by a single central arc outside their
common blue frame, with

    h²=(2Z-3Y)/3,        q²=(3Y-Z)/12.

These are ordinary chord neighbors, not a claim that the probes form
a native metagraph line. For one probe instead take B=TT_s with
2^s greater than h²; then h=sqrt(H(W) mod 2^s). Taking 2^s>(|X|!)²
gives a source-independent budget s=O(|X| log |X|).

A limitation of a decoder is visible in these same formulas. For
B=TT_s, the represented values are h²+2^s q², with q even for
nonempty X (swap the ordered slots freely). In particular they
all equal 1 mod 4. Changing the blue middle block B changes the
actual arithmetic restriction, although the two-coordinate read
remains exact. The companion mirror-apex construction supplies
another independent probe and a different obstruction.

Run `python check.py` for a direct Hamilton subset-DP replay across
all labeled sources through order four, with four middles.
