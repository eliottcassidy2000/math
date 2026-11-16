# Exact obstruction to a symmetric product odometer

The postulated radix sum is arithmetically valid. With k+1
four-state digits a_j in {0,1,2,3}, it produces every even
number from 0 to 3(3^(k+1)−1). In half-units, the old interval
through stage j−1 ends at B=3(3^j−1)/2; the four translated
intervals by 0,3^j,2·3^j,3·3^j touch or overlap because
3^j≤B+1. This proves the claimed range by induction.

But if the step at seam j is literally 2 times the product of
**all other seams' fixed factors** h_i≥1, this cannot give
1,3,9,… as half-step weights. As soon as one seam has weight 1,
all other h_i are 1. All the other weights must therefore be
the same number, namely the one remaining h_i. Even more
generally, with arbitrary weight ordering, k distinct
weights require k distinct h_i and then the minimum weight
cannot be one for k≥3. A physical odometer needs an asymmetric
or hierarchical read, e.g. each new rank scales a prior
prefix but not vice versa. This is analogous to a marked
Pythagorean parent/child transform retaining the ordered
parameters and cannot be replaced by a symmetric product
of side lengths.

## Strong obstruction for a particular rooted separator model

Use the earlier [polarized two-block head formula](../chanpak_signed_seam_bank_20261116/PRODUCT.md).
Two independently replaceable marked blocks D1,D2 with fixed
local Hamilton counts h1,h2 are separated by a nonempty
transitive spacer of length t. Fix outer transitive guards
with X=2^a,Y=2^b, a,b≥1. In terms of their two-slot profiles
h_i,q_i,E_i,L_i,J_i, the number of Hamilton paths is

    H = (q1 X−L1)(q2 Y+E2) 2^t
        + h1 h2 + h2(X E1−J1)−h1(Y L2+J2).

Take fixed finite digit banks in each position, require every
Cartesian pair to be permitted, and require their differences
to stay bounded for an unbounded sequence of spacer lengths.
Their high products must coincide, so positivity of both
factors makes each of q1 X−L1 and q2 Y+E2 individually constant
on its own bank. The remaining step on either bank is a
multiple of the opposite local count. All physical counts
are odd; since that local count is odd, the step is in fact
a multiple of twice the opposite count.

Suppose one bank provides an exact +2 step and the other is
also nontrivial. The opposite count must then be 1, so its
block is a transitive tournament (any directed triangle
already contributes two additional Hamilton paths by
the positive odd-cycle-packing formula). But a marked
transitive tournament cannot move on its own protected factor.
For order n, q=2^n, and L is the ordinary binary mask address:
a marked vertex i in forward indexing 0,...,n−1 contributes
2^(n−1−i) because its own slot can include any later sites.
Similarly E is the reverse mask address. Consequently

    q X−L ∈ ((X−1)2^n, X·2^n],
    q Y+E ∈ [Y·2^n, (Y+1)2^n).

Within each order, the mask is determined uniquely. The
displayed intervals are disjoint between distinct orders
(X,Y≥2), so allowing variable transitive size does not help.
The second bank is physically constant up to relabeling.

Thus this particular Cartesian two-bank separator cannot
carry even one +2 window beside another independently
varying search register with persistent bounded steps. The
proof does not apply to unmarked padding next to one bank,
to t=0 (the live endpoint term changes), or to a joint
state machine in which the two digits may move together.

## Even with no spacer, the two factor counts 1 and 3 obstruct an additive width-four pair

Keep outer guards fixed with X,Y≥2 powers of two and use
the adjacent marked-block formula (spacer removed):

    H=h1 h2+q1q2XY + X(q1E2+h2E1) − Y(h1L2+q2L1)
      − h2J1−h1J2−E2L1−T1S2.

Here T counts Hamilton paths ending in the block's marked set,
and S counts those starting in it. Let the first block have
Hamilton count 1 and the second have Hamilton count 3. Each
block's order is fixed, although the variants on the vertices
may differ. The first block is transitive; the second consists
in strong-component order of one directed triangle, with
transitive prefixes and suffixes (possibly empty). This last
claim follows from multiplication over strong components:
only one nontrivial component can occur, and a strong
tournament with more than three vertices has at least two
triangles (start with a triangle; if every exterior vertex were
uniform to it and no other triangle existed, its condensation
would make the tournament nonstrong). Alternatively it follows
directly from Busch's lower bound H≥5^((n−1)/3).

In the count-3 block, q,E,L are multiples of 3, regardless
of triangle position and masks. Thus for a fixed second-block
setting, every term depending on the transitive first-block
setting is a multiple of 3 except −T1S2. The quantity T1 is
binary. Hence the transitive side gives at most two residues
modulo 3 and cannot itself provide four consecutive odd counts.

Now let the transitive side have ≥3 distinct mask states,
and suppose exact additive independence, meaning the
mixed differences of H in the two steps vanish in the
Cartesian grid. A change of the first block varies its
L coordinate (an injective address) and its binary T.
These variations span two dimensions unless T is fixed;
in that exceptional case they still span the L-direction.
If T is fixed, the mixed difference forces E2 fixed,
which is all the shadow argument below uses except that
the term T1 S2 appears in the second block readout.
With T variable, mixed differences force both E2 and S2 fixed.
Thus in either case E2 is fixed. If T1 is fixed, it can be
absorbed exactly into the prefix as one fixed mark bit;
the proof below writes that bit as b in the formula.
If T1 varies, the S2 term is constant and can be ignored
(set b=0). This corrects the tempting but false claim that
three distinct TT masks automatically give two-dimensional span.

### Exact shadow for a moving triangle

For any fixed marked count-3 block order N+1, encode
the transitive vertices as ordinary bits (forward
address E, reverse address L), and replace its triangle
by two bits u,v with s=u+2v, where s is the number
of marked triangle vertices. The resulting N-bit shadow
Z has ordinary forward address E2/3. Thus fixing E2 fixes
the whole shadow; only the position of the two-bit triangle
encoding can vary. Let r be its number of following bits,
let E_R be their forward address, and let T_L,S_R be the
adjacent bits on the left and right, taken as 0 if absent.

Write K=Y L2+J2+b S2 with b∈{0,1}. Its shadow is 3 times
the corresponding functional on the transitive tape Z.
The exact difference is

    K−3K_shadow
      = 2 [u=0,v=1] + (v−u)(3Y·2^r+3E_R+2T_L−S_R)
        + b(S2−3 S_shadow).

To prove this, compose the marked profiles:
h=h1 h2, q=q1 q2,
E=q1 E2+h2E1, L=q2 L1+h1L2,
T=h1T2, S=h2S1,
J=h2 J1+h1J2+E2L1+T1S2.
The triangle has (h,q,E,L,T,S,J)
=(3,12,3s,3s,s,s,s(s−1)), and the corresponding TT2
block is (1,4,s,2u+v,v,u,2uv).
Substitution gives the displayed formula. The last term
is zero except when the triangle touches the first site:
otherwise S2=3S_shadow since its first marked singleton,
if present, starts each of the three whole Hamilton paths.
If the triangle starts the tape, that term is
b(s−3u)=2b(v−u). It simply increases the number D_r
in parentheses by 0 or 2 at its largest possible r.

There is one triangle position per r. For r=0, its D lies in
[3Y,3Y+2], before the optional increase at the largest r.
For r≥1, D lies in [3Y·2^r−1,(3Y+3)·2^r−1].
These bands have successive gaps at least 3; increasing
the largest-r value preserves that separation. Thus the
corrections 0, −D, D+2 across all positions cannot differ
by 2. Therefore the count-3 block at fixed E2 gives no
exact +2 step, even allowing its triangle to move.

Combining the residue and shadow tests: in this precise
adjacent h1=1,h2=3 architecture with an additive Cartesian
grid of two nontrivial banks and ≥3 states on the TT side,
neither side supplies the proposed width-four unit window.
Converse reverses the order and swaps E with L, giving
the same result when the two factor counts are swapped.
This does not preclude joint non-Cartesian changes, changing
core orders, or additional asymmetric ports.
