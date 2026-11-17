# A repeatable +2, a compensated ×7, and an exact Stern–Brocot inverse

**Status.** Three separate all-order results: (i) a single-arc +2
switch continuing the actual order-eight 289/291 pair into the
7p gap channel; (ii) a strong square family for every exponent,
with an obstruction to obtaining that entire family by one-site
cloning; (iii) an autonomous ×7 switch in a two-slot cyclic
construction, plus an exact inverse for a different one-apex
channel. The +2 switch does not stay at the Fermat squares;
a universal absolute selector and Chan–Pak remain open.

## 1. Which part of the Fermat-square flip repeats?

The literal [flip list](flip_edges.json), replayed against a full
archived isomorphism census through order eight, gives 20, 3, 2, 6
lower classes with a +2 reversal at p=1,2,3,4. These are numbers of
source classes, not edge counts. The [replay](replay.py) verifies
each edge independently by literal subset path counting. This
table is bounded; the following constructions are not.

Let X be a marked tournament vertex. Replace it by a module B,
where outsiders see all of B just as they saw X. For r possible
visits, let a_r count quotient words with every unmarked vertex
once and X exactly r times. Then

    H(G[B]) = sum_r a_r c_r(B),

where c_r counts ordered nonempty r-path covers. Read the runs in
appearance order for a bijection. This is the inherited
[walk-content transfer](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
In particular c_1=u=H(B), c_2=R−2u for R=P_B(2) with empty
slots allowed. For a transitive B of m vertices, c_r=r! S(m,r)
because the vertices in any one run have only their increasing order.

**An unbounded physical +2.** Take outgoing rows
(8,13,25,176,99,71,143,55), flip arc {4,7}, and replace vertex 0
by TT_m in both graphs. The complete two slot coefficient lists
are (289,151) and (291,151). For a short certificate, permute the
seven nonmodule vertices, classify each gap for X as impossible,
optional, or mandatory, then take the coefficient of z^r in
z^(number mandatory)(1+z)^(number optional). The replay performs
that enumeration. So the count pair for every m≥1 is

    151·2^m−13, 151·2^m−11.

The difference stays exactly 2 under repeated extension of the
cloned spine; the changed arc is external to it. For general B
the difference is 2u, so it stays 2 precisely when B is transitive.
The upper count is divisible by 7 iff m=0 mod 3; for m=0 mod 6 it
is also coprime to 3 and lies in I_(m+3)=[18·2^(m+3),28·2^(m+3)].
The lower count is never divisible by 7. At m=18 the upper count
is 7·5,654,819 with a verified prime quotient. This is a genuinely
new infinite repair family seeded by the p=4 flip, not an
infinitude claim for its prime quotients.

**An unbounded strong square.** The seven-vertex rows
(118,124,120,49,96,64,8) with vertices 1 and 4 replaced by modules
A,B give two-visit coefficient matrix (25,10;10,4). Flipping
{2,3} gives (27,12;9,4); rows/columns index one or two visits in
A and B respectively. For vertex 1 there is only one ingress,
vertex 0. For vertex 4 its only exits are 5 and 6; using exit 5
before a return consumes 6 immediately, since 5 only beats 6.
If exit 6 is used first, exit 5 cannot enable a second return.
This caps both visits at two in either scaffold. Literal small
inflations solve the four coefficients by triangular inversion.
In optional-cover notation the lower H and the difference are

    H=(u+2R)(v+2S),       Δ=2uS−Rv.

For transitive modules of lengths i,j, Δ=2·2^j−2^i. Thus i=j=p−1
gives every square (2^p+1)^2 for p≥2, and i=j+1 gives an exact
neutral reversal on neighboring Fermat rectangles, but equal
lengths yield +2 only at p=2 for this marked edge. If both
modules have at least two vertices, both R and S are multiples
of 4 (pair each subset with its complement, then sum an even
number of odd terms), so this particular edge cannot add 2
with larger nontransitive substitutions either.

## 2. Two all-order obstructions that narrow the square problem

**An intercomponent arc cannot give +2 once H>3.** Put the arc
x→y across an ordinal cut A→B. Reversing it loses
the unweighted end-at-x count times the unweighted start-at-y
count. The new paths using y→x have one optional A-run before
that crossing and one optional B-run after it. Consequently
the gain is the product of the *weighted* start menu at x
and the weighted end menu at y; each menu weights a partial
path by the Hamiltonian count of its complement. No gained
path needs the removed x→y arc since it already uses y→x.
Write the menus as even/odd totals E_A,O_A and E_B,O_B.
The signed start-menu identity O−E equals the unweighted
end count at that vertex; converse gives the dual identity.
Thus the exact change is

    Δ = 2(E_A O_B + O_A E_B).

If both even totals are positive, the odd totals are at least 1
and Δ≥4. If E_A=0, x is a sink in A, O_A=H(A), and Δ=2
requires H(A)=1 and E_B=1. The dual case is symmetric.
A weighted even end total 1 forces a unique in-neighbor u
of y and a transitive complement of {u,y}. In that complement
u can lose only to the first vertex: otherwise an earlier
vertex gives a four-vertex path ending x'→u→y, a second
even contribution. Hence B is either transitive, or a
triangle followed by transitive vertices, and H(B)≤3.
This proves the claimed exclusion. An internal reversal
in one SCC also multiplies its difference by the counts
of all other SCCs. Thus a +2 square reversal must occur
inside the *unique* nonsingleton SCC. The ordinary ordinal
square of two Fermat blocks cannot serve as the source.

For a self-contained proof of the signed menu identity,
work with square-zero variables x_i, X=diag(x_i), and adjacency
matrix A. The signed start series at vertex v is
x_v[(I+AX)^(-1)1]_v; multiply by the total (with empty path)
to weight complements. To compute the end series use
A^T=11^T−I−A. The diagonal X term in I−A^T X contributes
only repeated variables and can be dropped after extracting
x_v times the resolvent row. Rank-one inversion then gives
the same signed start series divided by
1−1^T X(I+AX)^(-1)1. Summing the end series and adding the
empty path identifies this inverse denominator with the
total path series, by reversal. Extract the full squarefree
monomial. This is the inherited O−E relation with an
independent algebraic derivation.

**No fixed skeleton with just one transitive growing site can
hit infinitely many Fermat squares**, even at unrelated values
of m and p. Let r be the maximum quotient multiplicity of X.
Permuting its r−1 internal return runs gives a free action,
so (r−1)! divides a_r. Expanding the Stirling read in powers gives
H(G[TT_m])=a_r r^m+O((r−1)^m) if r≥2. If infinitely many
pairs m,p hit (2^p+1)^2, then p is proportional to m and

    log(a_r) + m log(r) − p log(4)

is exponentially small in m. A nonzero linear form in these
fixed positive rational logarithms is bounded below by an
inverse power of m, by [Matveev's bound](https://www.mathnet.ru/eng/im314).
Thus for all sufficiently large hits a_r r^m=4^p. If r=3,
3 divides r; if r≥4, 3 divides (r−1)! and hence a_r. Both
are impossible. If r=2, H=a_2·2^m+a_1−2a_2, so the
constant would have to equal 2^(p+1)+1, impossible infinitely
often. r=1 is constant. In the synchronized case m=p+c the
obstruction is elementary, directly from comparing the leading
bases and (r−1)! divisibility. Changing the skeleton, adding a
second growing site, or using nontransitive growth is essential.

## 3. The prescribed target chooses its Stern direction and scale

For the one-apex over a transitive spine, split a word into
left A and right B. Let x be the binary integer A, and y the
binary integer of the reverse-complement of B. The inherited
inversion rule for one apex gives the exact seam read

    H=xy+F(x)+F(y)−1+[x,y odd],
    F(0)=1, F(2z+1)=F(z), F(2z)=F(z)+z+(z mod 2) for z≥1.

The interior count of B equals F(y) by converse. The cross
inversion weight factors as xy, and only the adjacent seam
adds the parity indicator. In particular 1≤F(z)≤z+1,
and odd z has the stronger bound F(z)≤(z+1)/2.

**Exact rational-direction inverse.** Write (x,y)=d(a,b)
with gcd(a,b)=1 and positive a,b. Follow the ordinary
Stern–Brocot address or Euclidean quotients of a/b, but
keep the common scale d. For any prescribed H=N a direction
can support at most one scale, forced by

    ab d²+1 ≤ N ≤ ab(d+1)²,
    d = floor(sqrt((N−1)/(ab))).

For the upper bound use (a+b)d+2≤2abd+ab when ab≥2.
For a=b=1, the stronger odd bound on F or the vanishing
even adjacency term gives the same inequality. The last
step is an equality check, not a congruence search. This
is the constructive version of [THM-4057's
direction/scale distinction](../../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md);
Stern depth alone never retained the scale.

**Complete one-apex inverse.** After stripping ordinal
padding, split the active word in half. Its first bit is
1 and last bit 0. Up to swapping the halves by converse,
the binary weights satisfy 1≤x≤y<4x. Thus
4x²+5x+2≥N≥x²+1, and for each x one need only test

    ceil((N−F(x)−1)/(x+1)) ≤ y ≤ floor((N−F(x))/x).

Here x is comparable to sqrt(N), so there are O(1) candidates
per x. This is a complete selector for one-apex existence
using O(sqrt(N) log N) arithmetic steps, or a scoped certificate
of absence. It is not polynomial in bit length. For the
predetermined p=30013, N=210091 in a dyadic gap, the decoder
selects (x,y)=(275,762). The continued fraction of 275/762
is [0;2,1,3,2,1,2,1,5], scale 1. The literal word
1000100111010000010 gives 275·762+61+481−1=210091.
The scale-restoration example (a,b)=(1,1), N=91, forces
d=9 and realizes 91. This is a new target-first physical
selection, not an assertion of an unrecorded count.

**One-apex contact flips are more rigid.** In a lower word
A0B, use the same x,y. Flipping the middle 0 to 1 changes H by
2(ceil(y/2)−ceil(x/2)); the lower H is
2xy+x+(x mod 2)+F(x)+F(y)−1. Up to converse this describes
all flips at the apex. For +2 put U=ceil(x/2); then x is
2U or 2U−1, and y is 2U+2 or 2U+1. For U≥1 the lower H
lies strictly between 8U² and 8(U+1)² by the bounds on F.
The U=0 possibilities below H=9 are 1 and 3.
Thus U=floor(sqrt(N/8)) and only four signatures need
checking at any height. Replay through p=500 finds p=1,2
only in this restricted lens. It does not exclude higher
p or reversals not touching the apex. An adjacent
transitive-spine swap has difference
2+2 floor(x/2)+2 floor(y/2), where x,y are now the
two external halves around the swapped 10/01. It is a unit
only when both sides have weight at most one.

## 4. From a geometric factor 7 to a physical ×7 switch

For cyclic gluing C3[A,B,1] the read is
(R−u)(S−v)+2uv. Set a=R−u, b=S−v. Then the target equation

    ab+2uv = 7p = det((a,2u),(-v,b))

has a literal Euclidean inverse. Put g=gcd(a,2u).
A reduced Farey parent z of (a/g,2u/g) of determinant 1
can be found by Euclid; all candidate (-v,b) are (7p/g) z
plus an integer multiple of (a/g,2u/g). Retain b,v positive
and the physical profile (v,b+v). The earlier [proper-module
obstruction](../seven_gap_atlas_20261117/README.md) forces
g∈{1,7} for prime p≠7, and at most one arm can have common
factor 7. If 7 divides neither u nor v, the only unordered
pairs (R/u,S/v) mod 7 that hit the target channel are
(0,3), (2,6), (4,5), from (R/u−1)(S/v−1)≡5 mod 7.
Ordinal padding rotates ratios in the square class since
2 has order 3 mod 7. This classifies the *local phase*,
not the integer endpoints of the determinant line.

**A physical all-order factor-seven move** is in
[FACTOR_SEVEN.md](FACTOR_SEVEN.md) with [independent
replay](factor_seven.py). Change three ordinal constituents
of the second arm and add a singleton; this scales its
actual pair (v,S) by seven, from (56551,3390208000) to
(395857,23731456000). Two unchanged-H swaps cancel unwanted
odd prime factors in S, and the singleton resets the dyadic
factor. Because the cyclic read is linear in both
coordinates of that arm, it multiplies the *whole output*
by 7 for every A. This also works in the other arm by the
symmetric expression. It does **not** hide a p-path
module inside the larger tournament: the smaller p was
already a three-block output. For transitive A of orders
1,4,6 the smaller outputs are verified primes
3,390,264,551, 50,852,384,837 and 213,579,654,389,
and all three sevenfold outputs lie in actual dyadic gaps.

The fixed reader reduces a prescribed p to a smaller
physical profile (u,R) with
3390151449 R − 3390038347 u = p. Its coefficients are
coprime: gcd(S−v,3v−S)=gcd(S−v,2v)=1 for this arm.
But no arithmetic argument here forces a physical point
on that line. In fact the finite-mask arithmetic cannot
be universal on its own: [the density bound](FACTOR_SEVEN.md)
excludes the entire two-arm grammar generated by arbitrarily
many ordinal one-apex SCCs. The next target is a bank with
target-dependent primitive directions and richer physical
profiles, not yet another congruence lift.

## 5. 35 square-sum graphs and the real recurrence behind 768

See [ARCHIVE_MAP.md](ARCHIVE_MAP.md) for the history pass.
The first 35 square-sum graphs give a useful hostile model
for a scalar count: Q_27 has **35 unordered Hamilton paths**
but only **16 endpoint pairs**, all sharing endpoint 18.
The multiplicities are 4 pairs of multiplicity 1, 6 of
multiplicity 2, 5 of multiplicity 3, and 1 of multiplicity 4.
Vertex 18 has only neighbor 7 through Q_30, so the endpoint
is structural. It is allowed a second neighbor 31 via
18+31=49. The count-35 coincidence with the gap atlas is
not a pairing between problems. What transfers is that
count and continuation freedom are distinct. Compare the
other early forced ends: Q_15 has the unique 8-to-9 path;
Q_16 and Q_17 extend it, and Q_18 cannot extend it because
the new leaf attaches at interior vertex 7. See
[THM-4505](../../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md).

**Exact square-clock intersection.** The inherited
square-sum new-edge clock at n counts squares in [n+1,2n−1].
At a square n, a lost label is exactly replaced by a
new square when 2n−1 or 2n−2 is square. For
n=(2^p+1)^2 this occurs **exactly at p=1,2,4**
(n=9,25,289), not at 81. To see it, put x=2^p+1.
If 2x²−1=l², the consecutive integers (l−1)/2,(l+1)/2
multiply to 2^p(2^(p−1)+1). For p≥3 both are less than
2^p, yet the even member must be divisible by 2^p.
The remaining odd-branch check gives p=2. If
2x²−2=l², then (l/2)²=2^p(2^(p−1)+1).
For p≥2 the odd factor must be square; its root minus
and plus one are even powers of two differing by two,
so they are 2,4 and p=4. p=1 is a separate even-branch
hit. This is a genuine Pell/dyadic intersection, not a
criterion for tournament +2 flips.

The two 768·2^k−5 encodings from the gap atlas give the
closest structural control: one has a rooted doubling
clock, the other has a two-visit coefficient 3·256.
Equal output across all heights does not identify the
geometry. The compensated ×7 move exploits the complementary
fact: equal *full interface coordinates after scaling*, not
equal scalar counts, guarantee behavior in every exterior.

Run `python replay.py`, `python factor_seven.py`. Use
`python replay.py --square-sum` to replay all prefix counts
through 35 rather than only the Q_27 control. Add `--census` to
rebuild all isomorphism classes through order eight independently.
