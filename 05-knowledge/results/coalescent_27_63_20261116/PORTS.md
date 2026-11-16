# Scalar synchrony versus physical square-sum routes

Write Q(N) for the graph on 1 through N with an edge between distinct
labels whose sum is a square. The following are endpoint-resolved
routing statements, not maps between the invariants of tournaments,
Collatz, and square-sum graphs.

## One Collatz meeting, two different routing regimes

For the supplied Collatz shortcut fiber, let
```
x=27+1024j, y=63+2304j, D=9+320j, j>=0.
```
The two indicated base parity itineraries recur with an exact
meeting at 182+6561j.
The replay checks the 10- and 8-step itineraries and their affine slopes;
it is not a claim about unbounded Collatz dynamics.

Since y−x=4D, the two ports have a square-sum common neighbor
```
z=(D−1)^2−x=37+4096j(1+25j).
```
The two edge sums are (D−1)^2 and (D+1)^2. At j=0 this is the
physical route 27–37–63 with sums 8²,10².

If z' is any positive common neighbor, set u²=x+z', v²=y+z',
with 0<u<v. The equation v²−u²=4D implies u and v have the
same parity, since distinct parity gives an odd difference.
Thus a=(v−u)/2 and b=(v+u)/2 are positive integers with
a<b and ab=D, and conversely every such factorization gives
u=b−a, v=b+a and a candidate z'=(b−a)^2−x. Discard it if it
is nonpositive or repeats an endpoint (loops are not allowed).
So **the factorization selects the physical activation height**.

If D is prime, a=1,b=D is the only admissible factor pair: no
other common neighbor is possible, and the explicit z is valid.
It costs quadratic height in j, rather than the affine scalar
meeting height. There are infinitely many such j by [Dirichlet's
theorem, original proof in translation](https://arxiv.org/abs/0808.1408),
since gcd(9,320)=1. This is a two-edge restriction, not a lower bound
for arbitrary longer paths.

Conversely choose t=1+160h, j=9(t²−1)/320, h>=0. Then j is an
integer, D=9t², and a=t,b=9t give a different common neighbor with
square sums (8t)^2,(10t)^2:
```
z'=(176t²+9)/5.
```
It is strictly between x and y, even at h=0. Its activation
height is exactly y: there is **no overhead** beyond the endpoint.
Thus a single explicit synchronization family contains both
infinite quadratic-cost and zero-overhead two-edge regimes.

## Ports 8 and 9

In the inherited first spanning path on Q(15), the two leaves 8 and
9 are forced ends (THM-4505). Their port pair cannot have an edge
(sum 17), or a common positive neighbor at any height: positive
consecutive squares cannot differ by one. Their first actual
three-edge route is
```
8 – 41 – 40 – 9; edge colors (49,81,49).
```
For exact minimality, in Q(41) the two neighbor sets are respectively
{1,17,28,41} and {7,16,27,40}; their only square-sum cross pair
is 41+40=81. Removing vertex 41 therefore removes all three-edge
routes. This onset is about shortest port closure, not Hamilton
paths at smaller heights.

The route is a reflection-conjugacy identity: applying x↦49−x
to the consecutive port pair turns sum 17 into 2·49−17=81.
More generally an equal-outer-color three-edge route between
ports x,y requires the
necessary identity 2u²−v²=x+y, where u² is the outer color
and v² the central color. Such a representation is physically
valid only if the two intermediate labels u²−x,u²−y are positive
and all four labels distinct.

An exact bridge to the inherited golden-square switch: the pair
(8,9) encodes the upper edge of the first (7,13,17) four-cycle,
since 2·8·9=144, 8²+9²=145, and 144+145=(8+9)²=17².
The full switch (24,25,144,145) is inherited; its hypotenuse
17 is also that of the earlier (8,15,17) reflection graph.
This identifies an edge through arithmetic, not a spanning-count map.

The middle pair 40,41 is also part of the right triangle 9,40,41:
41²−40²=(41−40)(41+40)=9². That is exactly the triangle
behind the inherited repeated-square template of roots (31,41,49),
which contributes one of the two switches at 30,013.
One can keep this geometric relation in an infinite family.
Set (v+1)+u sqrt(2)=(2+sqrt(2))(3+2sqrt(2))^t for t>=1.
Its norm is 2, its integer part is even, and therefore v is odd
and u²=v+(v²−1)/2. The path
v−1, (v²+1)/2, (v²−1)/2, v
has colors u²,v²,u². At t=1 this is 8–41–40–9; at t=2
it is 56–1625–1624–57. No general onset claim is made.

## Ramanujan, taxicab, and the effect of releasing a color constraint

The 1472 reflection obstruction from the prior graph is
−tau(4)=2·32²−24². The taxicab identity is
1729=1³+12³=9³+10³=3·24²+1. Hence
```
1729−1472 = 257 = 1+2304−2·1024.
```
The same two affine periods in the Collatz fiber give the Fermat-prime
gap after doubling the shorter period, plus one. This is an exact
identity; it supplies no dynamics on its own.

The gap 257 has only one representation as a difference of positive
squares: 129²−128². Indeed factoring v²−u²=(v−u)(v+u)=257
forces v−u=1,v+u=257. Thus 1472 and 1729 have exactly one common
neighbor, 14912. But allowing distinct outer colors gives the
three-edge route
```
1472 – 49 – 120 – 1729; edge colors (39²,13²,43²).
```
It already exists at the minimum possible graph height 1729.
At that height the distance is exactly three: sum 3201 lies strictly
between 56² and 57², and the sole common neighbor is unavailable.
The two-edge distance first becomes attainable at 14912.

In fact **no** three-edge route with equal outer colors exists
for these ports at any height. It would give 2u²−v²=3201.
Modulo 3 this forces 3|u and 3|v, whereas 9 does not divide 3201.
Thus releasing a color-equality constraint here is necessary, not
just one way to improve a finite search.

## An irrational ratio from exact physical closures

For k>=5 with k≡2 (mod 3), let u=F(k+1), v=F(k), with
F(0)=0,F(1)=1. Fibonacci parity gives u even and v odd. Put
```
p=(v²−1)/2, q=p+1, x=u²−q, y=u²−p=x+1.
```
These are four distinct positive labels. The path x–q–p–y has
colors (u²,v²,u²) and lies within height y. It is a shortest
route at every height >=y: the ports are consecutive, excluding
a common neighbor, and x+y=2u²−v²≡7 (mod 8), excluding an edge.
Cassini gives u²−uv−v²=(−1)^k, and therefore
```
x+y = v F(k+3) + 2(−1)^k;
(x+y)/(p+q) = F(k+3)/F(k) + 2(−1)^k/F(k)² → phi³=2+sqrt(5).
```
The irrational limiting ratio has an exact finite edge certificate.
It does not yield integer densities or show that asymmetric
tournament ports with comparable scalar values physically embed.

The two decompositions of 1729 are the divisor pairs 13·133 and 19·91:
the Eisenstein cofactors are 7·19 and 7·13, by the exact
[THM-463 split-axis criterion](../../../01-canon/theorems/THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis.md).
That earlier theorem separates this cubic/Moser structure from any
unproved tournament split-prime law. See also
[Ono and Trebat-Leder, The 1729 K3 Surface](https://arxiv.org/abs/1510.00735).
The natural repair analogy is to switch an entire factor/port profile,
not to assume equality of scalars gives a continuation. The tournament
experiments in TOURNAMENT.md make that distinction literally observable.
