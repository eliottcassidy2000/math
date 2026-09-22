# Actual Collatz edges as a thin, angle-restricted family of primitive triangles

**2026-09-21. Status: PROVED scoped elementary results; FINITE-EXACT
controls. No global Collatz convergence, cycle completeness, or novelty
claim.** Positive states and odd acceleration are essential below.

## Inheritance and concept board

The closest proved mechanism is the coprime odd-square chart of
[THM-3756, odd-square ordinal Berggren descent](../../01-canon/theorems/THM-3756-odd-square-ordinal-berggren-affine-descent.md),
together with the exact affine guard in
[arithmetic braids](arithmetic_braids_20260917_collatz.md).
The canonical hostile is a legal edge whose Berggren parent is not legal.
The corrected near miss is transferring descent from the full triangle
tree to a subset it does not preserve. The least-used coordinate is the
ordered pair of square roots, with source direction retained.

Anchor: legal Collatz edges and a geometric encoding. Niche: exact angle
exclusion. Wildcard: compare the induced edge measure with the ambient
primitive-triangle measure.

| Object | Operation / invariant | Retained coordinate or test |
|---|---|---|
| Odd affine edge | 2^k y=3x+b | b, actual k, source and target |
| Primitive triangle | roots of c+/-even leg | gcd and direction |
| Berggren descent | strictly smaller parent | Does the parent remain a legal edge? |
| Shape | angle, projections and altitude | Lost scale and parity labels |
| Height population | count below c<=X | Edge sampling differs from lattice sampling |

## 1. A genuine edge encoding

For b=+1 or -1, let x,y be positive odd integers with
`2^k y=3x+b`, k>=1, and x!=y. Since powers of two are units modulo x,

```text
gcd(x,y)=gcd(x,3x+b)=gcd(x,b)=1.
```

The odd-square chart therefore gives the primitive triangle

```text
E(x,y)=(a,b_even,c)=(xy, |x^2-y^2|/2, (x^2+y^2)/2).        (E1)
```

Here `c+b_even=max(x,y)^2`, `c-b_even=min(x,y)^2`.
The two roots and their source/target order recover the edge exactly;
the unmarked triangle forgets the direction and the dynamical parameter.
With those roots recovered, legality is the exact power-of-two test
`(3x+b)/y=2^k`. The exponent is unique when the test succeeds.
Fixed edges at 1 would produce a zero leg and are excluded.

The same identity for general odd parameter b is
`gcd(x,y)=gcd(x,b)`. Dividing both roots by that content gives a primitive
triangle and divides its three sides by the content squared. The
[inverse/content lane](odd_square_20260921_inverse.md) states the precise
parameter transport, including b=-5.

This is an encoding of edges, not a claim that the Berggren map is
Collatz. On an actual sequence x,y,z, the encoded hypotenuses satisfy

```text
c(y,z)-c(x,y)=(z^2-x^2)/2.                              (E2)
```

Thus decrease of this geometric height is exactly two-step descent of
the original orbit. It need not occur: 7->11->17 under U_1 gives
85->205. The U_(-1) cycle 5<->7 gives constant height 37.

## 2. A sharp forbidden near-isosceles region

Let theta be the smaller acute angle of E(x,y), equivalently the
altitude-to-leg angle specified in the user's smaller subtriangle. Put
`t=min(x,y)/max(x,y)` and `alpha=arctan(t)`. Then

```text
theta=min(2alpha, pi/2-2alpha),
tan(theta)=min(2xy, |x^2-y^2|)/max(2xy, |x^2-y^2|).       (E3)
```

**PROVED sharp bounds.**

```text
b=+1: tan(theta)<=65/72,
       equality only at x=13,y=5, giving (65,72,97);
b=-1: tan(theta)<48/55,
       with supremum approached along x=3 mod16, x->infinity.  (E4)
```

The corresponding angles are approximately 42.0750 and 41.1121 degrees,
strictly below the full shape endpoint 45 degrees.

Proof by the actual halving exponent:

| Parameter | k | Ratio restriction | Bound on tan(theta) |
|---|---:|---|---|
| +1 | 1 | x>=3; 3/2<y/x<=5/3 | <=8/15 |
| +1 | 2 | x>=9; 3/4<y/x<=7/9 | <7/24 |
| +1 | 3 | x=13 mod16; 3/8<y/x<=5/13 | <=65/72 |
| +1 | >=4 | 0<y/x<=1/4 | <=8/15 |
| -1 | 1 | x>=5; 7/5<=y/x<3/2 | <5/12 |
| -1 | 2 | x>=7; 5/7<=y/x<3/4 | <=12/35 |
| -1 | 3 | x=3 mod16; 1/3<=y/x<3/8 | <48/55 |
| -1 | >=4 | 0<y/x<3/16 | <96/247 |

Each bound follows by the monotonicity of `2t/(1-t^2)` below
`sqrt(2)-1`, and `(1-t^2)/(2t)` above it. The growing k=1 case can
equivalently use `(r^2-1)/(2r)` with r=y/x. At k=3 the residue guard
is exact, so x=13 really attains the positive maximum and x=3+16j
approaches the negative supremum. Every other row has a strictly smaller
upper bound. These are all-exponent proofs, not census inferences.

This gives a concrete limitation to sweeping the full epsilon-to-45-degree
triangle family as a model for consecutive Collatz states: an entire
terminal shape interval never occurs in this encoding.

## 3. The full tree's descent and the shell permutation do not preserve legality

Under U_1, the edge 7->11 encodes (77,36,85), with odd roots (11,7).
Its unique Berggren parent has roots (7,3), hence triangle (21,20,29).
Neither 3->7 nor 7->3 is a U_1 edge: their actual images are 5 and 11.
Therefore universal descent in the ambient primitive-triangle tree does
not imply descent along this selected set of Collatz edges.

The genuine fixed-root permutation in the
[triangle lane](odd_square_20260921_triangles.md), `v -> |u-2v|`, also
fails the legality test. The U_(-1) edge 5->7 has roots (7,5); the shell
permutation sends it to (7,3), which is not legal in either orientation
for U_(-1). The failure persists despite exact preservation of primitivity.
This locates the missing implication, rather than denying the useful
triangle correspondence.

## 4. Exact edge-height asymptotics and an atomic shape law

Count directed nonfixed positive U_b edges whose encoded hypotenuse is
at most X; call this E_b(X), b=+/-1. Then

```text
E_b(X)=C_edge sqrt(X)+O(log(2X)),
C_edge=(1/sqrt(2)) sum_(k>=1) 1/sqrt(4^k+9)
      =0.5078191557... .                               (E5)
```

Proof. For each k>=1 the exact valuation condition is one odd residue
class modulo `2^(k+1)`:

```text
x = (2^k-b)*3^(-1) mod 2^(k+1).
```

The height condition is
`(x^2+((3x+b)/2^k)^2)/2<=X`. Its positive endpoint is

```text
R_(k,b)(X)=(-3b+2^k sqrt(2X(4^k+9)-b^2))/(4^k+9)
         =sqrt(2X/(1+9/4^k))+O(1),                     (E6)
```

uniformly in k for fixed b=+/-1 and X>=1. Counting a residue class in
this interval yields `sqrt(X)/sqrt(2(4^k+9))+O(1)`.
Every admitted x is at most sqrt(2X), and `2^k<=3x+1`, so only O(log X)
branches can occur. The omitted tail of the series contributes O(1),
since each term is at most `2^(-k)/sqrt(2)`. Removing the unique fixed
edge changes the count by one. This proves (E5).

The ambient count is `P(X)=X/(2pi)+O(sqrt(X)log(2X))`, proved with its
parity and height conventions in the [shape lane](odd_square_20260921_shapes.md).
Consequently the unmarked triangles produced by legal edges have density
zero among all primitive triangles. No classification of 2-cycles is
needed: the number of unmarked images is at most E_b(X).

There is a stronger sampling distinction. In uniform directed-edge
sampling under this height, the limiting shape measure is the countable
atomic measure

```text
sum_(k>=1) w_k delta_(theta_k),
w_k = 1/(C_edge sqrt(2(4^k+9))),
t_k=min(3/2^k,2^k/3),
theta_k=min(2 arctan(t_k),pi/2-2 arctan(t_k)).            (E7)
```

For fixed k, x tends to infinity in density inside its residue class,
so y/x tends to 3/2^k and its branch shape approaches theta_k.
The normalized branch count tends to w_k by (E5)-(E6).
The uniform tail estimate is O(2^(-K))+O(log(X)/sqrt(X)) for k>K,
which justifies passage from finitely many branches to the full measure.
This proves weak convergence against continuous functions on [0,pi/4].

The ambient hypotenuse-bounded primitive triangles instead have a uniform
limiting theta law. Neither spatial law is automatically the time law
of one Collatz orbit. A primitive-point density such as 6/pi^2 cannot
be assigned to these edges: primitivity here holds identically by gcd.

## 5. Reproduction and stopping boundary

Run `python 04-computation/experiments/odd_square_20260921_edges.py` and
the same command with `-O`. The certificate exhausts positive odd starts
through 100001 for b=+/-1; independently counts bounded edges by starts
and by valuation residue classes; checks the sharp witness and both
legality failures; and compares primitive-triangle populations at explicit
heights. Decimal constants are illustrations, not proof gates.

The exact new questions are about operations that preserve the legal
edge condition or combine a controlled number of legal steps with a
decreasing integer quantity. The two natural triangle operations tested
here do not do that. Equations (E1)-(E7) classify this bridge's strengths
and limits; none eliminates a nontrivial Collatz cycle or an infinite orbit.
