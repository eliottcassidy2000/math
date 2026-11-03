# Exact arithmetic carriers for the error-loop language

Status: **PROVED** for the six explicit macros below. The observed integer is
the absolute error multiplier, NOT the tournament's H or its automorphism order.

## A1. Complete valuation classification in the selected macro language

Use L2,L6,L45,L13 from THEORY and the neutral-magnitude macros B=b and Q=q2.
Their gains are-2,6,45,13,-1,1 and vertex costs8,9,19,17,3,4. These six
words form a prefix code, so parsing a macro sequence after its supplied
seed is unambiguous. Repeated matrices do not identify distinct literal words.

If their four nonneutral counts are a,b,c,d, the magnitude is

    g=2^a 6^b 45^c 13^d=2^u 3^v 5^w 13^x,
    (u,v,w,x)=(a+b,b+2c,c,d).

Hence the exact generated valuation cone is

    w,x>=0,  v>=2w,  u>=v-2w.                         (A1)

Conversely a=u-v+2w,b=v-2w,c=w,d=x constructs every such magnitude.
Thus3 and5 occur only with the stated other valuations; they are not free
generators individually. This classifies this chosen loop language only,
not all possible return gains of the65-symbol semigroup.

The error sign is `epsilon=(-1)^(a+#B)`, and the vertex cost is
`8a+9b+19c+17d+3#B+4#Q`. Both remain explicit clocks. Macro order still
changes H, as the noncommuting control in THEORY shows.

## A2. One exact multigraded generating function

Let z mark added vertices, x_p mark v_p(g), and s mark error sign, with s^2=1.
The presentation series is

    1 / [1-s z^3-z^4-s x2 z^8-x2 x3 z^9
                 -x3^2 x5 z^19-x13 z^17].             (A2)

This follows by uniquely choosing the next macro. It counts source words;
it is neither a generating function for distinct H values nor A000568.
For weighted Hamilton sums replace each scalar macro weight by its actual
3x3 matrix and use `v_seed (I-sum z^cost P_macro)^(-1) e_h`. This retains
order in the matrix products; the scalar error-character alone does not.

## A3. Liouville transports repeated factors; Mobius detects them

Write lambda(g)=(-1)^Omega(g), and mu(g)=0 when g has a square factor,
otherwise(-1)^omega(g). These definitions and the distinction between their
multiplication responses are inherited classical arithmetic, recovered in
[liouville_mobius_11_13/ARITHMETIC A1,A4](../liouville_mobius_11_13_20261026/ARITHMETIC.md).
Here we apply them to a new, explicitly realized error-multiplier language.

From(A1), the squarefree magnitudes are EXACTLY

    {1,2,6,13,26,78}.                                  (A3)

Indeed c must vanish, d<=1, and (a,b) is(0,0),(1,0),or(0,1).
Let N=1-s z^3-z^4, w=s z^8+z^9, t=z^17. The series counting squarefree
words, retaining error sign and vertex clock, is

    1/N + (w+t)/N^2 + 2wt/N^3.                        (A4)

There can be zero, one, or two nonneutral macros. With two, one is L13 and
the other L2 or L6, in either order; neutral words fill the intervening gaps.

For mu-weighted words replace w by `-s z^8+z^9` and t by `-z^17` in(A4).
For lambda-weighted words complete multiplicativity instead gives

    1 / [1-s z^3-z^4+s z^8-z^9+z^19+z^17].            (A5)

The two operations are different: a repeated L2 contributes lambda(4)=1
but mu(4)=0. B reverses error sign while leaving both arithmetic values
unchanged. Error orientation, Liouville parity and squarefreeness must not
be merged into one bit. Exact dynamic coefficients through added cost48
agree with(A4-A5), independently for s=1 and s=-1.

## A4. Smallest arithmetic response carriers

Two reachable magnitudes are equivalent if every future macro continuation
has the same requested arithmetic observation. Then the exact minima are:

| Observation | Minimum states |
|---|---:|
| lambda(g) |2|
| mu(g) |7|
| (lambda(g),mu(g)) |8|
| (lambda(g),mu(g),error sign) |16|

For mu, retain the six squarefree values(A3) and one absorbing nonsquarefree
state. The three available contexts with magnitudes1,2,13 distinguish all
seven states; their mu-response triples are

    1:(1,-1,-1), 2:(-1,0,1), 6:(1,0,-1),
    13:(-1,1,0), 26:(1,0,0), 78:(-1,0,0), dead:(0,0,0).

For the joint observer split dead into lambda=+1 and-1, represented by4 and45.
On squarefree states lambda=mu, so no other split is needed. All states are
reachable. B toggles error sign without altering g, proving the extra factor2
and its necessity when sign is reported. A clock modulo m can be retained
by a product with Z/m; no minimality claim is made for that product or for
the unbounded vertex count. The original source word remains a separate record.

## A5. The Hamilton observer still needs three linear coordinates

Restricting to these macros does not reduce the rational linear H observer.
For actual prefixes U, U B, U Q and contexts empty,B,Q, the response matrix is

    [[31,187,417],[187,1575,3325],[417,3425,7279]],

with determinant309024. Thus its Hankel rank is at least3; inherited3-state
closure supplies equality. Literal graph frontier counts check all nine entries.

There is also an exact ambient integral observation lattice. The four column
observations for empty,B,Q,L2 generate `Z e_h + 2Z e_x + 2Z e_y`; the Hermite
normal form is diag(1,2,2). Every macro preserves this lattice, by the parity
of its first column. Therefore for arbitrary integer row states the minimal
future-H quotient modulo m is

    Z/m x (Z/(m/gcd(m,2)))^2.                          (A6)

This is an integer-state response theorem, not an assertion that every quotient
class has an actual source of a prescribed size. It also does not make the
small arithmetic carriers in A4 into Hamilton-path or graph-reconstruction
carriers: those observe a different quantity.
