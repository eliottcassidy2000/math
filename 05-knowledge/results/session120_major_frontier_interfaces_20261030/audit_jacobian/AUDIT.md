# Independent audit: affine-in-x Hamiltonian image ideal

**Verdict: the main classification is correct.** This is an elementary restricted
sector theorem, not a solution or reduction of the general Jacobian conjecture.
The proof works over **every characteristic-zero field K**, not only an
algebraically closed one. No novelty claim is made.

## 1. Independent coefficient proof: no concealed cancellation

Let P=A(y)x+B(y), with A nonzero. The substitution T=P identifies K(y)[x]
with K(y)[T], and the Hamiltonian derivation D={P,-} becomes A partial_y
with T fixed. Write Q=sum q_j(y)T^j and h(T)=sum h_j T^j. Then

    DQ=h(P)  iff  q_j'=h_j/A for every j.

If h is nonzero, choose j with h_j nonzero and put R=q_j/h_j. This proves
R'=1/A directly; it does not assume integration inside the polynomial ring.
For every i, (q_i-h_i R)'=0, so q_i-h_i R belongs to K. (The constants of
K(y) under d/dy are K in characteristic zero.) Consequently

    Q=h(P)R(y)+F(P),  F in K[T].

In particular F cannot be a rational function with a denominator: the
coefficient equations already prove it is polynomial. Because F(P) is itself
polynomial, it cannot conceal a nonpolynomial h(P)R(y). If h=0, the same
argument gives ker D=K[P], independent of whether a primitive exists.

## 2. Independent rational-primitive classification

Extend scalars to an algebraic closure just for the pole argument. Suppose
R'=1/A. A polynomial A of degree one would give a simple pole in R', which
no rational derivative has. If deg A>=2, subtract the finite limit of R
at infinity: its polynomial part was constant because its derivative tends
to zero there. Now write R=N/D with coprime polynomials, n=deg N<d=deg D.
Let r be the number of distinct roots of D. At a root a of multiplicity k,
N'D-ND' has order exactly k-1: the leading coefficient is
-k N(a) times the nonzero local unit of D. Characteristic zero is essential.
Moreover

    deg(N'D-ND')=n+d-1,
    deg gcd(N'D-ND',D²)=d-r.

The first formula also holds for n=0. The reduced numerator of R' therefore
has degree n+r-1. Since 1/A has constant numerator, n+r=1. Necessarily
n=0 and r=1, so R is a nonzero scalar multiple of (y-b)^(-d), and
A=c(y-b)^(d+1). Conversely these A have the advertised primitive. A nonzero
constant A=c has primitive y/c.

Descent to K causes no additional case: a polynomial over K with its unique
geometric root b has b in K, since b=-a_(m-1)/(m a_m). Thus the complete list
over K is: nonzero constants, and c(y-b)^m with m>=2 and b,c in K, c nonzero.

## 3. Exact exponent and all companions

For A=c(y-b)^m write z=y-b, beta=B(b), and

    R=-1/[c(m-1)z^(m-1)].

Let t=ord_b(B-beta), with t=infinity for constant B, and v=min(m,t).
The z-valuation of P-beta is **exactly v** in K[x][z]. If t=m the leading
coefficient is cx plus a scalar, never the zero polynomial. If h has root
multiplicity k at beta, writing h(T)=(T-beta)^k h0(T) with h0(beta) nonzero
gives

    ord_z h(P)=k v.

Thus h(P)R is polynomial precisely if kv>=m-1. The ideal is

    I_P=(T-beta)^e K[T],  e=ceil((m-1)/v).

For finite t this equals ceil((m-1)/t), including t>m; for constant B it
is 1. All companions are h(P)R+F(P), F arbitrary in K[T]. In the constant-A
case all h are allowed, with companions (y/c)h(P)+F(P). In all other cases
only h=0 is allowed. This also proves the claimed ideal property rather
than presupposing it.

## 4. Boundary checks and scope

* Linear A is excluded by its simple pole. Multiple distinct roots are
  excluded even when all root multiplicities are at least two.
* Constant B is not an exception to the valuation calculation: v=m and e=1.
* A nonconstant irreducible A over K is handled by extending scalars; no
  hidden real/rational-field exception exists.
* P=y²x+y and Q=-xy-1 give {P,Q}=P, but cannot give a nonzero constant.
* A nonzero constant bracket implies A constant in this sector. The inverse
  of the resulting triangular map is explicit. This gives no coordinate
  reduction for arbitrary planar maps.
* A=0 is outside the above theorem. The producer's separate extension to
  P=B(y) nonconstant is also valid: D=-B' partial_x, and h(P) lies in its
  image iff B' divides h(B). At a critical point a with ord_a B'=r,
  ord_a(B-B(a))=r+1, so this means precisely that h vanishes at every
  distinct critical value. All companions are -x h(B)/B'+H(y); kernel K[y].
  Constant P must remain excluded from identifying K[T] with K[P].

## 5. Independent exact verification

`check_independent.py` uses only Python standard-library Fraction arithmetic,
with no SymPy and no imports from producer files. All failures use explicit
RuntimeError checks and remain active under optimization.

Finite universe:

* m=2,...,7; t=1,...,8 or infinity; three tuples (b,c,beta), namely
  (-2,3,5), (0,-2,-1), (3,1,0); all k=0,...,e+2. Exact division by y-b tests
  both sides of the sharp threshold; admitted companions receive a kernel
  term P² and are checked by direct differentiation and multiplication.
* Constant A=3, three B choices 0,y,y³+1, and h=P^k for k=0,...,4.
* Fifteen full matrix intersection cases with polynomial companion domain
  0<=deg_x<=3, 0<=deg_y<=10 and target span(1,P,P²,P³). The derivation matrix
  is constructed directly from monomial derivatives, not from integration.
  Rational rank gives the complete intersection dimension in this finite
  universe, equal to the predicted values. The domain is sufficiently large
  to contain a companion for every admitted target basis element in these
  selected cases; this is not an unrestricted nonexistence proof. Four cases
  have A=0, including B=(y²-1)² with repeated critical value 0, and
  B=y³+y whose critical points are not rational.
* For every d1,d2=1,...,4, D=y^d1(y-1)^d2 and n=0,...,d1+d2-1,
  N=1 when n=0 and N=y^n+1 otherwise: exact Euclidean gcd reduction of
  (N'D-ND')/D² verifies reduced numerator degree n+1, obstructing a
  reciprocal-polynomial derivative for these two-pole samples.

Normal Python and Python -O outputs are required to match byte-for-byte.
The accompanying receipt records the counts and exact matrix dimensions.
Finite tests support implementation and examples; Sections 1–3 establish
all-parameter claims.

## 6. Provenance correction and parameter-family extension

The producer subsequently recovered repository THM-3326/3348, which already
contain the A nonzero classification. Accordingly this audit certifies that
classification and its implementation, **not its novelty within the repo**.
The following parameter-family claims were then independently checked.

Over K[t], set P_t=t y^5 x+y^3 and h(P_t)=P_t². For every fixed nonzero
parameter, the unique companion modulo K[P_t] is

    Q0=-t y^6 x²/4-y^4 x/2-y²/(4t).

At t=0 the target is y^6 and the companion -xy^4/3 works. Nevertheless
there is no polynomial Q in K[t,x,y] realizing this identity globally.
Over K(t), all companions are Q0+H_t(P_t). Expand the finitely many
coefficients of H_t(T) at t=0. If H has pole order r>1, the leading pole
coefficient after substituting P_t is the nonzero polynomial H_{-r}(y³);
substitution T->y³ is injective, so it cannot vanish. This cannot cancel
Q0's pole of order one. Hence r<=1. Cancellation of the coefficient of
1/t would require H_{-1}(y³)=y²/4, impossible since every exponent on the
left is divisible by three. This argument excludes every polynomial lift,
not just the displayed choice of companion. It also excludes lifts whose
coefficients are rational functions of t regular at zero.

For P_t=t y^5 x+(y²-1)², the generic image ideal is ((T-1)²), since
ord_0((y²-1)²-1)=2 and ceil(4/2)=2. The special ideal at t=0 is
(T(T-1)), since the distinct critical values are 0 and 1. Neither ideal
contains the other. This is an exact example, not a universal assertion
about arbitrary specialization. In particular fibrewise existence need
not imply existence of a companion varying regularly with the parameter.
