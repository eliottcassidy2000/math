# The characteristic sequence (9,6,2) is polynomially realizable

## Main result

Over every characteristic-zero field k, set

- f(t)=t^6+4t^2;
- q(t)=t^9+6t^5+6t.

Then the exact identity

q(t)^2=f(t)^3+8f(t)+4t^2

holds. Consequently, if z=t^2 and p(z)=z^4+6z^2+6, then

k[f,q]=k[z,t p(z)].

Indeed, the identity places z=(q^2-f^3-8f)/4 in k[f,q]; conversely f=z^3+4z and q=t p(z). Thus the terminal degree pair (6,9), and the characteristic degree sequence (9,6,2), are genuinely realized by polynomial generators. No Keller hypothesis is used or asserted.

## Semigroup, codimension, and conductor

Write A=k[z,t p(z)] and B=k[t]. The even/odd decomposition gives

A=k[z] direct-sum t p(z) k[z],
B=k[z] direct-sum t k[z].

The two summands of A have disjoint degree parity, so its degree semigroup is exactly

Gamma(A)={0,2,4,...} union {9,11,13,...}=<2,9>.

Its four gaps are 1,3,5,7, and dim_k(B/A)=deg_z p=4. Also Frac(A)=k(t), since t=q/p(z).

The actual normalization conductor is p(t^2) B=(t^8+6t^4+6)B. To check this directly, write h=r(z)+t s(z). The conditions h in A and th in A respectively say p divides s and p divides r. Thus hB is contained in A precisely when h belongs to p(z)B. In particular, the source-conductor degree is 8.

The first nontrivial degree cancellation is

q^2-f^3-8f=4t^2,

which supplies the characteristic generator of degree 2 after degrees 9 and 6. Abstract characteristic admissibility here does not hide a realization obstruction.

## Immersion and nodal geometry

The polynomial p has no root at zero and is squarefree. In fact p(0)=6, and p'(z)=4z(z^2+3). If z^2=-3, then p(z)=-3, so gcd(p,p')=1.

The normalization parametrization t maps to (z,t p(z)) is immersed: its first derivative is 2t, which is nonzero away from zero, while the second derivative at zero is p(0)=6. Since f and q generate this same algebra, their parametrization is immersed too. Directly, differentiating 4t^2=q^2-f^3-8f shows that any common zero of f' and q' would have t=0, but q'(0)=6.

Over an algebraic closure, the model w^2=z p(z)^2 has exactly four ordinary nodes. Each root alpha of p is nonzero and simple and identifies the pair t=+sqrt(alpha), -sqrt(alpha); the tangent slopes in the (z,w) model are +sqrt(alpha)p'(alpha) and -sqrt(alpha)p'(alpha), hence distinct. The point over z=0 is smooth, and the partial derivative criterion yields no other affine singularities.

## A one-parameter family

For every v in k with v nonzero, define

f_v(t)=t^6+v t^2,
p_v(z)=z^4+(3v/2)z^2+3v^2/8,
q_v(t)=t p_v(t^2).

Then

q_v^2=f_v^3+(v^3/8)f_v+(v^4/64)t^2.

Thus every member realizes the same characteristic sequence, degree semigroup, delta=4, and conductor degree 8. Moreover p_v is squarefree and p_v(0) is nonzero. Its derivative is z(4z^2+3v); at a nonzero critical point z^2=-3v/4, the value of p_v is -3v^2/16, which is nonzero. The numerical example above is v=4.

## Optional elimination equation

Let U,V be plane coordinates and W=V^2-U^3-8U. A defining equation for the displayed (f,q) image is

W^3+64W-64U=0.

Setting z=W/4 converts this to U=z^3+4z; the identity then converts V^2=U^3+8U+4z into V^2=z(z^4+6z^2+6)^2. Conversely these substitutions recover W=4z, proving the two coordinate rings are isomorphic, rather than merely establishing a vanishing relation.

## Scope

This explicitly disposes of a proposed universal polynomial-realizability obstruction at (9,6,2). It does not construct a constant-Jacobian plane map, does not supply a relative C[x] realization with Keller compatibility, and does not resolve any Jacobian-conjecture frontier. It does show that a conductor-eight exclusion must use extra structure beyond abstract characteristic admissibility, polynomial plane generation, immersion, or nodality.
