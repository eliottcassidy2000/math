# Three-node two-jet Smith form from a minimal minor certificate

**PROVISIONAL PROOF CANDIDATE / UNDER INDEPENDENT AUDIT.** No theorem ID is
claimed yet. This is not a proved dependency until promotion. The finite
probe and algebra below extend the dyadic special family in
[THM-4419, two-jet prime-wall precision](../../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md).

## 1. Candidate exact formula

Let x0<x1<x2 be arbitrary integers. Set a=x1-x0, b=x2-x0,
g=gcd(a,b), P=ab(b-a)>0, u=a/g, v=b/g, and

```text
epsilon=gcd(3,g,u+v).
```

For the map from integer polynomials of degree below six to their values
and first Hasse derivatives at the three nodes, the Smith factors are

```text
(1,1,D1,D2/D1,D3/D2,D4/D3),
D1=g*gcd(g,2),
D2=gcd(g^4,6P),
D3=2*g^4*P*epsilon,
D4=P^4.                                              (1)
```

These are global integer factors, not just a dyadic list. Translation and
node permutation preserve the observer lattice. The quantities g and |P|
are visibly invariant; u+v modulo three is unchanged when the origin node
is changed, up to sign, so epsilon is invariant as well.

## 2. Mechanism and proof candidate

Translation is an integral unitriangular column change. Clearing the value
and derivative at zero leaves the residual matrix

```text
A = [a^2, a^3, a^4, a^5;
     2a,  3a^2,4a^3,5a^4;
     b^2, b^3, b^4, b^5;
     2b,  3b^2,4b^3,5b^4].                           (2)
```

Write D_h for the gcd of its h-minors, D_0=1. D1 is the gcd of
a^2,b^2,2a,2b: every remaining entry is already divisible by it.

For D2, every minor except the first two columns of the two derivative
rows has total degree at least four, and hence is divisible by g^4.
The exceptional minor is +/-6P. The same-node minors on those two columns
are a^4,b^4, whose gcd is g^4. This proves the D2 formula without a
generic-position assumption.

For D3, all sixteen minors have the common factor P*g^4. Up to exchanging
a and b (and changing signs), the two row types are

```text
I:  value(a), derivative(a), value(b);
II: value(a), derivative(a), derivative(b).
```

In increasing column triples 234,235,245,345, their factors before removing
P*g^4 are respectively

```text
I:  a^4*b^2*(a-b)^2 times
    (1, 2a+b, a(a+2b), a^2*b);
II: P times
    (2a^3(a-2b),
     a^3(4a^2-5ab-5b^2),
     2a^4(a^2+ab-5b^2),
     a^5*b(3a-5b)).                                 (3)
```

Signs do not affect the gcd. The first type after division contains
g*u^3*v*(u-v), which is even since uv(u-v) is even. In the second type,
the first and third expressions are visibly even. For the other two,
reduction modulo two gives a multiple of uv(u+v), also always even.
Thus every normalized minor is divisible by two.

The two degree-seven minors of type II, after normalization, are

```text
2u^3(u-2v), 2v^3(2u-v).                              (4)
```

Their gcd is 2*gcd(3,u+v). Indeed gcd(u,v)=1 excludes primes dividing
u or v from the common factor; any remaining common prime divides both
u-2v and 2u-v, hence divides 3. When three divides these linear forms,
their common three-adic valuation is exactly one, since their linear
combinations include 3u and 3v and neither u nor v is divisible by three.

If three divides u+v, all u,v,u-v are units modulo three. The first
type-I normalized minor g*u^3*v*(u-v) therefore removes the extra factor
three precisely when 3 does not divide g. If 3 divides g, every other
minor has at least one additional factor g (its total degree is at least
eight), so the factor three remains. This proves D3 in (1).

Finally the confluent Vandermonde determinant is P^4, giving D4. Successive
determinantal-divisor quotients give the asserted Smith form. Every witness
comes from an actual minor; no determinant-to-Smith shortcut is used.

## 3. Metric form and sharp precision

Fix a prime p. The sorted three pairwise difference valuations have the
form (e,e,f), f>=e. For p=2, necessarily f>e. Put d_h=v_p(D_h). Formula (1)
is equivalent to

```text
d1=min(2e,e+v_p(2)),
d2=min(4e,2e+f+v_p(6)),
d3=6e+f+v_p(2)+[p=3 and e=f>0],
d4=8e+4f.                                           (5)
```

Thus the complete p-Smith list depends only on the pairwise distance tree,
not the higher unit digits. The largest exponent, equivalently the sharp
uniform p-adic precision loss for recovering all six source coefficients,
is

```text
L=2e+3f-v_p(2)-[p=3 and e=f>0].                      (6)
```

This uses the same elementary Smith-coordinate argument as THM-4419; the
fixed degree box and the full two-jet observer are essential. Higher jets,
four or more nodes, and arbitrary moving source modules remain OPEN.

## 4. Inheritance, hostile controls, and status

Closest mechanism: THM-4419's all-minor proof for one dyadic triple.
Corrected near miss: THM-4010 separates the determinant from the Smith
list. Canonical hostile: adding 16 to nodes (0,8) lowers an earlier
exponent. Least-used sidecar: the ideal of all intermediate minors.

The paper-inspired move is to replace the entire search by a small
target-bearing witness set, while keeping the literal full enumeration as
an independent check. There is no map from point-set hexagons to this
integer matrix. Equations (3)--(4) are native witness minors, not geometric
triangles. The parity and exceptional-three factors are the information a
generic rational-rank argument would discard.

Initial FINITE-EXACT probe: every (0,a,b) with 1<=a<b<=64, at primes
2,3,5,7, has no same-distance/different-Smith pair. This is evidence only;
the formula rests on the gcd argument. Stronger independent exact audits,
source/output manifests, and promotion are still owed.
