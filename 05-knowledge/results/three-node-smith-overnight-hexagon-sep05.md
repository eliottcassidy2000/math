# Three-node two-jet Smith form from a minimal minor certificate

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** Canonical statement:
[THM-4429, arbitrary three-node Smith form](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md).
The finite probe and algebra below extend the dyadic special family in
[THM-4419, two-jet prime-wall precision](../../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md).

## 1. Exact formula

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

## 2. Mechanism and proof

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

Independent derivation of type II: differentiate the corresponding type-I
minor with respect to b. Only its value(b) row varies. This checks all eight
displayed formulas without importing a symbolic minor compiler.

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

The formula rests on the all-minor gcd proof, independently reviewed by
two other research lanes. FINITE-EXACT corroboration is deliberately broader
than its discovery family:

- Primary: all 2,016 triples (0,a,b), 1<=a<b<=64, their 8,064 profiles at
  p=2,3,5,7, all 69 residual minors, and 18,146 explicit gates.
- Independent: literal full six-by-six integer Smith forms on 2,999 triples:
  all -8<=x0<x1<x2<=16; 300 seeded translated/nonunit-dilated cases;
  and 399 depth/lift controls at p=2,3,5,7,11 through depth six. Thirty node
  permutations preserve the answer. Two mutations remove the exceptional
  factors two/three and are detected.
- Positive/hostile controls: (0,1,2) gives (1,1,1,1,4,4), (0,2,4) gives
  (1,1,4,4,32,128), (0,3,6) gives (1,1,3,27,324,324), and (0,8,16) gives
  (1,1,16,128,4096,131072). No metric-only conclusion for four nodes or
  higher jets is inferred from a finite census.

Reproduce from the repository root, also with `python3 -O`:

```bash
python3 04-computation/three_node_smith_probe_overnight_hexagon_sep05.py --height 64 --show-minors
python3 04-computation/three_node_smith_independent_overnight_hexagon_sep05.py
```

Normal and optimized runs agree byte-for-byte for each script. Semantic
SHA-256: primary `b573663ea9a0430170859e8577d4b8cb24365f9d484048e7fab12ffd82fcc3fb`;
independent `fefd67c54dea9188f05e53c2ea806064b321d057717d0a10ebed40237c38bb7f`.

Raw-LF SHA-256 manifest (paths under `04-computation/` and
`05-knowledge/results/`, respectively):

```text
three_node_smith_probe_overnight_hexagon_sep05.py
473f41a632f596a4cf37036474795e93b69203e353eea0478dc8f5370e932c2b
three_node_smith_independent_overnight_hexagon_sep05.py
6394cc11b70552d28f1f10b0a666b8e27e47eaa13653afd08e04c6560f5c3fab
three_node_smith_probe_overnight_hexagon_sep05.out
0e6d9484742e3d0029e56bc383f2be8f8df1a3e6c173bd6fae138855b33932c7
three_node_smith_independent_overnight_hexagon_sep05.out
a928080d2c96b34fcb4942c8dcd9bced6311fb46b2ad107c2e15c70a99728511
```

These are in-repository continuations, not an external priority claim.
