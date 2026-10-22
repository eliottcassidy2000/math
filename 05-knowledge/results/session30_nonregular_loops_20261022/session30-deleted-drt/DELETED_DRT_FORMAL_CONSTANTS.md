# Exact first and second formal constants for a deleted doubly regular tournament

## Scope and theorem

Let S be obtained by deleting a vertex from a doubly regular tournament of order q=n+1. Set M=S/n and d=S1. Assume the independently proved completion identities

S²=J−(n+1)I+ddᵀ, S1=d, Sd=−1, dᵀ1=0, d_i∈{−1,1}.

A1 and A2 are the complete five-core and 26-core finite resummed log-forest expressions, respectively. They are NOT the exact Hamiltonian count. The predictors κ1 and κ2 use the Euler/falling-factorial corrections defined below.

Along any sequence of such tournaments with n→∞, uniformly for z in any fixed bounded compact set,

A1(z) = −(z²+7z⁴/2)/n + (12z⁴−8z⁶/3)/n² + O(n⁻³),
A2(z) = (7z⁴+26z⁶)/n² + O(n⁻⁵ᐟ²).

In particular, with y_j=n(1−κ_j/ρ),

y1 = 3/2 − 13/(3n) + O(n⁻²),
y2 = 3/2 − 139/(24n) + O(n⁻³ᐟ²),
κ1/e = 1 − 5/(2n) + 20/(3n²) + O(n⁻³),
κ2/e = 1 − 5/(2n) + 65/(8n²) + O(n⁻⁵ᐟ²).

No assertion about the corresponding uniform asymptotic expansion of the EXACT Hamiltonian count is made. That requires a separate uniform bound on the omitted forest gas.

## 1. Exact scalar evaluation of A1

Write u=1/n, a=z²(n+1)/n², c=z²/n², and

g0=a/(1−a), g1=c/(1−c), h0=−a/(1+a), h1=−c/(1+c), o0=(z/n)/(1+a).

Let E=(J+ddᵀ)/n. This is the orthogonal projection onto span{1,d}. The propagators are

G=g0 I+(g1−g0)E,
H=h0 I+(h1−h0)E,
O=o0 S+(o1−o0)SE, where o1=(z/n)/(1+c).

Set gd=g0+2(g1−g0)/n, hd=h0+2(h1−h0)/n, gs=2(g1−g0)/n, hs=2(h1−h0)/n. The G and H off-diagonal entries vanish when d_i and d_j have opposite signs. When they have the same sign, E_ij=2/n, (SE)_ij=0 and O_ij=o0 S_ij. There are exactly n(n−2)/2 ordered unequal same-sign pairs. Therefore the full vacuum two-loop expression is EXACTLY

X1 = n(2gd hd−3gd²/2)
   + n gd(hd²−2gd²/3)
   + n(n−2)gs(hs²−o0²−2gs²/3)/2
   − n(hd−gd)²g1.

The path tree and lollipop are EXACTLY

T=−n g1,  L=−2n(hd−gd)g1.

Expanding these rational functions yields

X1=−7z⁴/(2n)+(8z⁴−8z⁶/3)/n²+O(n⁻³),
T=−z²/n+O(n⁻³), L=4z⁴/n²+O(n⁻³).

This proves the asserted A1 formula. The executable symbolic derivation is derive_a1.py.

## 2. Leading four-basis propagators and terminal elimination

Use the exact four-basis kernel convention R+, R−, Q++, Q+− from the certified engine, with Q=zS/n:

R+=Q(I−Q)⁻¹, R−=−Q(I+Q)⁻¹,
Q_st=−R_s R_t (I−Q²)(I+Q²)⁻¹.

Entrywise, uniformly in indices,

R± = n⁻¹(±zS−z²I)+O(n⁻²),
Q++ = n⁻¹z²I+O(n⁻²),
Q+− = −n⁻¹z²I+O(n⁻²).

These are entrywise, not operator-norm, assertions. They follow directly by reducing powers with S²=−(n+1)I+J+ddᵀ and S(J+ddᵀ)=d1ᵀ−1dᵀ. The latter rank-two matrices have bounded entries; all four scalar resolvents have denominators bounded away from zero for fixed bounded z and sufficiently large n.

A terminal index must be summed BEFORE substituting entrywise errors. On span{1,d}, S²=−I, so exactly

R±1 = ±(z/n)/(1+c) d − c/(1+c)1.

Consequently a leaf row sum contributes ±zd/n+O(n⁻²), a column sum has the opposite leading sign, and a Q_st leaf sum is O(n⁻²). These errors are uniform entrywise.

A weight-two core has E=V+2−ℓ edges, where ℓ is its terminal count. Eliminate all terminal indices. There remain v=V−ℓ indices and E factors of order n⁻¹, hence the trivial total scale is n^(v−E)=n⁻². The terminals of every weight-two core attach to nonterminal vertices, so these eliminations are independent. Replacing any factor by its O(n⁻²) error gives O(n⁻³). Choosing an internal off-diagonal-edge identity contribution forces two free indices to coincide and likewise gives O(n⁻³). A loop's identity contribution does not identify distinct indices and must be retained.

Thus each leading contraction reduces to a constant times n⁻E times a sum over v free indices of a product of S entries on the nonloop internal edges and one d factor at each terminal's attachment vertex.

## 3. Explicit bound for discarded colored diagrams

Every repeated internal edge with even multiplicity can be replaced by 1, since S_ij²=1−δ_ij. The equality cases discarded in this replacement contribute O(n^(v−1)), hence O(n⁻³) after normalization. Odd multiplicities reduce to one S factor. The remaining graph is simple and vertex weights are powers of d, hence bounded by one.

If this remaining signed graph contains an edge (a,b), fix all other vertex indices. Every factor involving a but not b becomes a bounded function f(a), and every factor involving b but not a becomes a bounded function g(b). All other factors become a bounded scalar. There is no second (a,b) edge because the graph is simple. Therefore

|Σ_ab S_ab f(a)g(b)| ≤ ||S||op ||f||2 ||g||2 ≤ √(n+1) n.

Summing the other v−2 indices gives at most n^(v−1)√(n+1)=O(n^(v−1/2)). After n⁻E normalization this is O(n⁻⁵ᐟ²). This argument permits arbitrary bounded vertex colors, including d, and needs no independence or unproved tournament counting claim.

If the remaining graph is empty, the sum factorizes into products of Σ_i d_i^m. Each odd m gives zero because Σd_i=0, while each even m gives n. Thus a leading diagram survives precisely when every internal nonloop multiplicity is even AND every residual d exponent is even.

This supplies the explicit quantitative bound on every discarded diagram, rather than assuming a blanket quasirandomness statement.

## 4. Exact finite integer-kernel contraction

For each surviving leading term, substitute into its certified four-basis integer tensor:

- Loop: coefficients (−z²,−z²,+z²,−z²).
- Internal nonloop edge: coefficients (+z,−z,0,0), with the S entry retained until its even multiplicity is checked.
- Leaf at the second endpoint: coefficients (+z,−z,0,0), multiplying d at the first endpoint.
- Leaf at the first endpoint: negate those leaf coefficients.

The exact rational calculation in derive_a2_leading.py reads symbolic-kernels.json, retaining each core's certified automorphism denominator. It produces a2-leading-certificate.json. Exactly FOUR of 26 cores survive:

1. One vertex with three loops: 12z⁶.
2. Two vertices, one loop at each and two parallel connecting edges: 14z⁶.
3. Two vertices with four parallel edges: 3z⁴.
4. One loop vertex with two attached terminal leaves: 4z⁴.

All other cores vanish at leading order by the criterion above. In particular, ALL seven one-leaf cores and the Y tree vanish. Hence

lim n² X2_vacuum(z)=3z⁴+26z⁶,
lim n² A2(z)=7z⁴+26z⁶.

The vacuum constant 29 at z=1 is DERIVED here from the exact finite kernels; it is not imported from numerical extrapolation. The deleted-tournament total is 33. The additional four is the two-leaf loop core, not a fitted correction.

## 5. Falling-factorial transform and the constant 65/8

Let D=z∂z, c_j=D^jX0, and b_j be the first four complete exponential Bell polynomials in these cumulants. Define

p1=(b2−b1)/2,
p2=(3b4−2b3−3b2+2b1)/24,
shift=(D²A1+(2c1−1)DA1)/2.

The formal predictions at z=1 are

Y1=A1+p1/n,
Y2=A2+A1²/2+(A1 p1+shift)/n+p2/n²,
κ1=ρ(1+Y1), κ2=ρ(1+Y1+Y2).

The exact spectrum gives

ρ=((n²+1)/(n²−1)) ((n²+n+1)/(n²−n−1))^((n−2)/2),
X0(z)=z²(1−1/n)+z⁶/(3n²)+O(n⁻³).

Thus ρ/e=1−1/n+5/(6n²)+O(n⁻³), p1(1)=3−5/n+O(n⁻²), p2(1)=155/6+O(n⁻¹), and shift(1)=−54/n+O(n⁻²). Since A1(1)=−9/(2n)+28/(3n²)+O(n⁻³),

Y1=−3/(2n)+13/(3n²)+O(n⁻³),
Y2=[33+81/8−27/2−54+155/6]/n²+O(n⁻⁵ᐟ²)
  =35/(24n²)+O(n⁻⁵ᐟ²).

Combining the factors proves the theorem. derive_a1.py verifies the rational arithmetic symbolically. Its formal SymPy O(n⁻³) output refers to the substituted coefficient polynomial; for the actual 26-core expression we conservatively retain the proved O(n⁻⁵ᐟ²) remainder.

## Validation and limitations

No finite Paley sample is an input to either coefficient proof or certificate. The first proof uses exact scalar resolvent algebra; the second uses an already independently audited finite kernel list and an explicit signed-density bound. Numeric samples were used only as optional orientation during exploration and are not evidence in the proof.

A rooted five-vertex census can affect finite-n values or higher-order terms; it cannot change these leading constants under the bound above. Neither this theorem nor its certificates prove a uniform asymptotic formula for exact Hamiltonian path counts, Collatz coverage, or any unrelated open conjecture.
