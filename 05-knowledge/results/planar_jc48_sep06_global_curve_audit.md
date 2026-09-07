# Independent audit: the concrete higher-cusp curve and Nori's theorem

**Status: PASS — full analytic/text/source audit and normal/optimized/frozen
replay. Nori's theorem is CITED. This excludes the declared curve as a
whole nonproperness set; it does not exclude all higher-cusp supports or
prove JC(2).** September 6, 2026.

Audited [primary proof](planar_jc48_sep06_global_curve.md),
[source](../../04-computation/planar_jc48_sep06_global_curve.py), and
[output](planar_jc48_sep06_global_curve.out). The exact declared curve is

```text
C = image(t -> (t^4+t^2, t^6+t^5+t^2)).
```

No primary source, output, or proof file was edited by this referee.

## 1. Primary theorem and application type

I opened the actual Numdam primary PDF of Madhav V. Nori,
[*Zariski's conjecture and related problems*](https://www.numdam.org/article/ASENS_1983_4_16_2_305_0.pdf),
Ann. Sci. École Norm. Sup. (4) 16 (1983), 305–344, and read
Proposition 3.27 on journal page **331** (PDF page 28, zero-based index
27), together with Introduction II on page 306. Its hypotheses are a
smooth projective surface, curves `D,E` meeting transversely, `D` nodal,
and `D_i^2>2r(D_i)` for each irreducible component of `D`. The kernel of
the map from `pi_1(X\(D union E))` to `pi_1(X\E)` is abelian. Its further
finite-index centralizer statement is not needed here.

This is the applicable theorem. A projective nodal-complement result
alone would not prove the affine conclusion, and equality in its
numerical hypothesis would not suffice. The proof correctly retains
all exceptional curves and the original line at infinity in `E`.

## 2. Complete normalization and singularity audit

The monic equation for `t` makes the parametrization finite onto its
image. The first divided difference factors as
`(s+t)(s^2+t^2+1)`. At `s+t=0`, the second difference is `s^4`, giving
only the diagonal critical pair at zero. On the other branch, the
pair-sum relation `q=(p^2+1)/2` reduces the second difference to the
displayed degree-five `H(p)`.

I checked the interpretation of every elimination gate. Its nonzero
discriminant gives five distinct nonzero pair sums; the resultant with
`p^2+2` keeps the two parameters of each pair distinct. The independent
parameter resultant has ten simple off-diagonal roots. The pair locus
is finite, so generic degree greater than one is impossible. Thus this
finite birational parametrization by `A1` is the normalization.

The common-fibre remainder has degree exactly three. If a common fibre
had three distinct parameters, its negative monic cubic would divide
`U-a`; comparison forces both

```text
A^2+A-1=0,   (A-2)(A^2+1)=0.
```

The primary's displayed Bezout combination equals five, proving
impossibility. A four-parameter fibre is also impossible from the cubic
remainder. This is a complete fibre argument, not an assumption that
different pair sums automatically yield different image points.

The derivative gcd is `t`; the critical image has no second parameter
because `gcd(U,V)=t^2`. The local target change gives
`x=t^2+t^4`, `y=t^5+3t^6+t^8`, hence an analytic `(2,5)` cusp.
The exact tangent ideal `(N,M,T)=(s+t,t^4)` puts every off-diagonal
collision outside the tangent-degeneracy locus. Therefore there are
exactly five ordinary nodes and no other finite singularities. The
pair sum `p=1` supplies the node `(-1,1)`, preserving the corrected fifth
node rather than the discarded four-node count.

The degree-six homogeneous parametrization has no base point and is
birational, so its image really is a sextic, with just one normalization
point at infinity. The infinity chart is the original projective
embedding; the finite triangular coordinate change is not silently
applied to that embedding.

## 3. Independent reconstruction of the two resolutions

At the finite cusp, direct successive coordinate divisions yield orders

```text
(2,5) -> (2,3) -> (2,1) -> (1,1) -> (1,1).
```

After the second blowup the strict curve is smooth but tangent to the
current exceptional curve. After the third it passes through the
intersection of two exceptional curves, with tangent ratio one. The
fourth blowup moves it to a point away from both old exceptional curves;
it is transverse to the new exceptional because its second coordinate
has order one. Hence all four centres are necessary for the stated
normal-crossing configuration, and their curve multiplicities are
`(2,2,1,1)`, with square sum ten.

At infinity, I independently used

```text
X=(w^2+w^4)/(1+w+w^4),   Z=w^6/(1+w+w^4),
rho=Z/X^3-1=2w+O(w^2).
```

The simultaneous curve/line resolution has successive orders

```text
(2,6) -> (2,4) -> (2,2) -> (2,1) -> (1,1) -> (1,1).
```

Through the first two blowups, the old line at infinity is one of the
coordinate axes through the strict curve. The third blowup follows the
diagonal tangent and places that line at `rho=-1`, while the curve
passes through `rho=0`. The fourth produces an exceptional crossing;
the fifth uses the tangent ratio `1/4` and separates the curve from
both old exceptional branches. Its final second coordinate has order
one and leading coefficient two. The correct multiplicities are
`(2,2,2,1,1)`, with square sum fourteen.

Thus this is a resolution of `Cbar union L_infinity` at infinity, not
just a resolution of the branch that forgets the line. The final strict
curve meets the reduced exceptional union in two smooth points, both
transversely and away from exceptional crossings. Its original five
nodes are untouched. The resulting inequality is exactly

```text
D^2 = 36-10-14 = 12 > 10 = 2r(D).
```

The separate genus check `2+3+5=10` agrees with the complete inventory,
but is not being used to infer transversality or count the node images.

## 4. Complement identities and the actual Keller conclusion

Let `E` be the reduced union of every exceptional curve and the strict
line at infinity. All centres lie over either the finite cusp origin
or that line. Consequently the blowdown identifies

```text
X\E = A2\{0},
X\(D union E) = A2\C.
```

The finite cusp belongs to `C`, so deleting its full exceptional fibre
has not deleted any additional point from the second complement.
Complex `A2\{0}` retracts onto `S^3` and is simply connected. Nori's
abelian kernel is therefore the entire affine complement group.

An abelian transitive permutation action is free after passing to its
faithful image: an element fixing one label fixes every label by
commutation and transitivity. Thus a single three-cycle fixing three
labels cannot occur in a connected degree-six cover of this complement.
This excludes the precise passport supplied by the odd-cusp theorem.

The primary also correctly gives an independent all-degree finish for
a Keller map with this whole support. Its proper-locus cover is
connected, each generic smooth-support meridian has an actual fixed
sheet, and the ambient affine plane is simply connected. The meridians
normally generate by transverse punctured-disk filling. Their images
in an abelian transitive group must all be trivial, so the entire
monodromy is trivial. This contradicts a nonautomorphic Keller map.
It does not require the odd-cusp degree spectrum as a dependency.

The local abstract cusp passport and the semilocal `A6` control in
[the alternating note](planar_jc48_sep06_alternating.md) remain valid.
The new exclusion is the missing actual-curve-complement obstruction,
not a refutation of those weaker objects.

## 5. Source replay and final acceptance

I read the complete producer and independently ran

```sh
python3 -B 04-computation/planar_jc48_sep06_global_curve.py
python3 -B -O 04-computation/planar_jc48_sep06_global_curve.py
```

Both modes pass **56 always-active gates**. A separate subprocess byte
comparison confirms normal stdout, optimized stdout, and the complete
frozen output are identical. The recomputed hashes are

```text
source SHA-256:
6a5ef57de9491090c013ced57355d888e5fb7959cd45c8f78fb7af7fabf0478d
output SHA-256:
c7016a7d0d45b28f777da16c784373712507128f5e7a1894361a262c0fe768cc
```

Final acceptance: **PASS**, without mathematical or source correction.
The exact producer verifies the one declared curve and chart identities;
the primary topological theorem and the complement transport supply the
analytic conclusion. No general intrinsic `(4,6)` classification,
higher-cusp exclusion, or full Jacobian-conjecture claim is approved.
This audit is frozen for integration.
