# Independent audit of the linear resolution criterion and equality sextic

**Status: independent analytic/source/replay PASS.** September 6, 2026.
Auditor: `orthogonal_returns`. No mathematical correction was needed.
The strict inequality gives an exclusion; equality supplies neither
abelianity nor a Keller realization.

Audited [primary proof](planar_jc48_sep06_resolution_budget.md),
[source](../../04-computation/planar_jc48_sep06_resolution_budget.py),
and [frozen output](planar_jc48_sep06_resolution_budget.out). This note
changes no primary artifact and reserves no theorem ID.

## 1. Actual theorem, intersection calculation, and complement transport

I independently opened the primary Numdam PDF of Madhav V. Nori,
[*Zariski's conjecture and related problems*](https://www.numdam.org/article/ASENS_1983_4_16_2_305_0.pdf),
and read Introduction II, journal p.306, and Proposition3.27, p.331
(PDF indices2 and27). The needed hypotheses are a smooth projective
surface, a nodal divisor `D`, transverse intersection with `E`, and
strict inequalities `D_i^2>2r(D_i)` for its irreducible components.
The conclusion used is abelianness of the kernel between the two
complement fundamental groups. The additional centralizer assertion
is not used. This is the applicable affine-complement route after the
specified blowups, not a projective-complement theorem used without
the line at infinity.

At a blowup centre where the strict curve has multiplicity `m`,
`K'.D'=K.D+m`, because the new exceptional has square minus one and
is orthogonal to pullbacks. This also holds with contribution zero
for a centre off the strict curve. Starting from `K_P2.Cbar=-3e`
therefore gives `K_X.D=-3e+sum m_j`. Since `D` is irreducible nodal,
`p_a(D)=g+N`. Adjunction proves exactly

```text
D^2-2N=3e-2+2g-sum m_j.
```

Genus has the correct positive sign. The cancellation of the node count
is justified by adjunction, not by discarding nodes. A further blowup
at a retained node decreases this margin by two and adds two to the
cost; a further blowup at a smooth curve point decreases it by one.
Thus the chosen-resolution dependence and the quantifier over every
admissible resolution are correct. No intrinsic minimum is computed.

With every finite exceptional included and every point of `T subset C`
blown up, removing `E` gives precisely `A2\T`. Removing `D` as well
gives precisely `A2\C`. The first space is simply connected: a disk
bounding a loop in real dimension four can avoid a finite collection
of points. Under the strict inequality the whole affine complement
group is therefore the abelian kernel supplied by Nori. The hypotheses
that `D` meets only smooth points of `E` transversely and has no
singularities besides the retained nodes have not been weakened.

The all-degree Keller consequence is sound. The inverse image of the
proper locus is connected, so its finite monodromy image is transitive.
Pulling back a defining nonconstant polynomial of the irreducible support
gives a nonempty curve in the source; quasi-finiteness prevents its
components from mapping to points. Hence a generic smooth-support
meridian fixes an actual sheet. In an abelian transitive image it fixes
every sheet. Normal generation by these meridians then makes the whole
image trivial. The resulting degree-one Keller map is an automorphism;
this is incompatible with a nonempty nonproperness support. Neither the
curve degree `e` nor a local cusp exponent has been substituted for the
mapping degree in this argument.

## 2. Complete affine curve inventory

The source uses the exact curve

```text
U=t^4+t^3+t^2,
V=2t^6+3t^5+2t^3+2t^2.
```

The monic equation for `t` over `U` proves finiteness. Independently,
the complete divided difference of `U` is
`p(p^2+p+1)-(2p+1)q`, where `p=s+t,q=st`.
The complete divided difference of `V` is
`2h5+3h4+2h2+2h1`, with `hj=(s^(j+1)-t^(j+1))/(s-t)`.
Their substitution gives the stated `-p^2 h(p)/(2p+1)`.
The denominator cannot vanish on the pair locus. Its zero solution
`p=q=0` is diagonal; the four simple nonzero roots of `h` each supply
one unordered pair of distinct parameters. Thus the pair locus is
finite, proving generic degree one.

The derivative gcd is exactly `t`, and the critical image has no other
parameter because `gcd(U,V)=t^2`. The local shear has second coordinate
`7t^5+8t^6+4t^7+2t^8`, while `U` has order two. Hence there is a single
analytic `(2,5)` critical branch. The tangent-degeneracy ideal has
support only at the diagonal cusp pair, so every off-diagonal collision
has two distinct tangent directions.

I independently reduced powers using `t^4=a-t^3-t^2`. This gives
`t^5=at-a+t^2`, `t^6=at^2-at+t^3`, and hence the displayed cubic
remainder `4t^3+(2a+5)t^2+at-3a-b`. Three distinct preimages force
this nonzero cubic to divide `U-a`. The three remainder equations,
the elimination to `j=4a^3+4a^2-63a`, and the Bezout constant `33957`
exclude every such fibre. Four distinct preimages are already impossible
for a nonzero cubic. Therefore the four pair sums give four **distinct
ordinary node images**, not a hidden triple or higher multiple point.

Finiteness and birationality identify `A1` with the normalization. Outside
the critical and multiple fibres a single immersed plane branch is smooth.
This completes the affine inventory; the genus check alone is not being
used to infer that inventory or transversality.

## 3. Projective embedding and both actual resolution chains

The degree-six homogeneous map has no base point and is birational.
Its image is therefore a sextic, not merely a parametrization whose
coordinate maximum happens to be six. There is one normalization point
over infinity and its target point is `[0:1:0]`.

I independently reconstructed the local rational charts using only
standard-library Fraction polynomial arithmetic, without importing
the producer or using its SymPy expressions. At infinity their inputs are

```text
X=z^2(1+z+z^2)/(2+3z+2z^3+2z^4),
Z=z^6/(2+3z+2z^3+2z^4).
```

Direct numerator/denominator valuations give the orders

```text
finite:  (2,5),(2,3),(2,1),(1,1),(1,1),
infinity:(2,6),(2,4),(2,2),(2,2),(2,1),(1,1),(1,1).
```

The terminal leading coefficients, first and second coordinates, are
respectively `(19/343,7)` at the finite cusp and `(-403/171500,35)`
at infinity. The independent first odd infinity residual has order nine
and coefficient `35/16`. In particular every asserted terminal smoothness
and transversality test is exact rather than based on decimal root fitting.

The divisor tracking was also checked separately from the valuations.
At infinity the first two blowups retain the strict line and new
exceptional as axes. At the third, the slope-four centering places the
line at `rho=-4`, away from the strict branch at zero. The fourth
blowup gives `(X,tau)` of orders `(2,1)`, with the curve tangent to the
new exceptional. The fifth gives two exceptional axes. Their tangent
ratio is `1/2450`; the sixth moves the curve away from both old axes,
leaving transverse intersection with the final exceptional.

The finite chain has the same final tangency/corner separation, with
ratio `1/49` and final second coordinate of leading coefficient seven.
It requires four blowups. Thus the exact centre multiplicities are
`(2,2,1,1)` and `(2,2,2,2,1,1)`. Their sum is16, their square sum is28,
and their singularity delta contribution is6. Consequently

```text
D^2=36-28=8,     N=4,     D^2-2N=0,
delta_finite+delta_infinity+N=2+4+4=10.
```

All ordinary nodes are untouched and all final curve/boundary intersections
are transverse away from exceptional crossings. This is a simultaneous
resolution with the original line at infinity, not just an abstract
resolution of the infinity cusp.

## 4. Strictness and the exact stopping boundary

The ordinary affine cusp `(t^2,t^3)` is a valid hostile to replacing
the strict Nori inequality by equality. Its finite resolution uses
`(2,1,1)`. Its projective infinity branch is smooth, with coordinates
`X=z,Z=z^3`; resolving its contact with `Z=0` requires three successive
multiplicity-one blowups, the last separating an exceptional corner.
Thus the full cost is7, equal to `3*3-2`, and its final square is zero.

Weighted radial scaling preserves the cusp and identifies the affine
complement with its link complement times a contractible radial factor.
The marked cusp relation admits the explicitly checked noncommuting
transposition images in `S3`. Hence abelianity at equality is false.
This control is not a Keller map. The new sextic's zero margin likewise
means only that **this strict criterion does not decide its complement**.
It is not evidence of a global `A6` representation or an affine source.

## 5. Full source read, replays, and pins

I read the complete source. It checks the exact declared two sextic
controls and ordinary cusp, not an enumeration of all curves. In particular
the pair denominator, separate tangent ideal, complete cubic divisibility
and Bezout obstruction, actual rational blowup charts, genus, and strict
equality distinction each have their own gates. The symbolic all-degree
identity retains genus explicitly. Its explicit exceptions remain active
under optimization.

I ran

```sh
python3 -B 04-computation/planar_jc48_sep06_resolution_budget.py
python3 -B -O 04-computation/planar_jc48_sep06_resolution_budget.py
```

Both pass **73 always-active exact gates**; normal, optimized, and frozen
stdout agree byte for byte. Recomputed pins match the primary note:

```text
source SHA-256:
47f1e5775e8447ddd199fe6c4f7f358ed3091031cb54cd1f54dd9a0480e2a387
output SHA-256:
d371ee9ccd8e4de711e5f243817bf428da0521e97921c216d8d14fd118e7087b
```

Final acceptance is **PASS**, with no mathematical or source repair.
The current audit is frozen. The new equality target remains outside the
strict criterion, and the wider Jacobian conjecture remains open.
