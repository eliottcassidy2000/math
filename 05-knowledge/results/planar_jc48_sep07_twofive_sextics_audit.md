# Independent full audit of the one-(2,5)-cusp (4,6) class

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS after the explicit
construction-order and degree-one repairs below.** This accepts the full
claim in [the primary proof](planar_jc48_sep07_twofive_sextics.md), including
the connected-family transport of the actual affine complement. It does
not infer a deformation theorem from the finite permutation census.
The primary source and output were not changed by this audit.

## 1. Exhaustiveness and marked local geometry

The allowed translations and nonzero affine scalings really reduce every
declared polynomial normalization to
`U=t^4+a t^3+b t^2`, `V=c t^6+e t^5+f t^3+g t^2`, `c!=0`.
They neither change the degree pair nor the whole-support question.
For `b!=0`, the order-three term after subtracting `(g/b)U` must vanish.
With `d=g/b`, the fifth coefficient after the local analytic shear
`V-dU+(d/b^2)U^2` is `e+2ad/b`; its nonvanishing is exactly the remaining
characteristic exponent five. That local shear is not being used as a
degree-preserving normalization of the global curve.

If `b=0,a!=0`, multiplicity two forces `g!=0`, and the independent leading
orders three and two give a (2,3) cusp. If `a=b=0`, the exact (2,5)
condition is `gf!=0`, as checked by the fifth term of `V^2-g^2U`.
Thus the exceptional normal form is retained rather than discarded by
division by `b`.

At infinity, the rational chart gives `ord(X)=2`, `ord(Z)=6` and the
seventh coefficient `(2e-3ca)/c^2`. If it vanishes in the `b!=0` case,
then `a=0` would make both polynomials even and violate birationality.
Consequently the rescaling to `beta=b/a^2` and `delta=2d/(ca^2)` is lawful.
The successive ninth, eleventh and thirteenth coefficients in the source
check the stated four types. The sextic degree and genus follow from the
birational, basepoint-free degree-six homogenization, so the genus ledger
is an actual curve calculation, not an inference about the degree of a
hypothetical Keller map.

I separately checked the proximity argument behind zero even coefficients.
After the first three infinity blowups, `rho=Z/X^3-c^2`; the line is at
`rho=-c^2`, and the previous exceptional intersection is at infinity.
A subsequent chart `(X,Y/X-k)` places the old exceptional component at
infinity, even when `k=0`. Its next finite branch point is therefore a
free centre. Once the orders are `(2,1)`, one free blowup and one satellite
blowup give the final transverse arrow. Hence vanishing `k4` or `k5`
does not alter the marked tree. In particular the infinity-nine result
applies to all of its stated normal forms, including the exceptional
`U=t^4` form. The separate
[geometry audit](planar_jc48_sep07_twofive_sextics_geometry_audit.md)
supplies independent explicit zero-coefficient controls and reconstructs
the infinity-thirteen tree and its marked homology values.

## 2. The good infinity-eleven locus is Zariski open

Here the base is `S=A1_C\{0,3/20}`, and

```text
U=t^4+t^3+beta t^2,
V=2t^6+3t^5+delta t^3+delta beta t^2,
delta=-(7+12beta)/8.
```

The load-bearing repair made during this audit is to **exclude second
preimages of the cusp before lifting the normalization through blowups**.
Without that order, the centre ideal could pull back to a noninvertible
ideal at an additional isolated preimage in the total parameter surface.
This is a gap in an unqualified global lifting argument, not a counterexample
to the class theorem.

The exclusion is algebraic and closed. Near the fixed parameter section,
`U=t^2(beta+t+t^2)` and the second factor is a unit. Thus no other preimage
of the cusp can approach `t=0` while the base stays in `S`. The point
`t=infinity` maps to infinity. The extra-preimage incidence is consequently
closed in `P1 x S`, so properness makes its image closed in `S`.

After this removal, each successive prescribed cusp-centre ideal pulls
back to a power of the parameter `t` times a unit. The characteristic
coefficient `-7/(4beta)` is a unit on this base. At infinity the analogous
coefficient `(3-20beta)/4` is a unit. The normalization therefore lifts
globally through these section blowups by invertibility of the pullback
ideals. The resulting ambient family is smooth and proper, and the source
of its normalization map is still the proper family `P1 x base`.

The nonimmersion locus is algebraically closed and has closed parameter
image by properness. Once removed, the map of the total normalization
family to the resolved ambient family is unramified: its vertical tangent
line injects, and on the base its differential is the identity. It is
also separated. Thus its diagonal in the fibre product is both open and
closed. The off-diagonal pair incidence is closed and proper over the
base; the pairwise-distinct triple incidence is likewise closed and
proper, using the three diagonal conditions. Nontransversality of the
two tangent lines is a closed algebraic condition on the pair incidence.
Its image, and the triple image, are Zariski closed by properness.

These arguments directly establish algebraic openness, without silently
equating analytic closure and a proper algebraic incidence. Nonbirational
normalizations are excluded: a generic pair over the same image has the
same tangent line and hence lies in the removed nontransverse incidence.
Conversely, a curve with the declared good singularities passes every
test. The genus ledger forces exactly three nodes. Therefore the resulting
set is exactly the good locus `B`. The audited `beta=1` representative
makes `B` nonempty. A nonempty Zariski open subset of the complex affine
line is the complement of finitely many points and is path connected.

## 3. Exact denominator hostile: beta=1/4 is actually good

I independently computed this control from its literal polynomials;
no source producer was imported. At `beta=1/4`, `delta=-5/4`. For
unordered preimages write `p=s+t`, `q=st`. The first divided equation is

```text
(2p+1)(2p^2+p-4q)/4=0.
```

The usual `2p+1` division therefore fails. At `p=-1/2`, however, the
second divided equation is the constant `-1/32`, so there is no missing
collision there. For all other pairs,

```text
q=p(2p+1)/4,
divided V=-p^2 h(p)/16,
h(p)=8p^3+20p^2+18p+7,
disc(h)=-2816,
Res_p(h,-p(p+1))=-7.
```

The `p=0,q=0` solution is the cusp diagonal. The three simple roots of
`h` give distinct unordered off-diagonal pairs. Their ordered parameter
projection is

```text
Res_p(h,t^2-pt+p(2p+1)/4)
 =(256t^6+640t^5+640t^4+320t^3+96t^2+28t+7)/4,
```

which has gcd one with its derivative. Thus these six preimages are all
distinct; no triple point can be concealed among the three pairs.
Also `gcd(U',V')=t`. Reducing the tangent determinant modulo the pair
quadratic and `h` gives

```text
(208p^2t+116p^2+288pt+178p+112t+91)/128.
```

Its norm along the quadratic, followed by resultant with `h`, is
`99648703/16777216 != 0`. Hence all three pairs are transverse.
The local characteristic coefficients are `-7` at the finite cusp and
`-1/2` at infinity; they give the required (2,5) and (2,11) branches.
Finiteness of the entire off-diagonal incidence proves birationality.
Therefore `beta=1/4` is a genuine member of `B`, not merely a parameter
which might be good. This concretely rejects excluding a parameter
because an auxiliary elimination denominator vanishes.

For reproduction of this additional control, use the displayed literal
`U,V`, form divided differences before specializing `p=-1/2`, then use
`sympy.rem` modulo `t^2-pt+q`. All later certificates above are one-variable
discriminants, gcds and resultants of the explicit displayed polynomials.
The control is supplementary; it does not replace the universal openness
argument of Section 2.

## 4. Simultaneous resolution transports the actual meridians

Over `B`, transversality makes each double-pair incidence locally etale
over the base, and properness makes it finite. The free swap of its two
preimages gives the finite etale multisection of the three actual node
centres. There are no triple coincidences. Blowing up that multisection is
therefore a simultaneous operation, without choosing global node labels.
The prescribed cusp and infinity resolutions and these node blowups give
a smooth proper ambient family. Its reduced total transform of the curve
and the marked infinity line is a relative normal-crossing divisor.
The strict curve is smooth after these resolutions; intersections of
components are smooth over the base, even if node branches permute around
a loop. No global numbering is required for relative normal crossings.

The local differential triviality proof in the primary is sufficient.
Along a compact base path, lift its real tangent field in relative divisor
charts by holding the local divisor coordinates fixed. A partition of
unity preserves both projection to the base field and tangency to every
divisor stratum. Properness ensures the lifted flow exists along the
whole compact path. This gives a diffeomorphism of the fibre pairs and
therefore of their open complements. Every centre belongs to the removed
divisor, so the open complement is precisely the original affine-plane
curve complement, not a projective quotient or a boundary link.

The distinguished strict curve remains marked. Complex normal orientation
is preserved by this continuous transport starting at the identity, even
though the patched real flow need not be holomorphic. Thus positive curve
meridians transport to positive curve meridians, with the expected access
conjugacies. In particular the two actual meridian generators of the
audited `beta=1` curve transport as generators of each complement group.
The nonzero scaling of its V-coordinate by eight is an affine complex
automorphism and does not alter this statement.

For a transitive action on `d>1` labels, a generator moving at most
`delta` labels contributes at most `delta-1` edges after replacing each
of its nontrivial cycles by a tree. Their union has the same orbits as
the group, hence needs at least `d-1` edges. The primary's general
inequality was qualified by `d>1` during this audit: the degree-one
identity action would otherwise be a literal exception. Nontrivial
transitivity forces `delta>=2`, so identity generators cause no issue.
At the actual passport `d=6,a=3,delta=3`, two generators contradict the
bound. Retained sheet counts are conjugacy invariant, so no marked-sheet
data are lost in transporting the meridians.

## 5. Infinity-thirteen argument and full exact replay

I checked the original vertex relations and all eliminations, retaining
the commutator `[c,d]=1`. A common square/cube in `A6` has type identity,
double transposition, or five-cycle. The double-transposition branch gives
an order-four `f`, impossible for a square in `A6`. At `c=1`, the orders
`1,2,4,5` of `d` cannot give the required single three-cycle meridian;
order three gives `e=d^2` and `mu=d^2`. This uses absence of order six
in `A6`. The remaining subgroup is a quotient of the (2,3,3) triangle
group, and an actual six-label action containing a single three-cycle
has support at most four. The primary's orbit argument is correctly
about the permutation action, not just the abstract name `A4`.

For five-cycle `c`, let `j` be the fixed label of `h=c^3`. If `mu` moves
`j`, its two transpositions first join the five-cycle and singleton,
then split the resulting six-cycle, giving exactly two cycles in `mu h`.
An even involution has four or six cycles, so `mu h=e^{-1}` is impossible.
Thus every fibre fixes `j`. This proves the group obstruction analytically.
The marked infinity-complement surjection used afterward is the actual
map already proved and source-checked in the cited plumbing supplier.

I read the entire source and replayed it normally and under `-O`.
Both outputs are byte-identical to the frozen output: **3,385 active
gates, 1,890 bytes**. The source checks the actual blowup weights and
marked homology, and ranges over every compatible common square/cube
and every terminal square root. Every original relation is rechecked
before counting. Its 1,120 assignments split into 40,360,720 with
respectively three,two,one fixed labels. Negative exponents are valid
because the exponent of `A6` divides sixty.

```text
python3 -B 04-computation/planar_jc48_sep07_twofive_sextics.py
python3 -B -O 04-computation/planar_jc48_sep07_twofive_sextics.py
source SHA256 02ef96f5a39fee7382fa81d0970fd06285c7bb3e180362bff234d4e2abad1e75
output SHA256 2713e10ba1715d4ae3b9d5d0095be91fa901ea3f63a16c9d53778e15b072e472
assignment SHA256 b35120fec36b3ded751aec092972d2e870466fb4d7defe0a193d41cc8412373f
```

The separate geometry audit additionally reconstructs the unsimplified
marked census and gives a nonvacuous infinity-thirteen curve control.
The full class theorem is accepted after the incorporated repairs.
Neither audit asserts that all rational sextics have this normal form,
that arbitrary infinity data determine an affine group, or that planar
JC is closed.
