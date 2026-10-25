# Six and seven are adjacent canonical thresholds, not isolated special integers

**Status:** PROVED here under the stated hypotheses; classical constructions
and adjunction mechanisms are recovered, not claimed as new. The all-incidence
non-torification dichotomy below is a deduction in this packet. The numerical
and matrix controls are FINITE-EXACT, not substitutes for the proofs.

## 0. Inheritance and the connection being made

- **Closest proved mechanism:** Session59, `toric-shared-roots/TORIFICATION_OBSTRUCTION.md`,
  Theorem T1: square roots of `n` distinct spanning plane lines, with no four
  concurrent, give a complete intersection with only A1 singularities and
  `K_Y=(n-6)H`; Theorem T2 uses its positive geometric genus against finite
  torification.
- **Hostile example:** `n-1` concurrent branch lines plus one transverse line.
  It has the same complete-intersection dualizing formula but its resolution
  is ruled, not general type. This invalidates an unrestricted canonical-class
  argument, **not** the non-torification conclusion.
- **Corrected near miss:** negative canonical class in dimensions at least
  three means Fano, not automatically rational; a singular ample dualizing
  sheaf need not give a general-type resolution.
- **Least-used sidecar:** both the dimension of the base and the full incidence
  arrangement; a root count alone omits both.

The live concept board is: canonical balance; coordinate sections; branched
covers versus their quotients; singularity discrepancies; birational torus
obstructions. The map from root labels to the cover preserves the function
field and the root equations. Passing only to its Hilbert series forgets
incidence and resolution discrepancies. Coordinate restriction preserves the
number of defining equations but changes dimension and canonical index.

## 1. The uniform smooth theorem

Work over `C`. Let `d>=1`, `r>=2`, `n>=d+1`, and let `L_0,...,L_(n-1)` be
linear forms on `C^(d+1)` in **linear general position**: every `d+1` of them
are independent. Thus their hyperplanes form a simple-normal-crossings
arrangement in `P^d`. Let `Lambda=P^d -> P^(n-1)` be their linear embedding.
Define `X(d,n;r)` as the inverse image of `Lambda` under

    [x_0:...:x_(n-1)] |-> [x_0^r:...:x_(n-1)^r].

Write `H=O_X(1)`, `c=n-d-1`, and `a=(r-1)n-r(d+1)=rc-n`.

### Theorem G1 (generalized Fermat root cover)

1. `X` is a smooth, geometrically integral projective `d`-fold and a complete
   intersection of `c` degree-`r` hypersurfaces in `P^(n-1)`.
2. The map `f:X -> P^d` is finite flat of degree `r^(n-1)` with deck group
   `(Z/r)^(n-1)`. Its ambient degree is `H^d=r^c`, and
   `f^*O_(P^d)(1)=O_X(r)`.
3. `omega_X=O_X(a)`. Consequently `a<0` gives a Fano variety, `a=0` gives
   trivial canonical bundle, and `a>0` gives an ample canonical bundle and
   general type. For `d>=2`, `H^i(X,O_X)=0` for `0<i<d`; thus the `a=0`
   case is Calabi–Yau, meaning here trivial canonical bundle with this
   intermediate structure-sheaf cohomology vanishing. For `d>=2` smooth
   projective complete intersections are also simply connected by Lefschetz. For `d=1`
   it is an elliptic curve. For `d=2` it is a K3 surface.
4. Its geometric genus is

       p_g(X)=[z^a] (1-z^r)^c/(1-z)^n,

   where a negative coefficient index means zero. If `0<=a<r`, this is
   simply `binom(n-1+a,n-1)`.
5. If `a>=0`, neither `X` nor any nonempty Zariski open subset of `X` has
   a finite decomposition into locally closed split tori.

**Proof.** Choose `d+1` independent forms as a basis. Every other root
coordinate satisfies one monic equation `x_j^r=sum b_ij x_i^r`. The quotient
ring is a free module over the polynomial ring in the first `d+1` roots,
with basis the monomials whose other exponents lie in `{0,...,r-1}`. The
monic equations form a regular sequence. This proves the complete-intersection
and degree assertions.

On each target coordinate chart the power map is the product of `n-1`
monic root extensions, hence finite flat of rank `r^(n-1)`. Its base change
is `f`. The classes of `L_j/L_0`, `j!=0`, are independent in
`C(P^d)^*/C(P^d)^(r)`: a relation with exponents in `{0,...,r-1}` has
valuation equal to its `j`-th exponent along `L_j=0`, forcing every exponent
to vanish. Kummer theory therefore gives a field of degree `r^(n-1)` at the
generic point. Flatness embeds each affine algebra into that field, proving
integrality. The coordinate-root actions give the deck group, with their
simultaneous scalar action removed projectively.

At any base point, at most `d` branch hyperplanes vanish, and their local
equations form part of a regular parameter system. Adjoining their `r`-th
roots simply replaces these parameters by their roots. Roots of the remaining
units give étale extensions. Thus `X` is smooth. This is the exact place
where linear general position enters.

Adjunction gives `omega_X=O_X(rc-n)`. The same answer follows independently
from ramification:

    K_X=f^*(K_(P^d)+(1-1/r) sum_j {L_j=0}),
    f^*O_(P^d)(1)=rH.

The Koszul resolution gives the Hilbert series and the cohomology vanishing.
For `a=0`, `H^0(omega_X)=C`; for `a>0`, the positive multiple of `H` is ample.
No rationality claim is made from `a<0` in higher dimensions.

Finally, a finite torification of an irreducible `d`-fold must contain a
`d`-dimensional torus stratum. A full-dimensional locally closed stratum is
dense and open, so the variety is rational. But for `a>=0`, the Hilbert
formula gives `H^0(omega_X)!=0`, a birational invariant of smooth projective
varieties that vanishes for projective space. The same function-field
argument applies to every nonempty open subset. QED.

The construction and canonical calculation are classical. In the literature's
notation `(d;k,N)`, translate **`k=r` and `N=n-1`**. Hidalgo–Hughes–Leyton-Álvarez,
*Automorphisms of Generalized Fermat manifolds*, arXiv:2010.04628v8,
Proposition 1 and §2.4 Proposition 2, give the smooth model and the canonical
balance. We use the proof's separated `a<0`, `a=0`, `a>0` cases, not an
unqualified positive-growth asymptotic outside its positivity range.
[Primary text](https://arxiv.org/html/2010.04628v8).

## 2. Why six and seven occur for square roots over a plane

For `r=2`, the coefficient is `a=n-2d-2`. Hence

| Root count | Canonical class | Ambient degree | Geometric genus |
|---|---|---:|---:|
| `n=2d+1` | `K=-H` | `2^d` | 0 |
| `n=2d+2` | `K=0` | `2^(d+1)` | 1 |
| `n=2d+3` | `K=H` | `2^(d+2)` | `2d+3` |

In the last row the root-coordinate embedding is the complete canonical
embedding, because there are no linear equations among the root coordinates.
Thus `p_g=2d+3` counts the independent canonical sections; it is not an
unrelated matching integer. Also `K^d=2^(d+2)`.

| Base dimension | Trivial-canonical root count | Next, canonically embedded member |
|---:|---|---|
| 1 | 4 roots: elliptic quartic in `P^3` | 5 roots: genus-5 canonical curve, degree 8 in `P^4` |
| 2 | 6 roots: K3, degree 8 in `P^5` | 7 roots: canonical surface, `p_g=7`, `K^2=16` in `P^6` |
| 3 | 8 roots: Calabi–Yau threefold, degree 16 in `P^7` | 9 roots: canonical threefold, `p_g=9`, `K^3=32` in `P^8` |
| 4 | 10 roots: Calabi–Yau fourfold, degree 32 in `P^9` | 11 roots: canonical fourfold, `p_g=11`, `K^4=64` in `P^10` |

The six-root K3 has `H^2=8`, hence **sectional genus 5**, since
`p_a(H)=1+H^2/2=5`. It is not a genus-6 polarized K3. The seven-root surface's
`p_g=7` is the dimension of its canonical sections, not a claim that it is a
genus-7 K3. On the nodal magic models use the pullback of `H` as a nef and big
quasi-polarization and invoke Session59's separate crepant resolution theorem.

Catanese already records the plane square-root sequence: four lines give a
quadric; five a degree-four del Pezzo; six a K3. He explicitly distinguishes
the maximal abelian cover from a single double cover of the same six lines.
*Catanese, Kodaira fibrations and beyond*, arXiv:1611.06617v1, §4.3, p.47.
[Author-hosted primary text](https://www.mathe8.uni-bayreuth.de/de/team/prof-fabrizio-catanese-old/pdf/156.pdf).

## 3. Varying the exponent: a complete small surface dictionary

For fixed `d,r`, a trivial-canonical member exists exactly when `r-1` divides
`d+1`, in which case

    n = r(d+1)/(r-1),      c=(d+1)/(r-1).

The first general-type root count is `floor(r(d+1)/(r-1))+1`. Its canonical
coefficient need not be one: adding one root increases `a` by `r-1`.

For surfaces, `a=(r-1)(n-3)-3`. Thus the complete **uniform-degree** lists are:

| Condition | `(r,n)` | Model | `(p_g,K^2)` |
|---|---|---|---|
| `a=0` | `(2,6)` | 3 quadrics in `P^5` | `(1,0)` |
| `a=0` | `(4,4)` | quartic in `P^3` | `(1,0)` |
| `a=1` | `(2,7)` | 4 quadrics in `P^6` | `(7,16)` |
| `a=1` | `(3,5)` | 2 cubics in `P^4` | `(5,9)` |
| `a=1` | `(5,4)` | quintic in `P^3` | `(4,5)` |

These are obtained by factoring `3` or `4`, respectively. A mixed-degree
`(2,3)` complete intersection in `P^4` is another K3 of degree six, but is
**not** in this simultaneous uniform-root model. The familiar threefold
uniform Calabi–Yau cases are `(r,n)=(2,8),(3,6),(5,5)`.

A recent primary preprint explicitly studies these Calabi–Yau generalized
Fermat varieties and their quotients: Billi–Palomino–Sarti,
arXiv:2610.21752v1 (2026-10-15), Proposition 3.6. Its Calabi–Yau condition
is exactly G1 after the notation shift above. Its quotient classification
is **not imported here** and does not automatically apply to singular magic
arrangements. [Primary text](https://arxiv.org/html/2610.21752v1).

### Cover-degree guard

The exponent-`r` cover above has total degree `r^(n-1)`, **not** `r`.
A single cyclic cover `w^r=product_j L_j` is a different object. If every
branch weight is one and it is constructed globally over `P^d`, one requires
`r|n`; its canonical divisor is

    K_Z = pi^*O_(P^d)(-d-1+(r-1)n/r).

This expresses the same ramification balance but gives neither the same
cover degree nor the same genus or projective embedding. Its normalization
and resolution need their own singularity check.

## 4. A real dimension-changing bridge: coordinate hyperplane sections

### Proposition G2

Under G1, for `d>=2`, setting one root coordinate `x_j=0` gives a smooth
root variety `X(d-1,n-1;r)`. It preserves `c=n-d-1` and raises the canonical
coefficient by one:

    a(d-1,n-1;r)=a(d,n;r)+1.

**Proof.** The base is the hyperplane `L_j=0`, isomorphic to `P^(d-1)`.
The other restricted forms are in linear general position: a dependency of
`d` restrictions would yield a dependency with `L_j` among `d+1` original
forms. Their inverse-image equations are exactly the remaining root equations.
The canonical increment also follows from adjunction to `H`. QED.

For three quadrics this gives the exact **genus-five ladder**:

    7-root Fano threefold, K=-H, H^3=8
       -> 6-root degree-eight K3, K=0, H^2=8
       -> 5-root canonical genus-five curve, K=H, deg H=8.

For four quadrics it gives:

    8-root Calabi–Yau threefold, degree 16
       -> 7-root canonical surface, p_g=7, K^2=16
       -> 6-root genus-17 curve, K=2H, deg H=16.

Here each arrow is the specified coordinate section, not deleting a root
while keeping the base dimension fixed. The magic arrangements are not in
linear general position: at triple points the restrictions of two other lines
can coincide on `L_j`. Their coordinate sections can therefore be singular or
reducible. G2 is a comparison-family bridge, not an unsupported smoothness
claim about those special sections.

## 5. The local singularity sidecar and an exact hostile example

For this section let `d=2`, retain distinct spanning lines, but drop general
position. Suppose exactly `m>=3` lines pass through a point. Normalize a
nonvanishing line to one. After separating the roots of unit radicands by an
étale extension, the local cover is the affine cone cut out in `A^m` by
`m-2` homogeneous degree-`r` equations expressing the dependent linear forms
in terms of two independent ones.

The projectivized cone is the smooth curve `C=X(1,m;r)`, because distinct
lines in the pencil correspond to distinct points of `P^1`. Write `D=r^(m-2)`.
Then

    deg O_C(1)=D,
    g(C)=1+(r(m-2)-m)D/2.

Blowing up the vertex is the total space of `O_C(-1)` and resolves the cone.
Its exceptional curve `E=C` has `E^2=-D`. Ambient blowup adjunction gives

    K_res = pi^*K_cone + b(m,r)E,
    b(m,r) = m-1-r(m-2).

Indeed the blowup of `A^m` contributes `(m-1)E`, and the `m-2` degree-`r`
equations each subtract `rE`. There are `r^(n-m-1)` such vertices above the
original `m`-fold point: choose the `n-m` nonzero roots modulo common scalar.
Away from these vertices the cover is smooth. The complete intersection is
Cohen–Macaulay and smooth in codimension one, hence normal.

| Exponent | Concurrent lines | Exceptional curve | Discrepancy |
|---:|---:|---|---:|
| 2 | 3 | conic, genus 0, self-intersection -2 | 0 |
| 2 | 4 | degree-4 elliptic curve | -1 |
| 3 | 3 | plane cubic, genus 1 | -1 |
| 4 | 3 | plane quartic, genus 3 | -2 |

Thus the square-root triple case is exactly crepant A1. Higher concurrence
or higher exponent can introduce nonrational exceptional curves, which an
adjunction-only proof would miss.

### Hostile cone H1

Take `L_0=z` and `L_i=x+t_i y` for `1<=i<=n-1`, with distinct `t_i`.
The latter `n-1` lines concur at `[z:x:y]=[1:0:0]`. The root surface is the
projective cone, with vertex the `x_0` coordinate point, over `X(1,n-1;2)`.
Its normality follows from the preceding local argument. Its minimal
resolution of singularities is a ruled surface over a curve of genus

    g = 1+(n-5)2^(n-4).

For `n>=7`, the singular surface has ample dualizing sheaf `O(n-6)`, but
its smooth resolution has

    kappa=-infinity,   p_g=0,   q=g,   K^2=8(1-g).

This is an explicit counterexample to transferring positive dualizing class
to general type without a discrepancy check. For `n=7`, it is a cone over a
genus-17 curve, with discrepancy `-3` and resolution `K^2=-128`, not `16`.

The ruled assertion can be checked either by the projectivized bundle
resolution of a cone or directly from the function field: divide the root
coordinates by one nonzero concurrent root. The remaining ratios give the
curve's function field and the chosen root supplies one independent rational
parameter.

## 6. Stronger survivor: no incidence hypothesis is needed for non-torification

### Theorem G3 (all-incidence plane obstruction)

Let `n>=4`, `r>=2`, and let the `n` branch lines be distinct and span the
plane dual space, over `C`. If

    a=(r-1)n-3r >= 0,

then the normal root surface and every nonempty open subset admit no finite
torification by split tori. In particular **every arrangement of at least
six distinct spanning square-root lines is nonrational and non-torifiable**,
not only arrangements without fourfold points.

**Proof by dichotomy.** If there is an `m`-fold point whose local exceptional
curve `C=X(1,m;r)` has genus at least one, project onto just the `m` root
coordinates belonging to those concurrent lines. Their forms span a
2-dimensional vector space, so their root equations define this curve `C`.
The rational map from the surface onto `C` is dominant: after taking powers,
it maps onto the pencil `P^1`, and its image in the irreducible curve must
therefore have dimension one.

A rational surface cannot dominate a positive-genus smooth curve. Indeed a
nonconstant rational map from `P^2` to that curve restricts nonconstantly to
a suitable line; it extends to a nonconstant morphism `P^1 -> C`, impossible
by Riemann–Hurwitz. Thus the root surface is nonrational in this case.

Otherwise all exceptional curves have genus zero. The genus formula shows
that the only singular case is `r=2,m=3`; these are the crepant A1 nodes.
When `r>=3` there are no triple points, so the surface is smooth. In either
case adjunction and, if needed, the crepant resolution give a smooth
projective model with `p_g=h^0(O_X(a))>=1`, so it is nonrational.

The finite-torification-to-rationality argument of G1 completes the proof,
including for all nonempty opens. QED.

For squares, the dichotomy is especially transparent:

    no fourfold point -> p_g>=1 by the crepant complete intersection;
    some fourfold-or-higher point -> positive-genus pencil-root quotient.

The hypothesis `n>=6` is a universal sufficient threshold, not an iff based
solely on root count. Five general-position square-root lines give a rational
del Pezzo; four concurrent square-root lines plus one transverse line already
give an elliptic ruled resolution and are nonrational. The distinction also
shows why a rationality criterion cannot be read off a Hilbert series alone.

This theorem blocks finite torifications and split toric models. It does not
block blueprints or every proposed framework called geometry over `F_1`.

## 7. Optional quantitative extension: resolved invariants of all plane arrangements

Let `t_m` count the points on exactly `m` lines, and let `Y` be the resolution
that blows up each cone vertex as above. This is the minimal resolution of
singularities, not necessarily the surface's minimal model. The exceptional
curves from different vertices are disjoint. With

    a=(r-1)n-3r,  b_m=m-1-r(m-2),  D_global=r^(n-3),

one obtains

    K_Y^2 = D_global * [a^2 - sum_(m>=3) t_m*b_m^2].

This follows because `pi^*H.E=0`, each `E^2=-r^(m-2)`, and the number of
such curves over an `m`-fold point is `r^(n-m-1)`.

For square roots the formulas simplify to

    K_Y^2 = 2^(n-3) * [(n-6)^2 - sum_(m>=4)t_m*(m-3)^2],
    e(Y)  = 2^(n-3) * [12-4n+t_2+sum_(m>=3)m*t_m].

For the Euler formula, the branch complement and the union of open line
pieces have Euler characteristics

    3-2n+sum_(m>=2)(m-1)t_m,
    2n-sum_(m>=2)m*t_m.

Their inverse-image covering degrees are `2^(n-1)` and `2^(n-2)`. Double
points have `2^(n-3)` reduced preimages each. At higher points replace each
preimage vertex by its exceptional curve, whose Euler characteristic is
`-(m-4)2^(m-2)`. Adding these strata gives the stated formula. This counts
reduced topological strata, not scheme-theoretic branch fibers.

Since `sum binom(m,2)t_m=binom(n,2)`, write `e_smooth` and `chi_smooth` for
the general-position values with the same `n`. Then

    e(Y)=e_smooth - 2^(n-4) sum_(m>=4)t_m*m*(m-3),
    chi(O_Y)=chi_smooth - 2^(n-6) sum_(m>=4)t_m*(m-3)*(m-2).

The last line uses Noether's formula. It computes `chi=1-q+p_g`, **not**
`p_g` alone: nonrational singularities can change `q`. Triple points cancel
from all these resolved numerical invariants, exactly explaining Session59's
incidence-independent values. Fourfold points do not cancel.

These are instances of the classical Hirzebruch–Kummer arrangement method;
we give the derivation to expose its hypotheses and the local mechanism.
Hirzebruch, *Arrangements of Lines and Algebraic Surfaces* (1983), is the
primary historical source for these covers and their surface invariants.
[Primary scan](https://hirzebruch.mpim-bonn.mpg.de/244/1/hirzebrucharr.pdf).

## 8. Verification and boundaries

Run:

    python3 check_thresholds.py

The standard-library test performs independent Koszul/free-basis Hilbert
calculations; exact Vandermonde general-position and Jacobian-support checks;
coordinate-section controls; threshold tables; curve-genus and local-discrepancy
checks; and exact incidence, Euler-stratum, Noether, and hostile-cone controls. `threshold_results.json` records the bounded
universe. The formulas and cone dichotomy above are proofs, not extrapolations
from those tests.

No statement here identifies the magic arrangements with a generic member,
proves an integer-square solution exists, classifies their Picard groups,
identifies their K3 surfaces as Kummer surfaces of abelian varieties, or
connects the root count six/seven to root systems `E_6/E_7`. None of those
claims follows from matching integers. In particular “Kummer cover” is the
name of an abelian branch-cover construction; it does not by itself mean
“Kummer surface,” the resolution of an abelian-surface sign quotient.
