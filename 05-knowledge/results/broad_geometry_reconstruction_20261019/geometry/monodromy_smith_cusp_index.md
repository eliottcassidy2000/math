# Exterior monodromy can miss the toric cusp index

Status: PROVED standalone integral linear algebra and mapping-torus
homology. The toric conclusion is conditional on the explicit translation
quotient hypotheses of THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.
No analytic S6 existence claim is imported.

## Inheritance and the missing sidecar

The earlier Hopf/S6 audit verifies the manuscript's unimodular rank-two
cusp shear and the exterior invariant ranks `(1,2,4,2,1)`. THM-3991 gives
`chi(W)=2d` and `d` irreducible components for an index-`d` translation
quotient of the A2 fan. The new result identifies exactly which integral
monodromy degrees see that index and which erase it.

## Theorem C: complete exterior Smith profile

Let `L` be free of rank four and `N:L->L` satisfy `N^2=0`, `rank N=2`.
Set `M=I+N` and `K=ker N`. Then `K` is primitive of rank two, and `N`
induces an injective rank-two map

    B: L/K -> K.

Let its Smith factors be `s|t`, with `s,t>0`. The nonzero Smith factors of
`wedge^q(M)-I` are

    q=0: none; q=1: (s,t); q=2: (s,s); q=3: (s,t); q=4: none.

Consequently the invariant ranks are always `(1,2,4,2,1)`, and

    coker(M-I) = Z^2 + Z/s + Z/t,
    coker(wedge^2(M)-I) = Z^4 + (Z/s)^2,
    coker(wedge^3(M)-I) = Z^2 + Z/s + Z/t.

Here `Z/1=0`. In particular saturation of the degree-two image means only
`s=1`; it does not imply `B` is unimodular. Saturation in degree one or
three is equivalent to unimodularity.

### Proof

The kernel of a homomorphism of free abelian groups is primitive. Since
`N^2=0`, its rank-two image lies in its rank-two kernel and has finite index
there. The quotient `L/K` is free, so choose an integral splitting. Change
bases in the two rank-two factors to put the induced map in Smith form.
There is thus a basis `e1,e2,f1,f2` of `L` with

    N e1=N e2=0, N f1=s e1, N f2=t e2.

The degree-one statement follows immediately. In degree two use the basis

    e1e2, e1f1, e1f2, e2f1, e2f2, f1f2,

where juxtaposition denotes wedge. The only nonzero images under
`wedge^2(M)-I` are

    e1f2 -> t e1e2,
    e2f1 -> -s e1e2,
    f1f2 -> st e1e2+s e1f2-t e2f1.                       (1)

The operator has rank two. Its entries have gcd `s`, while its two-by-two
minors have gcd `s^2`: all are divisible by `s^2`, and the minor on rows
`e1e2,e1f2` and columns `e2f1,f1f2` is `-s^2`.
Thus its two nonzero Smith factors are `(s,s)`.

Since `det M=1`, integral exterior duality identifies `wedge^3(M)` with
`M^(-T)`. The difference `M^(-T)-I=-M^(-T)(M^T-I)` has the same Smith
factors as `M-I`, proving degree three. Degrees zero and four are fixed.
Subtracting ranks from the exterior dimensions proves the invariant ranks
and cokernel descriptions.

## Corollary D: the collar's middle homology also loses the index

Let `Y_M` be the mapping torus of the orientation-preserving automorphism
of the real four-torus whose action on `H1` is `M`. Then, as abstract
abelian groups,

    H0(Y_M)=Z,
    H1(Y_M)=Z^3 + Z/s + Z/t,
    H2(Y_M)=Z^6 + (Z/s)^2,
    H3(Y_M)=Z^6 + Z/s + Z/t,
    H4(Y_M)=Z^3,
    H5(Y_M)=Z.

Proof: the integral Wang sequence gives a short exact sequence

    0 -> coker(wedge^k M-I) -> Hk(Y_M)
      -> ker(wedge^(k-1) M-I) -> 0.

The final term is free, so this sequence splits as abstract abelian groups.
Use Theorem C in each degree. The splitting is not claimed canonical, nor
is a cohomology-ring classification asserted. The torsion pairing pattern
also passes the necessary Poincare-duality check: torsion in degrees one
and three agrees, while degree-two torsion comes in two equal factors.

## Corollary E: recover the cusp index in degree one, not degree two

Assume the toric setting of THM-3991 in rank two: the vertex/cocharacter
lattice is `K`, the locally finite periodic triangulation is unimodular,
the deck translations have image `Gamma=B(L/K) subset K`, and the toric
central fibre quotient `W` is compact, with the stated free proper action
on its neighborhood. Then

    d=[K:Gamma]=s t,
    number of irreducible components of W=s t,
    chi(W)=2s t.                                         (2)

The first equality is Smith normal form; the remaining equalities are
THM-3991. Equivalently,

    chi(W)=2 * |Tors coker(M-I)|.

This last formula is valid in the stated full-rank translation-fan grammar,
not for arbitrary degenerations with this rational monodromy type.

If such a cusp is the sole nonzero-Euler fibre of a compact torus-fibred
complex threefold, sphere homology requires `s=t=1`. This isolates the
unimodular boundary from the whole integral shear family, but does not
verify that any global filling exists.

## An infinite hostile family

Take `B_d=diag(1,d)` for any positive integer `d`. All the following remain
constant with `d`:

- all exterior invariant ranks `(1,2,4,2,1)`;
- the complete degree-two Smith normal form `(1,1,0,0,0,0)`;
- `H2(Y_M)=Z^6` as an abstract abelian group.

But `H1(Y_M)=H3(Y_M)` at the torsion level contains `Z/d`, while the
conditional toric cusp has `d` components and Euler characteristic `2d`.
Thus degree-two saturation, rational monodromy invariants, and middle
collar homology do not determine component count or the Euler number.
The exact lost coordinate is the second Smith factor `t`; retaining
`coker B`, or the degree-one integral monodromy cokernel, restores it.

This does not refute the manuscript's actual calculation: its displayed
`B0` is unimodular. Instead it explains why that precise integrality input
is indispensable and gives a controlled family of alternative cusp data.

## A top-degree specialization obstruction

Under the same toric quotient hypotheses, the top homology of `W` has rank
`d`, since its `d` irreducible compact complex-surface components supply
independent fundamental classes, with intersections of real dimension at
most two. Equivalently, remove the conductor and use top-dimensional
relative chains; identifying lower-dimensional strata cannot create a
relation between the component fundamental classes.

Hence the rational top cohomology has dimension `d`. A nearby smooth
complex two-torus has one-dimensional `H^4`. If the toric space is realized
as a proper degeneration with the usual cohomological specialization map

    sp^4:H^4(W;Q)->H^4(F;Q),

then injectivity of that map is impossible for `d>1`. Thus the all-degree
injective invariant-cycle conclusion used in the manuscript cannot extend
unchanged to these nonunimodular translation quotients. This argument only
uses ranks; it does not assert the exact integral kernel, its image, or an
analytic realization of the degeneration.

## Connection contract

Source: an integral square-zero rank-two monodromy operator, with the
saturated kernel lattice retained.
Target 1: its exterior degree-two operator and middle mapping-torus
homology.
Map: exterior power and the Wang exact sequence.
Preserved predicate: degree-two image saturation and the abstract middle
homology group.
Lost information: second Smith factor `t` and the cusp lattice index `st`.
Sidecar: `coker B`, equivalently the degree-one integral monodromy cokernel.
Target 2: the periodic toric cusp with deck image `B(L/K)`.
Map: identify the shear image with the translation-period sublattice.
Preserved predicate: finite index, hence component and Euler counts under
the explicit quotient hypotheses.
Decisive hostile: `B=diag(1,d)` with unbounded `d`; positive control `d=1`
recovers the manuscript's local lattice profile and Euler value two.

## Literature direction and boundary

The source's infinite A2 fan is a specialized Mumford degeneration. The
primary background is David Mumford, *An analytic construction of
degenerating abelian varieties over complete rings*, Compositio Mathematica
24 (1972), 239-272:
https://www.numdam.org/item/CM_1972__24_3_239_0/

A particularly relevant modern primary survey is Engel, de Gaay Fortman,
and Schreieder, *Combinatorics and Hodge theory of degenerations of abelian
varieties: A survey of the Mumford construction*, version January 28, 2026:
https://philip-engel.github.io/SurveyMumford.pdf
and https://arxiv.org/abs/2507.15695 . Its Sections 3 and 5 develop toric
models and collapse maps; Figure 20 describes how a four-valent crossing
is replaced by two three-valent crossings under a small resolution.

A useful next construction is to compare the square fan of a product of
nodal elliptic degenerations with its A2 diagonal subdivision: locally
`xy=uv=t` is the ordinary-double-point model and the two diagonal choices
are its small toric resolutions. This is a concrete geometric route to
realizing the `B=diag(1,d)` lattice controls, but global compatibility and
analytic quotient existence have not been proved in this note. It remains
a proposed construction, not a proved corollary.

## Reproducibility

`verify_geometry_extensions.py` checks all five exterior degrees for 100
shears `1<=s<=10, t=s*m, 1<=m<=10`, plus four non-diagonal integral shear
controls and the unbounded-family sample `d=1,2,3,5,12,101`. All arithmetic
uses exact integer Smith forms. The unrestricted theorem is proved above;
the finite checks are independent controls on signs, ordering, and minors.
