# Finite dual-field closure reconstructs the inverse; finite coefficients alone do not

2026-10-13. **PROVED:** the elementary finite-closure certificate and the
characteristic-p boundary below. **FINITE-EXACT:** the declared rational and
modular controls. **CITED, supplied accepted inputs:** the two paper interfaces
in section 1. **OPEN:** obtaining closure for every planar Keller map. This
packet neither proves that assertion nor reopens the excluded cotangent route.

## 1. Inheritance and the paper mechanisms

The current mathematical route is
[THM-4620, Keller dual-field module/kernel/de Rham dictionary](../../01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md)
and its strengthened sequel
[THM-4626, generic line fibres imply local nilpotence](../../01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md).
Generic rank zero, zero cokernel, surjectivity and local nilpotence of a dual
field are already equivalent. A nonzero torsion cokernel is not an outstanding
escape route. In the plane this criterion for invertibility is inherited,
not a new Jacobian theorem here.

The canonical hostile is the inherited map of
[THM-4621, coordinate-pair fibrations of THM-1300](../../01-canon/theorems/THM-4621-thm1300-coordinate-pair-fibrations-are-rational-functions-minus-critical-points.md):
all its dual fields fail local nilpotence, though their kernel ladders have
the same displayed form as the invertible case. The corrected near miss is
[the cusp-length note](keller_module_length_cusp_20261012.md). Its replacement
of the guessed length four by three is stated **on the fibre stratification**,
as a sketch with a finer-stratum caveat. We use neither a global length-three
claim nor the false inference that each fibre-count drop needs a new simple
constituent. Generic cover rank, module length and a bounded vector-space
rank are different objects.

The less-used sidecar is the complete action on actual polynomial generators,
including the source point at which their joint values are evaluated. The
anchor is planar inverse reconstruction; the niche is an exact finite
closure consumer; the wildcard is the failure of prime-characteristic
nilpotence as evidence for characteristic-zero closure. The concept board
is: marked inclusion, bilinear products, derivative closure, integral lifting,
and generic sheet rank.

The supplied paper **Reconstruction of Function Fields from Mod-ell Milnor
K-Theory**, OpenAI, October 5, 2026, Theorem 1.1, is accepted as a cited input.
Its datum consists of the full degree-one and degree-two mod-ell Milnor
groups with their bilinear product. Section 5 retains a common incidence
compositum, rather than taking a norm which could kill classes; Lemma 7.1
retains simultaneous values of finite tuples. A finite coefficient field is
not a finite-dimensional datum. Its conclusion reconstructs the abstract
perfect function field with constants and the stated scalar/Frobenius
ambiguities; it does not supply a polynomial affine model or a marked
Keller inclusion. Original supplied attachment:
[mod-ell-bogomolov-pop.pdf](C:/Users/Eliott/Downloads/mod-ell-bogomolov-pop.pdf).
Extracted local text:
`C:/Users/Eliott/.codex/scratch/math-papers-20261013/mod-ell-bogomolov-pop.pdf.txt`,
pages 2--4, 19--24 and 26--28.

The supplied **The Deligne--Drinfeld conjecture**, OpenAI, September 23, 2026,
Theorem 1.1, is also accepted, not independently audited here. Its sections
2, 4 and 8 distinguish an integral deletion detector, saturated lattices,
and the subsequent dimension comparison. The useful imported mechanism is
that reduction and lifting need an actual integral compatibility argument.
Original supplied attachment:
[paper-28.pdf](C:/Users/Eliott/Downloads/paper-28.pdf).
Extracted local text:
`C:/Users/Eliott/.codex/scratch/math-papers-20261013/paper-28.pdf.txt`, pages 3--5
and section 4. Neither supplied headline is a premise of our elementary
closure proof. The transfer is the design of the retained data and the
test for a lost coordinate.

## 2. Exact finite closure, with all quantifiers retained

Let k have characteristic zero, A=k[x_1,...,x_n], and let F=(F_1,...,F_n)
be a polynomial map with nonzero constant Jacobian. Its polynomial dual
derivations D_j satisfy

    D_j F_i = delta_ij,       [D_i,D_j]=0.

Suppose a supplied finite-dimensional k-subspace W of A contains
1,x_1,...,x_n and satisfies **D_j W subset W for every j**. Put M=dim W.
Then F has a polynomial inverse, of total degree at most M-1.
Conversely, every polynomial automorphism has such a finite W.

Here closure means an exact identity between source polynomials. Replacing
each derivative by its projection to W does not satisfy the hypothesis.
No closure under multiplication is required: the actual derivation law
on A transports the finite generator information to their products.

**Proof of the finite bound.** THM-4620 gives no nonzero eigenvalue for any
D_j on k-bar(x). For completeness, if D_j h=lambda h with lambda nonzero,
differentiate the minimal polynomial of h over k-bar(F). Its coefficients
satisfy rational one-variable equations a'=c a with c nonzero. Such an
equation has no nonzero rational solution: a denominator would have to
divide its lower-degree derivative, and a polynomial fails by degree.
The minimal polynomial would consequently be a power of T, forcing h=0.
Thus each restriction D_j|W is nilpotent.

Commuting nilpotent matrices have a common strict triangular flag, already
over k. Indeed, successively restrict to nonzero kernels of the commuting
operators to find a common killed vector, then induct on the quotient.
Every product of M matrices is therefore zero. In particular

    D^alpha x_i = 0 whenever |alpha| >= M.                 (1)

Leibniz now gives finite simultaneous Taylor expansions on every polynomial
in the x_i. Equivalently, the commuting exponentials define a polynomial
additive-group action, taking F to F+t. This deduction uses the full
polynomial-ring multiplication, not unrelated finite vector-space matrices.

Fix any retained k-point a of the original affine source. The inverse is

    G_i(T) = sum_(|alpha|<M)
       (D^alpha x_i)(a) (T-F(a))^alpha / alpha!.           (2)

Exponentiating first and then evaluating at a gives F(G(T))=T. The coordinate
functions of G are algebraically independent because T lies in their
generated algebra. Thus G has dense image, and G(F(G(T)))=G(T) implies
G(F(x))=x by density. This proves both inverse identities and the degree
bound. In the converse direction, write x_i=G_i(F); the span of 1 and all
mixed target derivatives of these finitely many G_i, pulled back by F,
is finite and is stable under every D_j. It contains all x_i. QED.

For the plane, one can instead use a finite subspace stable under a single
dual field, then the inherited planar kernel/slice theorem. The implemented
consumer keeps all dual fields so that the explicit inverse formula (2)
works in every dimension without choosing an extra kernel coordinate.

## 3. A complete success certificate and an honest bounded search

Define W_r as the span of 1 and the mixed derivatives D^alpha x_i with
|alpha|<=r. Exact equality W_(r+1)=W_r certifies closure and is enough for
section 2. A strict inclusion merely says that this truncation has not
closed. It is not a certificate of noninvertibility.

The script's search instead adds independent derivative columns in batches;
the same criterion applies. It has explicit dimension and round caps, and
returns `None` if they are exhausted. For example, the invertible map
(x+y^3,y) fails the dimension-three cap, while (x+y^2,y) already has the
closed space span(1,x,y). This is the cheapest test against mistaking a cap
for a theorem.

A `Closure` packet retains the variables, the entire source map, an
independent polynomial basis, every rational action matrix, and the source
anchor. `audit` recomputes the Jacobian, all literal derivative identities,
coordinate membership, matrix commutation and nilpotence. `reconstruct`
uses (2) and checks **both** exact polynomial compositions, with simultaneous
substitution. A source-map change with the old matrices is rejected.
Boolean, floating-point, string-expression and malformed matrix aliases
are rejected before use. The implementation is over QQ; the proof is over
any characteristic-zero field.

For the declared positive family

    F_(a,b)=(x+y^a, y+(x+y^a)^b),       1<=a,b<=3,

the independent inverse is

    (X-(Y-X^b)^a, Y-X^b).

The complete closure dimensions are

    (a,b): (1,1)3 (1,2)3 (1,3)4 (2,1)3 (2,2)6
           (2,3)8 (3,1)4 (3,2)10 (3,3)14.

Both anchors (0,0) and (1,-1) reconstruct the same inverse. This checks
that the retained anchor is an evaluation choice, not a new sheet or a
replacement source.

## 4. Prime-characteristic nilpotence is automatic, not a closure test

Let F be **any** Keller map over a field of characteristic p>0. For each
dual field,

    D_j^p = 0 on the entire rational function field.       (3)

The p-th iterate of a derivation is a derivation, by the binomial Leibniz
formula in characteristic p. It kills each F_i and hence k(F). The extension
k(x)/k(F) is finite separable because the Jacobian is invertible. A derivation
vanishing on the base vanishes on a separable algebraic extension. This
proves (3) without invertibility of F.

Consequently reducing a characteristic-zero Keller map at arbitrarily many
good primes and observing nilpotence at **the varying depth p** proves
nothing about its characteristic-zero local nilpotence. The inherited
THM-1300 map is an explicit hostile: all nine coordinate/dual pairs vanish
at depth p for p=3,5, while its characteristic-zero V_3^k(x) is nonzero for
all k by THM-4621. The checker independently computes k=1,...,5, with total
degrees 11,21,31,41,51. These finite calculations illustrate the cited
all-depth statement; they do not replace it.

**A lawful lifting repair.** For an integral derivation
delta=sum f_i partial_i, put d=max(1,max deg f_i) and
C=max_i ||f_i||_1. Then, on any coordinate,

    ||delta^N x_i||_1 <= C^N product_(j=0)^(N-1)(1+j(d-1)). (4)

Indeed deg(delta h)<=deg h+d-1 and
||delta h||_1<=C deg(h)||h||_1. For a rational dual field first multiply it
by a common integer denominator. At a **fixed** N, a prime larger than the
right side of (4), with a zero reduction, proves literal zero over the
integers. Alternatively, zero reduction at infinitely many good primes at
the same fixed N proves zero. Neither repair follows from (3): its depth
changes with p and its coefficient bill is unpaid.

The elementary non-Keller derivation x^2 partial_x makes the price exact:
delta^N x=N!x^(N+1). Its third derivative vanishes modulo 3 but is 6x^4
over QQ. The right side of (4) is exactly N!. This auxiliary example only
tests the lifting implication; the preceding Keller example supplies the
actual characteristic-p hostile in the intended category.

## 5. Bilinear reconstruction needs the marked inclusion and affine model

For a generic finite separable extension L/E, the full trace pairing
Tr_(L/E)(uv) is nondegenerate. This is inherited in the repository's
[THM-3784, rational Keller different/codifferent trace duality](../../01-canon/theorems/THM-3784-rational-keller-tower-different-codifferent-trace-duality.md),
not a new trace theorem. In L=k(t)[z]/(z^3-t), the marked basis (1,z,z^2)
has Gram matrix

    [ 3   0   0 ]
    [ 0   0  3t ]       determinant = -27t^2.
    [ 0  3t   0 ]

Both Tr(z) and Tr(z^2) vanish, while the full pairing has rank three.
Scalar trace cancellations cannot identify the sheet degree. Nor is this
generic degree the composition length of a holonomic module.

Even full abstract field reconstruction requires another sidecar for JC:
both k(x,y) and its proper subfield k(x^3,y) are abstract rational function
fields. The marked inclusion has degree three. The latter map is not
Keller, which is precisely why this example tests the field-information
claim rather than purporting to be a planar counterexample. Even a marked
birational field map can lose the affine model: (x,y)->(x,y/x) is a rational
isomorphism but is not regular on x=0. Formula (2), in contrast, checks the
entire original polynomial source and both polynomial inverse identities.

| Source | Retained map or predicate | Missing data / decisive test |
|---|---|---|
| Supplied mod-ell field datum | Full bilinear product and common tuple incidence | Finite coefficient field does not make the datum finite-dimensional |
| Abstract reconstructed field | Field operations | Marked inclusion and regular affine coordinate ring |
| Finite derivative span | Every dual action on actual coordinate generators | Exact closure, not projected closure or a rank census |
| Closed characteristic-zero span | One common flag and finite Taylor action | Retained evaluation point and polynomial multiplication |
| Prime reduction | D_j^p=0 automatically | Fixed-depth integral lifting with a coefficient bill |
| Generic trace scalars | Some averaged sheet data | Full marked bilinear pairing |

## 6. Precise next obligation and reproduction

The positive instrument is a finite inverse compiler. Its next mathematical
obligation is a reason that the **actual** derivative spans close for a
specified unresolved planar Keller class, or a degree/height estimate that
makes a fixed-depth modular lifting test succeed. No such uniform estimate
is claimed. Adding a dimension cap, a changing prime, a new abstract field
description, or a repeated finite kernel identity cannot discharge it.

This uses the repository's patterns “retain local incidence through
counting” and “separate descent, ambient scale, and regularity debt.” It
adds a concrete counterindication: prime-depth nilpotence is universal in
the reduction and so does not detect the target distinction.

[Script](../../04-computation/experiments/jc_bilinear_reconstruction_20261013a.py)
and [saved output](jc_bilinear_reconstruction_20261013a.out):

```text
python -B -X utf8 04-computation/experiments/jc_bilinear_reconstruction_20261013a.py
python -B -O -X utf8 04-computation/experiments/jc_bilinear_reconstruction_20261013a.py
```

The finite universe is nine two-shear maps at two source anchors, 25 integer
source points per receipt, exact inverse compositions, the cubic trace
packet, all three THM-1300 dual fields on all three coordinates at p=3,5,
five rational derivative depths, eight integral coefficient bills, and
malformed/false packet controls. All tests use exact arithmetic and remain
active under `-O`. Normal and optimized runs match the saved output at
**565 exact checks**. LF SHA256:
`1ee212755c3f79fa824f006086fec016169303805f663e6cc0550369f9a93ae8`.
