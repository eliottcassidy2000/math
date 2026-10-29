# What the symplectic F1 skeleton keeps, and which paired sectors it loses

2026-10-29. **PROVED in finite symplectic linear algebra**, with exact finite
controls. The underlying counting identities and Weyl-group interpretation
are classical. This packet's contribution is an explicit, source-preserving
connection contract for the requested doubling/squaring analogy, together
with the coupling-rank filtration and hostile controls. It proves no
finiteness theorem or realization theorem for a Shafarevich–Tate group.
Source snapshot: `d92d6c0d742c7bab872d0e62769e70e9f0595d2c`.

## 1. Inherited F1 boundary and the precise literature input

The repository's [marked root-algebra repair](../../normseven_f1_20261025/KUMMER_AND_F1.md)
retains root order, base radicands and additive branch equations. Its
[twist-average calculation](../../normseven_f1_20261025/TWIST_AVERAGE.md)
finds the formal value seven at q=1 but a negative coefficient in q-1,
blocking a positive split-torus interpretation. The later
[determinant obstruction](../../normseven_arithmetic_f1_20261026/f1/RESULT.md)
separates root-fibre recovery, an ambient Laurent lift and extension over an
affine boundary. None licenses forgetting additive equations or a marked
pairing merely because a numerical specialization is simple.

There is a genuine but limited F1 precedent here. Lorscheid's
[*Tits–Weyl models of algebraic groups*](https://arxiv.org/abs/1201.1324)
constructs a Weyl extension for models including symplectic groups.
His [blueprint foundations](https://arxiv.org/abs/1103.1745) incorporate
additive relations rather than treating every arithmetic object as a bare
monoid. We use the Weyl-group/counting shadow below; we do not construct an
F1 model of an abelian variety, its Selmer groups or Sha.

For the arithmetic boundary, [Poonen–Stoll](https://math.mit.edu/~poonen/papers/sha.pdf)
distinguishes alternating from merely antisymmetric Cassels–Tate pairings.
The square-or-twice-a-square conclusion, under its hypotheses, is not a
claim that every such pairing is symplectic. The elementary hostile in
section6 makes the distinction visible already at group order four.

## 2. The two basic exact recurrences

Let V be a 2n-dimensional vector space over F_q with a nondegenerate
**alternating** bilinear form. Let L_n(q) count its Lagrangian n-dimensional
subspaces, and let F_n(q) count complete isotropic flags
0<L_1<...<L_n. Write [j]_q=1+q+...+q^(j-1), with empty products equal to one.

Every line is isotropic. Choosing its orthogonal quotient of dimension
2n-2 yields

    F_0(q)=1,  F_n(q)=[2n]_q F_(n-1)(q),
    F_n(q)=product_(i=1)^n [2i]_q.                         (1)

Counting incidences between a Lagrangian and a line in it gives

    L_n(q)[n]_q=[2n]_q L_(n-1)(q),
    L_0(q)=1,  L_n(q)=(q^n+1)L_(n-1)(q),
    L_n(q)=product_(i=1)^n (q^i+1).                       (2)

The number I_(n,k)(q) of isotropic k-dimensional subspaces is

    I_(n,k)(q)= GaussianBinomial(n,k;q)
                       product_(i=n-k+1)^n(q^i+1).       (3)

For a direct proof, count ordered isotropic independent k-tuples. After
j vectors the next lies in their common orthogonal space but not their
span: q^(2n-j)-q^j choices. Dividing by the number of ordered bases of a
k-space gives product_(j=0)^(k-1)(q^(2(n-j))-1)/(q^(k-j)-1), equal to (3).
These arguments work in characteristic two because the form is explicitly
alternating, not just antisymmetric.

At q=1 the polynomial specializations are

    F_n(1)=2^n n!,     L_n(1)=2^n,
    I_(n,k)(1)=2^k binom(n,k).                           (4)

A fixed symplectic basis divides its 2n coordinate vectors into n paired
lines. A coordinate Lagrangian chooses one line from each pair, giving
2^n choices. Ordering those choices gives 2^n n! coordinate isotropic flags.
The signed permutation group C2 wr S_n acts simply transitively on these
ordered choices; it is the Weyl group of type C_n. This describes the thin
coordinate object. It is not a claim that literal matrices over a one-element
field form Sp_(2n).

For comparison, choosing an actual symplectic basis gives

    |Sp_(2n)(F_q)|=q^(n²) product_(i=1)^n(q^(2i)-1).       (5)

The basis recurrence first chooses a nonzero vector and then a partner
with pairing one, giving (q^(2n)-1)q^(2n-1) times the lower-dimensional
count. The raw polynomial in (5) is zero at q=1 for n>0. It is the
normalized value after dividing by (q-1)^n, or the flag polynomial (1),
that gives 2^n n!. Forgetting this normalization confuses group and flag
counts.

## 3. Direct sums contain mixed Lagrangians

Let U,W have symplectic dimensions 2m,2n. A Lagrangian L in U orthogonal-sum W
need not be the direct sum of a Lagrangian in U and one in W. The following
is an exact closed composition formula:

    L_(m+n)(q)=sum_(r=0)^min(m,n)
                I_(m,m-r)(q) I_(n,n-r)(q) |Sp_(2r)(F_q)|. (6)

Here r is a **coupling rank**, with a direct geometric meaning. Set
A=L intersection U and B=L intersection W, of dimensions a,b. Projection of
L to U lies in A^perp. Its dimension is m+n-b, so n-b<=m-a. The symmetric
inequality forces equality. Hence m-a=n-b=r.

Pass to A^perp/A and B^perp/B, each symplectic of dimension2r. The quotient
L/(A+B) intersects each summand trivially and projects isomorphically onto
both. It is the graph of an anti-symplectic isomorphism f, characterized by

    <f(x),f(y)>_W = -<x,y>_U.

Conversely such a graph lifts to a unique L with those intersections. There
are I_(m,m-r) choices for A and I_(n,n-r) for B. Once one anti-symplectic
isomorphism is chosen, all others form a torsor under Sp_(2r)(F_q), whose
cardinality is (5). This proves (6), including characteristic two where
minus equals plus. It is a counting torsor, not a canonical identification
between the two quotient spaces.

The r=0 term is exactly the separable product L_m(q)L_n(q). No normalization,
unmarked source averaging or orbit division is hidden in (6).

## 4. The F1 limit erases a measured hierarchy of coupling

Set t=q-1 and let C_(m,n,r)(q) denote the rth summand of (6). The isotropic
counts have nonnegative coefficients in q and hence in t. Formula (5) has
exactly r factors divisible by t. Therefore C_(m,n,r)(1+t) has nonnegative
integer coefficients and **exact vanishing order r**, with

    [t^r] C_(m,n,r)(1+t)
       = 2^(m+n-r) r! binom(m,r) binom(n,r).             (7)

Indeed I_(m,m-r)(1)=2^(m-r)binom(m,r), similarly for n, and the leading
coefficient of Sp_(2r)(1+t) is product_(i=1)^r2i=2^r r!.

At t=0 every genuinely coupled sector disappears, leaving

    L_(m+n)(1)=L_m(1)L_n(1)=2^(m+n).                     (8)

The first discarded term is exactly

    L_(m+n)(1+t)-L_m(1+t)L_n(1+t)
       = mn 2^(m+n-1)t + O(t²),                         (9)

for m,n>=1. A jet through t^s cannot observe the contribution of any sector
r>s. This is a filtration of these explicit counting polynomials, not a
lower bound on every possible representation of the finite-field geometry.

The smallest hostile is m=n=1:

    L_2(q)=(q+1)(q²+1),
    L_1(q)^2=(q+1)^2,
    L_2(q)-L_1(q)^2=q(q²-1).                            (10)

At q=2 there are15 Lagrangians, only9 separable; six mixed graphs are lost by
the product model. At q=1 both counts are four. Thus exact equality in the
thin limit does not certify exact compositional factorization over fields.
The positive t coefficients here contrast with the old seven-line negative
coefficient obstruction: one has a positive but lossy counting limit, while
the other blocks that particular positive torus decomposition altogether.
Neither conclusion asserts a global torification for the relevant Sha data.

## 5. What this does and does not mean for doubling and squaring

Within this model, adjoining one hyperbolic plane changes dimension
2n ->2n+2 and multiplies vector-space cardinality by q². Its thin Lagrangian
count doubles:2^n ->2^(n+1). Taking a second copy changes dimension2n ->4n
and squares vector cardinality and the thin coordinate-Lagrangian count.
These are exact typed analogues of additive-size and doubling-size recursions.
At finite q, however, the Lagrangian count does not simply square, because
(6) includes mixed sectors. The coupling map is precisely what the thin
product forgets.

This is a comparison of operations, not an identification with a particular
tournament recursion. Such an identification would require the tournament
construction to carry an actual nondegenerate alternating form, specify its
integral lattice or coefficient field, and intertwine the relevant direct-sum
operation. A match between x+2 and x*2 alone supplies none of those maps.

Likewise Sha is not supplied with a variable q whose value one can be
substituted into (2). Even a finite p-primary paired group need not be an
F_p vector space. The elementary symplectic layer above is one controlled
model of pairing data, not a universal description of Sha.

## 6. Square cardinality does not identify the marked pairing

On G=(Z/2)^2, take Q/Z-valued bilinear pairings represented by half of

    H=[0 1]             I=[1 0].
      [1 0]               [0 1]

Both are nondegenerate and antisymmetric, since1/2=-1/2 modulo Z. Their
underlying groups have the same square order four. H is alternating and all
four elements have self-pairing zero. I is not alternating; only0 and(1,1)
have self-pairing zero. Its characteristic element c, defined by
b(x,x)=b(x,c), is(1,1), and b(c,c)=0. Thus even a zero characteristic
self-pairing does not imply that the original marked pairing is alternating.

This is an abstract finite-pairing hostile, not a claim that both are
realized by an abelian variety with specified polarization. It shows exactly
why a square-cardinality or matching-only F1 shadow cannot decide the
2-primary marked pairing. The separate pairing lane gives the arithmetic
hypotheses; this note does not replace them.

## Verification and loss ledger

`check_symplectic_skeleton.py` uses only the Python standard library and
explicit exceptions. It verifies (6) as polynomial identities for0<=m,n<=6,
all sector vanishing orders and leading coefficients, and literal finite-field
isotropic-subspace counts through dimension6 over F2,F3 and dimension4 over
F5. RREF subspaces and alternating products are enumerated directly.
Intersection dimensions independently classify mixed sectors. Both pairings
in section6 are checked element by element.

Source: finite alternating spaces with specified summands and their additive
forms. Target: coordinate pair choices / polynomial point counts / thin
q=1 specialization. Preserved: pair ordering and the displayed recurrence
identities. Lost at q=1: every positive-rank mixed graph sector and all
q-dependent additive data. Required restoration: coupling-rank sectors,
quotient spaces, anti-symplectic gluing maps and, for Sha comparisons, the
actual 2-primary characteristic element and local-global arithmetic.
