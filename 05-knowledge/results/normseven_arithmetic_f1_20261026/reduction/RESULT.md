# Reduction at three: lost freeness and the exact missing-sheet law

2026-10-26. **PROVED** formulas R1-R3 below; **FINITE-EXACT** the finite
Fano census and executable controls. All point counts mean reduced rational
points, not scheme lengths. Characteristic two is excluded.

## Inheritance and scope

[BOUNDARY B1-B3 and section 5](../../normseven_boundary_20261025/BOUNDARY.md)
already prove normalization, conductor, and the missing 2^r sheets when
the rank-r inertia injects into both quotient decks. The characteristic-three
extra intersections are already in [TWIST_AVERAGE](../../normseven_f1_20261025/TWIST_AVERAGE.md).
[Session63 twist mass](../../session63_integral_observers_root_arrangements_20261025/twists/RESULT.md)
already proves general twist mass, fixed-extension parity and open-cover
moments. We combine these mechanisms to compute precisely where the
inertia-injectivity hypothesis fails, and how the normalization discrepancy
changes. No new assertion about Campedelli classification in positive
characteristic is made.

## R1. Arbitrary inertia, without assuming free quotient kernels

Let G=F_2^n=H direct_sum K and let I be an inertia subspace. Put

    r=dim I, a=dim(I intersect H), b=dim(I intersect K), d=r-a-b.

Then the natural map

    G/I -> G/(H+I) x G/(K+I)

is surjective and has exactly 2^d elements in every fibre. The source has
2^(n-r) elements and the target 2^(n-r-d).

**Proof.** (H+I)+(K+I)=G, so the map is surjective. The target dimension is
(dim K-r+a)+(dim H-r+b)=n-2r+a+b=n-r-d. Its kernel therefore has dimension
d. Equivalently it is ((H+I) intersect (K+I))/I. Nonnegativity follows from
the direct sum (I intersect H) direct_sum (I intersect K) inside I. QED.

This is the reduced geometric-fibre map of a constant elementary abelian
cover and its two coarse quotients, whenever its geometric fibre is G/I.
The old 2^r formula is exactly the a=b=0 case. In the nonfree case the
quotient-stack identity [X/H] x_[X/G] [X/K] = X remains valid as a
2-fibre product. Its base is [X/G], not the coarse plane; the quotient
stacks cannot be silently replaced by the coarse quotient schemes.

## R2. The two chiral norm-seven pairs at the exceptional prime

Keep the seven lines and their fixed branch order

    (E,E+x,E-x,E+z,E-z,E+x+z,E-x-z)

and the four inherited binary label rows

    U2=(1,2,4,5,6,3,7), D2=(1,2,4,6,3,7,5),
    U4=(1,2,4,5,7,3,6), D4=(1,2,4,7,3,6,5).

Each row maps G=F_2^7/<ones> to F_2^3; H is its kernel. The rows and
complementarity are combinatorial and do not change on reducing the base.
The actual inertia subgroups do change: in characteristic three the triples
136 and 245 join 012,034,056. A triple has a one-dimensional intersection
with a kernel exactly when the xor of its three labels is zero.

| Pair | Triple | a for U | b for D | d | X fibre | Coarse product fibre | Sheets per product point |
|---|---|---:|---:|---:|---:|---:|---:|
| U2,D2 | 136 | 1 | 0 | 2 | 8 | 2 | 4 |
| U2,D2 | 245 | 0 | 1 | 2 | 8 | 2 | 4 |
| U4,D4 | 136 | 0 | 1 | 2 | 8 | 2 | 4 |
| U4,D4 | 245 | 1 | 0 | 2 | 8 | 2 | 4 |

For example U2(136)=2 xor 5 xor 7=0. Each displayed triple has three
distinct nonzero labels, so its rank is two when the xor vanishes and
three otherwise. At the original triples a=b=0, so d=3 and the old eight
missing sheets remain. At double and single points d=r. Thus all four
named kernels cease to act freely at characteristic three, each at one
new triple. The norm-seven arithmetic prime 7 is *not* a bad-incidence
prime for this arrangement: its triples and fibre law are the same as
over Q. The roles of 3 and 7 come from different maps.

**Normalization survives.** In each chart the radicands remain distinct
linear primes in every odd characteristic. The root algebra is normal:
single/double crossings are smooth after etale unit extensions, while a
triple has an ordinary quadratic node (three lines with pairwise independent
normals, characteristic not two). The character-basis proof of BOUNDARY B1
uses prime divisibility and works verbatim with 2 invertible. Consequently
the coarse product is still integral with normalization X and the same
conductor ideal O_X(-7H). Freeness and rational point lifting are separate
properties; their failure does not invalidate normalization.

**FINITE-EXACT census.** The 5,040 label permutations give 30 embedded
Fano kernels and 120 unordered complementary pairs. Against the original
three triples there are 16 free kernels and 32 both-free pairs; against
all five characteristic-three triples these become 12 and 18. This is a
finite combinatorial census, not a count of isomorphism classes of surfaces.

## R3. Equal twist mass hides a measurable normalization discrepancy

Fix a rational base point over F_q, q odd, with a geometrically transitive
fibre G/I. Take all 2^n constant G-twists over this very field. Frobenius
acts on G/I by translation, so twisting makes its translation uniform on
G. Write X_g and Y_g for the rational point counts in this fibre of the
root cover and the coarse product, respectively. Put

    F=2^(n-r), M=2^(n-r-d).

Their joint distribution is exactly

| (X_g,Y_g) | Number of twists |
|---|---:|
| (F,M) | 2^r |
| (0,M) | 2^(r+d)-2^r |
| (0,0) | 2^n-2^(r+d) |

**Proof.** The root fibre is rational precisely when the translated
Frobenius is in I; then all F points are rational. Both quotient fibres
are rational precisely when it lies in J=(H+I) intersect (K+I), of dimension
r+d; then all M product points are rational. Since I is contained in J,
the displayed partition follows. A different initial Frobenius only
translates the partition. QED.

In particular

    E[X_g]=E[Y_g]=1,
    Prob(X_g>0 | Y_g>0)=2^(-d),
    E[|X_g-Y_g|]=2(1-2^(-d)),
    E[(X_g-Y_g)^2]=2^(n-r)-2^(n-r-d).

At an original triple the distribution is (8,1) eight times and (0,1)
56 times. At a new characteristic-three triple it is (8,2) eight times,
(0,2) 24 times and (0,0) 32 times. Thus the rational lift rate among
successful coarse-product twists rises from 1/8 to 1/4; this describes
fresh twists at the specified source, not a prescribed twist over Q.

For either named pair, sum these pointwise moments over all base points.
The inherited stratum counts (r=1,2,3) are (7q-26,12,3) in characteristic
greater than three and (7q-20,6,5) in characteristic three. The last five
triples split into three original and two new ones. Hence

    sum_P E[|X_g(P)-Y_g(P)|] = 7q-11/4                  (all odd q),
    sum_P E[(X_g(P)-Y_g(P))^2] = 112q-251             (char >3),
                                  = 112q-215         (char =3).

The first statistic happens to hide the exceptional prime, whereas the
second detects it. These are sums of *pointwise* discrepancies. They are
not E[|#X_g-#Y_g|] or the variance of total counts: source correlations
would have to be retained for those quantities. The shared means do not
supply a section of normalization or a torification.

## Reproduction and connection contract

Run from the repository root:

```sh
python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/reduction/check_reduction.py > /tmp/reduction.normal.json
python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/reduction/check_reduction.py > /tmp/reduction.optimized.json
```

Compare with the committed receipts. The checker enumerates F_p projective
points at p=3,5,7,11,13, all Fano markings, all complementary pairs, and
literal coset actions for all 64 translations at every actual stratum in
characteristics three and five. Its orbit implementation does not use the
rank formulas to construct the counts. Normal and optimized runs must agree.

Source: marked root cover, inertia and complementary kernels. Target:
coarse quotient-product fibre and twist distribution. Map: the two coset
projections. Preserved: common base point, exact rational counts and twist
labels. Lost: the d-dimensional comparison torsor. Necessary sidecar:
I intersect H, I intersect K and the actual Frobenius translation, not
just the number of vanishing roots. Decisive hostile: triple 136 in
characteristic three. Positive control: all characteristic-seven strata.
