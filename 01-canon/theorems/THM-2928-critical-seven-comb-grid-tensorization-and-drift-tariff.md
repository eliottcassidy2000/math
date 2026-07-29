---
id: THM-2928
title: Critical seven-comb grid tensorization and drift tariff
status: >
  PROVED.  On a nonempty 1/L-grid carrier, aligned danger combs factor
  through the normalized multiplier circle.  Seven aligned combs leave at
  least 15/154 of the carrier.  After k aligned combs, every pointwise
  (7-k)-comb drift cover pays a strict reciprocal tariff controlled by the
  normalized safe mass and tooth-component count.  For six aligned combs,
  the remaining speed obeys the 39/61 cone and an exact finite
  phase-address intersection.  On a literal six-body carrier, the two
  six-comb safe floors collapse that intersection to finitely many clocks,
  all of which fail; hence the fully aligned and one-drift branches are
  empty.  This does not close the branches with at least two drifts or
  prove LRC(14).
source: root-2026-07-29 with independent hostile audits by seven-wall-tensor-audit and critical-residue-tree
depends_on:
  - THM-594-pair-overlap-law-mirsky-newman-floor
  - THM-1094-exact-two-comb-component-theorem
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-1234-sharp-five-comb-compatibility-floor
  - THM-2182-endpoint-grid-product-and-tail-overlap-sidecar
related:
  - THM-2184-two-scale-tail-continuation-profile
  - THM-1132-sharp-measure-horn-constant-dissolves-r6-wall
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-2893-complement-cap-finite-core-flag-lemma
  - HYP-7678
  - HYP-7870
script: 04-computation/lrc14_critical_one_drift_clock_thm2928.py
output: 05-knowledge/results/lrc14_critical_one_drift_clock_thm2928.out
---

# THM-2928 -- critical seven-comb grid tensorization and drift tariff

## 1. Grid carrier and the inherited drift chart

Write

```text
T=R/Z,                 D_w={t in T: ||wt||<1/14}.       (1)
```

Fix a positive integer `L`, a nonempty set of cell addresses
`J subset Z/LZ`, and the carrier

```text
C_J=union_(j in J) (j/L,(j+1)/L),       h=mu(C_J)=|J|/L. (2)
```

Grid endpoints never affect the measure statements.  Once at least one
aligned comb has been removed, they do not affect the pointwise statements
either: every grid endpoint belongs to every `D_(La)`.

For an arbitrary positive speed, write

```text
w=La+b,                 0<=b<L.
```

On cell `j`, put `t=(j+u)/L`.  Then

```text
||wt||=||(a+b/L)u+bj/L||.                               (3)
```

This is the `Q=L` specialization of the exact midpoint-cell coordinate in
[THM-2184, two-scale tail continuation profile][thm2184].  It identifies the
mixed-residue datum precisely:

```text
slope detuning:       b/L;
cell phase word:      bj mod L.                         (4)
```

When `b=0`, both disappear.  The aligned locus is therefore a genuine fixed
locus of the cell transport, not a scalar approximation to a general tail.

[thm2184]: THM-2184-two-scale-tail-continuation-profile.md

## 2. Boolean endpoint-grid tensorization

Let `A` be any finite set of positive integers.  For every bounded function
`F` of the Boolean danger word,

```text
integral_(C_J) F((1_(D_(La))(t))_(a in A)) dt
 =h integral_T F((1_(D_a)(u))_(a in A)) du.             (5)
```

Indeed, on each selected cell,

```text
La(j+u)/L=aj+au=au mod 1,
```

and `dt=du/L`.  Every selected cell therefore carries the same complete
Boolean word, and summing its `1/L` copy over `|J|` cells proves `(5)`.

For a product of safe indicators, `(5)` is exactly the endpoint-grid product
mechanism of
[THM-2182, endpoint-grid product][thm2182].  Equation `(5)` records its
Boolean extension; it is not a claim that the base product mechanism is new.

[thm2182]: THM-2182-endpoint-grid-product-and-tail-overlap-sidecar.md

Put

```text
R_A=T minus union_(a in A)D_a,       u_A=mu(R_A).        (6)
```

Then

```text
mu(C_J minus union_(a in A)D_(La))=h u_A,                (7)
```

and, inside every selected grid cell, the residual is a copy of `R_A`
scaled by `1/L`.  This is the exact toothpick self-similarity used below.

## 3. Tooth count and its scale law

Suppose `A={a_1,...,a_k}` has distinct entries and `1<=k<=7`.  Let

```text
N_A=sum_i a_i
```

and let `nu_A` be the number of **all** connected components of `R_A`,
including possible singleton components.  In every case used below `u_A>0`:
for `k<=6` this follows from the union bound, and for `k=7` it follows from
Section 4.

Each `D_(a_i)` has `a_i` open teeth.  The `k` teeth containing zero all
overlap, so those `k` components contribute only one component to their
union.  Further mergers can only lower the count.  Hence

```text
1<=nu_A<=N_A-k+1.                                      (8)
```

The common origin tooth also means every component of `R_A` has a
nonwrapping representative inside `(0,1)`.

For every positive integer `d`,

```text
R_(dA)={t:dt mod 1 belongs to R_A},
u_(dA)=u_A,
nu_(dA)=d nu_A.                                        (9)
```

The first identity is immediate from `D_(da)={t:dt in D_a}`; the degree-`d`
circle covering gives the other two.  Thus the tooth functional

```text
H(A)=7k u_A/(6nu_A)                                    (10)
```

has the exact homogeneity

```text
H(dA)=H(A)/d.                                          (11)
```

For one multiplier, `u_{a}=6/7`, `nu_{a}=a`, and therefore

```text
H({a})=1/a.                                            (12)
```

This is the recursive functional form behind the ordinary slow-gap harmonic
invoice.

## 4. The critical seven-comb fixed locus is quantitatively safe

For every seven distinct positive integers `A`, the global strict-spectrum
Hunter theorem
[THM-1221][thm1221] gives

```text
u_A>=15/154.                                           (13)
```

Combining `(7)` and `(13)`,

```text
mu(C_J minus union_(a in A)D_(La))>=15h/154>0.          (14)
```

Thus seven distinct aligned combs cannot cover any nonempty grid carrier.
More locally, every selected cell contains an aligned-safe interval of
length at least

```text
15/(154L nu_A)>=15/(154L(N_A-6)).                      (15)
```

The independent divisor-minimal Fourier argument of
[THM-594][thm594] gives the weaker but purely spectral floor

```text
u_A>=2sin^2(pi/7)/(7pi^2).                             (16)
```

The stronger value `(13)`, not `(16)`, is used in `(14)`--`(15)`.

[thm1221]: THM-1221-seven-wall-strict-spectrum-hunter-floor.md
[thm594]: THM-594-pair-overlap-law-mirsky-newman-floor.md

There is an important equality boundary.  For every positive integer `a`,
`(5)` gives

```text
mu(C_J intersect D_(La))=h/7.                          (17)
```

Consequently the sum of seven aligned singleton coverages is exactly `h`.
The strict scalar first-apex gate at the `p=7` wall really does stop at
equality, as anticipated by
[THM-2893][thm2893].  What repairs the aligned branch is the joint Boolean
law and its positive safe defect `(13)`, not a stricter singleton estimate.

[thm2893]: THM-2893-complement-cap-finite-core-flag-lemma.md

## 5. The strict reciprocal tariff after aligned deletions

Let `1<=k<=6`, let `1<=r<=6`, and suppose positive speeds

```text
z_1,...,z_r
```

cover the aligned residual pointwise:

```text
C_J minus union_(a in A)D_(La)
 subset union_(q=1)^r D_(z_q).                         (18)
```

Then the general tariff is

```text
sum_(q=1)^r 1/z_q
 >7(7-r)u_A/(6Lnu_A).                                 (19)
```

At the critical seven-wall cardinality `r=7-k`, this becomes

```text
sum_(q=1)^(7-k) 1/z_q
 >H(A)/L
 =7k u_A/(6Lnu_A)                                     (20)
 >=k(7-k)/(6Lnu_A)
 >=k(7-k)/(6L(N_A-k+1)).                              (21)
```

### Proof

Choose a largest positive-length component `K` of `R_A`.  It is a compact
interval with

```text
|K|>=u_A/nu_A.                                        (22a)
```

In any selected grid cell, its copy `I` is a compact interval of length

```text
ell=|K|/L>=u_A/(Lnu_A).                               (22b)
```

The sharp periodic discrepancy estimate in
[THM-1094, equation (10)][thm1094] says, for every interval `I`,

```text
mu(I intersect D_z)<=ell/7+6/(49z).                    (23)
```

The union of the `r` danger combs is open and contains the compact interval
`I`.  Hence a slightly longer interval `I'` of length `ell'>ell` is still
covered.  Applying `(23)` to `I'` gives

```text
ell'<=r ell'/7+(6/49)sum_q 1/z_q.
```

Rearranging and using `(22b)` proves the strict inequality
`(19)`.  For `r=7-k`, the union bound
`u_A>=1-k/7=(7-k)/7` and `(8)` give `(20)`--`(21)`.  QED.

For an almost-everywhere cover, the compact enlargement is unavailable and
the same proof gives the corresponding non-strict inequality.

[thm1094]: THM-1094-exact-two-comb-component-theorem.md

At the critical cardinality `r=7-k`, at least one drift speed satisfies

```text
min_q z_q
 <6rLnu_A/(7ku_A)
 <=6Lnu_A/k
 <=6L(N_A-k+1)/k.                                     (24)
```

For `k=1`, equations `(12)` and `(19)` recover exactly the strict slow-gap
pressure

```text
(La) sum_(q=1)^6 1/z_q>1,                              (25)
```

the cellwise starting rung of
[THM-1176, slow-gap harmonic crowding][thm1176].  For larger `k`, `(19)` is
its recursive residual version; the state is the pair `(u_A,nu_A)`, not
mass alone.

[thm1176]: THM-1176-seven-wall-slow-gap-harmonic-crowding.md

### Current all-cardinality safe floors

The crude union floor in `(21)` is not the strongest proved input.  For
distinct multiplier sets, the current uniform safe floors are

```text
k             1       2        3         4          5          6
u_A >=       6/7    66/91    55/91    558/1183   478/1365    61/273. (25a)
```

The `k=2,3,4` entries follow from the pair/triple and spanning-tree
consequences of [THM-1166][thm1166].  For `k=5`, the
`R_5>=44/273` floor of THM-1234 and the uniform spanning-tree average give
`478/1365`; its six-comb quadratic consequence gives `61/273`.
Substituting `(25a)` into `(20)` strengthens the tariff coefficient
`7ku_A/6` to

```text
c_k = 1, 22/13, 55/26, 372/169, 239/117, 61/39,       (25b)

sum 1/z_q > c_k/(Lnu_A)
             >=c_k/(L(N_A-k+1)).                      (25c)
```

[thm1166]: THM-1166-seven-wall-fano-gcd-discrepancy.md

## 6. Six aligned combs and one drift: the `39/61` cone

Now take `k=6`, `r=1`, and write the aligned speeds as

```text
w_i=La_i.
```

The six-comb consequence of
[THM-1234][thm1234] gives

```text
u_A>=61/273.                                           (26)
```

Equations `(19)` and `(26)` imply

```text
z<Lnu_A/(7u_A)
 <=(39/61)Lnu_A
 <=(39/61)L(N_A-5)                                    (27)
 =(39/61)(sum_i w_i-5L).
```

In integral form,

```text
z<=ceil(39L(N_A-5)/61)-1.                              (28)
```

This is a strict additive cone for the one-drift branch.  It is only a
uniform corollary: using the exact `(u_A,nu_A)` in the first inequality of
`(27)` can be substantially stronger.

[thm1234]: THM-1234-sharp-five-comb-compatibility-floor.md

## 7. Exact one-drift phase and gcd sidecar

Continue with `k=6`.  List **all** connected components of `R_A`, including
any singleton components, as

```text
K_s=[m_s-rho_s,m_s+rho_s] subset (0,1),       rho_s>=0. (29)
```

Then the one drift comb covers the entire aligned residual pointwise if and
only if

```text
||z(j+m_s)/L||+z rho_s/L<1/14
             for every j in J and every s.             (30)
```

Indeed, the copy `(j+K_s)/L` is connected and compact.  It lies in `D_z`
exactly when it lies strictly inside one open `z`-tooth; comparing its
centre and half-length gives `(30)`.

Equation `(30)` has a finite address form.  On the circle `R/LZ`, define

```text
A_z=intersection_s
 {x: ||(x+z m_s)/L||<1/14-z rho_s/L}.                  (31)
```

An individual set in `(31)` is empty when its displayed radius is
nonpositive.  Otherwise it is an open arc.  The exact condition is

```text
zJ mod L subset A_z.                                   (32)
```

This intersection retains the relative positions of *all* residual teeth;
replacing it by one width is a lossy relaxation.

Let `rho_max` be the largest half-length of a positive component and let
`span_L` denote the length of the shortest circular arc containing a finite
residue set.  Since

```text
2rho_max>=u_A/nu_A,
```

condition `(32)` forces

```text
span_L(zJ mod L)
 <L/7-2zrho_max
 <=L/7-z u_A/nu_A.                                    (33)
```

In particular the right side must be positive.  If `g=gcd(z,L)`, the image
`zJ` lies in the subgroup `gZ/LZ` and multiplication by `z` has fibers of
size `g`.  Hence the exact and width-only cardinality screens are

```text
|J|<=g # (A_z intersect (gZ/LZ)),                      (34)

|J|<=g ceil((L/7-z u_A/nu_A)/g).                       (35)
```

Thus the metric cone `(27)` does not leave an unstructured finite box.  It
leaves a simultaneous rational tooth-address intersection with an explicit
gcd compression.

There is a useful density form of this compression.  Put

```text
g=gcd(z,L),       d=L/g,       c=z/g,       P=J mod d. (35a)
```

Then `gcd(c,d)=1`, every residue modulo `d` has exactly `g` lifts modulo
`L`, and `(30)` is equivalent to

```text
P subset B_(c,d)(A)
 ={r in Z/dZ:
   ||c(r+m_s)/d||+c rho_s/d<1/14 for every s}.         (35b)
```

The largest-component argument proving `(35)` now gives

```text
h<=|P|/d
 <=ceil(d/7-c u_A/nu_A)/d.                            (35c)
```

In particular, since `ceil(x)<x+1`,

```text
g>L(h-1/7).                                            (35d)
```

Thus every carrier with `h>1/7` forces a bounded quotient clock.  If
`h>=61/273`, then

```text
g>22L/273,                 d=L/g<273/22,              (35e)
d<=12.
```

For fixed `A`, positivity in `(35c)` also gives

```text
c/d=z/L<nu_A/(7u_A),                                  (35f)
```

so only finitely many coprime phase clocks `(c,d)` remain.

## 8. Literal six-body one-drift closure

For a fixed finite speed set `F`, put

```text
G_F=T minus union_(v in F)D_v,
L_F=lcm(14v:v in F).                                   (36)
```

Every boundary of `G_F` has the form

```text
(14q+/-1)/(14v),
```

so, modulo its finitely many endpoints, `G_F` is a union of `1/L_F` cells.
Whenever `mu(G_F)>0`, the preceding theorem applies with `L=L_F` (or with
any smaller denominator that exactly resolves the same carrier).

Now suppose `F` consists of six distinct body speeds, six of the seven tail
speeds are the distinct aligned speeds `La_i`, and the last tail speed is
`z`.  Take `J` to be the exact set of open `1/L`-cells contained in `G_F`.
The finitely many grid endpoints are all in every aligned danger comb and
therefore do not occur in the residual cover problem.

Apply the six-comb quadratic consequence of THM-1234 twice, first to `F`
and then to `A`.  It gives

```text
h=mu(G_F)>=61/273,             u_A>=61/273.            (36a)
```

Assume for contradiction that `D_z` covers the aligned residual.  Equations
`(35c)` and `(35e)` give `d<=12`, while

```text
h<=ceil(d/7-c u_A/nu_A)/d
  <=ceil(d/7)/d.                                      (36b)
```

Among `1<=d<=12`, comparison with `61/273` leaves only

```text
d in {1,2,3,4,8}.                                     (36c)
```

The carrier address set is invariant under circle reflection:

```text
j in J  iff  L-1-j in J.                              (36d)
```

For `d=2,4`, equation `(36b)` permits only one occupied residue, but
reflection sends it to `-1-r mod d` and has no fixed residue.  Hence these
two clocks are impossible.  For `d=1,3`, the unique occupied residue is,
respectively, `0,1`.  For `d=8`, density forces exactly two occupied
residues, necessarily one of the reflected pairs

```text
{r,7-r},                  0<=r<=3.                    (36e)
```

It remains to test the phase windows, not the original speed box.  For an
occupied residue `r`, `(30)` says

```text
R_A subset E_(c,d,r)
 ={u in [0,1]: ||c(u+r)/d||<1/14}.                    (36f)
```

After the change of variable `y=(u+r)/d`, THM-1094's exact interval
discrepancy gives

```text
u_A<=mu(E_(c,d,r))<=1/7+6d/(49c).                     (36g)
```

Together with `(36a)`, this forces

```text
c<=819d/539.                                           (36h)
```

Since `gcd(c,d)=1`, only the following exact endpoint table remains.  The
displayed `d=8` masses are the maxima over the four reflected pairs in
`(36e)`.

```text
d  occupied phase(s)  possible c          simultaneous window masses
1       {0}               1               1/7
3       {1}             1,2,4             0, 3/14, 3/28
8     {r,7-r}       1,3,5,7,9,11          1/7,1/21,1/35,1/49,2/63,2/77.
                                                               (36i)
```

The table follows by subdividing `[0,1]` at the rational endpoints

```text
u=(d/c)(q+/-1/14)-r.
```

Its largest entry is `3/14<61/273`, contradicting `(36a)` and `(36f)`.
Therefore:

> **Literal one-drift closure.**  Six aligned tail combs and one arbitrary
> tail comb cannot cover the safe carrier of six distinct body combs.

The fully aligned case has the stronger quantitative conclusion

```text
mu(G_F minus union_(i=1)^7 D_(La_i))
 >=(61/273)(15/154)=305/14014.                         (36j)
```

The exact rational table `(36i)` is reproduced by two independent interval
integrations in
`04-computation/lrc14_critical_one_drift_clock_thm2928.py`.  Ordinary and
`python3 -O` output must agree byte for byte with
`05-knowledge/results/lrc14_critical_one_drift_clock_thm2928.out`.
Their frozen SHA-256 hashes are

```text
script  0d0d4acd71b5987995c41b8a19b3c0de7ee1fde6a01d1b6c9d0e0cb4a32ddd6a
output  f46211d9f441f72753059c01a22cc6a5765d200b61bcd4541ccfbe9578708a3b.
```

## 9. Scope and information audit

The proved transport is

```text
source:       labelled carrier cells and danger-tooth incidences;
map:          t=(j+u)/L;
preserved:    the complete Boolean word when every residue b is zero;
destroyed off the fixed locus:
              slope detuning b/L and cell phase bj/L;
necessary mixed-residue sidecar:
              J, bJ mod L, component centres/radii, and gcd fibers. (37)
```

For an arbitrary grid carrier, equations `(27)`--`(35)` remain necessary
rather than sufficient and need not eliminate every one-drift row.  The
double use of the six-comb floor and reflection is what closes the literal
six-body one-drift branch.  The theorem does **not** close an arbitrary
mixed-residue seven-wall, the branches with two or more drift speeds, or
LRC(14).  Their next exact object is a multi-slope version of the finite
address clock `(35b)`.  A tournament on the seven labels forgets the metric
widths, cell phases, and gcd fibers in that clock and is therefore not an
equivalent quotient.
