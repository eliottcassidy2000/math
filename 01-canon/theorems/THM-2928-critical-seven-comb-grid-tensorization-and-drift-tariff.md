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
  empty.  A component-free phase-load tariff bounds one slope in every
  k>=2 mixed branch, and every fixed five-aligned/two-drift chart reduces
  to an exact finite pair-clock quotient with the sharper second-slope
  bound c_2/d_2<13 max A.  At the 2/25 deletion threshold, either that
  bound improves to c_2/d_2<(25/3) max A or all five aligned multipliers
  lie in an explicit body-uniform finite staircase.  An independent body-
  projection sieve kills 240,560 of 251,536 body/divisor rows.  Its fixed-
  phase address capacity reduces 3,066,274 denominator-pair occurrences to
  23,755.  An exact relaxed arithmetic-progression mask kills all 22,813
  survivors having a denominator at most 1,000; an arbitrary-class load
  relaxation kills 940 of the remaining 942; and a quotient-fiber
  transversal obstruction kills the last two diagonal pairs.  Hence the
  literal five-aligned/two-drift branch is empty uniformly.  The critical
  two-torus carrier integral is positive on every fixed seven-tail affine
  ray;
  uniformly, seven tails in one common canonical-ruler quotient block are
  safe once the quotient is at least 315,586.  Branches with three or more
  drifts remain open; this does not prove LRC(14).
source: root-2026-07-29 with independent hostile audits by seven-wall-tensor-audit and critical-residue-tree
depends_on:
  - THM-594-pair-overlap-law-mirsky-newman-floor
  - THM-1094-exact-two-comb-component-theorem
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-1234-sharp-five-comb-compatibility-floor
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-2182-endpoint-grid-product-and-tail-overlap-sidecar
  - THM-2184-two-scale-tail-continuation-profile
  - LRC(<=13)
related:
  - THM-1135-r6-harmonic-tail-finite-box
  - THM-2184-two-scale-tail-continuation-profile
  - THM-1132-sharp-measure-horn-constant-dissolves-r6-wall
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-2893-complement-cap-finite-core-flag-lemma
  - HYP-7678
  - HYP-7870
script: 04-computation/lrc14_critical_one_drift_clock_thm2928.py
output: 05-knowledge/results/lrc14_critical_one_drift_clock_thm2928.out
support_script: 04-computation/lrc14_two_drift_body_projection_support_thm2928.py
support_output: 05-knowledge/results/lrc14_two_drift_body_projection_support_thm2928.out
address_script: 04-computation/lrc14_two_drift_relaxed_address_mask_thm2928.py
address_output: 05-knowledge/results/lrc14_two_drift_relaxed_address_mask_thm2928.out
address_completion_script: 04-computation/lrc14_two_drift_arbitrary_class_fiber_thm2928.py
address_completion_output: 05-knowledge/results/lrc14_two_drift_arbitrary_class_fiber_thm2928.out
affine_profile_script: 04-computation/lrc14_j7_affine_profile_min_frequency_thm2928.py
affine_profile_output: 05-knowledge/results/lrc14_j7_affine_profile_min_frequency_thm2928.out
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

### The component-free phase-load tariff

The largest-tooth tariff `(19)` is strongest when `nu_A` is controlled.
There is a complementary estimate that forgets components but retains the
complete load on one grid cell.  For `j in J`, put

```text
E_q(j)={u in [0,1]: ||z_q(j+u)/L||<1/14}.             (25d)
```

Apply THM-1094's interval discrepancy to the physical cell
`[j/L,(j+1)/L]` and rescale by `L`.  It gives

```text
mu(E_q(j))<=1/7+6L/(49z_q).                           (25e)
```

The proper compact set `R_A` lies in the open union of the `E_q(j)`.
Consequently that union has measure strictly larger than `u_A`, and

```text
L sum_(q=1)^r 1/z_q
 >(49/6)(u_A-r/7).                                   (25f)
```

For an almost-everywhere cover the same inequality is non-strict.  At the
critical cardinality `r=7-k`, the safe floors `(25a)` give

```text
k                 1      2      3        4         5        6
b_k              0     7/78   7/26    119/338   308/585   77/117

L sum 1/z_q > b_k.                                    (25g)
```

For `k>=2`, some drift therefore satisfies

```text
min_q z_q/L < (7-k)/b_k.                              (25h)
```

In particular, five aligned combs and two drifts force

```text
min(z_1,z_2)/L<585/154.                               (25i)
```

This is a slope bound independent of `nu_A`; it does not control the
denominator of that rational slope uniformly as `L` varies.  For `k=1`,
the coefficient vanishes and the component-sensitive tariff `(19)` remains
the live estimate.

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
script  16e7a630eb4f721836dfc95c0adf97bf9a1a252518b119450b16c71e08305dcd
output  6eb2ad0d031a94d7e85db37a30cb48a2e0abb03e4caa080f05f870a23d449746.
```

## 9. Five aligned combs and two drifts are empty

Keep the literal six-body carrier `G_F`, take its canonical resolving
denominator `L=L_F` from `(36)`, and let `J` be its full selected-cell set.
Now let `A` contain five distinct aligned
multipliers and suppose two drift combs cover the residual.  Relabel them so
that `z_1<=z_2`.  Equation `(25i)` gives

```text
z_1/L<585/154.                                        (37a)
```

Write each drift uniquely as

```text
z_q=(L/d_q)c_q,       d_q=L/gcd(z_q,L),       gcd(c_q,d_q)=1. (37b)
```

Thus `d_q` divides `L`, and `(37a)` already leaves finitely many first
clocks:

```text
d_1 divides L,             c_1/d_1<585/154.           (37c)
```

For a residue `r in S_(d_1):=J mod d_1`, define

```text
E_(c_1,d_1,r)
 ={u in [0,1]: ||c_1(u+r)/d_1||<1/14},

B_r=R_A minus E_(c_1,d_1,r).                          (37d)
```

At least one `B_r` has positive measure.  One proof uses the cited LRC
through thirteen: the partial family

```text
F union LA union {z_1}
```

has twelve nonzero speeds, so it has a time at which all distances are at
least `1/13`.  The strict `1/13-1/14` margin supplies a positive safe
interval.  That interval cannot meet a grid endpoint, where every aligned
speed vanishes, and therefore normalizes into one positive part of some
`B_r`.  Independently, the exact address/reflection audit in the companion
script proves that one drift cannot cover `R_A` even almost everywhere:
the only surviving five-aligned clocks are `(d,c)=(3,1),(8,1),(8,3)`,
whose simultaneous phase masses are respectively `0`, at most `1/7`, and
at most `1/21`, all below `478/1365`.

There is a uniform second-clock bound.  Put

```text
a_*=max A,                 N_A=sum_(a in A)a.          (37e)
```

Because five distinct positive multipliers have `a_*>=5`, equation `(37a)`
gives `z_1<5L<=La_*`.  Also `v<=L/14` for every `v in F`.  Thus `La_*` is
the largest speed in the twelve-speed partial family.  The `1/13` witness
is therefore safe at level `1/14` on a closed physical circle arc `I` of
length

```text
|I|=1/(91La_*).                                       (37f)
```

Indeed, the distance-to-integer function for a speed `v` is
`v`-Lipschitz, so the closed radius-`1/(182La_*)` arc around the witness
stays safe at level `1/14` for every member of the partial family.  Under a
pointwise countercover, every point of `I`, including both endpoints, must
therefore lie in `D_(z_2)`.  The connected arc `I` lies in one component of
that open set, whose teeth have length `1/(7z_2)`.  Closed containment in an
open tooth is strict, and hence

```text
1/(91La_*)<1/(7z_2),

c_2/d_2=z_2/L<13a_*.                                  (37g)
```

There is also a useful, but weaker, BV sidecar.  On one normalized cell,
the first drift is the integer `c_1` comb on an interval of length `1/d_1`.
THM-1094's tooth count gives

```text
number of meeting teeth
 <=c_1/d_1+8/7
 <585/154+8/7
 =5327/1078<5.                                        (37h)
```

Thus at most four first-drift teeth cut `R_A`.  Since `(8)` gives
`nu_A<=N_A-4`, every positive `B_r` has at most `N_A` positive-length BV
components.  The arc `I` cannot meet a body-grid endpoint, because every
aligned comb is dangerous there.  It therefore lies in one selected cell
and normalizes to a subinterval of some `B_r` of length `1/(91a_*)`.
Writing the mass and positive-length BV component count of that `B_r` as
`mu_r` and `K_r`, we have `mu_r>=1/(91a_*)` and `K_r<=N_A`.  Applying
THM-2162 to `(r_2+B_r)/d_2` gives

```text
mu_r/d_2<=mu_r/(7d_2)+6K_r/(49c_2),
```

and recovers the older almost-everywhere estimate

```text
c_2/d_2<=13N_Aa_*,
```

but the connected closed-arc argument `(37f)`--`(37g)` is stronger for the
pointwise problem.  Hence both clocks lie in the explicit finite box

```text
d_1,d_2 divide L,
gcd(c_q,d_q)=1,
c_1/d_1<585/154,
c_2/d_2<13 max A.                                     (37i)
```

Row-specific exact values of `K_r/mu_r` retain additional BV information,
although the universal pointwise face `(37g)` is usually stronger.

On the direct exactly-six-in-window LRC(14) boundary,
`F subset {1,...,14}`.  In that case

```text
L_F divides L_14:=14 lcm(1,...,14)=5,045,040.          (37ia)
```

Thus all possible clock denominators, and the first numerator through
`(37c)`, already lie in one body-uniform finite set.  The remaining
unbounded scale is `max A`; its projective high slope
`c_2/(d_2 max A)` lies in the bounded interval `(0,13)` and remains
coupled to the five-multiplier shape.  After the fully aligned and
one-drift branches have been discharged, a genuinely two-drift row also
has `d_1,d_2>1`; a clock with `d_q=1` is aligned and belongs to an earlier
branch.

### The `2/25` deletion dichotomy and reciprocal staircase

Let

```text
P=F union LA union {z_1}
```

be the twelve-speed family obtained by deleting `z_2`, and write

```text
M(P)=max_t min_(v in P)||vt||.
```

There are two branches.

If `M(P)>=2/25`, the same Lipschitz argument as `(37f)` gives a closed arc
safe for `P` at level `1/14` of length at least

```text
2(2/25-1/14)/(La_*)=3/(175La_*).
```

It must lie strictly inside one open `z_2`-tooth.  Therefore

```text
z_2/L=c_2/d_2<(25/3)a_*.                              (37ib)
```

Suppose instead that `M(P)<2/25`.  Order

```text
a_1<a_2<a_3<a_4<a_5
```

and put

```text
p=4/25,          eta=p(1-p)=84/625,
C_0=585/154,
C_s=max(C_0,a_s)                 (1<=s<=4).            (37ic)
```

For `0<=s<=4`, let

```text
P_s=F union {z_1} union {La_1,...,La_s}.
```

This has `7+s` speeds and maximum at most `LC_s`.  Cited LRC through
thirteen supplies a point with clearance at least `1/(8+s)`, and
Lipschitz fattening supplies a closed `2/25`-safe arc of length at least

```text
ell_s=
 2(1/(8+s)-2/25)/(LC_s).                              (37id)
```

The remaining `r=5-s` aligned combs must cover that arc.  The exact
interval discrepancy at radius `2/25` is

```text
mu(I intersect D_w^(2/25))
 <=p|I|+eta/w.                                        (37ie)
```

Indeed, the centered period-one indicator of an interval of length `p`
has a primitive of range exactly `p(1-p)=eta`; scaling by `w` proves
`(37ie)`.  Summing it over the remaining aligned speeds and using
`(37id)` gives the reciprocal invoice

```text
sum_(i=s+1)^5 1/a_i >= k_s/C_s,

k_s=
 2(1/(8+s)-2/25)(1-(5-s)p)/eta.                       (37if)
```

The five exact tariffs are

```text
s             0        1        2        3        4
k_s         15/112    1/6     13/84    17/154    1/24
(5-s)/k_s  112/3      24      252/13   308/17     24. (37ig)
```

Since the remaining reciprocals are each at most `1/a_(s+1)`,
`(37if)` recursively yields

```text
a_1 <=       141,
a_2 <=     3,384,
a_3 <=    65,597,
a_4 <= 1,188,463,
a_5 <=28,523,112.                                     (37ih)
```

The harmonic inequalities `(37if)`, rather than only their rectangular
consequence `(37ih)`, are the sharper finite search object.  Together with
`L|5,045,040`, `(37a)`, and `(37g)`, they make the entire
`M(P)<2/25` subbranch body-uniformly finite.  No emptiness census of that
large box is claimed.  The unbounded branch is confined instead by the
strict projective cone `(37ib)`.

### The body projection and relaxed address-mask sieve

The quotient has a second exact projection that forgets the aligned shape
but retains the complete body address support.  Put

```text
Y_D(A)=union_(r in S_D)(r+R_A)/D.                    (37ii)
```

The copies in `(37ii)` lie in distinct `1/D` cells.  Consequently

```text
mu(Y_D(A))=(|S_D|/D)u_A.                             (37iii)
```

If `(37k)` holds, then

```text
Y_D(A) subset D_(a_1) union D_(a_2).
```

The two quotient speeds are distinct.  The sharp global pair floor from
THM-594/THM-1166 therefore gives

```text
mu(D_(a_1) union D_(a_2))
 =2/7-rho(a_1,a_2)
 <=2/7-1/91
 =25/91.
```

Together with `u_A>=478/1365`, this proves the support obstruction

```text
|S_D|/D<=375/478.                                    (37iv)
```

The exact all-body census ranges over

```text
F in C({1,...,14},6),             D divides L_F.
```

There are `251,536` such body/divisor rows.  Condition `(37iv)` eliminates
`240,560`, leaving `10,976`.  Every surviving divisor satisfies `D>=42`,
and each body has at most six surviving divisors.  The unique row at the
minimum divisor is

```text
F=(1,2,3,4,6,12),       L_F=168,       D=42,
|S_D|/D=32/42=16/21.                                 (37v)
```

There is a numerator-free denominator screen on every remaining row.  Fix
one `u in R_A` and define

```text
M_i(u)={r in Z/DZ: ||c_i(r+u)/d_i||<1/14}.
```

Multiplication by `c_i` permutes `Z/d_i Z`.  An open arc of length `1/7`
contains at most `ceil(d_i/7)` points of that equally spaced grid, and each
such residue has `D/d_i` lifts.  Hence

```text
|M_i(u)|<=(D/d_i)ceil(d_i/7),

|S_D|/D
 <=ceil(d_1/7)/d_1+ceil(d_2/7)/d_2.                 (37vi)
```

Among the `3,066,274` divisor pairs with

```text
d_1,d_2>1,                  lcm(d_1,d_2)=D,
```

on the support-hard rows, `(37vi)` leaves exactly `23,755` pair
occurrences on `6,292` rows.

Reflection supplies a useful visible subcase.  The involution

```text
sigma(r)=D-1-r
```

preserves `S_D`.  When `d_1=2`, it swaps the two parity classes, so each
contains exactly half of `S_D`; the first clock can meet at most one class.
The other denominator `d_2` must therefore satisfy

```text
|S_D|<=2(D/d_2)ceil(d_2/7).                          (37vii)
```

This kills `6,754` of the `6,756` cardinality survivors with denominator
two.  The two remaining parity rows also fail the stronger address test
below.

For a general denominator `d`, the actual fixed-`u` mask is the lift of at
most `ceil(d/7)` consecutive residues after multiplication by a unit modulo
`d`.  Enlarge it to exactly that many residues and allow the two clocks to
choose their phases independently.  This is an upper relaxation of the
literal mask problem.  For every pair with

```text
min(d_1,d_2)<=1000,
```

the exact verifier maximizes the first enlarged arithmetic-progression load
on `S_D` and grants the second its full cardinality capacity.  That already
kills `22,811` of the `22,813` pairs.  On each of the two survivors, for
both possible first parity masks, project the residual to
`P subset Z/d_2 Z`.  Containment in a `k_2=ceil(d_2/7)` term cyclic
arithmetic progression would necessarily imply

```text
|P|<=k_2,                    |P-P|<=2k_2-1.           (37viii)
```

Both rows violate the difference-set inequality.  Thus no denominator pair
with a denominator at most `1,000` survives even this relaxed necessary
test.  The entire remaining denominator ledger consists of `942`
occurrences, both denominators greater than `1,000`, on exactly two rows:

```text
F                         D=L_F     |S_D|/D       pairs
(1,4,5,7,9,11)           194040     2308/8085      371
(1,5,7,8,9,11)           388080     3029/10780     571. (37ix)
```

There is also an actual-slope consequence near the top of the support
window.  If

```text
x=|S_D|/D>2535/3346,
delta=x(478/1365)-13/49>0,
```

then a cover would force

```text
rho(a_1,a_2)<=1/49-delta.
```

Write `a_i=g alpha_i` with coprime `alpha_1<alpha_2`.  LEM-043 gives

```text
rho(a_1,a_2)>=1/49-1/(7alpha_2),
alpha_2<=1/(7delta).                                 (37ixa)
```

Moreover `g<=a_1<(585/154)D`, so `(37ixa)` bounds the actual second slope
`a_2=g alpha_2`, not only its projective ratio.  The exact census places
`1,150` quotient rows in this body-uniform finite-slope sector; its largest
displayed integer cap is

```text
a_2<=26,927,449,547.
```

The separate identity

```text
gcd(g,D)
 =gcd(D/d_1,D/d_2)
 =D/lcm(d_1,d_2)
 =1
```

is a useful clock sidecar but is not itself the size bound.

For this symmetric denominator ledger, reorder the two clocks so that
`d_1<=d_2`; this relabelling is independent of the earlier ordering by
speed.  The remaining `942` pairs admit a still coarser, and therefore
cheaper, completion.  Put

```text
k_d=ceil(d/7),
lambda_d(s)=#{r in S_D:r=s mod d},
Lambda_d=sum of the k_d largest values of lambda_d.  (37ixb)
```

The first actual mask is contained in the lift of some `k_(d_1)` residue
classes, even if all arithmetic-progression structure is forgotten.
Therefore it meets at most `Lambda_(d_1)` points of `S_D`.  Granting the
second mask its full ambient capacity gives the necessary inequality

```text
|S_D|<=Lambda_(d_1)+(D/d_2)k_(d_2).                  (37ixc)
```

The exact load ledger violates `(37ixc)` for `940` of the `942` pairs.  Its
only survivors are the two diagonal pairs

```text
F                         D          (d_1,d_2)       two capacities
(1,4,5,7,9,11)           194040     (D,D)           27720+27720
(1,5,7,8,9,11)           388080     (D,D)           55440+55440. (37ixd)
```

They fail for a structural reason hidden by cardinality.  Write `D=7k`.
Every enlarged diagonal mask has the form

```text
B(s,h)={s+jh mod D:0<=j<k},             gcd(h,D)=1.  (37ixe)
```

Because `h` is also a unit modulo `k`, projection of `(37ixe)` to
`Z/kZ` is a bijection.  Thus one mask contains exactly one point in each
modulo-`k` fibre, and two masks contain at most two.  The first row in
`(37ixd)` has four support points

```text
29701, 57421, 112861, 168301 = 1981 mod 27720,
```

while the second has three

```text
59401, 114841, 225721 = 3961 mod 55440.
```

So neither support word can be covered by two diagonal masks.  As a hostile
control, the complete fibre-multiplicity histograms are

```text
row 1: N_0=3960, N_1=0, N_2=17472, N_3=4704, N_4=1584;
row 2: N_0=7920, N_1=0, N_2=33516, N_3=14004.       (37ixf)
```

This finishes the necessary-mask ledger: all `23,755` denominator-capacity
survivors are impossible.

The two support implementations (cyclic bitsets and merged integer arcs),
the denominator and parity ledgers, and every mask relaxation are replayed
with checks active under optimized Python.  Ordinary and optimized outputs
are byte-identical.  The frozen hashes are

```text
support script  778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a
support output  648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f
address script  870498c4f0a2d97a2d42bce593c44283c77a141fb08669b4a91133e39db5c276
address output  74f7c270034dc40b4de3d33b9abf67481435d1e97eb6e52f2448d0a152cb68d7
completion script c0a07747c300c7e82d3da27b4f498425ffb72dc950f797381d1f0d7e9096655c
completion output 3ff1c54a2818b9a0f061912758584200ffa4dd5b549ea91ba2c6c7e6a92f3638.
```

To reconnect the relaxed ledger to the exact pointwise cover, put

```text
D=lcm(d_1,d_2),          S_D=J mod D,
a_q=c_qD/d_q,            X_r=(r+R_A)/D.               (37j)
```

Then `D` divides `L`, the `a_q` are integers, and the original residual is
covered pointwise if and only if

```text
X_r subset D_(a_1) union D_(a_2)       for every r in S_D. (37k)
```

This is a finite rational interval-arrangement test.  Isolated points pay
zero in `(37h)` but must remain in the final pointwise check `(37k)`.  If
`(37k)` held, then for any fixed `u in R_A` its values at
`(r+u)/D`, `r in S_D`, would give

```text
S_D subset M_1(u) union M_2(u).
```

The carrier `R_A` has positive mass, while the completed necessary-mask
ledger proves that no such inclusion exists for any admissible denominator
pair.  This contradiction proves:

> **Literal two-drift closure.**  Five aligned tail combs and two arbitrary
> tail combs cannot cover the safe carrier of six distinct literal body
> combs.

For reuse in higher-drift branches, the exact pair sidecar is
carrier-local.  If
`a_1=g alpha`, `a_2=g beta`, with `gcd(alpha,beta)=1`, and `P_(alpha,beta)`
is a periodic primitive of

```text
1_(D_alpha)1_(D_beta)-rho(alpha,beta),
```

then for interval components `[L_s,R_s]` of `X_r`,

```text
mu(X_r intersect D_(a_1) intersect D_(a_2))
 =mu(X_r)rho(alpha,beta)
  +(1/g)sum_s[P_(alpha,beta)(gR_s)-P_(alpha,beta)(gL_s)]. (37l)
```

The endpoint term depends on the primitive pair and the located carrier.
It cannot be replaced by the global overlap `rho(a_1,a_2)`: THM-1166's
hostile interval `[1/7,6/7]` has zero local overlap for every consecutive
pair `(N,N+1)`, although its global overlap is positive.

The local-current identity remains valid, but the two-drift branch no longer
needs its growing numerator box: the coarser address quotient is already
empty uniformly over all literal bodies and aligned shapes.

## 10. The critical affine-profile residual is empty

THM-2184 leaves one honest critical possibility: a fixed seven-tail
two-torus profile could in principle have zero safe mass.  The
divisor-minimal Fourier coordinate rules it out on every positive grid
carrier.

Let `C` be a nonempty union of `1/L` cells and fix

```text
c_i in Z_(>0),             r_i in Z,                 1<=i<=7,
W_i(N)=NLc_i+r_i.                                      (39)
```

For all sufficiently large `N`, the displayed speeds are positive; an LRC
packet also requires them to be distinct.  Put

```text
d(y)=1_(||y||<1/14),

Phi_(c,r)(t)
 =integral_T product_(i=1)^7(1-d(c_i x+r_i t)) dx,

P_(C;c,r)=integral_C Phi_(c,r)(t)dt.                  (40)
```

THM-2184's grid transfer gives

```text
|mu(C intersection intersection_i D_(W_i(N))^c)
       -P_(C;c,r)|
 <=5||r||_1/(2NL).                                    (41)
```

The key point is that its integral on every positive carrier is always
positive, even though the pointwise profile can vanish at isolated slow
phases.  Fix `t` and let

```text
m_t(x)=sum_(i=1)^7 d(c_i x+r_i t).
```

Every summand has mass `1/7`, so `integral m_t=1`.  The negative part of
`m_t-1` is precisely the uncovered set and has integral `Phi_(c,r)(t)`;
the mean-zero identity makes the positive part have the same integral.
Therefore

```text
||m_t-1||_(L1)=2Phi_(c,r)(t).                         (42)
```

Let

```text
c_0=min_i c_i,              I_0={i:c_i=c_0}.
```

At `x`-Fourier frequency `c_0`, no comb with larger slope contributes.
Since

```text
hat d(1)=sin(pi/7)/pi>0,
```

equation `(42)` gives the quantitative divisor-minimal obstruction

```text
Phi_(c,r)(t)
 >=sin(pi/7)/(2pi)
   |sum_(i in I_0) exp(2pi i r_i t)|.                 (43)
```

If `P_(C;c,r)=0`, then `(43)` makes the exponential polynomial on the
right vanish almost everywhere on one open carrier cell.  It must therefore
vanish identically.  After equal exponents `r_i` are grouped, however, all
its coefficients are positive integers.  It is not the zero polynomial.
This contradiction proves

```text
P_(C;c,r)>0                                           (44)
```

for every fixed positive slope vector and integer residue vector.  Equations
`(41)` and `(44)` give the effective terminal

```text
N>5||r||_1/(2L P_(C;c,r)).                            (45)
```

In particular, if the minimum slope `c_0` is unique, `(43)` has the
pointwise floor

```text
Phi_(c,r)(t)>=sin(pi/7)/(2pi).                        (46)
```

Thus a profile approaching zero must first repeat its divisor-minimal
frequency.  This is the shifted two-scale analogue of the
Mirsky--Newman divisor-minimal obstruction, now with the slow phase retained.

There is a rational quantitative form when all seven leading slopes agree.
Absorb their common value into `N`, and suppose the fixed residues
`r_1,...,r_7` are distinct.  If

```text
g_1(t),...,g_7(t)
```

are the cyclic gaps between the seven centers `-r_i t`, then

```text
Phi_(1,r)(t)
 =sum_j max(g_j(t)-1/7,0)
 =(1/2)sum_j |g_j(t)-1/7|.                            (47)
```

Consequently the profile vanishes exactly when the centers are a coset of
`(1/7)Z/Z`.  For any nonzero residue difference `d=r_i-r_j`, summing the
gap errors along the arc between those two centers gives

```text
Phi_(1,r)(t)
 >=(1/2)dist(dt,(1/7)Z/Z)
 =||7dt||/14.                                         (48)
```

On an interval of length `ell`, put

```text
n=7|d|,       n ell=q+s,       q=floor(n ell),       0<=s<1.
```

The exact translated-interval minimum is

```text
inf_(|I|=ell) integral_I ||nt||dt
 =(q+s^2)/(4n).                                       (49)
```

When `ell<=1/7`, the facts `n=7|d|` and `n ell<=|d|`
give from `(49)`

```text
integral_I Phi_(1,r)(t)dt>=ell^2/8.                  (50)
```

An independent exact census of all `3,003` literal six-body carriers finds

```text
min_F longest_component(G_F)=23/1092
```

at `F=(1,6,7,8,10,13)`.  Taking an interval of that length in `(50)` gives
the body-uniform profile floor

```text
P_(G_F;1,r)>=529/9539712.                             (51)
```

This yields a concrete uniform finite sector.  Let

```text
R_14=14 lcm(1,...,14)=5,045,040
```

and suppose seven distinct tails lie in one common quotient block,

```text
w_i=kR_14+b_i,             0<=b_i<R_14.
```

Here `sum_i b_i<7R_14`, so `(41)` and `(51)` give

```text
safe mass
 >529/9539712-35/(2k)>0
```

for every

```text
k>=315,586.                                            (52)
```

Thus the common-ruler block index is uniformly finite.  More generally,
`(44)` closes every fixed seven-tail affine ray, with its exact rational
profile providing the ray-specific terminal.

The dependency-free referee checks the equal-slope gap and coset laws,
`57,015` pair-distance inequalities, `12,480` exact interval formulas,
the all-body component minimum, the arithmetic in `(51)`--`(52)`, and
general minimum-frequency hostile controls.  Its checks remain active under
optimized Python; ordinary and optimized outputs are byte-identical.  The
frozen hashes are

```text
script  b12f0927312c1ea56b2e7fce1937b82cd49e7978c0d820c9c7621a9be838f13a
output  418e871e5a26b805ce3e161d85c5b500753e3a43090ae111a494caa4a9e93ddc.
```

This section empties the fixed affine-ray zero-profile residual, not the
union over changing slopes and residues.  A genuinely escaping packet must
therefore change its normalized affine data, or move between scale charts,
rather than travel along one fixed two-torus ray.

## 11. Scope and information audit

The proved transport is

```text
source:       labelled carrier cells and danger-tooth incidences;
map:          t=(j+u)/L;
preserved:    the complete Boolean word when every residue b is zero;
destroyed off the fixed locus:
              slope detuning b/L and cell phase bj/L;
necessary mixed-residue sidecar:
              J, bJ mod L, component centres/radii, and gcd fibers. (38)
```

For an arbitrary grid carrier, equations `(27)`--`(35)` remain necessary
rather than sufficient and need not eliminate every one-drift row.  The
double use of the six-comb floor and reflection is what closes the literal
six-body one-drift branch.  Section 9 now closes the literal six-body,
five-aligned/two-drift branch uniformly: its decisive quotient retains the
body address multiplicities but deliberately forgets the aligned shape,
clock phases, and numerator sizes.  The theorem does **not** close an
arbitrary mixed-residue seven-wall, the branches with three or more drift
speeds, or LRC(14).  The next exact object is a higher-mask address-balance
or carrier-local multi-endpoint current, not a global pair overlap.  A
tournament on the seven labels forgets the metric widths, cell phases, gcd
fibres, and located endpoint current, and is therefore not an equivalent
quotient.
