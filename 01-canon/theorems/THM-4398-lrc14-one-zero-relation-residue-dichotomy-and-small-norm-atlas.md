---
id: THM-4398
title: "LRC14 one-zero relation residue dichotomy and small-norm atlas"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4386/4393 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Primitive relations on ternary-unit speed triples
  have either three unit coefficients or exactly one zero coefficient modulo
  three. The latter have nonzero defect and one live affine class. All 14
  one-zero presentation sectors of l1 norm at most 14 have sharp global
  maxima at or below 6/77. Seam entry, arbitrary nonresonance, synchronization,
  and LRC(14) remain OPEN.
source: root + jc4385_elliptic independent referee / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
  - THM-4393-lrc14-minimal-ternary-unit-norm-eighteen-shell
related:
  - THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality
  - THM-4394-lrc14-minimal-ternary-unit-norm-twenty-shell
primary_script: 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398.py
primary_output: 05-knowledge/results/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398.out
primary_script_sha256: ac37abac2e52e2800548854bc69ace568a10a21bb18482876d9bd16d23eb5cbb
primary_output_sha256: a2a29363ccd3673f0982e7f936755da54f8b2f773dde7cde36ed9d0d59f4f866
independent_audit_script: 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.out
independent_audit_script_sha256: 587434c319e4c7fac3ad9db933f55334ce75005cdee3598e0d5373874211dbad
independent_audit_output_sha256: bdd5352801ea0b47ada2ac82df6221b7f89de898d2e7e73ee72184ff504007b0
hash_basis: raw LF bytes
audit: >
  PASS. The primary imports the audited THM-4393 exact carrier primitives;
  the independent verifier imports no repository mathematics, reconstructs
  the finite-field gate, two independent finite universes, cube slices,
  affine carriers, physical circle intersections, and unit-sector controls,
  and rejects an intentional tripwire. Normal, optimized, and hash-seeded
  runs agree. No mathematical or status discrepancy remains.
---

# THM-4398 -- one-zero-mod-three relation master formula and complete `l1<=14` atlas

**Status: PROVED ELEMENTARY RELATIVE TO THM-4386/4393 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**  The residue dichotomy, affine carrier formula,
cube-slice tail bound, and all fourteen sharp presentation-sector maxima are
derived below.  The primary reuses audited carrier primitives; the independent
verifier imports no repository mathematics and compares every quotient
presentation in every finite proof box with a definition-level physical
circle construction.

These are **presentation sectors**, not minimal-relation sectors: a triple is
included whenever it has the stated relation, even if it also has a shorter
unit relation.  This distinction is load-bearing because every sharp
maximizer below lies in the already proved unit `(1,1,2)` sector.  Nothing
here supplies seam entry, synchronization, a bound for arbitrary triples, or
an LRC proof.  **LRC(14) remains OPEN.**

## 1. The mod-three relation-defect dichotomy

Let `w=(w1,w2,w3)` have each coordinate nonzero modulo three, and let
`c in Z^3` be a primitive nonzero relation `c.w=0`.  Put

```text
a_i=c_i w_i mod 3.
```

Primitivity forbids all three `a_i` from vanishing.  Exactly two cannot
vanish, since the remaining nonzero term would contradict `sum a_i=0`.
Therefore there are precisely two residue types:

```text
UNIT:      every c_i is nonzero mod 3;
ONE-ZERO:  exactly one c_i is zero mod 3.                (1)
```

At a physical scale-three failure phase write

```text
n_i=nint(w_i y),       e_i=w_i y-n_i,       |e_i|<3/14,
o_i=-w_i^(-1)n_i mod 3,       {o_1,o_2,o_3}=F_3,
delta=c.n=-c.e.
```

In the unit case, three nonzero elements of `F_3` summing to zero are equal:
`a_1=a_2=a_3=A`.  Hence

```text
delta=-A(o_1+o_2+o_3)=0 mod 3.                          (2)
```

In the one-zero case, suppose `a_3=0`.  Then
`a_1=A`, `a_2=-A`, and

```text
delta=A(o_2-o_1) != 0 mod 3,                            (3)
```

because the owners are distinct.  This proves the defect half of the
dichotomy without a size bound.

The carrier half is equally sharp.  Let

```text
C=w cross n,       Lambda_w={C:C.w=0}.
```

Modulo three, the all-nonzero elements of `Lambda_w` are exactly

```text
+/-u,       u=(w1^(-1),w2^(-1),w3^(-1)).                (4)
```

Indeed, the three nonzero values `w_i C_i` sum to zero and therefore are all
equal.  For a primitive relation choose `v.c=1` and put

```text
C_delta=w cross (delta v).
```

The complete integral defect fibre is

```text
C_delta+k c,       k in Z.                              (5)
```

To see integrality directly, choose `u.w=1`; then `n=C cross u` is an
integral inverse to `C=w cross n` modulo `Zw`.  Differences of two lifts at
fixed defect lie with `w` in `c^perp`, so their carriers differ by an integer
multiple of primitive `c`; conversely surjectivity supplies `h` with
`w cross h=c` and `c.h=0`.

In the unit case `c` is parallel to `u` modulo three.  The zero-residue defect
line in (5) contains the three classes `0,+u,-u`, so exactly **two** values
of `k mod 3` are owner-live.  In the one-zero case `c` is not parallel to the
full-support vector `u`; equations (3)--(4) place `+u` and `-u` in the two
opposite nonzero defect residues.  Each physically allowed fixed-defect line
therefore has exactly **one** live value of `k mod 3`.

An exhaustive finite-field check independently verifies this mechanism for
all `64` admissible `(w,c)` residue pairs: eight speed-residue triples and the
eight nonzero kernel relations for each:

```text
16 unit cases:      live-fibre counts (delta=0,1,2) = (2,0,0);
48 one-zero cases:  live-fibre counts (delta=0,1,2) = (0,1,1).
```

This is a structural dichotomy, not a numerical pattern observed in the
sector scan.

## 2. The fourteen coefficient shapes and their defect states

For odd speeds, `c.w=0 mod 2` implies `||c||_1` is even.  Enumerating sorted
positive primitive magnitude triples with exactly one coordinate divisible
by three gives precisely:

```text
l1=6:   (1,2,3)
l1=8:   (1,1,6), (1,3,4)
l1=10:  (2,3,5)
l1=12:  (1,2,9), (1,3,8), (1,5,6), (2,3,7), (3,4,5)
l1=14:  (1,1,12), (1,3,10), (1,4,9), (1,6,7), (3,4,7). (6)
```

Strict eligibility and (3) give

```text
|delta| < (3/14)||c||_1,       delta nonzero mod 3.
```

Thus norms six and eight have defects `{-1,+1}`, while norms ten, twelve,
and fourteen have defects `{-2,-1,+1,+2}`.  No endpoint equality is admitted.

## 3. Exact roof, cube slices, and tail bound

For any carrier define the six-term intersection roof

```text
L_w(C)=max(0,min(
  2r/w1, 2r/w2, 2r/w3,
  r/w1+r/w2-|C3|/(w1w2),
  r/w1+r/w3-|C2|/(w1w3),
  r/w2+r/w3-|C1|/(w2w3))),       r=3/14.                (7)
```

It is exactly the common length of the three nearest-integer intervals.  The
one-zero relation formula is therefore

```text
mu(F_w)=sum_(allowed delta)
        sum_(k in the unique owner-live class mod 3)
          L_w(C_delta+k c).                             (8)
```

The three pair terms in (7) give strict rational bounds on `k`, making (8)
an exact finite sum.

For real `t`, `f_delta(t)=L_w(C_delta+t c)` has interval superlevel sets, so
it is nonnegative, compactly supported, and unimodal.  Its integral depends
only on the magnitude shape `p=(a,b,c)`:

```text
I_delta(p)=area([-r,r]^3 intersect {p.e=-delta})/||p||

 =1/(2abc) sum_(S subset {1,2,3}) (-1)^|S|
   (|delta|+r(a+b+c)-2r sum_(i in S)a_i)_+^2.           (9)
```

The verifier reconstructs every value twice: by (9), and by exact rational
polygon clipping of `[-r,r]^2` against
`|delta+a e1+b e2|<=c r` followed by division by `c`.

One residue class of a step-three sample obeys the unimodal rectangle bound

```text
sum_j f(alpha+3j) <= integral(f)/3 + sup(f),
sup(f)<=2r/W=3/(7W),       W=max(w).
```

If `D_p` is the allowed defect set, this gives the global envelope

```text
mu(F_w) <= B_p + 3|D_p|/(7W),
B_p=(1/3) sum_(delta in D_p) I_delta(p).                (10)
```

## 4. Complete sharp atlas

The exact bulk constants, analytic thresholds, and maxima are:

| shape | positive slice bulks `(I_1,I_2)` | `B_p` | equality threshold | proof height | sharp mass | unique unordered winner |
|---|---|---:|---:|---:|---:|---|
| `123` | `(1/147,-)` | `2/441` | `945/31` | 30 | `8/245` | `{1,5,7}` |
| `116` | `(17/588,-)` | `17/882` | `8316/569` | 14 | `6/77` | `{1,5,11}` |
| `134` | `(1/56,-)` | `1/84` | `792/61` | 12 | `6/77` | `{1,5,11}` |
| `235` | `(1/49,1/2940)` | `61/4410` | `83160/3109` | 26 | `6/77` | `{1,5,11}` |
| `129` | `(1/49,5/588)` | `17/882` | `143640/4181` | 34 | `46/665` | `{5,7,19}` |
| `138` | `(53/2352,5/784)` | `17/882` | `25704/755` | 34 | `58/833` | `{5,7,17}` |
| `156` | `(19/980,1/196)` | `4/245` | `1190/37` | 32 | `58/833` | `{5,7,17}` |
| `237` | `(23/1029,4/1029)` | `6/343` | `539/19` | 28 | `6/77` | `{1,5,11}` |
| `345` | `(39/1960,2/735)` | `19/1260` | `28080/833` | 33 | `6/91` | `{1,7,13}` |
| `1112` | `(3/196,3/196)` | `1/49` | `1932/61` | 31 | `12/161` | `{1,11,23}` |
| `1310` | `(9/490,11/980)` | `29/1470` | `57960/1853` | 31 | `12/161` | `{1,11,23}` |
| `149` | `(71/3528,11/1176)` | `26/1323` | `6237/212` | 29 | `6/77` | `{1,5,11}` |
| `167` | `(25/1372,11/1372)` | `6/343` | `13965/421` | 33 | `46/665` | `{5,7,19}` |
| `347` | `(167/8232,1/168)` | `6/343` | `4998/145` | 34 | `8/119` | `{5,11,17}` |

Here “equality threshold” is the exact real number

```text
[3|D_p|/7] / [M_p-B_p].
```

For every integer height strictly above the displayed proof height, (10) is
strictly below `M_p`.  Exhaustion through the displayed height therefore
proves the global maximum in that presentation sector.

The winner defect masses are:

| shapes | masses |
|---|---|
| `123` | `mu_-1=mu_1=4/245` |
| `116,134` | `mu_-1=mu_1=3/77` |
| `235,237` | `mu_-1=mu_1=3/77`, `mu_-2=mu_2=0` |
| `129,167` | `mu_-1=mu_1=3/133`, `mu_-2=mu_2=8/665` |
| `138,156` | `mu_-1=mu_1=3/119`, `mu_-2=mu_2=8/833` |
| `345` | `mu_-1=mu_1=3/91`, `mu_-2=mu_2=0` |
| `1112,1310` | all four masses `3/161` |
| `347` | `mu_-1=mu_1=3/119`, `mu_-2=mu_2=1/119` |

The `149` winner has two charts, and their defect split is not intrinsic:

```text
c=(1,-9,4):  mu_-1=mu_1=3/77,  mu_-2=mu_2=0;
c=(9,-4,1):  mu_-2=mu_2=3/77,  mu_-1=mu_1=0.           (11)
```

Both compute the same carrier dictionary and physical mass.  Equation (11)
is a useful warning not to treat `delta` as chart-independent.

## 5. Finite universes and definition-level controls

Each proof universe consists of sorted primitive distinct positive odd
three-unit speed triples with at least one signed placement of the named
coefficient shape.  No shorter-relation filter is applied.

| shape | height | triple incidences | relation rays | positive combs | carrier components |
|---|---:|---:|---:|---:|---:|
| `123` | 30 | 51 | 55 | 33 | 70 |
| `116` | 14 | 2 | 2 | 2 | 4 |
| `134` | 12 | 3 | 3 | 2 | 4 |
| `235` | 26 | 20 | 22 | 15 | 30 |
| `129` | 34 | 22 | 22 | 20 | 46 |
| `138` | 34 | 27 | 27 | 24 | 60 |
| `156` | 32 | 27 | 28 | 24 | 66 |
| `237` | 28 | 18 | 19 | 14 | 34 |
| `345` | 33 | 38 | 41 | 32 | 80 |
| `1112` | 31 | 8 | 8 | 7 | 18 |
| `1310` | 31 | 20 | 20 | 17 | 44 |
| `149` | 29 | 17 | 19 | 16 | 38 |
| `167` | 33 | 23 | 23 | 21 | 56 |
| `347` | 34 | 24 | 24 | 21 | 54 |

The row sums are `300` sector incidences, `313` relation rays, and `604`
carrier-component incidences; they represent `136` distinct speed triples,
of which `114` have positive combs.

Universe generation is checked by two independent finite routes: solving
each labelled relation for the third sorted speed, and brute force over every
speed triple in the box.  For every one of the `313` relation rays, (8) is
compared as an exact dictionary `C -> Fraction` with a separate construction
which intersects the three lists

```text
|w_i y-n_i|<3/14
```

on the physical circle and retains precisely the cells whose owners are
`{0,1,2}`.  Every carrier and length agrees, not only the total measure.

## 6. Equality and overlap with the proved unit sectors

Every sharp winner above has a primitive unit relation of shape `(1,1,2)`:

| physical winner | unit relation | nonunit shapes attaining their maxima there |
|---|---|---|
| `{1,5,7}` | `(2,1,-1)` | `123` |
| `{1,5,11}` | `(1,2,-1)` | `116,134,235,237,149` |
| `{5,7,19}` | `(1,2,-1)` | `129,167` |
| `{5,7,17}` | `(2,1,-1)` | `138,156` |
| `{1,7,13}` | `(1,-2,1)` | `345` |
| `{1,11,23}` | `(1,2,-1)` | `1112,1310` |
| `{5,11,17}` | `(1,-2,1)` | `347` |

The verifier evaluates each winner a second time through the unit zero-defect
chart from
[`THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md`](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md)
and checks exact carrier-dictionary equality with every nonunit chart and the
physical comb.  Thus repeated values in the atlas are genuine identities of
one physical object, not masses to add.

Across the fourteen presentation sectors, the largest value is `6/77`.  It
appears in five rows but at the single underlying speed triple `{1,5,11}`.
This agrees with the existing unit-sector maximum and does not create five
new extremizers.

As a deliberately non-minimal common-height sidecar, at `H=34` there are

```text
399 nonunit sector incidences,
418 nonunit relation rays,
150 distinct speed triples,
124 triples also carrying a unit l1<=14 relation,
123 triples carrying more than one nonunit sector label.
```

The high overlap rate is why this theorem is stated as a presentation atlas.
It is not a partition of speed triples or phase space.

## 7. Verification and scope firewall

The independent verifier uses only Python's standard library, exact
`Fraction` arithmetic, and explicit checks that survive `python -O`.  It
performs `67,529` live checks, including the full finite-field dichotomy,
two-route pattern and speed universes, two-route cube slices, all analytic
cutoffs, every quotient/physical dictionary, every sharp winner, and every
unit-sector equality.

The primary verifier independently replays all fourteen sharp rows and their
individual tail cutoffs using THM-4393's audited exact primitives.  It checks
every one of the same `313` relation rays against definition-level physical
circle dictionaries.  Its `628` relation-chart carrier incidences deliberately
count a component again when one triple has two relations in the same sector;
the referee's `604` physical carrier incidences count it once per sector
triple.  This is a documented counting convention, not a discrepancy.

Reproduction:

```powershell
python -B 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398.py
python -B -O 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398.py
python -B 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.py
python -B -O 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.py
$env:PYTHONHASHSEED='314'
python -B -O 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.py
python -B -O 04-computation/lrc14_one_zero_relation_residue_dichotomy_small_norm_atlas_thm4398_independent_audit.py --tripwire
```

Normal, optimized, and hash-seeded streams are byte-identical to their frozen
canonical outputs; all four raw-byte hashes are in the theorem metadata, and
the tripwire run fails as required.  The result closes exactly the fourteen
one-zero-mod-three full-support primitive relation presentation sectors of
`l1<=14`.  It neither asserts minimality nor controls relation-free triples.
**LRC(14) remains OPEN.**
