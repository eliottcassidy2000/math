---
id: THM-4422
title: "LRC14 projection deficit and Beatty-row reduction"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. The raw projection
  target is exactly a maximum boundary-deficit inequality; every coordinate
  projection has an exact one-dimensional mod-three Beatty-row compiler; all
  three signed norm-four relation families satisfy the 6/77 target at every
  height; and two exact triples rule out every speed-independent convex
  average. The universal projection inequality and LRC(14) remain open.
source: root + projection_inequality / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality
  - THM-4413-lrc14-owner-transversality-gap-and-complete-norm-eighteen-empty-atlas
  - THM-4420-lrc14-near-doubling-ray-network-closure
script: 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422.py
output: 05-knowledge/results/lrc14_projection_deficit_beatty_row_reduction_thm4422.out
script_sha256: 3628feb6a5b5cdd0b8b543f7c26ec417fe6d52e06d789d8f784c7790f1fa3c86
output_sha256: 9e6fab64ebd4548848e1f4387e128574bb0fe161e339a0cb5b854455d738bb64
independent_script: 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422_independent_referee.py
independent_output: 05-knowledge/results/lrc14_projection_deficit_beatty_row_reduction_thm4422_independent_referee.out
independent_script_sha256: eb267fb0fa9de30a2c1e53cee604b3614988a6d854f806b178d25cdb3e69e8ff
independent_output_sha256: a07cce08a4a5c4ce1f52c09db8ca861ae1dcb529cb99c30d6b358127d744b970
hash_basis: raw LF bytes
audit: >
  PASS. The primary exact standard-library verifier rebuilds
  every carrier and every row for all 2,910 eligible triples through height
  79, audits the three norm-four families through height 499, checks both
  finite bases and the two averaging hostiles, and contains no optimizable
  assertions. Its 280,598 explicit gates pass, and ordinary and optimized
  outputs byte-match. A no-import referee independently rebuilds all 8,730
  low row compilers, the complete signed norm-four derivation and 105-row
  base, the layer-cake integral, and the two averaging hostiles; normal and
  optimized runs byte-match after 874,137 further exact gates.
---

# THM-4422 -- LRC14 projection deficit and Beatty-row reduction

**PROVED ELEMENTARY RELATIVE TO THM-4414 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** This theorem
replaces most of the remaining degree-zero capacity search by exact boundary
and one-dimensional discrepancy problems. It does not prove the universal
projection inequality, chart entry, synchronization, or `LRC(14)`, which
remains **OPEN**.

## 1. Raw projections and the exact deficit dual

**Audited continuation (2026-09-06):** the
[full-cap carrier theorem](../../05-knowledge/results/overnight_20260906_lrc_cap_carriers.md)
proves the sharp stronger bound `min_i S_i<=204/5957`, uniquely at
`(23,29,37)`, when the complete dictionary is noncollinear and has no
three collinear points. Owner-colored convex midpoint closure gives at
most eight carriers; independently replayed finite bases and a strict
all-height tail close this class. The forced progression count for larger
dictionaries does not alone pay the weighted deficits below. The independent
[one-ray](../../05-knowledge/results/lrc14_one_ray_overnight_hexagon_sep05.md)
and [two-ray](../../05-knowledge/results/lrc14_two_ray_overnight_hexagon_sep05.md)
closures further leave at least three primitive directions and an affine
collinear triple in any failure. Physical entry and LRC(14) remain open.

The independently audited [midpoint payment continuation](../../05-knowledge/results/overnight_20260906_lrc_midpoint_deficit.md)
proves `D_i=B_i+R_i` with optimal disjoint same-bucket curvature payment
`B_i`. A complete infinite two-line family has unbounded compulsory
progressions but `B_i=0`; its positive deficit is affine. Retaining `R_i`
is therefore necessary to reconstruct the deficit. This does not refute
curvature sufficiency after the already count-safe and cap classes are removed.

Let `w=(a,b,c)` be sorted, primitive, distinct, positive, odd, and nonzero
modulo three. Write

```text
q=3/(7c),
p_i(C)=3(w_j+w_k)-14|C_i|,

Lambda(w)={C in Z^3:
  C dot w=0, each C_i nonzero modulo 3, each p_i(C)>0}.  (1)
```

Put `N=|Lambda(w)|` and rename the three THM-4414 raw projection sums

```text
S_i=sum_(C in Lambda(w))
  min(q,p_i(C)/(14w_jw_k)).                              (2)
```

Define the normalized boundary deficits

```text
D_i=sum_(C in Lambda(w))
  [1-c p_i(C)/(6w_jw_k)]_+.                             (3)
```

Then, exactly and at every height,

```text
S_i=q(N-D_i),
min_i S_i<=6/77  iff  max_i D_i>=N-2c/11.               (4)
```

Indeed, for one summand `h_i=p_i/(14w_jw_k)`,

```text
min(q,h_i)=q min(1,h_i/q)
          =q(1-[1-cp_i/(6w_jw_k)]_+).                   (5)
```

Summing proves the first identity, while `6/77=q(2c/11)` proves the second.
In particular,

```text
N<=2c/11                                                (6)
```

is already a complete certificate: no projection selection is required.

The deficits are literal boundary-strip masses in the carrier lattice:

```text
D_1=sum_C [14|C_1|-3(c-b)]_+/(6b),
D_2=sum_C [14|C_2|-3(c-a)]_+/(6a),
D_3=sum_C [14c|C_3|-3c(a+b)+6ab]_+/(6ab).              (7)
```

Thus the hard regime is no longer an undifferentiated carrier count. It asks
why the surplus `N-2c/11` must enter at least one weighted boundary strip.

## 2. Exact one-dimensional row compiler

Fix the first coordinate and set

```text
d=gcd(b,c),      b=db_0,      c=dc_0.                   (8)
```

Primitivity gives `gcd(a,d)=1`. The kernel equation forces `C_1=dn`; choose
Bezout coefficients `b_0u+c_0v=1`. Every integer solution in this row is

```text
C_1=dn,
C_2=-anu+c_0t,
C_3=-anv-b_0t,             t in Z.                     (9)
```

Because `d` and every speed are units modulo three, liveness implies
`3` does not divide `n`. Moreover, three nonzero elements of `F_3` sum to zero
only when they are equal, so

```text
b_0C_2=an=c_0C_3  (mod 3).                              (10)
```

Equation `(10)` selects exactly one residue class of `t` modulo three.

The two remaining roof conditions cut this progression by the open interval

```text
max((-A_2+anu)/c_0,(-anv-A_3)/b_0)
  <t<
min(( A_2+anu)/c_0,(-anv+A_3)/b_0),                    (11)

A_2=3(a+c)/14,      A_3=3(a+b)/14.
```

Let `M_1(n)` count the members of the selected residue class in `(11)`. Then

```text
S_1=sum_n M_1(n)
 min(3/(7c),[3(b+c)-14d|n|]/(14bc)),                   (12)
```

where `3` does not divide `n` and `d|n|<3(b+c)/14`. This is an exact compiler,
not an estimate.

The interval inherited from the `C_2` roof has length

```text
3d(a+c)/(7c).
```

Allowed values of `t` are three apart, and the `C_3` roof can only shorten
the intersection. Therefore

```text
M_1(n)<=ceil(d(a+c)/(7c)).                               (13)
```

Cyclic relabelling gives the other two projections. When `gcd(b,c)=1`, the
right side of `(13)` is one, so the row word is binary. Its endpoints move
affinely with the natural address `n`, while `(10)` deletes two phase classes;
the original two-dimensional count has become a finite mod-three Beatty-type
word without losing the projection consumer.

## 3. Complete signed norm-four relation families

Suppose the speeds satisfy one of the three possible sorted signed relations
with coefficient multiset `(1,1,2)`:

```text
c=2a+b,         c=2b-a,         c=a+2b.                 (14)
```

Then every live carrier lies on exactly the corresponding primitive ray

```text
(2,1,-1),       (1,-2,1),       (1,2,-1),               (15)
```

up to integer multiplication and global sign.

For `c=2a+b`, write the kernel equation as

```text
a(C_1+2C_3)+b(C_2+C_3)=0.
```

Since `gcd(a,b)=1`, for some integer `ell`,

```text
C_1+2C_3=b ell,       C_2+C_3=-a ell.                   (16)
```

Viewed as two allowed intervals for `C_3`, their centers are
`|ell|c/2` apart, while their radii sum to `3c/7`. Since `c/2>3c/7`, the
intervals cannot meet for nonzero `ell`. Thus `ell=0`, giving the first ray.
The case `c=a+2b` is identical after using

```text
a(C_1+C_3)+b(C_2+2C_3)=0,                               (17)
```

and gives the third ray. In the arithmetic-progression case `c=2b-a`,

```text
a(C_1-C_3)+b(C_2+2C_3)=0,
```

so `C_1-C_3=b ell`; the roofs give

```text
|C_1-C_3|<3(b+c)/14+3(a+b)/14=6b/7<b,                  (18)
```

again forcing `ell=0` and the middle ray. This proves completeness, not merely
the existence of the displayed short relations.

Select the projection whose carrier coordinate has absolute coefficient two.
For a positive carrier `kv`, its cap-transition point `A`, roof endpoint `B`,
and selected projection are

| speed relation | ray `v` | projection | `A` | `B` |
|---|---|---:|---:|---:|
| `c=2a+b` | `(2,1,-1)` | `S_1` | `3a/14` | `3(a+b)/14` |
| `c=2b-a` | `(1,-2,1)` | `S_2` | `3(b-a)/14` | `3b/14` |
| `c=a+2b` | `(1,2,-1)` | `S_2` | `3b/14` | `3(a+b)/14` |

In each row `A+B=3c/14`, and the exact selected sum is

```text
S_r=2qT,
T=sum_(1<=k<B, 3 does not divide k)
  min(1,(B-k)/(B-A)).                                   (19)
```

Let `R(t)` count positive integers `k<=t` not divisible by three. Layer cake
gives

```text
T=(1/(B-A)) integral_A^B R(t)dt.                        (20)
```

The three-term block count is

```text
R(t)<=(2t+2)/3.                                         (21)
```

Using `A+B=3c/14` in `(20)--(21)` yields

```text
T<=c/14+2/3,
S_r<=3/49+4/(7c)<=6/77        for c>=35.                (22)
```

The exact finite base has only 25 admissible triples with `c<35` and closes
all of them. For a sharp leader audit, `(22)` is below `24/343` once `c>=67`;
the 105 exact rows with `c<67` have the unique first three values

```text
(1,5,11)  : 6/77,
(1,11,23) : 12/161,
(1,25,49) : 24/343.                                    (23)
```

Therefore every signed norm-four family satisfies the required projection
ceiling at every height. Only the first triple attains equality.

## 4. No speed-independent convex selector

There do not exist constants `alpha_i>=0`, independent of `w`, with
`sum_i alpha_i=1` and

```text
sum_i alpha_i S_i(w)<=6/77                              (24)
```

for every eligible triple. Two exact carrier sets suffice:

```text
Lambda(1,5,7)   ={+/- (2,1,-1)},
Lambda(5,11,17) ={+/- (1,-2,1), +/- (2,-4,2)},          (25)

S(1,5,7)   =(8/245,6/49,4/35),
S(5,11,17)=(122/1309,8/119,12/119).                     (26)
```

For the first row, subtracting `6/77` gives

```text
(-122/2695,24/539,2/55).
```

The two positive entries are at least `98/2695`, so `(24)` forces

```text
alpha_1>=49/110.                                       (27)
```

For the second row the deviations are `(20,-14,30)/1309`, forcing

```text
alpha_2>=10/17.                                        (28)
```

But `49/110+10/17=1933/1870>1`, impossible for two coordinates of a convex
weight vector. This refutes every constant weighted-pigeonhole proof, while
leaving speed-dependent or carrier-dependent selection open.

## 5. Finite inheritance and the next hostile

The all-height statements above are analytic. Separately, exact enumeration
through height 79 gives

```text
eligible triples                                      2910
automatic rows N<=2c/11                               2796
dense rows                                               114
dense rows with pairwise gcd profile (1,1,1)            114
dense rows on one norm-four direction                   113.              (29)
```

The sole multi-direction dense row is

```text
w=(19,23,29),
Lambda(w)=+/-{(1,8,-7),(10,-7,-1),(11,1,-8)},
(11,1,-8)=(1,8,-7)+(10,-7,-1),                         (30)

(S_1,S_2,S_3)=(156/4669,192/3857,3840/88711).          (31)
```

Thus the finite remainder first leaves a one-ray ordinal model at an additive
`A_2` carrier circuit. The next sharp task is to prove a projection bound for
general live hexagons `+/-{u,v,u+v}`, or find the first larger carrier
configuration needing a different sidecar. The universal inequality remains
**OPEN**.

Run

```powershell
python -B 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422.py
python -B -O 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422.py
python -B 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422_independent_referee.py
python -B -O 04-computation/lrc14_projection_deficit_beatty_row_reduction_thm4422_independent_referee.py
```
