---
id: THM-2257
title: "Depth-three common-core 169-image sieve exclusion"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. For every odd positive H and
  positive q with 13 not dividing Hq, the set image under multiplication by
  169 of the guard-hole/danger overlap E_H intersect D_q has measure at least
  164775/426496, strictly greater than 457/1183. The proof combines an exact
  169-fibre block-intersection table, THM-2080's unequal-comb overlap law,
  seven exact small reduced pairs, and direct mesh-strip arguments in the
  four high-multiplicity residue classes. In THM-2250's all-equal normalized
  core, every such image would have to lie in a three-time danger union of
  measure exactly 457/1183. This contradiction closes the last (3,4,5)
  depth-three branch. The 165 first-depth-one rows remain open, so LRC(14)
  is not claimed.
source: codex-2026-07-25-common-core-169-image-sieve
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2250-depth-three-pair-incidence-partition-reduction
related:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2252-depth-three-cross-block-spectral-resonance-gate
script: 04-computation/lrc14_depth_three_common_core_169_image_sieve_thm2257.py
output: 05-knowledge/results/lrc14_depth_three_common_core_169_image_sieve_thm2257.out
script_sha256: a0a6a8485be9e65ee54a3ec4bb52e17fe171a5bd26a4afe0e7580c30e7853d23
output_sha256: d541db3e71ee9a634015f5e6e3c935e592205437446daab6b58869ff0b8ea103
hash_basis: working-tree bytes (LF)
---

# THM-2257 -- the common core cannot survive the 169-image sieve

Put

```text
E_h={x in R/Z:||hx||<1/7},
D_q={x in R/Z:||qx||<1/14},
Sx=13x mod 1.                                         (1)
```

The image estimate proved here is

```text
measure(S^2(E_H intersection D_q))
 >=164775/426496
 =457/1183+2879/72077824
 >457/1183                                             (2)
```

for every odd positive `H` and positive `q` satisfying
`13 does not divide Hq`. The left side means the measure of the set image,
not an image counted with multiplicity. The constant in (2) is a convenient
uniform certificate, not a claimed sharp value.

The estimate supplies exactly the missing sidecar in THM-2250. In its
all-equal normalized-core branch the negative support would force every set
in (2) into a set of measure `457/1183`, which is impossible. Thus the
profile `(3,4,5)` is empty.

## 1. Remove a common dilation without losing the image

Write

```text
g=gcd(H,q),             H=ga,             q=gb.       (3)
```

Since `H` is odd and both original coefficients are thirteen-units, `a` is
odd, `gcd(a,b)=1`, and `13 does not divide gab`.

For `d>=1`, write `T_d(x)=dx`. If `gcd(g,169)=1`, then every set `A` in the
circle satisfies

```text
T_169(T_g^(-1)A)=T_g^(-1)(T_169 A).                  (4)
```

The forward inclusion is commutativity. Conversely, suppose
`gy=169z` with `z in A`. Choose integers `alpha,beta` with

```text
169alpha+g beta=1
```

and put `x=alpha y+beta z`. Then `169x=y` and `gx=z`, proving the reverse
inclusion. Haar measure is invariant under integer preimages, so (4) gives

```text
measure(T_169(E_H intersection D_q))
 =measure(T_169(E_a intersection D_b)).              (5)
```

It therefore suffices to work with the coprime reduced pair `(a,b)`.

## 2. The exact fibre multiplicity object

Fix `y` and write the 169 roots of `y` as

```text
x_k=(y+k)/169,              k in Z/169Z.              (6)
```

Label a root by

```text
j=bk mod 169.
```

This is a permutation because `b` is a thirteen-unit. On the circle of
circumference 169, the condition `x_k in D_b` asks for integer residues in
an open arc of length

```text
2(169/14)=169/7=24+1/7.
```

Its root labels therefore form a translate of a consecutive block of at
most 25 residues.

Put

```text
r=b a^(-1) mod 169.                                  (7)
```

Since `ak=r^(-1)j mod 169`, the condition `x_k in E_a` says that
`r^(-1)j` lies in an open arc of length

```text
2(169/7)=338/7=48+2/7.
```

Its labels are consequently a translate of `r` times a consecutive block
of at most 49 residues. Define the exact worst fibre intersection

```text
B_25={-12,...,12},       B_49={-24,...,24},

K_r=max_(c in Z/169Z)
    |B_25 intersection (c+rB_49)|.                   (8)
```

Smaller 24- or 48-blocks embed in the displayed maximal blocks, so `K_r`
also bounds every endpoint state. If

```text
n(y)=#{x in E_a intersection D_b:169x=y},
```

then

```text
n(y)<=K_r 1_(T_169(E_a intersection D_b))(y).        (9)
```

Root disintegration gives

```text
integral n(y)dy=169 measure(E_a intersection D_b).
```

Integrating (9) yields the basic image inequality

```text
measure(T_169(E_a intersection D_b))
 >=169 measure(E_a intersection D_b)/K_r.            (10)
```

This is the useful reframing: the image-support problem is controlled by a
finite binary relation between a 25-block and a multiplier image of a
49-block. The orientation parameter is the genuine residue `r`; replacing
it by a cosmetic tournament would discard the interval geometry.

## 3. The 156-unit fibre table

Exact enumeration of (8) over the 156 units modulo 169 gives

```text
K_r       number of r

 8             42
 9             72
10             26
12              2       r in {75,94}
13              8       r in {2,4,42,57,112,127,165,167}
17              2       r in {56,113}
25              4       r in {1,84,85,168}.          (11)
```

There are no other values. The full finite universe is asserted and
replayed by the companion script.

THM-2080 gives both

```text
measure(E_a intersection D_b)>=1/42                  (12)
```

and, with `N=ab`,

```text
measure(E_a intersection D_b)
 =2/49+(2/N)F,
F>=-1/8,

measure(E_a intersection D_b)>=2/49-1/(4N).          (13)
```

For `K_r<=10`, equations (10) and (12) already give

```text
measure(T_169(E_a intersection D_b))
 >=169/420
 =457/1183+163/10140.                                (14)
```

For the next three rows, (10) and (13) cross the target at

```text
K_r    N at least    image lower       gap over 457/1183

12        19         17407/44688       24295/7552272
13        23          1755/4508         2287/761852
17       128         164775/426496       2879/72077824. (15)
```

Only a tiny finite suffix remains below those thresholds. Under

```text
a odd, gcd(a,b)=1, 13 does not divide ab,
r=b a^(-1) mod 169,
```

the complete list and exact image masses are

```text
K_r     (a,b)                         image mass

12      (9,1)                         1
13      (1,2), (1,4), (3,2)           1, 1, 1
17      (1,56), (1,113), (3,1)        267/392, 555/791, 1. (16)
```

Every entry exceeds the last lower bound in (15). This proves (2) outside
the four `K_r=25` classes.

## 4. The four exceptional classes are resonant strips

We use one elementary mesh fact. Every open circular interval of length
strictly greater than `1/169` contains a point of every translate of the
`1/169` grid.

First suppose `r=1` or `r=-1`. Then, for some integer `ell` and
`epsilon in {1,-1}`,

```text
b=epsilon a+169ell.                                  (17)
```

For fixed `y`, let `x` run through the roots `169x=y` and put `z=ax`.
Because `a` is a unit modulo 169, these `z` form a translate of the
`1/169` grid. Relation (17) gives

```text
bx=epsilon z+ell y.                                  (18)
```

The centered arcs of radii `1/7` and `1/14` in the variable `z` have an
overlap component longer than `1/169` whenever

```text
||ell y||<3/14-1/169.
```

The mesh fact and (18) therefore prove

```text
{y:||ell y||<3/14-1/169}
 subset T_169(E_a intersection D_b).                 (19)
```

If `ell` is nonzero, the left side has measure

```text
2(3/14-1/169)=493/1183
 =457/1183+36/1183.                                  (20)
```

If `ell=0`, positivity and coprimality force `(a,b)=(1,1)`; the same mesh
argument makes the image the whole circle.

It remains to treat `r=84,85`. Since

```text
2(84)=-1 mod 169,        2(85)=1 mod 169,
```

there are `epsilon in {1,-1}` and an integer `ell` such that

```text
2b=epsilon a+169ell.                                  (21)
```

Here `ell` cannot vanish: `a` is odd, while `2b` is even. For a root of
`y`, put `s=bx`. The `s` again form a translate of the `1/169` grid, and
(21) gives

```text
2s=epsilon ax+ell y.                                  (22)
```

Choose the representative `w` of `ell y` with `|w|=||ell y||`. Near
`s=0`, the two requirements

```text
||s||<1/14,             ||2s-w||<1/7
```

are intervals of radius `1/14` whose centers are separated by `|w|/2`.
Their overlap length is

```text
1/7-|w|/2.
```

It is longer than `1/169` whenever

```text
||ell y||<2/7-2/169.
```

Consequently

```text
{y:||ell y||<2/7-2/169}
 subset T_169(E_a intersection D_b),                 (23)
```

and this strip has measure

```text
2(2/7-2/169)=648/1183
 =457/1183+191/1183.                                 (24)
```

The factor `2/169` in (23), rather than `1/169`, is necessary because the
overlap is measured in `s` and has slope one half in `w`. Keeping that
factor is a hostile endpoint audit; (24) still has ample margin.

Equations (11)--(24) prove the uniform image estimate (2).

## 5. The all-equal `(3,4,5)` core is empty

Now use THM-2250's notation. In its only surviving equality partition put

```text
u_0=u_1=u_2=u.
```

For any one of the five thirteen-unit coefficients `q_i`, the set

```text
A_i=E_H intersection D_(q_i)                         (25)
```

lies in `{R<0}`: the guard bit vanishes there while the `q_i` danger bit is
one. THM-2250's two-odd-checkpoint inclusion gives

```text
A_i subset W_0 intersection W_2                     (26)
```

up to a null set. Since `W_2=W_0 composed S^2`, equation (26) implies

```text
S^2(A_i) subset W_0.                                 (27)
```

In the common-core branch,

```text
W_0=D_u union D_(13u) union D_(169u).                (28)
```

For the single danger process `X_t=1_(D_u) composed S^t`, the stationary
safe probability is `6/7`. The exact backward root law is

```text
P(X_t=0 | X_(t+1)=0)=11/13.
```

Thus

```text
measure(W_0)
 =1-P(X_0=X_1=X_2=0)
 =1-(6/7)(11/13)^2
 =457/1183.                                          (29)
```

But (2), (27), and (29) would give simultaneously

```text
measure(S^2(A_i))>457/1183
and
measure(S^2(A_i))<=457/1183,
```

a contradiction. Hence THM-2250's all-equal partition is empty. Since that
theorem already excludes the other four equality partitions, no scalar
depth-three profile `(3,4,5)` survives.

## 6. Exact verification and failure boundary

The companion performs three independent finite tasks with exact integers
and rational arithmetic.

1. It enumerates every unit `r mod 169` and every translate in (8), checks
   the complete distribution (11), and asserts the exceptional residue
   lists.
2. It independently enumerates every reduced pair below the three product
   thresholds in (15), obtaining exactly (16), then constructs the source
   comb intervals and their 169-image union with `Fraction` endpoints.
3. It replays hostile `r=-1` resonances

   ```text
   (a,b)       exact image mass       excess over 493/1183

   (159,10)       464/1113                 29/188097
   (135,34)       394/945                  31/159705
   (87,82)        254/609                   5/14703
   (1,168)        491/1176                155/198744. (30)
   ```

The first probe in (30) is the smallest value found in the exact
`a,b<=300` discovery scan, but neither that finite scan nor a sharpness
claim is used in the proof. The analytic exceptional bound (20) lies just
below these hostile values.

Reproduce the theorem ledger with

```bash
python3 04-computation/lrc14_depth_three_common_core_169_image_sieve_thm2257.py
```

This theorem closes the last depth-three profile only. THM-2239's 165
first-depth-one rows remain outside its scope, so LRC(14) remains open.
QED.
