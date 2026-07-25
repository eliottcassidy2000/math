---
id: THM-2178
title: "Mod-14 transverse-code dichotomy for rank-two relation harvesting"
status: >
  PROVED + VERIFIED-EXACT. Let L be a saturated rank-two lattice of integer
  relations on a nonzero thirteen-speed row. If the length-thirteen code
  L mod 14 has minimum Hamming distance at least four, its ambient
  codimension-two torus contains a denominator-14 point whose every
  coordinate has distance at least 1/7 from zero. A quantitative local-Haar
  and relative-Fejer argument then forces, for a zero-safe row, a third
  independent relation of explicit height at most
  100*616^22*||B||_infinity^44-1, where B is any saturated lattice basis.
  Hence THM-2164 gives an unconditional fork: either rank W_H is at least
  three at one explicit universal H, or the mod-14 relation code has a
  nonzero word of support at most three. The latter is a sparse radix digit,
  not necessarily a sparse integer relation, and the theorem does not prove
  LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2054-relative-fejer-whole-product-decorrelation
  - THM-2164-relative-packet-rank-harvesting
related:
  - THM-2163-radix-relation-carry-descent
  - THM-2171-radix-quotient-order-pumping
script: 04-computation/lrc14_mod14_transverse_code_rank_harvest_thm2178.py
output: 05-knowledge/results/lrc14_mod14_transverse_code_rank_harvest_thm2178.out
script_sha256: 1ab04acdb1862c383d570bca7f44d21e2ada6d2ae7a0973835f9b1e3208e3ca2
output_sha256: 5475a4fddd2418a772f46dd8d17fae7ae78799a60b9a8a669c06af740ec387b6
hash_basis: working-tree bytes (LF)
---

# THM-2178 -- mod-14 transverse-code rank harvesting

Put

```text
J=[1/14,13/14],
S(v)={t in R/Z:v_i t in J for every i}.               (1)
```

The new carrier is the reduction modulo `14` of a known relation plane.
Its Hamming distance gives a sharp finite alternative:

```text
large transverse distance -> a deep safe torsion point -> a new relation;
small transverse distance -> a support-at-most-three radix digit.          (2)
```

The first implication is quantitative. It is the first rank-two-known-lattice
instance of the relative packet principle in THM-2164.

## 1. The mod-14 relation code

Let `L` be a saturated rank-two sublattice of `Z^13`, and choose a lattice
basis as the rows of

```text
B in Z^(2 x 13).                                      (3)
```

Define the basis-independent relation code

```text
C_14(L)={yB mod 14:y in (Z/14Z)^2}.                   (4)
```

Changing the lattice basis left-multiplies `B` by an element of
`GL_2(Z/14Z)`, so it does not change this code or its Hamming weights.
Saturation says that the gcd of the `2 x 2` minors of `B` is one. Hence
`B` has rank two modulo both `2` and `7`, and the map in (4) is injective.
Thus `C_14(L)` has exactly `196` words.

Let

```text
A={2,3,...,12} subset Z/14Z.                          (5)
```

> **Transverse-code lemma.** If every nonzero word of `C_14(L)` has
> Hamming weight at least four, then there is an
>
> ```text
> a in A^13,                    Ba=0 mod 14.           (6)
> ```

### Fourier count

Let `zeta=exp(2 pi i/14)` and

```text
S(z)=sum_(a in A) zeta^(za).                          (7)
```

Character orthogonality counts the solutions of (6) exactly:

```text
N_B
 =1/14^2 sum_(y in (Z/14Z)^2)
          product_(i=1)^13 S((yB)_i).                 (8)
```

Here `S(0)=11`. For `z!=0 mod 14`,

```text
S(z)=-(1+zeta^z+zeta^(-z))
    =-(1+2 cos(pi z/7)).                              (9)
```

Among the nonzero residues,

```text
|1+2 cos(pi z/7)|<=1+2 cos(pi/7)<141/50.             (10)
```

The last rational bound is elementary. Since `pi>3`, cosine is decreasing
on the relevant interval and its alternating Taylor polynomial gives

```text
cos(pi/7)<cos(3/7)
 <1-(3/7)^2/2+(3/7)^4/24
 =17471/19208
 <91/100.                                             (11)
```

There are `195` nontrivial characters in (8). Under the distance-four
hypothesis, each corresponding product has at least four nonzero factors.
Therefore

```text
N_B
 >=11^9/196 [11^4-195(141/50)^4]                     (12)

 =6805833364678152211/245000000
 >0,
```

because, after multiplying the bracket by `50^4`,

```text
11^4*50^4-195*141^4=14431688605>0.                   (13)
```

This proves the lemma.

The point

```text
x_0=a/14 in (R/Z)^13                                 (14)
```

lies on the relation subtorus

```text
K_L={x:Bx=0 mod 1}.                                   (15)
```

Every coordinate of `x_0` has circular distance at least `1/7` from zero.
Thus this is not merely a weak boundary witness: it lies a uniform distance
inside the `1/14`-safe cube.

## 2. Quantitative local Haar mass

Assume now that

```text
L subset Lambda(v):={m in Z^13:m.v=0},               (16)
```

where every `v_i` is nonzero. Put

```text
R=||B||_infinity.                                    (17)
```

Choose a nonzero `2 x 2` minor `D` of `B` and write

```text
Delta=|det D|>=1.                                     (18)
```

Project `K_L` onto the other eleven coordinates. Since the integer torus
endomorphism defined by `D` is surjective, this projection is onto
`T^11`, with `Delta` local inverse branches of equal Haar weight.
Saturation makes `K_L` connected; it also ensures that every nontrivial
coordinate character used below pushes Haar measure to Haar measure.

Use the branch through `x_0`. If every free coordinate moves by less than

```text
c=Delta/(1232 R^2),                                  (19)
```

then each of the two dependent coordinates moves by at most

```text
(22 R^2/Delta)c=1/56.                                (20)
```

Indeed, every entry of `adj(D)` is at most `R`, every column of the
remaining matrix is at most `R`, and there are eleven free columns.
Also `Delta<=2R^2`, so `c<1/28`.

Let `I_i` be the interval of radius `1/28` about `(x_0)_i`. Equations
(19)--(20) give the explicit Haar lower bound

```text
mu_(K_L)(product_i I_i)
 >=(2c)^11/Delta
 =Delta^10/(616^11 R^22)
 >=1/(616^11 R^22).                                  (21)
```

Every `I_i` lies strictly inside `J`, so their product is a safe box.

## 3. Relative Fejer harvesting

Define

```text
A_R=616^11 R^22,
N_R=100 A_R^2,
H_R=N_R-1=100*616^22 R^44-1.                         (22)
```

> **Rank-harvesting theorem.** If `mu(S(v))=0` and the code in (4) has
> minimum distance at least four, then there is
>
> ```text
> m in Lambda(v)\L,             ||m||_infinity<=H_R. (23)
> ```

### Proof

Let `f_i=1_(I_i)` and let `p_i=F_(H_R)*f_i`, with the Fejer convention of
THM-2054. Each coordinate character is nonzero on `K_L`: otherwise
`e_i` would lie in `L subset Lambda(v)`, contradicting `v_i!=0`.
The same is true on the line `t|->vt`.

For either Haar space, product telescoping and THM-2054's BV estimate give

```text
|integral product_i f_i-integral product_i p_i|
 <=13 epsilon_N,

epsilon_N<=(1+2 log N_R)/(2N_R).                     (24)
```

Suppose (23) fails. A Fourier tuple in `[-H_R,H_R]^13` survives the line
average exactly when it lies in `Lambda(v)`, and survives the `K_L`
average exactly when it lies in the annihilator `L`. The failure of (23)
makes these two finite index sets equal. Hence

```text
integral_T product_i p_i(v_i t)dt
 =integral_(K_L) product_i p_i(x_i)dx.                (25)
```

The left unsmoothed product vanishes almost everywhere because its support
is contained in `S(v)`. The right unsmoothed product has mass at least
(21). Equations (24)--(25) would therefore imply

```text
1/A_R<=13(1+2 log N_R)/N_R.                          (26)
```

But `A_R>=2`, `log 100<5`, and `log A_R<=A_R-1`, so

```text
13(1+2 log(100 A_R^2))
 <13(4A_R+7)
 =52A_R+91
 <100A_R.                                             (27)
```

After division by `N_R=100A_R^2`, (27) contradicts (26). This proves
(23). Since `L` is saturated, `m notin L` is independent over `Q` of the
two rows of `B`. QED.

## 4. The unconditional LRC(14) fork

Let `v` be a distinct positive thirteen-speed row with `mu(S(v))=0`.
THM-2164 proves

```text
dim_Q W_105(v)>=2.                                   (28)
```

Choose independent `u,w` in that height-`105` relation box and let `L` be
the saturation of their span. A saturated basis of `L` may be chosen with

```text
||B||_infinity<=26*105^2=286650.                     (29)
```

For completeness, the Euclidean covolume of `L` is no larger than the
area of `u,w`, hence at most `13*105^2`. A Gauss-reduced lattice basis has
angle at least `60` degrees; its longer vector has norm at most twice the
covolume because the shorter nonzero integer vector has norm at least one.
This proves the slightly loose integer bound (29).

Set

```text
H_*=100*616^22*286650^44-1.                          (30)
```

Applying Sections 1--3 gives the unconditional dichotomy

```text
dim_Q W_(H_*)(v)>=3,                                 (31)
```

or

```text
C_14(L) contains a nonzero word of support <=3.      (32)
```

Thus failure of effective rank-three harvesting is no longer an unnamed
interference channel. It is a finite, basis-independent sparse digit in
the reduction modulo `14` of the known relation plane.

If `yB mod 14` has support `S` in (32), the actual relation `yB` has every
coefficient outside `S` divisible by `14`. This is the exact relation-carry
interpretation:

```text
known rank-two plane
 -> one base-14 digit supported on at most three coordinates
 -> higher digits retain the remaining relation current.                (33)
```

## 5. Boundaries and hostile controls

The support-three alternative in (32) is **not** an integer relation
supported on three coordinates. Coefficients outside its support vanish
only modulo `14`; dividing them out requires a carry and may restore all
thirteen coordinates.

Likewise, the point `x_0` lies on the ambient codimension-two torus
`K_L`, not necessarily on the original one-parameter orbit. It becomes a
new bounded relation through the relative-Fejer comparison, not directly
an LRC witness.

The distance condition is sufficient, not necessary. The saturated lattice
generated by `e_1-e_2,e_2-e_3` has code distance two but plainly admits
many points from `A^13`; the sparse-digit branch can coexist with a safe
ambient point. Conversely, the affine-column basis

```text
B_aff=
 [1 1 ... 1
  0 1 ... 12]                                        (34)
```

is saturated and has exact code distance six. The companion dynamic program
counts

```text
176136387828
```

solutions of (6) for this positive control.

Finally, the height (30) is deliberately crude. The theorem does not claim
that rank three closes LRC(14), that every sparse radix digit lifts to a
sparse relation, or that the support-three branch is empty. Its gain is a
rigorous new terminal interface:

```text
effective third relation
or
three-coordinate base-14 carry packet.               (35)
```

The companion verifies every rational inequality, the code cardinalities
and distances in both controls, and the exact finite solution count. QED.
