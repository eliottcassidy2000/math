---
id: THM-2393
title: "C3-safe double-fibre capacity and the sole common-core residual"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In
  THM-2392's ten oriented middle-depth-two cage types, delete the
  quotient high blocker as well as q_* and c_3. The resulting high-safe
  base has exact mass 396/637. Exact seven-root empty-fibre counts force
  positive clean-hole mass in eight types. A second thirteen-root peel,
  using the inherited 13-unit q_* and 13-divisible C_3/c_3 pair, forces
  it in the ninth. Combined with THM-2391, every no-clean survivor has
  M=1 and the sole orientation (a,b)=(1,1), equivalently
  (C_1,C_2,c_1,c_2)=h(1,13,13,169). All other packets have clean-hole
  mass at least 1/26754, hence a fixed q_*-labelled charged cell at least
  1/1391208 and a same-parent C_7 x C_13 cell at least 1/9042852. This is
  a structural reduction, not a row exclusion,
  ledger decrement, target landing, or proof of LRC(14).
source: codex-2026-07-26-double-fibre-capacity
depends_on:
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
related:
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
  - THM-2372-hard-septimal-signed-stalk-and-toothpick-divisibility
script: 04-computation/lrc14_c3_safe_double_fibre_capacity_thm2393.py
output: 05-knowledge/results/lrc14_c3_safe_double_fibre_capacity_thm2393.out
script_sha256: 516d6f5d7a839be77d8b3f9dd097e78ee81140b37ab7cc2d2f8021f599dba318
output_sha256: 74c473f73d67db4fdf702142e072ebac466764085f4b4bbfeb490cceb8564bb6
hash_basis: working-tree bytes (LF)
---

# THM-2393 -- C3-safe double-fibre capacity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2392 reduces the only unpriced middle-depth-two cage to ten oriented
ratios. THM-2391 independently says that only one of those ratios can
survive above the adjacent septimal layer. The present theorem changes
the object once more: it counts **empty low-cage fibres**, first over
seven roots and then over thirteen.

The two counts have different jobs:

```text
seven roots:
  six units of guard/lower-q capacity leave a K-hole;

thirteen roots:
  a nonempty empty-fibre word has at least ten roots,
  while the 13-unit q_* danger word has at most two.             (1)
```

After both peels the ten-ratio bank collapses to one literal nested
chain. The result does not make the remaining chain contradictory.

## 1. Setup and the compatible low cage

Retain THM-2388/2392's last-lane notation

```text
K=1_(E_H)+sum_(i=1)^5 1_(D_(q_i)),

X=D_(q_*)^c intersection D_(c_3)^c,

B=D_(C_1) union D_(C_2) union D_(C_3),

c_j=13C_j,

S={K=0} intersection X intersection B^c,

delta=mu(S).                                                   (2)
```

Here `q_*` is the unique top septimal `q` label. THM-2391 proves that
the guard, the four lower `q` labels, and `c_1,c_2` are all septimal
units, while

```text
nu_7(q_*)=M>=1,

nu_7(C_3)=nu_7(c_3)>M.                                        (3)
```

In the only profile not already given a uniform clean-hole floor by
THM-2392, orient the reduced low pair as

```text
(C_2,c_1)=13h(a,b),             gcd(h,91)=1,

(C_1,C_2,c_1,c_2)
 =h(b,13a,13b,169a).                                          (4)
```

The compatible low-low cage is one whole Boolean object:

```text
U_(a,b)
 =(D_b union D_(13a))
   intersection
   (D_(13b) union D_(169a)).                                  (5)
```

After the measure-preserving common dilation by `h`, the low event

```text
(D_(C_1) union D_(C_2))
 intersection
(D_(c_1) union D_(c_2))                                       (6)
```

has measure `mu(U_(a,b))`.

THM-2392 proves that `delta<1/26754` forces one of exactly ten oriented
types. Its exact whole-union table is

```text
(a,b)     mu(U_(a,b))

(1,1)       193/1183
(1,2)       114/1183
(2,1)       239/2366
(1,3)       263/3549
(3,1)        95/1183
(4,1)       331/4732
(2,3)        43/546
(3,2)        95/1183
(3,4)       491/7098
(4,3)       331/4732.                                        (7)
```

## 2. Exact mass of the high-safe base

Delete the quotient high blocker too:

```text
X_0
 =D_(q_*)^c intersection D_(c_3)^c intersection D_(C_3)^c.    (8)
```

This deletion is load-bearing. Deleting only `c_3` leaves `C_3`
available to absorb a `K`-hole.

The pair `C_3,c_3=13C_3` has exact danger overlap `1/91`, so

```text
mu(D_(C_3)^c intersection D_(c_3)^c)
 =1-2/7+1/91
 =66/91.                                                       (9)
```

Now average over the seven shifts

```text
x -> x+k/7^(M+1),                 k in F_7.                    (10)
```

The two high-blocker indicators are invariant by (3). Away from their
finitely many aligned endpoints, `D_(q_*)` occupies exactly one root.
Therefore its safe indicator contributes exactly `6/7`, not merely a
union-bound floor:

```text
mu(X_0)
 =(6/7)(66/91)
 =396/637.                                                     (11)
```

## 3. The seven-root empty-cage lemma

Disintegrate Haar measure over

```text
x_r=(y+r)/7,                         r in F_7.                 (12)
```

The base variable in (12) must use the divided high coefficients. Put

```text
Y_0
 =D_(q_*/7)^c intersection D_(c_3/7)^c
    intersection D_(C_3/7)^c.                                (12a)
```

All three coefficients in (8) are divisible by seven, and

```text
X_0=T_7^(-1)(Y_0),                       T_7(x)=7x.             (12b)
```

Consequently Haar invariance gives

```text
mu(Y_0)=mu(X_0)=396/637,                                  (12c)
```

and `y in Y_0` means that **every** root `x_r` in (12) is safe
for `q_*,c_3,C_3`. This base/point distinction is load-bearing; writing
`y in X_0` here would compare the wrong physical points. MISTAKE-264
records the post-promotion typing repair; no numerical or logical
consequence changed.

For a low-cage set `U`, write

```text
N_U(y)=#{r:x_r in U}.                                         (13)
```

The mean root count is exact:

```text
integral N_U(y)dy=7mu(U).                                    (14)
```

More importantly, if

```text
y in Y_0                  and                 N_U(y)=0,        (15)
```

then at least one of the seven roots belongs to `S`.

Indeed, (12a) makes all three high conditions safe at every root in
(12). At the seven roots, the guard word has size two and the four lower `q` words have
size one each. Their total incidence is six, so some root has `K=0`.
The scalar cover puts that root in `D_(c_1) union D_(c_2)`. Since
`N_U=0`, it is outside `D_(C_1) union D_(C_2)`; since `C_3` is safe,
it is outside `B`. It is therefore a clean hole.

Let

```text
p_0(a,b)=mu{y:N_(U_(a,b))(y)=0}.                    (16)
```

Counting clean roots and using only
`mu(Y_0 intersection {N_U=0})>=mu(Y_0)+p_0-1` gives

```text
delta
 >=(1/7) max(0,396/637+p_0(a,b)-1).                 (17)
```

The exact endpoint sweep on the common grids in (5) gives

```text
(a,b)      p_0(a,b)

(1,1)          0
(1,2)        66/169
(2,1)        62/169
(1,3)        88/169
(3,1)       248/507
(4,1)        93/169
(2,3)        84/169
(3,2)       248/507
(3,4)       281/507
(4,3)        93/169.                                  (18)
```

Thus (17) forces positive clean-hole mass in every row of (7) except

```text
(a,b)=(1,1), (2,1).                                  (19)
```

The smallest positive floor from this first peel occurs at `(1,2)`:

```text
delta>=101/57967.                                    (20)
```

## 4. A thirteen-root peel closes `(2,1)`

For `(a,b)=(2,1)`, retain on the seven-root base the empty-cage set

```text
A_0={y:N_(U_(2,1))(y)=0}.                            (21)
```

Disintegrate it over the thirteen roots

```text
y_s=(z+s)/13,                         s in F_13.      (22)
```

Exact endpoint counting gives the complete occupancy law

```text
#{s:y_s in A_0}       base-z mass

0                         7/13
10                        5/13
12                        1/13.                      (23)
```

Hence the set of bases with a nonempty `A_0` word has mass `6/13`,
and every nonempty word has at least ten roots.

In the `y` coordinate, the two high blockers have speeds `C_3/7` and
`c_3/7`. The strict profile has

```text
nu_13(C_3)>=4,                 nu_13(c_3)>=5,         (24)
```

so their joint-safe indicator is invariant on (22) and still has mass
`66/91`. Its base set meets the nonempty-`A_0` base set in mass at
least

```text
66/91+6/13-1
 =17/91.                                             (25)
```

Finally write `q_*=7u`. The inherited thirteen-unit typing gives

```text
13 does not divide u.                                (26)
```

On a thirteen-grid, the open `D_u` arc contains at most two roots.
Thus every base in (25) has at least

```text
10-2=8                                               (27)
```

roots which lie in `A_0` and are `q_*`-safe. Therefore

```text
mu_y(
 A_0 intersection D_u^c
     intersection D_(C_3/7)^c
     intersection D_(c_3/7)^c
)
 >=(8/13)(17/91)
 =136/1183.                                          (28)
```

Each such `y` forces at least one clean root in (12), exactly as in
Section 3. Consequently

```text
delta>=136/8281>0                                   (29)
```

for `(a,b)=(2,1)`.

## 5. Final residual and quantitative consequence

Only `(a,b)=(1,1)` remains. THM-2391 supplies an independent septimal
filter:

```text
M>=2:
  only (4,3) can satisfy the blocker-partition congruence modulo 49;

M>=3:
  none of the ten types can satisfy it modulo 343.   (30)
```

The type `(4,3)` is already closed by (17)--(18). Hence every no-clean
survivor has

```text
M=1,

(a,b)=(1,1),

(C_1,C_2,c_1,c_2)=h(1,13,13,169).                   (31)
```

This is a literal equality, not only a residue:

```text
C_2=c_1.                                             (32)
```

It is the common-core nested chain

```text
D_h, D_(13h), D_(169h).                              (33)
```

Every packet outside (31) has, by THM-2392 and Sections 3--4,

```text
delta>=1/26754.                                      (34)
```

THM-2392's distinguished-top deletion therefore yields a fixed
`q_*`-labelled singleton/adjacent charged cell of mass

```text
>=1/(52*26754)
 =1/1391208,                                         (35)
```

and its same-parent transverse tensor yields a fixed
`F_7 x F_13` cell of mass

```text
>=1/(338*26754)
 =1/9042852.                                         (36)
```

All nonzero septimal and target colours of that latter cell are live.
Neither (35) nor (36) supplies the terminal endpoint reference needed
to exclude a scalar row.

## 6. Sharp boundary and next target

The equality type explains why empty-fibre capacity stops. Its exact
seven-root occupancy law is

```text
N_(U_(1,1))=1 on mass 145/169,

N_(U_(1,1))=2 on mass 24/169.                        (37)
```

There is no empty fibre. Algebraically,

```text
U_(1,1)
 =D_13 union (D_1 intersection D_169),               (38)
```

so the middle tooth `D_13` supplies one low-low address on every
seven-root fibre. An abstract positive control places the two guard
roots and four lower-`q` singleton roots disjointly on the other six
addresses. Hence neither the six-unit capacity count nor a pointwise
seven-address tournament can eliminate (31).

The next exact target is the transition law of that mandatory middle
address under multiplication by thirteen. THM-2232 gives the exact
one-core Markov eigenline, and THM-2257 gives a sharp `169`-image sieve,
but neither currently identifies the `K`-hole with the mandatory
`D_(13h)` address in the present owner-typed word. The missing sidecar
is that labelled hole/address transition, not another unlabelled mass
bound.

No thirteen-adic profile is excluded, the scalar ledger remains `165`,
and LRC(14) remains open.

## 7. Exact companion

The dependency-free exact companion:

- reconstructs every whole-union mass in (7);
- computes the full seven-root occupancy distribution, including
  (18) and the equality boundary (37);
- verifies the high-safe factorization (9)--(11);
- checks every rational floor in (17), (20), and (28)--(36);
- computes the complete thirteen-root occupancy law (23);
- exhausts the thirteen-grid danger count and confirms the maximum two;
- replays THM-2391's modulo-49/modulo-343 filter on all ten oriented
  types; and
- retains explicit abstract positive controls at the common-core
  boundary.

Run

```bash
python3 04-computation/lrc14_c3_safe_double_fibre_capacity_thm2393.py
python3 -O 04-computation/lrc14_c3_safe_double_fibre_capacity_thm2393.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_c3_safe_double_fibre_capacity_thm2393.out
```

after LF normalization. Every executable check raises explicitly under
optimized Python.
