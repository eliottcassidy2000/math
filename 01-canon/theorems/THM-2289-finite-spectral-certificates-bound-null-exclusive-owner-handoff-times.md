---
id: THM-2289
title: "Finite spectral certificates bound null exclusive-owner handoff times"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT EXACT ARITHMETIC REFEREE. In every
  one of the 165 first-depth-one scalar profiles, choose the labelled
  exclusive owner supplied by THM-2263. The set of times k at which its
  ancestry-aware blocker-only return has zero Haar measure has cardinality
  at most M_strict=((115919^9-1)(115919^7-1))/2 on the 150 strict rows and
  M_repeat=((331095^9-1)(331095^7-1))/2 on the 15 repeated-first rows.
  Consequently a genuine blocker-to-other-blocker handoff occurs within
  the corresponding coefficient-independent number of steps after the
  owner's prescribed expiration. Its Haar mass is strictly larger than
  12331416859/792352055420125200 on a strict row and
  10450627633/246356440619713023600 on a repeated-first row. Every
  M_branch+1 distinct times contain such a quantitatively positive time.
  The mechanism is that any exceptional time has a genuinely crossing
  Jackson certificate A+13^k B=0 with A,B nonzero, while the same signed
  certificate cannot occur at two times. The theorem does not force the
  exact expiration time, excludes no valuation profile, and does not prove
  LRC(14).
source: codex-2026-07-25-bounded-handoff-probe
depends_on:
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2193-uniform-rank-six-safe-torus-floor
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
related:
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
script:
  - 04-computation/lrc14_finite_spectral_null_handoff_thm2289.py
  - 04-computation/lrc14_finite_spectral_null_handoff_referee_thm2289.py
output:
  - 05-knowledge/results/lrc14_finite_spectral_null_handoff_thm2289.out
  - 05-knowledge/results/lrc14_finite_spectral_null_handoff_referee_thm2289.out
script_sha256:
  - ec78100473cac073e33f8be2b3a36a167e304727f281650c001e4188b5bed648
  - b9436eb2decefb68b30f78568d9c29de7e5cc556171f1c5f55a2b42ea2ce2c56
output_sha256:
  - 1a448087b57d5479ac4d00bc0ad70effcd58a943765830e46cecc6a33bd31f0f
  - 766b5287a809b19c631c68b4a05c9dd4d6f06e239070976670059a49f18c7133
hash_basis: working-tree bytes (LF)
---

# THM-2289 -- null exclusive-owner handoff times have finite spectral type

Use the scalar five-unit/three-blocker notation

```text
T(x)=13x mod 1,

D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i),                 (1)
```

with the almost-everywhere cover

```text
A_0 subset D_(c_1) union D_(c_2) union D_(c_3)       (2)
```

and one of the `165` first-depth-one valuation profiles

```text
(lambda_1,lambda_2,lambda_3)=(1,b,c),

either 2<=b<c (150 strict rows),
or     b=1<c (15 repeated-first rows),     5<=c<=19, (3)

c_h=13^(lambda_h)u_h.
```

All coefficients are positive integers, `H` is odd, the usual scalar
distinctness hypotheses hold, and
`H,q_1,...,q_5,u_1,u_2,u_3` are thirteen-units.

THM-2263 selects one label `j in {1,2,3}` for which the exclusive source
and blocker-only target

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(h!=j)D_(c_h),

R_j=A_0 minus D_(c_j)                                (4)
```

satisfy the branch-dependent bounds

```text
strict:
  measure(E_j)>=alpha_strict:=15041431/593783190;

repeated-first:
  measure(E_j)>=alpha_repeat:=5229541/593783190;

both branches:
  measure(R_j)>=beta:=2593/90090.                    (5)
```

The two source floors are THM-2263's labelled one-third owner floors.
For the target, `A_0` has mass at least `961/6930`, while the uniform
depth-positive guard/blocker cap is at most `10/91`: THM-2273's
parity-sharp fold gives

```text
d odd:   5/49+5/(49*13^d),
d even:  5/49+5/(294*13^d),
```

whose unique global maximum for `d>=1` is the odd value at `d=1`,
namely `10/91`. Subtraction gives the displayed `beta`. This derivation is
independent of which label THM-2263 selects, including the deep label in a
repeated-first row.

Let

```text
Z_j={k in Z_(>=0):
       measure(E_j intersection T^(-k)R_j)=0}.       (6)
```

Define the quantitative branch floors

```text
gamma_strict
 =12331416859/792352055420125200,

gamma_repeat
 =10450627633/246356440619713023600.                 (6a)
```

The theorem proves the stronger coefficient-independent branch bounds

```text
#{k:measure(E_j intersection T^(-k)R_j)
       <=gamma_strict}
 <=M_strict
 =((115919^9-1)(115919^7-1))/2
 =531427494109850809274382490322199270940210549004984817720024492918545763422455682

on strict rows, and

#{k:measure(E_j intersection T^(-k)R_j)
       <=gamma_repeat}
 <=M_repeat
 =((331095^9-1)(331095^7-1))/2
 =10428262671922434393098039414106300324101457669423374402844539686260225238346754941835938

on repeated-first rows. In particular `#Z_j` obeys the same applicable
bound.                                                     (7)
```

Put `M_branch=M_strict` or `M_repeat` according to the row. Then some

```text
k in {lambda_j+1,...,lambda_j+1+M_branch}           (8)
```

has

```text
measure(E_j intersection T^(-k)R_j)>gamma_branch.   (9)
```

More strongly, every `M_branch+1` distinct nonnegative times contain a time
satisfying (9).

## 1. The two rectangular events

Up to null endpoints, each condition in (4) is one heptimal circular
interval condition on an integer character. Put

```text
C={z:||z||>1/7},
D={z:||z||<1/14},
G={z:||z||>=1/14}.                                  (10)
```

The source `E_j` is the pullback of the nine-coordinate rectangle

```text
C x G^5 x D x G^2                                  (11)
```

along

```text
x |-> (Hx,q_1x,...,q_5x,c_jx,c_hx,c_lx),            (12)
```

where `{j,h,l}={1,2,3}`. The target `R_j` is the pullback of the
seven-coordinate rectangle

```text
C x G^6                                             (13)
```

along

```text
x |-> (Hx,q_1x,...,q_5x,c_jx).                      (14)
```

Thus the source block has `r=9` coordinates and the target block has
`s=7`. This exact grouping is important: no inclusion-exclusion expansion
and no boundary-count estimate is used.

## 2. A mixed-rectangle crossing lemma

Let `I` be any one of the three circular intervals in (10). For `N>=2`,
let `Q_N` be the normalized squared-Fejer kernel of THM-2193 and put

```text
q_I=Q_N*1_I,                   H_N=2N-2.             (15)
```

THM-2193 proves

```text
0<=q_I<=1,
FourierSupport(q_I) subset [-H_N,H_N],

||q_I-1_I||_1<epsilon_N:=3/(2N).                    (16)
```

The last estimate is uniform over circular intervals: translating an
interval indicator by `t` costs at most `2||t||` in `L^1`.

Let `alpha` denote the applicable source floor in (5). Let `Q_E` and
`Q_R` be the products of the nine and seven smoothed factors in (11) and
(13), after their integer-character pullbacks. Product
telescoping, Haar invariance under every nonzero integer dilation, and (5)
give

```text
integral Q_E>alpha-9epsilon_N,

integral Q_R>beta-7epsilon_N.                       (17)
```

For every `k>=0`, a second product telescope gives

```text
|integral Q_E(x)Q_R(T^k x) dx
  -measure(E_j intersection T^(-k)R_j)|
 <16epsilon_N.                                      (18)
```

Put

```text
gamma_N
 =(alpha-9epsilon_N)(beta-7epsilon_N)-16epsilon_N.
```

Whenever

```text
alpha-9epsilon_N>0,
beta-7epsilon_N>0,

gamma_N>0,                                          (19)
```

every time satisfying

```text
measure(E_j intersection T^(-k)R_j)<=gamma_N
```

obeys

```text
(integral Q_E)(integral Q_R)
  >integral Q_E(x)Q_R(T^k x) dx.                    (20)
```

If the two nonconstant Fourier supports did not meet after reflection,
orthogonality would make the two sides of (20) equal. Hence there are
nonzero integer frequencies `A,B` with

```text
A+13^k B=0,                                         (21)
```

where `A` belongs to the Fourier support of `Q_E` and `B` belongs to the
Fourier support of `Q_R`.

Expand one nonzero convolution summand in each block. With

```text
w=(H,q_1,...,q_5,c_j,c_h,c_l),
v=(H,q_1,...,q_5,c_j),                              (22)
```

there are vectors

```text
a in [-H_N,H_N]^9 intersection Z^9,
b in [-H_N,H_N]^7 intersection Z^7                 (23)
```

such that

```text
A=a.w!=0,              B=b.v!=0,

a.w+13^k b.v=0.                                      (24)
```

This is a genuinely crossing certificate. A relation supported in only
one block cannot certify a null time.

## 3. Exact branch degrees and coefficient alphabets

For a strict row take

```text
N_strict=33810,
H_strict=2N_strict-2=67618.                         (25)
```

The two positive lower factors and exact crossing margin are

```text
alpha_strict-9epsilon_N
 =4766997229/191198187180,

beta-7epsilon_N
 =117991/4144140,

(alpha_strict-9epsilon_N)(beta-7epsilon_N)-16epsilon_N
 =12331416859/792352055420125200
 >0.                                                (26)
```

This `N_strict` is minimal for the displayed coarse error ledger. At
`N=33809` the same margin is

```text
-91689987875828/15286538167789662548775<0.          (27)
```

For a repeated-first row take

```text
N_repeat=96570,
H_repeat=2N_repeat-2=193138.                        (28)
```

Here the two positive lower factors and crossing margin are

```text
alpha_repeat-9epsilon_N
 =11044460029/1274258725740,

beta-7epsilon_N
 =5543557/193333140,

(alpha_repeat-9epsilon_N)(beta-7epsilon_N)-16epsilon_N
 =10450627633/246356440619713023600
 >0.                                                (29)
```

Again the degree is certificate-minimal for this ledger: at `N=96569`
the margin is

```text
-21997372328/8518227246964664775<0.                 (30)
```

These are certificate minimalities for (16)--(19), not optimal recurrence
constants.

There is a useful exact reduction in each coefficient alphabet. Every
nonzero Fourier coefficient of each interval in (10) vanishes when its
index is divisible by seven. Convolution with `Q_N` cannot create a
missing coordinate mode. Thus every coordinate of `a` and `b` lies in

```text
Xi_H={0} union {m in Z:0<|m|<=H, 7 does not divide m}.    (31)
```

For the strict branch,

```text
floor(67618/7)=9659,

#Xi_(H_strict)
 =2(67618)+1-2(9659)
 =115919=:Q_strict.                                 (32)
```

For the repeated branch,

```text
floor(193138/7)=27591,

#Xi_(H_repeat)
 =2(193138)+1-2(27591)
 =331095=:Q_repeat.                                 (33)
```

For either branch there are at most

```text
(Q_branch^9-1)(Q_branch^7-1)                        (34)
```

ordered coefficient pairs with both blocks nonzero. Simultaneous sign
reversal does not change a certificate, and it acts freely on these pairs.
The number of signed certificate classes is therefore at most the
corresponding integer `M_branch` in (7).

## 4. One certificate class belongs to at most one time

This is the decisive injection. Suppose one signed pair `(a,b)` certified
two distinct times `k<l`. Reversing both signs in one occurrence if
necessary, equations (24) would give

```text
a.w+13^k b.v=0,
a.w+13^l b.v=0.                                     (35)
```

Subtracting,

```text
(13^l-13^k)b.v=0.                                   (36)
```

The integer factor in parentheses is nonzero, while the crossing property
in (24) says `b.v=B!=0`. This is impossible. Thus the nonempty certificate
sets attached by Section 2 to distinct times of handoff mass at most
`gamma_branch` occupy disjoint signed classes. Equations (32)--(34) prove
the stronger bounds in (7).

Notice what the argument does not use: no upper bound on

```text
H+sum_i q_i+sum_h c_h
```

and no boundary complexity of `E_j` or `R_j`. This supplies the
coefficient-independent horizon missing from the BV route, at the price of
a much smaller quantitative mass floor and an enormous finite certificate
count.

## 5. The endpoint is a genuine blocker handoff

At a source point of `E_j`, the guard is active, all five unit masks are
unavailable, and `c_j` is the unique available blocker. At a target point
of `R_j`, the guard is active, all five unit masks and `c_j` are
unavailable. The global cover (2) therefore forces at least one of the
other two blockers.

Consequently every positive set in (9) is an ancestry-compatible

```text
c_j -> {c_h,c_l}.                                   (37)
```

handoff. Partitioning the target by a measurable first available blocker
shows that one named target label receives mass strictly larger than
`gamma_branch/2`. The theorem does not say that the orbit's first switch
occurs at `k`; it identifies the source and endpoint service labels.

Combining Sections 4--5 gives a coefficient-independent bounded
post-expiration occurrence on both branches. On the strict branch,
THM-2288 remains complementary: for its possibly different shallow label,
it gives the explicit mass

```text
14772292477/132461154025200
```

at every sufficiently large time. The two conclusions are complementary.
THM-2286 gives the analogous larger, coefficient-dependent delayed floor on
the repeated-first branch.

## 6. Sharp stopping boundaries and hostile controls

Three limitations are load-bearing.

1. **No exact-expiration conclusion.** The proof counts all possible
   crossing certificates. It gives no reason why the first one should
   occur at `lambda_j+1`.

2. **A bounded certificate can live arbitrarily late.** For any `K>=1`,
   the two thirteen-unit integers

   ```text
   h=1,                 q=13^K+1                    (38)
   ```

   obey the height-one crossing identity

   ```text
   (-h+q)+13^K(-h)=0,                               (39)
   ```

   with both partial sums nonzero. This is an arithmetic hostile control,
   not a scalar cover or LRC counterexample. It proves that the size of one
   certificate alone cannot bound its time uniformly; the finite-type
   injection across *all* exceptional times is essential.

3. **The uniform mass and horizon are certificate-scale.** The explicit
   floors in (6a) are positive but much smaller than the BV delayed floors,
   while the certificate counts in (7) are enormous. Equations (27) and
   (30) show only minimality for the coarse `3/(2N)` Jackson ledger, not
   optimality of either constant.

The theorem applies to all `150` strict profiles and all `15`
repeated-first profiles `(1,1,c)`. It excludes no profile and does not prove
LRC(14). Its new structural content is that no first-depth-one scalar
counterexample can suppress a fixed quantitative blocker-only ancestry
return at infinitely many times: outside at most `M_branch` uniformly
bounded spectral exceptions, every time has mass greater than
`gamma_branch`.

## 7. Exact verification

The primary companion verifies both labelled source floors, the common
target floor, both exact crossing inequalities, both minimal `N` values for
the coarse Jackson ledger, the degrees, the non-seven coefficient
alphabets, both signed certificate counts, and delayed arithmetic controls.
The independent referee uses cleared integer numerators and literal
enumerations of both coefficient alphabets.
Reproduce with

```bash
python3 04-computation/lrc14_finite_spectral_null_handoff_thm2289.py
python3 -O 04-computation/lrc14_finite_spectral_null_handoff_thm2289.py

python3 04-computation/lrc14_finite_spectral_null_handoff_referee_thm2289.py
python3 -O 04-computation/lrc14_finite_spectral_null_handoff_referee_thm2289.py
```

Normal and optimized transcripts are byte-identical to their stored
outputs. QED.
