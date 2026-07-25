---
id: THM-2275
title: "Mixed scalar relation and guard-blocker crossing"
status: >
  PROVED + VERIFIED-EXACT. Every scalar five-unit/three-blocker cover has
  a nonzero relation of coefficient height at most 20 on its nine scalar
  coefficients. It also has a genuine relation of height at most 462
  crossing the six guard/unit coordinates against the three blockers,
  with every nonzero coefficient prime to seven. Through THM-2203 these
  become original-row relations supported on the fixed nine-coordinate
  section, of heights at most 40 and 924. More strongly, every one of the
  120 interior first-depth-one profiles has scalar relation rank at least
  two by height 9841: THM-2266 supplies a support-two relation of height at
  most 9841, and a trigger-adapted safe-factor cut supplies an independent
  crossing relation of height at most 708 or 2116. The two relations lift
  to the fixed original-row section with uniform height at most 19682.
  No scalar profile is excluded and LRC(14) remains open.
source: codex-2026-07-25-mixed-scalar-relation
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-2085-explicit-height-57-rank-seven-selberg-gate
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2145-two-block-spectral-crossing-and-6-plus-7-carry
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
related:
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
script: 04-computation/lrc14_mixed_scalar_relation_crossing_thm2275.py
output: 05-knowledge/results/lrc14_mixed_scalar_relation_crossing_thm2275.out
script_sha256: 703b6daead45608a42783c713abfce3065c860cb44cea90cb6ccc6bcb375a4c4
output_sha256: 6d1699f9c2678892c92bc0906a99f60f36e6de634bda81ae48f6471fe206c6ed
hash_basis: working-tree bytes (LF)
---

# THM-2275 -- mixed scalar relations after the quotient

Put

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

Assume the live scalar hypotheses of THM-2198:

```text
H,q_1,...,q_5 are positive thirteen-units;
H is odd;
q_1,...,q_5 are pairwise distinct;
c_1,c_2,c_3 are distinct positive multiples of thirteen; (2)
```

and, outside a null set,

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j).            (3)
```

Write

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3) in Z_(>0)^9.         (4)
```

The theorem proves two complementary bounded relations on `w_*`.
The first is smaller but need not cross the natural cut:

```text
there is 0!=r in Z^9 with
r.w_*=0,                         ||r||_infinity<=20. (5)
```

The second has labelled crossing:

```text
a_H H+sum_(i=1)^5 a_i q_i
   =-sum_(j=1)^3 b_j c_j !=0,                        (6)

|a_H|,|a_i|,|b_j|<=462,                              (7)
```

where every displayed coefficient is either zero or prime to seven.

For every interior first-depth-one profile, THM-2266 supplies a primitive
support-two relation `p` of height at most `9841`. A safe-factor cut adapted
to the support of `p` supplies a second relation `s`, independent of `p`,
with

```text
||s||_infinity<=708       if supp(p)={H,c_1},
||s||_infinity<=2116      if supp(p)={q_i,c_1}.       (7a)
```

Consequently every such scalar row has relation rank at least two by
height `9841`.

The source is the scalar cover itself. The proof does not attempt to
project an arbitrary relation supplied on the original thirteen-speed row.
That distinction is what makes (5)--(6) survive the scalar quotient.

## 1. The exact original-row section

THM-2203 identifies the nine scalar labels with a fixed coordinate section
of the original thirteen-speed row. In its notation the three blocker
coefficients are `c_j=13s_j`, so the selected coordinates are

```text
S_I=(8H,16q_1,...,16q_5,16c_1,16c_2,16c_3).         (8)
```

The four forgotten original coordinates are

```text
(2h_0,4h_1,x,y).                                     (9)
```

Consequently a scalar relation

```text
r_H H+sum_i r_iq_i+sum_j r_*j c_j=0                 (10)
```

lifts to the original row, supported on `I`, by the integral coefficient
vector

```text
(2r_H,r_1,...,r_5,r_*1,r_*2,r_*3).                  (11)
```

Indeed the original-row sum in (11) is sixteen times (10). Thus (5)
immediately gives an original-row relation supported on `I` of height at
most

```text
40.                                                  (12)
```

Likewise (6) lifts to a supported original-row relation of height at most

```text
924,                                                 (13)
```

and its two partial sums are `+16A` and `-16A`, where the nonzero carry
`A` is the left side of (6). The diagonal section therefore preserves
crossing as well as linear independence.

The adaptive relations in (7a) lift with heights at most `1416` and
`4232`, respectively. The THM-2266 pair relation lifts with height at most
`19682` in the guard-owner case and at most `9841` in the owner-unit case.
Thus the supported original-row relation lattice has rank at least two by
the uniform height `19682`.

## 2. A mixed-length signed Selberg tensor

The joint safe event in (3) has nine interval coordinates. The guard-safe
coordinate is the interval

```text
{x:||x||>1/7},
```

of length `5/7`, while each terminal-safe coordinate is

```text
{x:||x||>=1/14},
```

of length `6/7`. Endpoint choices are immaterial for Haar measure.

THM-2085 proves the following signed tensor statement for arbitrary circle
intervals. For intervals of lengths `alpha_1,...,alpha_k`, take degree-`K`
Vaaler lower and upper polynomials with

```text
epsilon=1/(K+1).
```

If no nonzero relation on the orbit has all coordinate coefficients at
most `K`, the tensor minorant integrates exactly to

```text
B_K(alpha_1,...,alpha_k)
 =product_l(alpha_l+epsilon)
   *[1-sum_l 2epsilon/(alpha_l+epsilon)].             (14)
```

Apply (14) with

```text
K=20,                  epsilon=1/21,
alpha_0=5/7,           alpha_1=...=alpha_8=6/7.      (15)
```

The bracket is

```text
1-2/16-8(2/19)=5/152,                                (16)
```

so the exact constant is

```text
B_20
 =(16/21)(19/21)^8(5/152)
 =8938717390/794280046581
 >0.                                                  (17)
```

If (5) failed, THM-2085 would give the null joint safe event in (3)
measure at least (17), a contradiction. This proves (5).

The adjacent degree is a sharp boundary for this particular tensor
certificate. At `K=19`, `epsilon=1/20`, the defect budget and bracket are

```text
2epsilon/(5/7+epsilon)
 +8*2epsilon/(6/7+epsilon)
 =13762/13589,

1-13762/13589=-173/13589<0.                          (18)
```

Each budget term decreases when its coordinate degree increases. Hence no
anisotropic Vaaler box whose maximum degree is at most nineteen makes this
same signed tensor positive. This is certificate minimality, not a claim
that twenty is the optimal true relation height or the best possible
trigonometric construction.

## 3. The mixed-interval form of THM-2145

We isolate the exact extension of THM-2145 used for the crossing relation.
Let labelled coordinates be partitioned into blocks `F,E`. For coordinate
`l`, suppose

```text
0<=p_l<=1,
||p_l-1_(J_l)||_1<=eta_l,
supp Fourier(p_l) subset [-K_l,K_l].                 (19)
```

Let `G_F,G_E` be the corresponding exact block-safe events and put

```text
eta_F=sum_(l in F)eta_l,
eta_E=sum_(l in E)eta_l.                             (20)
```

If

```text
measure(G_F)>=alpha,       measure(G_E)>=beta,
measure(G_F intersection G_E)=0,                    (21)

alpha-eta_F>0,             beta-eta_E>0,

(alpha-eta_F)(beta-eta_E)>eta_F+eta_E,               (22)
```

then there is a relation crossing the two blocks, with coordinate `l`
bounded by `K_l`.

The proof is the proof of THM-2145 Section 1 with coordinate-dependent
intervals and errors. Put

```text
P_F=product_(l in F)p_l(v_l t),
P_E=product_(l in E)p_l(v_l t).                      (23)
```

Product telescoping gives

```text
integral P_F>=alpha-eta_F,
integral P_E>=beta-eta_E,
integral P_FP_E<=eta_F+eta_E.                        (24)
```

If the Fourier supports met after reflection only at zero,

```text
supp Fourier(P_F) intersection
       -supp Fourier(P_E)={0},                       (25)
```

then Fourier orthogonality would give

```text
integral P_FP_E=(integral P_F)(integral P_E),        (26)
```

contradicting (22)--(24). Therefore a nonzero frequency `k` belongs to
the first support and `-k` to the second. A nonzero convolution summand on
each side supplies a coefficient vector with both block sums equal to
`+k` and `-k`. This proves the mixed lemma. It asserts a common frequency,
not independence of the resulting relation.

## 4. The scalar `6+3` crossing

Use the normalized squared Fejer kernel

```text
F_N(t)=(1/N)(sin(pi Nt)/sin(pi t))^2,
J_N=F_N^2/integral F_N^2.                            (27)
```

For any circle interval `J`, put

```text
p_(N,J)=J_N*1_J.                                     (28)
```

Then

```text
0<=p_(N,J)<=1,
supp Fourier(p_(N,J)) subset [-2N+2,2N-2],
||p_(N,J)-1_J||_1<3/(2N).                           (29)
```

The last bound is independent of the interval length. Translation by `y`
changes the indicator of a circle interval in `L^1` by at most
`2||y||`; for a long or wrapping interval this is the same statement
applied to its complement. For completeness, Parseval and the elementary
pointwise estimate give

```text
integral F_N^2=(2N^2+1)/(3N)>2N/3,

F_N(y)<=min(N,1/(4N||y||^2)).                        (30)
```

Splitting the numerator integral at `1/(2N)` yields

```text
integral ||y||F_N(y)^2dy<1/2.
```

After normalizing by the first line of (30),

```text
integral ||y||J_N(y)dy<3/(4N).                       (30a)
```

and averaging the translation bound proves (29).

Take the six-coordinate block

```text
F=C_H minus union_(i=1)^5D_(q_i).                    (31)
```

THM-2137 gives its sharp uniform floor

```text
measure(F)>=delta_5:=961/6930.                       (32)
```

Take the three-coordinate blocker-safe block

```text
E=(R/Z) minus union_(j=1)^3D_(c_j).                  (33)
```

THM-1166's sharp three-comb union bound gives

```text
measure(E)>=1-36/91=55/91.                           (34)
```

The scalar cover (3) is exactly

```text
measure(F intersection E)=0.                         (35)
```

At

```text
N=232,                   2N-2=462,                   (36)
```

the six-coordinate error, three-coordinate error, and total joint error
from (29) are

```text
eta_F=9/232,            eta_E=9/464,
eta_F+eta_E=27/464.                                  (37)
```

The two lower factors are

```text
961/6930-9/232=80291/803880,
55/91-9/464=24701/42224,                             (38)
```

and the exact crossing margin is

```text
(961/6930-9/232)(55/91-9/464)-27/464
 =8134831/33943029120
 >0.                                                  (39)
```

The mixed lemma now proves (6)--(7).

There is also an exact coefficient congruence. The Fourier coefficient of
`p_(N,J)` at frequency `m` is the product of the corresponding squared-
Fejer coefficient and the interval coefficient. The first is positive
through the support. For an interval of length `d/7`, with `d=5` or `6`,
the second vanishes at a nonzero integer frequency precisely when

```text
7 divides m.                                         (40)
```

Thus the nonzero convolution summands in the proof of the mixed lemma use
only zero coefficients and seven-units. This proves the last assertion
following (7).

The symmetric bandwidth boundary is again certificate-relative. At
`N=231`, the same three lower/error terms are

```text
691/6930,                  1171/2002,       9/154,
```

and the resulting margin is

```text
-1649/13873860<0.                                    (41)
```

The first lower factor is nonpositive for `N<=64`. For `N>=65`, the margin
in this ledger is increasing with `N`: as a quadratic in `x=1/N`, its
derivative is negative throughout `0<x<=1/65`. Hence `232` is the first
successful common bandwidth for the safe-mass and `3/(2N)` error
certificate. Equation (41) does not say that a sharper Jackson moment,
another positive kernel, or the actual row fails at `N=231`.

More generally, bandwidth `N_U` on the six guard/unit coordinates and
`N_D` on the three blocker coordinates works whenever both lower factors
are positive and

```text
(961/6930-9/N_U)(55/91-9/(2N_D))
     >9/N_U+9/(2N_D).                                (42)
```

It gives coefficient bounds `2N_U-2` and `2N_D-2` on the two sides.

## 5. Why the ambient theorems do not already imply this

THM-2144 gives a height-twenty-nine relation somewhere on the original
thirteen coordinates. THM-2145 gives a crossing relation across any
chosen original `6+7` cut. Neither conclusion forces the four coefficients
on (9) to vanish. An ambient relation containing those coordinates has no
relation-valued image under the scalar projection.

Linear elimination makes the missing rank precise. Cancelling four
forgotten coordinates from arbitrary ambient relations requires at least
five independent rows before a nonzero supported combination is forced.
THM-2144 guarantees only one row at height twenty-nine, and its first
rank harvest guarantees only rank two at height 105 unless a subset-sum
relation occurs. THM-2145's every-cut use supplies no rank by itself,
because the relation chosen for different cuts need not change.

Here is an exact hostile control for that last implication. On

```text
(1,2,4,...,2^11,2^12-1),                             (43)
```

the single relation with coefficients `+1` on the first twelve entries and
`-1` on the last crosses every six-versus-seven cut. A six-set not
containing the last coordinate has positive partial sum. If it contains
the last coordinate, the other five powers sum to at most

```text
2^7+2^8+2^9+2^10+2^11=3968<4095,                    (44)
```

so its partial sum is negative. Thus one rank-one relation can witness
every cut.

The direct scalar proofs in Sections 2 and 4 avoid this quotient loss and
then lift their outputs back through the fixed diagonal section (8). They
do not claim that the ambient relation theorems were false or weaker in
their stated domains.

## 6. Trigger-adapted cuts give unconditional rank two

Consider any of the 120 interior first-depth-one profiles. Write

```text
c_1=13u_1,
(lambda_1,lambda_2,lambda_3)=(1,b,c),
3<=b<=c-2.                                           (45)
```

THM-2266 proves that one of the six reduced pairs

```text
(H,u_1),             (u_1,q_i), 1<=i<=5             (46)
```

has coprime product at most `757`.

For the first pair, put

```text
g=gcd(H,u_1),       a=H/g,       d=u_1/g.            (47)
```

The corresponding relation on the actual scalar row (4) is

```text
13d H-a c_1=0.                                      (48)
```

For the pair `(u_1,q_i)`, put

```text
g=gcd(u_1,q_i),     a=u_1/g,     d=q_i/g.            (49)
```

Its actual scalar relation is

```text
d c_1-13a q_i=0.                                    (50)
```

In both cases

```text
gcd(a,d)=1,          13 does not divide ad,          (51)
```

because all coefficients in (46) are thirteen-units. Hence (48)--(50) are
primitive integer relation vectors. Their heights are respectively

```text
max(13d,a),                max(d,13a),               (52)
```

and are at most

```text
13*757=9841.                                         (53)
```

Let `p` denote this primitive support-two relation. We choose the
factorization of the null joint safe event according to `supp(p)`.

### 6.1 The guard-owner trigger

Suppose `supp(p)={H,c_1}`. Put those two safe coordinates in the first
block:

```text
A={H,c_1},
B={q_1,...,q_5,c_2,c_3}.                             (54)
```

The first block-safe event is the intersection of an interval of mass
`5/7` and one of mass `6/7`, so the union bound gives

```text
measure(G_A)>=4/7.                                   (55)
```

The seven speeds in `B` are distinct positive integers: the five `q_i`
are pairwise distinct thirteen-units, whereas `c_2,c_3` are distinct
multiples of thirteen. THM-1221 therefore gives

```text
measure(G_B)>=15/154.                                (56)
```

With the squared-Fejer approximants of Section 4, the two block errors and
their sum are

```text
eta_A=3/N,       eta_B=21/(2N),
eta_A+eta_B=27/(2N).                                 (57)
```

At `N=355`, the exact margin in (22) is

```text
(4/7-3/355)(15/154-21/710)-27/710
 =21177/135854950
 >0.                                                  (58)
```

The adjacent value fails:

```text
margin at N=354 =-3/15010072<0.                      (59)
```

An exact scan of the rational margin shows that `355` is the first
bandwidth accepted by this ledger. The mixed crossing lemma therefore
gives a scalar relation `s` of height

```text
2N-2=708                                             (60)
```

whose `A`-partial sum is nonzero. But `p` is supported inside `A` and its
`A`-partial sum is `p.w_*=0`. Hence `s` cannot be a rational multiple of
`p`.

### 6.2 The owner-unit trigger

Suppose instead that `supp(p)={q_i,c_1}`. Use

```text
A={q_i,c_1},
B={H,q_1,...,q_(i-1),q_(i+1),...,q_5,c_2,c_3}.      (61)
```

THM-1166's sharp pair overlap gives

```text
measure(G_A)>=1-2/7+1/91=66/91.                      (62)
```

The second block consists of the odd guard and six distinct positive
terminal coefficients. THM-2137's six-plus-one scalar-tail estimate is
uniform in those six coefficients and gives

```text
measure(G_B)>=delta_6=191/6930.                      (63)
```

The error ledger is again (57). At `N=1059`,

```text
(66/91-3/1059)(191/6930-21/2118)-27/2118
 =44021/78582173670
 >0,                                                  (64)
```

whereas

```text
margin at N=1058 =-861505/47060301288<0.             (65)
```

The exact rational scan makes `1059` the first accepted bandwidth for
this ledger. It produces a relation `s` of height

```text
2N-2=2116                                            (66)
```

with nonzero `A`-partial sum. As in Section 6.1, the support-two relation
`p` is internal to `A` and has zero `A`-partial sum. The two relations are
therefore independent.

Together the two cases prove

```text
dim_Q span{m in Z^9:m.w_*=0, ||m||_infinity<=9841}
 >=2.                                                 (67)
```

for every interior profile. The fixed diagonal lift in Section 1 preserves
the two-dimensional span and gives supported original-row rank at least
two by height `19682`.

### 6.3 The old dependency rays are diagnostic only

For comparison, let `r` be the height-twenty relation from (5). If `r`
is proportional to the primitive `p`, Bezout's identity forces `r=np`
with nonzero integer `n`. Equations (48)--(52) then give exactly the nine
guard-owner rays

```text
H=a u_1,
a in {1,3,5,7,9,11,15,17,19},                       (68)
```

and the nineteen owner-unit rays

```text
q_i=d u_1,
d in {1,2,3,4,5,6,7,8,9,10,11,12,
      14,15,16,17,18,19,20}.                         (69)
```

These rays diagnose when the two *smallest previously selected* relations
can coincide. They are not exceptions to (67): the adaptive crossing
relation has nonzero partial sum on the support block and bypasses all of
them.

## 7. Exact verification and scope

The companion uses only integer and `Fraction` arithmetic. It checks:

```text
the exact degree-19 and degree-20 tensor budgets;
the positive tensor constant (17);
the N=231 and N=232 crossing margins;
the first common-bandwidth boundary under this ledger;
the N=354/355 guard-owner adaptive margins and first boundary;
the N=1058/1059 owner-unit adaptive margins and first boundary;
all Fourier residues through height 462;
the 3,643/1,822 generic THM-2266 pair census;
the 3,279/1,640 thirteen-unit subcensus;
primitivity of every actual pair lift;
the diagnostic nine rays (68) and nineteen rays (69);
all scalar and fixed-section lift heights.            (70)
```

Reproduce with

```bash
python 04-computation/lrc14_mixed_scalar_relation_crossing_thm2275.py
python -O 04-computation/lrc14_mixed_scalar_relation_crossing_thm2275.py
```

The two modes produce identical stdout. After the platform newline is
normalized to `LF`, both reproduce the stored transcript byte for byte.

The theorem proves bounded relation structure, not an exclusion. In
particular:

1. the height-twenty relation need not cross the unit/blocker cut;
2. the height-462 crossing relation may be a support-two
   commensurability and need not raise relation rank;
3. the adaptive relation proves rank two but does not select a blocker
   owner or retain a root digit,
   endpoint current, or post-expiration target;
4. THM-2261's expiration-surjectivity obstruction and THM-2263's missing
   labelled target inclusion remain;
5. none of the 165 scalar profiles is removed.

The gain is a finite, quotient-faithful relation sidecar: one small
supported relation on every scalar row, one bounded nonzero unit-to-blocker
carry, and unconditional bounded relation rank two on all 120 interior
profiles. LRC(14) remains open. QED.
