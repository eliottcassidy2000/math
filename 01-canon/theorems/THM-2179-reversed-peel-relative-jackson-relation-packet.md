---
id: THM-2179
title: "Reversed-peel relative Jackson relation packet at defect at least seven"
status: >
  PROVED + VERIFIED-EXACT. At the stronger radius 3/41, every thirteen-speed
  row of AP defect at least seven either has safe measure greater than
  478970390236831/39525379884148950000 or has a body-touching integer
  relation of coefficient height at most 140, with every nonzero coefficient
  a 41-unit. The proof lifts all small-core factors onto one torus coordinate
  and every noncore factor onto its own coordinate, thereby retaining every
  internal small-core relation and all higher overlap cancellation. On the
  named hostile row all six scalar peel covariances are positive, so replacing
  their absolute values by their signed sum gives exactly the same negative
  level-one bound; the higher-overlap packet is load-bearing. The relation
  alternative is not a closure of defect at least seven or of LRC(14).
source: codex-2026-07-24-reversed-peel-relation-packet
depends_on:
  - THM-2145
  - THM-2146
related:
  - THM-2054
  - THM-2086
  - THM-735
  - THM-2171
script: 04-computation/lrc14_reversed_peel_relation_packet_thm2179.py
output: 05-knowledge/results/lrc14_reversed_peel_relation_packet_thm2179.out
script_sha256: 846b94e89a661367c98ff03b156cdf5af32a09e15a359723b914b48fe1c35d66
output_sha256: 7863caa0824646b18cf4d894b59e9c342f7ab49bba23f6f172d58e476d4db2b6
hash_basis: working-tree bytes (LF)
---

# THM-2179 -- reversed-peel relative Jackson relation packet

Put

```text
h=3/41,
G={x in R/Z:||x||>=h},
G_h(A)={t:at in G for every a in A}.                   (1)
```

The radius in (1) is strictly larger than `1/14`. A positive-measure exit
here therefore proves the ordinary lonely-runner target for that row, but a
failure to obtain this stronger exit is not a counterexample to LRC(14).

## 1. The defect-seven relation-packet alternative

Let `V` be a set of thirteen distinct positive integers. Split it as

```text
F=V intersection {1,...,13},
E=V\F,
j=|F|,                  d=|E|=13-j.                   (2)
```

Thus the AP defect of `V` is `d`. Suppose `d>=7`, equivalently `j<=6`.

> **Theorem.** At least one of the following holds:
>
> 1. the stronger safe set has the explicit positive-measure exit
>
>    ```text
>    measure(G_h(V))
>      >478970390236831/39525379884148950000
>      >0;                                             (3)
>    ```
>
> 2. there are integers `(a_f)_(f in F)` and `(b_e)_(e in E)` such that
>
>    ```text
>    sum_(f in F) a_f f+sum_(e in E)b_e e=0,
>    some b_e!=0,
>    |a_f|,|b_e|<=140,                                 (4)
>    ```
>
>    and every nonzero coefficient in (4) is not divisible by `41`.

We call (4) **body-touching** or **large-involving**. It need not genuinely
cross the cut. If

```text
sum_(e in E)b_e e=0,
```

then the body part alone is a nonzero bounded relation. Otherwise both block
partial sums are nonzero and (4) is genuinely crossing. This distinction is
load-bearing: the theorem does not claim the second case always occurs.

## 2. A positive Jackson approximant at radius `3/41`

Use THM-2145's normalized Jackson kernel

```text
J_N=F_N^2/integral F_N^2
```

and define, now for the interval `G` in (1),

```text
q_N=J_N*1_G.                                          (5)
```

The kernel is nonnegative and has integral one, so

```text
0<=q_N<=1,
integral q_N=35/41,
degree(q_N)<=2N-2.                                    (6)
```

Translation of one circular interval changes its indicator in `L1` by at
most `2||x||`. Consequently the proof of THM-2145 gives

```text
||q_N-1_G||_1<=eta_N,
eta_N=2 integral ||x||J_N(x)dx.                       (7)
```

Take

```text
N=71,                    H=2N-2=140.                  (8)
```

For the integer Jackson coefficients `C_k` of THM-2145, exact arithmetic
gives

```text
C_0=238631,
sum_(1<=k<=139, k odd) C_k/k^2>290903.                (9)
```

The circle-distance Fourier series and `pi<355/113` therefore give

```text
eta_71
 <357148519/60146943550
 <297/50000.                                          (10)
```

The nonzero Fourier coefficients are also explicit. For `k!=0`,

```text
Fourier(1_G,k)=-sin(6 pi k/41)/(pi k).                (11)
```

Every Jackson multiplier is positive throughout `[-140,140]`. Since `41`
is coprime to `6`,

```text
Fourier(q_71,k)!=0
iff
k=0
or
0<|k|<=140 and 41 does not divide k.                  (12)
```

This is the source of the `41`-unit assertion in (4).

## 3. The relative lift keeps the entire small block

Label `E=(e_1,...,e_d)` and work on `T^(d+1)` with coordinates
`(x_0,x_1,...,x_d)`. Assign the characters

```text
chi_f=(f,0,...,0)              for f in F,
chi_(e_i)=the i-th body basis vector,
lambda=(1,e_1,...,e_d).                               (13)
```

Then

```text
chi_v dot lambda=v
```

for every `v in V`. The actual line product and lifted product are

```text
Q_line(t)=product_(v in V)q_71(vt),

Q_lift(x)=product_(f in F)q_71(fx_0)
          product_(i=1)^d q_71(x_i).                  (14)
```

The crucial point is that all small factors still share `x_0`. Thus the lift
retains their complete relation lattice, exact union geometry, and every
higher overlap. Only the noncore factors are decorrelated.

Put

```text
I=integral_T Q_line,
A_F=integral_T product_(f in F)q_71(ft),
J=integral_(T^(d+1))Q_lift
  =A_F(35/41)^d.                                      (15)
```

The product telescope, (7), and Haar invariance give

```text
|I-measure(G_h(V))|<=13 eta_71,                       (16)
|A_F-measure(G_h(F))|<=j eta_71.                      (17)
```

Now compare the finite Fourier expansions of `I` and `J`. A scalar
constant-term tuple satisfies

```text
sum_(f in F)a_f f+sum_(i=1)^d b_i e_i=0,              (18)
```

whereas a lifted vector constant-term tuple satisfies

```text
sum_(f in F)a_f f=0,
b_1=...=b_d=0.                                        (19)
```

Every tuple in (19) is a tuple in (18). If there is no body-touching relation
of the coefficient shape in (4), (12) says that the two finite index sets
are identical. Hence

```text
I=J.                                                   (20)
```

Equivalently, if `I!=J`, some scalar-but-not-vector Fourier tuple supplies
exactly (4). Notice that zero modes are allowed. This is why arbitrary
internal relations in the small block cause no alias and no hypothesis.

## 4. The exact defect-sensitive margin

THM-2146 gives the exact minimum

```text
B_j=min_({F subset {1,...,13}:|F|=j}) measure(G_h(F))
```

for `0<=j<=6`:

| `j` | `B_j` |
|---:|---:|
| 0 | `1` |
| 1 | `35/41` |
| 2 | `59/82` |
| 3 | `1615/2706` |
| 4 | `239/492` |
| 5 | `2729/7380` |
| 6 | `153101/568260` |

Equations (15)--(17), under the no-relation equality (20), imply

```text
measure(G_h(V))
 >=B_j(35/41)^d
   -(13+j(35/41)^d)eta_71.                            (21)
```

The seven exact comparisons `7<=d<=13`, `j=13-d`, give

```text
B_j(35/41)^d
 >=B_6(35/41)^7
 =281440305453125/3162030390731916.                   (22)
```

Also

```text
13+j(35/41)^d<=13+6(35/41)^7.                         (23)
```

Insert the strict cap (10) into (21)--(23). The resulting uniform margin is

```text
B_6(35/41)^7
 -(13+6(35/41)^7)(297/50000)

=478970390236831/39525379884148950000
>0.                                                   (24)
```

This proves (3) whenever (4) is absent, and therefore proves the theorem.
The adjacent value `N=70` fails the same exact global margin ledger. Thus
`N=71`, or coefficient height `140`, is certificate-minimal for this
particular Jackson/core-floor calculation; no optimality for the true
relation height is asserted. QED.

## 5. Why the six signed scalar discrepancies cannot repair the hostile row

The proposed reversed-peel repair was to replace six absolute discrepancies
by their signed sum. The named hostile row shows that this does not touch the
actual defect.

Take

```text
F=(1,2,3,4,6,8),
E=(95,163,187,206,208,214,332).                       (25)
```

For `f in F`, put

```text
epsilon_f
 =measure(G_h(E) intersection D_f)
  -(6/41)measure(G_h(E)),                             (26)

D_f={t:||ft||<3/41}.
```

Exact rational interval arithmetic gives

```text
measure(G_h(E))
 =7521335361151/22863470734060,

measure(G_h(F))=20/41,

measure(G_h(E union F))
 =470973614624713/3388366362787692.                   (27)
```

Every one of the six signed discrepancies is positive:

| `f` | exact `epsilon_f` |
|---:|---:|
| 1 | `484544917897629/115769184061912810` |
| 2 | `123423920266531/13619904007283860` |
| 3 | `14683462055431813/1389230208742953720` |
| 4 | `4106490555883167/463076736247651240` |
| 6 | `3427539824959769/231538368123825620` |
| 8 | `6524856171419367/463076736247651240` |

Consequently

```text
sum_f epsilon_f=sum_f |epsilon_f|
 =85546520069055739/1389230208742953720.              (28)
```

The separately charged base and its signed level-one bound are

```text
(5/41)measure(G_h(E))
 =7521335361151/187480460019292,

(5/41)measure(G_h(E))-sum_f epsilon_f
 =-727156708364069/33883663627876920
 <0.                                                   (29)
```

Thus even perfect use of the six scalar signs leaves exactly the same
negative number as the absolute-value estimate.

The missing coordinate is the higher-overlap packet. Let `p_k` be the mass
of points in `G_h(E)` lying in exactly `k` of the six danger combs. Then

```text
measure(G_h(E union F))
 =(5/41)measure(G_h(E))-sum_f epsilon_f
   +sum_(k=2)^6 (k-1)p_k,                             (30)
```

and here

```text
sum_(k=2)^6 (k-1)p_k
 =1812297618203733/11294554542625640
 >0.                                                   (31)
```

This correction changes the failed `-0.02146...` level-one estimate into the
true positive mass `0.13899...`. The aggregate safe-mask covariance is

```text
measure(G_h(E union F))
 -measure(G_h(E))measure(G_h(F))
 =-2983319810838331/138923020874295372.               (32)
```

So the repair is not a favorable universal covariance sign either. It is the
relative lift in Section 3: retain the whole small-block Boolean product and
all of its overlaps before taking any Fourier absolute value.

## 6. Consequence, transfer, and exact boundary

The theorem converts the entire defect-`>=7` branch into the exact dichotomy

```text
strong 3/41 positive-measure exit
or
height-140 body-touching relation packet with 41-unit coefficients. (33)
```

This is a substantial reduction of the unbounded-height branch, but it is
not its closure:

1. the relation may have support two or lie wholly inside the noncore body;
2. it need not be independent of relations already known;
3. bounded coefficients do not bound the speed magnitudes;
4. THM-2171 gives a finite algebraic representative for fixed relation
   systems, but its target-mixing example forbids replacing the original row
   without a phase/current argument; and
5. absence of the stronger `3/41` exit is not failure of the `1/14` LRC
   predicate.

The reusable mechanism is broader than the displayed constants. Put all
factors whose internal correlations must be preserved on common torus
coordinates, put the genuinely free factors on separate coordinates, and
compare line and lifted constant terms only after positive whole-product
smoothing. A discrepancy then produces a bounded relation that touches a
decorrelated coordinate. This is the relative-character idea of THM-2054
combined with THM-2145's sharper Jackson kernel and THM-2146's exact
small-core floor.

## 7. Exact referee

The companion performs the following exact checks:

1. direct convolution and the closed Jackson formula agree on all `141`
   coefficients at `N=71`, and every multiplier in the support is positive;
2. the odd Jackson sum, both rational eta caps, the seven defect margins, the
   displayed fraction (24), and the hostile adjacent `N=70` check all pass;
3. the height-`140` support has exactly `274` signed nonzero `41`-unit modes;
4. the hostile row is evaluated by two independent rational interval
   algorithms--complete boundary cells and merged danger intervals;
5. the six epsilon signs, their sum, the negative level-one bound, the full
   danger-depth distribution, the higher-overlap correction, and the
   aggregate covariance are checked exactly.

Normal and optimized Python executions agree, and after newline normalization
they match the stored LF transcript. An independent proof audit separately
checked the character lift, the line and lift telescope directions, the
body-touching scope, the `N=71` sharpening, and every displayed rational
margin.
