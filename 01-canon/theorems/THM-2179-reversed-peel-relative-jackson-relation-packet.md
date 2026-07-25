---
id: THM-2179
title: "Reversed-peel relative Jackson relation packets at defect at least six"
status: >
  PROVED + VERIFIED-EXACT. At the stronger radius 3/41, relative to any chosen
  dilated AP g*{1,...,13}, a thirteen-speed row of defect at least seven either
  has safe measure greater than
  478970390236831/39525379884148950000 or has an outside-core-touching
  41-unit integer relation of coefficient height at most 140. At defect six,
  the corresponding alternatives have safe margin
  287560991216713/1253243752424235000 and coefficient height 180. The new
  defect-six input is the independently exact seven-core floor
  B_7=39965/211068, uniquely attained by {1,5,7,8,9,11,13}; Jackson N=91
  closes its ledger while N=90 does not. The proof retains every internal
  core relation and all higher overlap cancellation. A two-speed
  down-conversion lemma proves that body endpoints can all have denominator
  divisible by 41 while fixed frequencies 1 through 13 retain nonzero
  limiting Fourier mass; denominator separation cannot repair scalar peeling.
  These body-touching relations need not genuinely cross the core/body cut,
  so neither packet closes LRC(14). Under the ordinary zero-safe premise and
  defect at least seven, combining the height-140 packet with THM-2169 gives
  an anchored rank-two base-41 carrier with at most 27200916 carry pairs;
  saturation makes the same plane universal-radix at a slightly larger count.
source: codex-2026-07-24-reversed-peel-relation-packet
depends_on:
  - THM-2145
  - THM-2146
  - THM-2169
related:
  - THM-2054
  - THM-2086
  - THM-2166
  - THM-735
  - THM-2171
  - THM-2187
script: 04-computation/lrc14_reversed_peel_relation_packet_thm2179.py
output: 05-knowledge/results/lrc14_reversed_peel_relation_packet_thm2179.out
script_sha256: 8d3e990db5d656e7667e402b1d94f09b6e39dcd737f1fb7c6d4c64fc9d219949
output_sha256: 7f5dfef83cad9ff167d521dd7d38f34bf1bad272e0d0c574edc0b32a7702ec7c
hash_basis: working-tree bytes (LF)
---

# THM-2179 -- reversed-peel relative Jackson relation packets

Put

```text
h=3/41,
G={x in R/Z:||x||>=h},
G_h(A)={t:at in G for every a in A}.                   (1)
```

The radius in (1) is strictly larger than `1/14`. A positive-measure exit
here therefore proves the ordinary lonely-runner target for that row, but a
failure to obtain this stronger exit is not a counterexample to LRC(14).

## 1. The dilated-AP defect-six relation-packet alternatives

Let `V` be a set of thirteen distinct positive integers and fix any integer
scale `g>=1`. Put

```text
A_g={g,2g,...,13g},
F=V intersection A_g,
E=V\F,
C={c in {1,...,13}:gc in F},
j=|F|,                  d=|E|=13-j.                   (2)
```

Every speed lying in the chosen dilated AP is assigned to `F`; `E` is the
disjoint remainder. Thus `d` is the defect of `V` relative to `A_g`. Suppose
`d>=6`, equivalently `j<=7`. Define

```text
M_d = 478970390236831/39525379884148950000  if d>=7,
      287560991216713/1253243752424235000   if d=6,

H_d = 140 if d>=7,
      180 if d=6.                                      (2a)
```

> **Theorem.** At least one of the following holds:
>
> 1. the stronger safe set has the explicit positive-measure exit
>
>    ```text
>    measure(G_h(V))>M_d>0;                            (3)
>    ```
>
> 2. there are integers `(a_f)_(f in F)` and `(b_e)_(e in E)` such that
>
>    ```text
>    sum_(f in F) a_f f+sum_(e in E)b_e e=0,
>    some b_e!=0,
>    |a_f|,|b_e|<=H_d,                                 (4)
>    ```
>
>    and every nonzero coefficient in (4) is not divisible by `41`.

We call (4) **body-touching** or **outside-core-touching**. The elements of
`E` need not be numerically larger than those of `F`. The relation need not
genuinely cross the cut. If

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

Choose the Jackson degree according to the defect:

| regime | `N` | `H=2N-2` | `C_0` | certified odd sum |
|---|---:|---:|---:|---:|
| `d>=7` | 71 | 140 | 238631 | `sum_(1<=k<=139, k odd) C_k/k^2>290903` |
| `d=6` | 91 | 180 | 502411 | `sum_(1<=k<=179, k odd) C_k/k^2>614083` |

Write these choices as `N_d` and `H_d`. For the integer Jackson coefficients
`C_k` of THM-2145, the circle-distance Fourier series and `pi<355/113` give
the strict caps

```text
eta_71
 <357148519/60146943550
 <297/50000,

eta_91
 <586539659/126632692550
 <93/20000.                                           (10)
```

The nonzero Fourier coefficients are also explicit. For `k!=0`,

```text
Fourier(1_G,k)=-sin(6 pi k/41)/(pi k).                (11)
```

Every Jackson multiplier is positive throughout `[-H_d,H_d]`. Since `41`
is coprime to `6`,

```text
Fourier(q_(N_d),k)!=0
iff
k=0, or 0<|k|<=H_d and 41 does not divide k.          (12)
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
Q_line(t)=product_(v in V)q_(N_d)(vt),

Q_lift(x)=product_(f in F)q_(N_d)(fx_0)
          product_(i=1)^d q_(N_d)(x_i).               (14)
```

The crucial point is that all retained-core factors still share `x_0`. Thus
the lift retains their complete relation lattice, exact union geometry, and
every higher overlap. Only the outside-core factors are decorrelated.

Put

```text
I=integral_T Q_line,
A_F=integral_T product_(f in F)q_(N_d)(ft),
J=integral_(T^(d+1))Q_lift
  =A_F(35/41)^d.                                      (15)
```

The product telescope, (7), and Haar invariance give

```text
|I-measure(G_h(V))|<=13 eta_(N_d),                    (16)
|A_F-measure(G_h(F))|<=j eta_(N_d).                   (17)
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

Integer multiplication on the time circle is surjective and Haar preserving,
so the definition of `C` in (2) gives

```text
measure(G_h(F))=measure(G_h(C)).                       (20a)
```

No coprimality or primitivity assumption on `g` is needed. Put

```text
B_j=min_({C subset {1,...,13}:|C|=j}) measure(G_h(C))
```

THM-2146 gives the rows through `j=6`. The exact companion to this theorem
extends the census through `j=7`:

| `j` | `B_j` |
|---:|---:|
| 0 | `1` |
| 1 | `35/41` |
| 2 | `59/82` |
| 3 | `1615/2706` |
| 4 | `239/492` |
| 5 | `2729/7380` |
| 6 | `153101/568260` |
| 7 | `39965/211068` |

For the last row, both a complete boundary-cell evaluator and an independently
merged danger-interval evaluator inspect all
`binomial(13,7)=1716` seven-cores. They agree term by term, and the minimizer
is unique:

```text
C={1,5,7,8,9,11,13}.                                  (20b)
```

Equations (15)--(17), under the no-relation equality (20), imply

```text
measure(G_h(V))
 >=B_j(35/41)^d
   -(13+j(35/41)^d)eta_(N_d).                         (21)
```

For `d>=7`, the seven exact comparisons `7<=d<=13`, `j=13-d`,
give

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

For `d=6`, `j=7`, the newly exact row instead gives

```text
B_7(35/41)^6
 =73466285703125/1002595001939388,

13+7(35/41)^6
 =74619214508/4750104241.                            (24a)
```

Insert the strict `eta_91<93/20000` cap from (10). The remaining margin is

```text
B_7(35/41)^6
 -(13+7(35/41)^6)(93/20000)

=287560991216713/1253243752424235000
>0.                                                   (24b)
```

Thus both rows of (2a) prove (3) whenever (4) is absent. The adjacent values
`N=70` for the `d>=7` uniform ledger and `N=90` for the `d=6` ledger fail
the same rational-cap tests. Thus `N=71` and `N=91` close while their
immediate predecessors `N=70` and `N=90` do not under the displayed
rational-cap ledgers. This does not prove that no nonadjacent smaller cutoff
could close through a different estimate, nor any optimality for the actual
smoothing error or true relation height. QED.

## 5. Why endpoint denominators and six scalar signs do not repair the route

### 5.1 Near-equal large speeds down-convert to fixed low frequencies

There is a tempting extra input in the reversed direction: every boundary
created by a body speed `v>=14` has a large denominator. This does **not**
make the body indicator have small low-frequency Fourier mass.

Let

```text
s(x)=1_{||x||>=3/41},
c_k=Fourier(s,k),
E_N={N,N+1}.                                         (25a)
```

Thus

```text
c_0=35/41,
c_k=-sin(6 pi k/41)/(pi k),              k!=0.       (25b)
```

Use the convention

```text
Fourier(phi,k)=integral_T phi(t)e(-kt)dt.
```

> **Large-denominator down-conversion lemma.** For every fixed
> `f in {1,...,13}` and every `N>f`,
>
> ```text
> Fourier(1_(G_h(E_N)),f)
>   =c_f^2+R_(N,f),
>
> |R_(N,f)|
>   <=(1/6)((N-f)^(-2)+N^(-2)).                      (25c)
> ```
>
> In particular,
>
> ```text
> limit_(N->infinity) Fourier(1_(G_h(E_N)),f)
>   =(sin(6 pi f/41)/(pi f))^2>0.                    (25d)
> ```

Indeed

```text
1_(G_h(E_N))(t)=s(Nt)s((N+1)t).
```

The frequency equation for a product Fourier term is

```text
aN+b(N+1)=f.
```

All its integer solutions are

```text
a=q(N+1)-f,                 b=f-qN,       q in Z.
```

Consequently

```text
Fourier(1_(G_h(E_N)),f)
 =sum_(q in Z)c_(q(N+1)-f)c_(f-qN).                 (25e)
```

This convolution is absolutely convergent: for `q>=1`, both nonzero
indices have magnitude at least `q(N-f)`, while for `q<=-1`, both have
magnitude at least `|q|N`. Since `|c_k|<=1/(pi|k|)` for `k!=0`, the terms
with `q!=0` have total absolute value at most

```text
(1/pi^2)sum_(q>=1)q^(-2)
  ((N-f)^(-2)+N^(-2))
 =(1/6)((N-f)^(-2)+N^(-2)).
```

The `q=0` term is `c_(-f)c_f=c_f^2`. This proves (25c)--(25d).
Equivalently, the coefficient pair `(-f,f)` uses the exact relation

```text
(-f)N+f(N+1)=f:                                     (25f)
```

two high carrier frequencies heterodyne to the fixed low frequency `f`.

On the other hand, every boundary point belonging to the `v`-comb has
the form

```text
(41k+/-3)/(41v).                                    (25g)
```

Its numerator is coprime to `41`, so its reduced denominator is still a
multiple of `41`. Hence every boundary of `G_h(E_N)` has reduced
denominator at least `41`, uniformly in `N`, while each of the thirteen
fixed low Fourier coefficients in (25d) stays away from zero.

This identifies the first failed implication in the proposed
endpoint-denominator repair:

```text
large boundary denominators
  -/-> small low-frequency mass of the body indicator. (25h)
```

What is missing is the relation lattice of the body boundaries. Nearly
equal large speeds can subtract before any frequency-decay estimate is
available. The scalar discrepancy against a small comb is a signed trace
of these coefficients on the lattice `fZ`; denominator size alone does
not control that trace. LEM-011 computes a different object--the Fourier
transform of the uncovered-measure function on independent phase
variables--and does not provide the missing one-dimensional bound for
`1_(G_h(E))`.

### 5.2 The six signed scalar discrepancies still fail on the hostile row

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

At every chosen scale `g`, the theorem converts the defect-`>=6` branch
relative to `g{1,...,13}` into the exact dichotomies

```text
d>=7: strong 3/41 positive-measure exit
       or height-140 outside-core-touching 41-unit relation;

d=6:  strong 3/41 positive-measure exit
       or height-180 outside-core-touching 41-unit relation. (33)
```

This is a substantial reduction of the unbounded-height branch, but it is
not its closure:

1. the relation may have support two or lie wholly inside the outside-core
   body;
2. it need not be independent of relations already known;
3. bounded coefficients do not bound the speed magnitudes, and an
   outside-core speed need not be numerically large;
4. THM-2171 gives a finite algebraic representative for fixed relation
   systems, but its target-mixing example forbids replacing the original row
   without a phase/current argument; and
5. absence of the stronger `3/41` exit is not failure of the `1/14` LRC
   predicate.

There is an important overlap but no dominance with THM-2166. For `g=1`,
defect six, and zero ordinary `1/14` safe measure, the stronger safe exit in
(33) is impossible, so THM-2179 forces a relation of global height `180`.
This is numerically below THM-2166's far-factor height `298`. But THM-2179's
relation may lie wholly in the outside-core body and supplies neither a
nonzero cut carry nor a sparse core representation. THM-2166 instead forces
a genuinely crossing relation, with far height `298`, core support at most
two, core height `57`, and cut carry at most `708`. Moreover, THM-2179 allows
an arbitrary chosen dilate `g{1,...,13}` and is an all-row
strong-exit/relation dichotomy, whereas THM-2166 is stated for an undilated
seven-core under the zero-safe premise. THM-2179 therefore has lower global
height and broader trigger/scale scope; THM-2166 has the stronger
crossing/carry/sparsity conclusion. Neither dominates.

The reusable mechanism is broader than the displayed constants. Put all
factors whose internal correlations must be preserved on common torus
coordinates, put the genuinely free factors on separate coordinates, and
compare line and lifted constant terms only after positive whole-product
smoothing. A discrepancy then produces a bounded relation that touches a
decorrelated coordinate. This is the relative-character idea of THM-2054
combined with THM-2145's sharper Jackson kernel, THM-2146's small-core
floors, and the exact seven-core extension above.

## 7. Zero-safe fixed-base-41 corollary

Assume now that the ordinary radius-`1/14` safe measure of `V` is zero and
that, at a chosen AP scale `g`, the defect in (2) satisfies `d>=7`. Since

```text
G_(3/41)(V) subset G_(1/14)(V),                       (34)
```

the positive exit in (3) is impossible. Take the height-`140` relation `m`
from the other branch and primitive-normalize it. Every nonzero coefficient
remains a `41`-unit because the normalizing gcd is itself a `41`-unit.
Choose an outside-core label `e` with

```text
m_e!=0 mod 41.                                        (35)
```

THM-2169 applied to deletion of `e` gives a primitive nonzero relation `u`
with

```text
u_e=0,                    ||u||_infinity<=1247.       (36)
```

Primitivity supplies some `j!=e` with `u_j!=0 mod 41`. Therefore

```text
det [[m_e,m_j],[u_e,u_j]]=m_e u_j!=0 mod 41.          (37)
```

The two relations are independent over both `F_41` and `Q`. Their
unrestricted simultaneous base-`41` digit fibres have size `41^11`.
Primitive coefficient bounds give

```text
||m||_1<=12*140+139=1819,
||u||_1<=11*1247+1246=14963.                          (38)
```

Thus there are at most

```text
1818*14962=27200916                                   (39)
```

positive-row carry pairs, and THM-2171's ordered algebraic cap is

```text
26*27200916=707223816.                                (40)
```

This raw pair retains the deletion anchor `u_e=0` and the `41`-unit pivot.
Applying THM-2187 to its rational-plane saturation instead retains `m`,
gives a universal-radix second basis vector of height at most `1247`, and
has the slightly larger universal bound

```text
1818*(13*1247-2)=29467962                             (41)
```

carry pairs. Saturation may lose the zero at coordinate `e`, so the raw
base-`41` and saturated universal-base forms are complementary rather than
identical.

This corollary is conditional on zero-safety and `d>=7` at the chosen scale.
It does not say that either pivot stays live at every level, retain phase
current, or close LRC(14). THM-2179's distinctive gain is the scale-local
outside-core touch; THM-2187 already makes arbitrary saturated planes
universal-radix.

## 8. Exact referee

The companion performs the following exact checks:

1. direct convolution and the closed Jackson formula agree on all `141`
   coefficients at `N=71` and all `181` coefficients at `N=91`, with every
   multiplier in both supports positive;
2. the two odd Jackson sums, all rational eta caps, the seven `d>=7`
   margins, the `d=6` margin, and the adjacent `N=70` and `N=90` controls
   all pass;
3. two independent rational interval algorithms agree on every one of the
   `1716` seven-cores and certify both `B_7` and its unique minimizer;
4. the height-`140` and height-`180` supports respectively have exactly
   `274` and `352` signed nonzero `41`-unit modes;
5. the hostile row is evaluated by two independent rational interval
   algorithms--complete boundary cells and merged danger intervals;
6. the six epsilon signs, their sum, the negative level-one bound, the full
   danger-depth distribution, the higher-overlap correction, and the
   aggregate covariance are checked exactly.

Normal and optimized Python executions agree, and after newline normalization
they match the stored LF transcript. An independent proof audit checked the
global arrangement, termwise union comparison, Jackson ledger, arbitrary
scale lift, relation scope, and THM-2166 comparison, and accepted the addendum.
