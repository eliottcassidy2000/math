---
id: THM-2209
title: "Sharp quadratic reversed peel and joint Fourier discrepancy ledger"
status: >
  PROVED + FINITE-EXACT + HOSTILE-AUDITED. For every body good set and
  j>=2 peel combs, the sharp quadratic pointwise minorant
  1_(K=0)>=1-K+(2/j)binom(K,2) gives a pair-retaining reversed-peel lower
  bound. It already certifies THM-2179's canonical six-peel hostile row with
  positive margin 733141701884261/8470915906969230. Independently, the entire
  signed level-one correction has one Cauchy--Schwarz bound whose body
  energy is an exact inclusion--exclusion ledger of at most 2^j-1 THM-732
  discrepancies. This does not prove that every reversed-peel row exits.
source: codex-2026-07-24-quadratic-reversed-peel
depends_on:
  - THM-2179-reversed-peel-relative-jackson-relation-packet
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
related:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-594-pair-overlap-law-mirsky-newman-floor
  - LEM-011-exact-fourier-transform-uncovered-measure
---

# THM-2209 -- the sharp quadratic reversed peel

Let `E` be a finite speed body, let

```text
G=G_h(E),                 g=1_G,       m=measure(G),
F={v_1,...,v_j},          j>=2,
D_v={t:||vt||<h},         p=measure(D_v)=2h,          (1)
```

and put

```text
X_v=1_(D_v),              K=sum_(v in F) X_v.         (2)
```

The safe measure after adjoining the peels is

```text
L_h(E union F)=integral g 1_(K=0).                    (3)
```

## 1. The sharp quadratic minorant

For every integer `0<=k<=j`,

```text
1_(k=0)>=1-k+(2/j)binom(k,2).                         (4)
```

At `k=0` and `k=1` this is equality. For `k>=1`, the right side is

```text
-(k-1)(j-k)/j<=0.                                    (5)
```

Moreover, `2/j` is the largest coefficient `c` for which

```text
1_(k=0)>=1-k+c binom(k,2)
```

holds throughout `{0,...,j}`: the endpoint `k=j` forces `c<=2/j`.
Thus (4) is the sharp quadratic minorant with the normalized constant and
linear terms fixed.

Multiply (4) by `g` and integrate. The exact pair-retaining reversed-peel
bound is

```text
L_h(E union F)
 >=m-sum_v measure(G intersection D_v)
      +(2/j)sum_(u<v) measure(G intersection D_u intersection D_v).
                                                               (6)
```

For `j=0`, equation (3) is simply `m`; for `j=1`, the exact first-order
identity is `m-measure(G intersection D_v)`. The coefficient `2/j` is used
only for `j>=2`.

## 2. Covariance form

Define

```text
A_2(F)=sum_(u<v) measure(D_u intersection D_v),

E_1=sum_v [measure(G intersection D_v)-mp],

E_2=sum_(u<v) [
      measure(G intersection D_u intersection D_v)
      -m measure(D_u intersection D_v)].              (7)
```

Then (6) is equivalently

```text
L_h(E union F)
 >=m[1-jp+(2/j)A_2(F)]-E_1+(2/j)E_2.                 (8)
```

Unlike the absolute THM-735 bound, (8) retains the first genuinely missing
higher-overlap coordinate. It is still only a lower bound: triple and
higher intersections can improve the true safe measure further.

## 3. The canonical hostile row is already pair-certifiable

Use THM-2179's exact radius `h=3/41` row

```text
F=(1,2,3,4,6,8),
E=(95,163,187,206,208,214,332).                       (9)
```

Let `p_k` be the Haar mass of points in `G_h(E)` lying in exactly `k` of
the six peel danger combs. THM-2179 gives the exact signed level-one bound

```text
B_1
=-727156708364069/33883663627876920                  (10)
```

and the full higher-overlap correction

```text
R_*
=sum_(k=2)^6 (k-1)p_k
=1812297618203733/11294554542625640.                 (11)
```

For `j=6`, the quadratic term in (6) is

```text
(1/3)S_2,          S_2=sum_(k=2)^6 binom(k,2)p_k.    (12)
```

The exact pair mass, independently recomputed from THM-2179's interval
partition, is

```text
S_2=3659723515901113/11294554542625640.
```

Consequently

```text
B_1+(1/3)S_2
 =733141701884261/8470915906969230
 >0.                                                   (13)
```

Thus the row which defeats every scalar-sign/level-one repair is certified
at the optimally scaled pair layer. Full Boolean depth is not needed for
this hostile instance. Even without computing `S_2` directly,
`binom(k,2)>=k-1` and (11) give the positive lower bound

```text
B_1+(1/3)S_2
 >=3308356432438/103303852524015.                     (13a)
```

## 4. One joint Fourier bound for the signed level-one term

Put

```text
S_F=sum_(v in F)(X_v-p),
V_F=integral S_F^2
   =jp(1-p)+2sum_(u<v)(measure(D_u intersection D_v)-p^2).
                                                               (14)
```

Then

```text
E_1=integral g S_F.                                  (15)
```

Use the Fourier convention

```text
hat(g)(n)=integral_T g(t) exp(-2 pi i n t)dt.
```

Every nonzero Fourier mode of `S_F` belongs to

```text
Omega_F=union_(v in F) vZ.                            (16)
```

The inclusion can be strict because overlapping comb modes may cancel. By
Parseval and Cauchy--Schwarz,

```text
|E_1|^2<=V_F Delta_F(E),                              (17)
```

where the exact union-energy envelope is

```text
Delta_F(E)
 =sum_(n in Omega_F\{0}) |hat(g)(n)|^2
 =sum_(emptyset!=A subset F) (-1)^(|A|+1)
      disc_(lcm(A))(G).                               (18)
```

Here

```text
disc_d(G)=sum_(r!=0)|hat(g)(rd)|^2.
```

Equation (18) is ordinary inclusion--exclusion on the frequency subgroups:

```text
intersection_(v in A) vZ=lcm(A)Z.
```

For `j<=6`, it uses at most `2^j-1<=63` exact THM-732 discrepancy terms.
The peel-side factor `V_F` uses only the `binom(j,2)<=15` pair overlaps.
Consequently the signed level-one certificate

```text
L_h(E union F)
 >=(1-jp)m-sqrt(V_F Delta_F(E))                       (19)
```

is one exact rational-ledger square-root test after the interval endpoints
of `G` are known. It keeps the six peel signs together until the final
Cauchy--Schwarz step.

## 5. Boundary and next exact test

The hostile row proves that no level-one estimate, even with perfect scalar
signs, can certify every reversed peel. The pair minorant (6) is the first
level that survives that control. It is also the cheapest next classifier:
for each bounded relation-template row, take the at-most-six small peels
against the large-speed body and evaluate, in order,

```text
signed level one (10),
joint Fourier ledger (19),
quadratic pair certificate (6).                       (20)
```

Only rows failing (6) should pay for the optimal cubic minorant or the full
Boolean packet. Those failures must retain endpoint current and the
THM-2197/2198 ownership word before any residue quotient.

The theorem does not show that (6) is positive on every LRC template,
control the cubic remainder uniformly, or prove LRC(14). Its exact advance
is that the canonical obstruction to signed reversed peeling is not a
pair-level obstruction.

An independent hostile audit checked the endpoint cases `j=0,1`, sharpness
at `K=j`, every covariance normalization, the union-envelope distinction
in (18), and all hostile-row fractions. It independently recomputed (13)
from exact interval masses and found no defect. QED.
