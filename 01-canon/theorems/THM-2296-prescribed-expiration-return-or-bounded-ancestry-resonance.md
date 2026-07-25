---
id: THM-2296
title: "Prescribed-expiration return or bounded ancestry resonance"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every one of the 165 first-depth-one scalar
  rows, choose the labelled exclusive owner supplied by THM-2263 and go to
  its prescribed expiration time. Either a genuine blocker-to-other-blocker
  return already has Haar mass at least
  39002430583/53493927587100 on a strict row or
  13560199813/53493927587100 on a repeated-first row, or the exact ancestry
  multiplicity and the same blocker-only residual share a nonzero integer
  Fourier atom n with 1<=n<=4SS_j-1<=4S^2-1. Equivalently, the source atoms at
  13^(lambda_j+1)n and the current blocker-residual atom at n are both
  nonzero. The bilinear endpoint-Prony mechanism gives a universal finite
  bank under THM-2199's primitive speed ceiling. Dividing a ramified common
  atom to a nonzero mod-13 residue requires pushing the target too, which
  destroys current blocker service; periodic hostile controls prove that
  loss is real. No scalar profile is excluded, the exact THM-2276 pair
  frequency is not selected, and LRC(14) remains open.
source: codex-2026-07-25-prescribed-expiration-bilinear-prony
depends_on:
  - THM-2199-effective-positive-subspace-rank-lift
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2289-finite-spectral-certificates-bound-null-exclusive-owner-handoff-times
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2271-expiration-support-forces-a-weighted-owner-absorber-cut
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
  - THM-2288-shallow-owner-bv-mixing-and-delayed-blocker-handoff
  - THM-2291-repeated-owner-bv-mixing-and-delayed-blocker-handoff
script: 04-computation/lrc14_prescribed_expiration_bilinear_prony_thm2296.py
output: 05-knowledge/results/lrc14_prescribed_expiration_bilinear_prony_thm2296.out
script_sha256: 6cd156bc849f192dd7fca245b16e5b69ba15d7316627965caf7ae8dbd1953959
output_sha256: 740f66a13e943f4b3406431c418da7c85271c3c8404cf046962c91f087b43bac
hash_basis: working-tree bytes (LF)
---

# THM-2296 -- fixed-time blocker return or bounded ancestry resonance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The delayed-handoff theorems prove that every remaining scalar row eventually
returns a positive exclusive-owner flow to another blocker. They deliberately
stop before the prescribed expiration time. The endpoint-Prony theorem proves
bounded atoms in marked residue classes, but does not make one atom interact
with current blocker service.

There is an exact fixed-time alternative between those two missing objects:

```text
large blocker-only return at prescribed expiration

or

one bounded exact Fourier atom shared by
  the ancestry multiplicity and the current blocker-only residual.       (1)
```

The mechanism is bilinear rather than residuewise. The product of two Fourier
coefficients of step functions is an exponential sum whose nodes are
**differences of endpoints**. A nonzero covariance therefore has a bounded
common atom.

## 1. The selected owner and its exact clock

Use

```text
T(x)=13x mod 1,

D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},

A_0=C_H minus union_(i=1)^5 D_(q_i).                 (2)
```

Assume the scalar cover

```text
A_0 subset D_(c_1) union D_(c_2) union D_(c_3)      (3)
```

almost everywhere and one of the `165` first-depth-one profiles

```text
c_h=13^(lambda_h)u_h,

(lambda_1,lambda_2,lambda_3)=(1,b,c),

either 2<=b<c  (150 strict rows),
or     b=1<c   (15 repeated-first rows),

5<=c<=19.                                            (4)
```

All normalized cores and all six guard/unit coefficients are thirteen-units,
`H` is odd, and the usual scalar distinctness and primitivity hypotheses hold.
Put

```text
S=H+sum_(i=1)^5 q_i+sum_(h=1)^3 c_h.                (5)
```

THM-2289, using THM-2263 and the parity-sharp guard/blocker cap, selects one
label `j` for which

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(h!=j)D_(c_h),

R_j=A_0 minus D_(c_j)                               (6)
```

satisfy

```text
strict:
  measure(E_j)>=alpha_s:=15041431/593783190;

repeated-first:
  measure(E_j)>=alpha_r:=5229541/593783190;

both:
  measure(R_j)>=beta:=2593/90090.                   (7)
```

Retain the selected depth and define its prescribed expiration clock

```text
k_j=lambda_j+1,                         2<=k_j<=20. (8)
```

Let `P` be the normalized Perron operator

```text
P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13),               (9)
```

and define

```text
g_j=P^(k_j)1_(E_j),             h_j=1_(R_j).        (10)
```

The function `g_j` is the exact ancestry **multiplicity density**, including
all sibling sheets and normalized to preserve mass. The transfer identity is

```text
rho_j
 :=measure(E_j intersection T^(-k_j)R_j)
  =integral g_j h_j.                                (11)
```

Thus `rho_j` is already the desired blocker-only return at the prescribed
clock. At almost every source point only `c_j` is available, and at almost
every target point neither `c_j`, the guard absorber, nor any of the five unit
masks is available. The global cover forces one of the other two blockers.

## 2. Bilinear endpoint-Prony lemma

We first prove the general finite-rank mechanism.

> **Bilinear endpoint-Prony lemma.** Let `g,h` be real step functions on the
> circle with respectively `J` and `K` nonzero jumps. If
>
> ```text
> g_hat(q) conjugate(h_hat(q))!=0                    (12)
> ```
>
> for some nonzero integer `q`, then (12) holds for some positive integer
>
> ```text
> 1<=q<=JK-1.                                       (13)
> ```

Write the jump of `g` at `x` as `Delta_g(x)`, and similarly for `h`.
Distributional differentiation gives, for every nonzero integer `n`,

```text
2pi i n g_hat(n)
 =sum_(x in Jump(g)) Delta_g(x)exp(-2pi i n x),      (14)

2pi i n h_hat(n)
 =sum_(y in Jump(h)) Delta_h(y)exp(-2pi i n y).      (15)
```

Consequently

```text
C_n
 :=(2pi n)^2 g_hat(n)conjugate(h_hat(n))

 =sum_(x,y) Delta_g(x)Delta_h(y)
             exp(-2pi i n(x-y)).                    (16)
```

After combining equal endpoint-difference nodes and deleting zero combined
coefficients, (16) is an exponential sum

```text
C_n=sum_(ell=1)^L gamma_ell z_ell^n,

1<=L<=JK,                                           (17)
```

with distinct nonzero nodes `z_ell`.

If a nonzero `L`-node exponential sum vanished at `L` consecutive integers,
the corresponding consecutive-power matrix would have determinant

```text
product_ell z_ell^(n_0)
product_(a<b)(z_b-z_a)!=0.                          (18)
```

Every coefficient `gamma_ell` would vanish, a contradiction. Assumption
(12) says that (17) is not the zero sequence. Periodicity of both step
functions gives

```text
C_0=(sum_x Delta_g(x))(sum_y Delta_h(y))=0.          (19a)
```

If `C_n` also vanished for `n=1,...,L-1`, the sequence would have `L`
consecutive zeros at `0,...,L-1`. It therefore cannot vanish at all of

```text
n=1,2,...,L-1.                                      (19b)
```

In particular `L>=2`, and some positive `n<=L-1<=JK-1` satisfies (12),
proving the lemma. QED.

The product `JK` is the honest endpoint-difference rank. Applying the
one-function endpoint bound separately would produce two atoms that need not
be the same; the bilinear sequence is what retains equality of frequency.

## 3. The prescribed-time dichotomy

The zero Fourier terms of (10) are

```text
g_j_hat(0)=measure(E_j),
h_j_hat(0)=measure(R_j).                            (20)
```

Parseval gives

```text
rho_j-measure(E_j)measure(R_j)
 =sum_(n!=0) g_j_hat(n)conjugate(h_j_hat(n)).        (21)
```

The series on the right is absolutely convergent by Cauchy--Schwarz, since
both Fourier coefficient sequences lie in `ell^2`.

Define the branch floors

```text
delta_s=alpha_s beta
       =39002430583/53493927587100,

delta_r=alpha_r beta
       =13560199813/53493927587100.                 (22)
```

If

```text
rho_j>=delta_s       on a strict row,
rho_j>=delta_r       on a repeated-first row,       (23)
```

then the first arm of (1) holds, quantitatively and at the prescribed clock.

Otherwise (7) implies

```text
rho_j<measure(E_j)measure(R_j).                      (24)
```

The left side of (21) is nonzero. At least one nonzero summand in (21) is
therefore nonzero, and the bilinear endpoint-Prony lemma supplies a common
positive atom bounded by the jump product.

This proves the exhaustive alternative:

```text
rho_j>=delta_branch

or

there exists n>=1 with
  g_j_hat(n)!=0
  and h_j_hat(n)!=0.                                (25)
```

The two arms need not be disjoint: a large return may also have common
nonzero atoms. What matters is that failure of the quantitative return forces
the spectral arm.

## 4. Jump counts and the exact finite bank

Every boundary point of `E_j` belongs to the combined endpoint bank of the
guard, five unit combs, and three blocker combs. A comb of coefficient `a`
has at most `2a` boundary points, so

```text
#Boundary(E_j)<=2S.                                 (26)
```

For a finite step function `f`, every jump of `P^k f` is the image under
`T^k` of a source jump. Indeed, away from those finitely many images, all
inverse branches in (9) stay in fixed source constancy cells. Colliding
images can only cancel jumps. Hence

```text
J_j:=#Jump(g_j)<=#Boundary(E_j)<=2S.                (27)
```

The target `R_j` uses only the guard, five units, and the selected blocker.
Put

```text
S_j=H+sum_(i=1)^5q_i+c_j<=S.                        (28)
```

Then

```text
K_j:=#Jump(h_j)<=2S_j.                              (29)
```

Equations (13), (27), and (29) sharpen (25) to

```text
1<=n<=J_jK_j-1<=4SS_j-1<=4S^2-1.                   (30)
```

The Perron normalization retains the source ancestry exactly:

```text
(P^k f)_hat(n)=f_hat(13^k n).                       (31)
```

Therefore the spectral arm in (25) is precisely

```text
1_(E_j)_hat(13^(lambda_j+1)n)!=0,
1_(R_j)_hat(n)!=0,                                  (32)

1<=n<=4SS_j-1.
```

This is an exact fixed-time ancestry resonance with the **same** current
blocker-only residual, not merely two atoms in the same residue class.

There is also a universal, although enormous, bank. Let `V_*` be
THM-2199's explicit `197`-digit primitive speed ceiling. THM-2203's scalar
section gives

```text
S<=V_*/8+8V_*/16=5V_*/8.                            (33)
```

Thus the terminal common atom may always be chosen with

```text
n<=25V_*^2/16-1,                                    (34)
```

the `393`-digit integer printed by the companion. Since `k_j<=20`, the
source ancestry atom in (32) is bounded by

```text
13^20(25V_*^2/16-1),                                (35)
```

the printed `415`-digit integer. Equations (34)--(35) make the alternative
a genuinely finite global bank, not an ineffective compactness statement.

## 5. The first information-losing map

Suppose the common atom in (32) is ramified:

```text
n=13^t n_0,             13 does not divide n_0.     (36)
```

Applying (31) to both functions gives

```text
(P^(k_j+t)1_(E_j))_hat(n_0)!=0,
(P^t1_(R_j))_hat(n_0)!=0.                           (37)
```

This lawfully reaches a nonzero mod-thirteen residue, but the target in
(37) is now `P^t1_(R_j)`, not the current-service indicator `1_(R_j)`.
Pushing the target averages its blocker eligibility over predecessor sheets.
It no longer says that the terminal point is serviced by another blocker.

This isolates the first information loss:

```text
exact prescribed-time pair:
  ancestry multiplicity g_j  <-> current blocker residual 1_(R_j)

divide a common atom by 13:
  push both sides by P

lost:
  current target service and its named blocker label.          (38)
```

The same projection explains the apparent mismatch between THM-2276 and
this theorem. THM-2276's normalized pair frequency is a thirteen-unit, hence
a nonzero root-character channel. The scalar Perron ancestry `P1_F` retains
only source frequencies divisible by thirteen. It is the zero root-character
projection. A proof cannot simultaneously use the pair's nonzero character
and terminal scalar service without retaining the rooted character and a
current-service sidecar.

THM-2269 retains the former as a marked root mask. The present theorem
retains the latter through `1_(R_j)`. The missing object is their
gauge-faithful signed pairing on the same sibling fibres.

## 6. Ramification is a real obstruction

The loss in (38) is not removable by a general step-function argument.
For any positive integer `d`, choose a relative length `a=1/4` and define
two adjacent periodic interval unions

```text
G_d=union_(r=0)^(d-1)
       [r/d,(r+a)/d),

H_d=union_(r=0)^(d-1)
       [(r+a)/d,(r+2a)/d).                          (39)
```

They are disjoint, positive, and each has exactly `2d` jumps. Both are
`1/d`-periodic, so

```text
G_d_hat(n)=H_d_hat(n)=0              unless d divides n. (40)
```

At `n=d`, both coefficients are nonzero because the normalized interval
length `a=1/4` is not integral. Their first positive common atom is therefore
exactly `d`.

Taking

```text
d=13^t                                                (41)
```

gives arbitrarily deep common-spectrum ramification. There need not be a
thirteen-unit common atom before both functions are pushed `t` times. This
is a universal harmonic-analysis control, not a scalar cover or LRC
counterexample, but it proves that target-service loss in (38) requires
genuinely global LRC structure to repair.

The example also shows why some boundary-complexity dependence is necessary:
the first common atom tends to infinity with the number of jumps. The
quadratic `JK-1` bound is uniform and explicit, not claimed sharp.

## 7. Frontier effect and exact stopping boundary

The theorem changes the fixed-time frontier on all remaining rows:

```text
165 selected exclusive-owner flows
  -> prescribed expiration, not a delayed time
  -> quantitative blocker-only return
     OR
     bounded exact ancestry/current-service resonance.          (42)
```

It therefore removes “nothing is known at prescribed expiration” as a
correct description. The unresolved problem is now a two-arm finite one.

It does **not** shrink THM-2276's `696` exact base-frequency cancellation
bank. The atom `n` in (32) is selected by endpoint-difference covariance and
need not equal the shallow pair multiplier `m u_1`, any multiple of it, or a
specified relation carry. Conversely, the local THM-2276 cancellation
witness does not satisfy the global cover and does not decide which arm of
(42) a genuine counterexample must enter.

The `45` non-interior rows are included in (42), unlike the anchored
`c_1` rank-three route. On those rows the selected owner may be one of the
two depth-one labels or the deep label; equation (8) retains its actual
clock.

The next decisive tests are:

1. on the return arm, charge the fixed positive blocker handoff against an
   owner-cut or gap-capacity inequality;
2. on the resonance arm, determine whether the common atom can be chosen in
   the shallow pair carry subgroup;
3. retain the nonzero root character while testing current blocker service,
   rather than applying the scalar projection in (38);
4. use THM-2284's anchored address bank only after proving that its selected
   address preserves the same atom in (32).

No row is excluded and LRC(14) remains open.

## 8. Exact companion

The companion uses explicit integer and `Fraction` arithmetic. It verifies:

1. the exact `150+15` profile census and both products in (22);
2. the Perron/Fourier normalization (31) through three depths;
3. forty exact consecutive-power Vandermonde determinants with negative and
   positive starting exponents;
4. the complete symbolic jump ledger in (26)--(30);
5. THM-2199's primitive ceiling, the scalar bound (33), and the exact
   `393`- and `415`-digit banks; and
6. the ramified hostile family through depth six.

Reproduce with

```bash
python3 04-computation/lrc14_prescribed_expiration_bilinear_prony_thm2296.py
python3 -O 04-computation/lrc14_prescribed_expiration_bilinear_prony_thm2296.py
```

Normal and optimized transcripts are byte-identical to the stored output.
The universal Parseval, distributional derivative, jump-image, and
Vandermonde arguments are proved above rather than inferred from the finite
checks. QED.
