---
id: THM-2830
title: "Disjoint and transport-ordered positive-cone factorial moment-three detection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  nonzero positive adjacent-difference cone V, the quotient
  R_n=L(d_n V^2)/L(d_n V) is strictly increasing, with a sharp universal
  one-step floor.  Integration by parts proves the mixed cubic orientation
  for two cones separated by a support cut, with only the adjacent
  singleton equality.  More generally, weighted stochastic dominance,
  and in particular a shifted coefficient monotone-likelihood-ratio
  order, proves the orientation for interlaced positive cones.  Factorial
  moments one through three detect every independent plane in these
  classes; the constant-free two-charge Gaussian lift is detected by
  moment at most six.
source: root/disjoint-adjacent-cone-factorial-orientation-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
related:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_disjoint_positive_cones_pascal_ratio_thm2830.py
output: 05-knowledge/results/gmc_disjoint_positive_cones_pascal_ratio_thm2830.out
script_sha256: f061612d06c8698f6913d5ad4624b8fe0ad38e40263146ee7e1abe3876c03b0e
output_sha256: 63d5ac47b58fa6892211261a0ed79bde7dae0d5d52dd59ca6eb81ad4cef6b738
hash_basis: LF-normalized bytes
---

# THM-2830 -- Pascal transport for positive factorial cones

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(s^n)=n!,                         f_n=s^n/n!,
d_i=f_(i+1)-f_i.                                      (1)
```

Fix `b>=1` and finite nonzero coefficient families

```text
U=sum_(i<b)lambda_i d_i,            lambda_i>=0,
V=sum_(j>=b)mu_j d_j,               mu_j>=0.           (2)
```

Then

```text
D(U,V):=2L(V^3)L(UV)-3L(UV^2)L(V^2)>=0,               (3)
```

with equality exactly when

```text
U is proportional to d_(b-1),
V is proportional to d_b.                             (4)
```

There is a broader transport theorem.  For arbitrary finite nonzero

```text
U=sum_(i>=0)lambda_i d_i,            lambda_i>=0,
V=sum_(j>=1)mu_j d_j,                mu_j>=0,           (5)
```

put `nu_i=mu_(i+1)`.  If

```text
nu_jlambda_i>=nu_ilambda_j                 for i<j,    (6)
```

then `(3)` holds even when the supports overlap or interlace.  Equality
under `(6)` occurs exactly when `lambda=c nu`, equivalently `U=cV'`.

## 1. Positive tensors and the Pascal quotient

Put

```text
H_n(j)=L(d_nd_j)=binom(n+j,n),
T_n(p,q)=L(d_nd_pd_q).                                 (7)
```

THM-2828 gives `T_n(p,q)>0`.  For any nonzero positive cone
`V=sum_jmu_jd_j`, define

```text
A_n=L(d_nV)=sum_jmu_jH_n(j)>0,
B_n=L(d_nV^2)=sum_(p,q)mu_pmu_qT_n(p,q)>0,
R_n=B_n/A_n.                                           (8)
```

The load-bearing theorem is

```text
R_(n+1)>R_n                         for every n>=0.     (9)
```

It follows from a strict symmetric three-index kernel.

## 2. A universal cyclic-kernel floor

Since

```text
H_(n+1)(r)/H_n(r)=1+r/(n+1),
```

define

```text
C_n(p,q)=T_n(p,q)/(H_n(p)H_n(q)),
m_n(p,q)=(n+1)[T_(n+1)(p,q)/T_n(p,q)-1],              (10)

J_n(p,q,r)=sum_cyc C_n(p,q)[m_n(p,q)-r].              (11)
```

Then

```text
K_n(p,q,r)
 :=sum_cyc[
    T_(n+1)(p,q)H_n(r)-T_n(p,q)H_(n+1)(r)
  ]
 =H_n(p)H_n(q)H_n(r)J_n(p,q,r)/(n+1).                 (12)
```

We prove

```text
J_n(p,q,r)>=6(n+1),
K_n(p,q,r)>=6H_n(p)H_n(q)H_n(r)>0,                    (13)
```

with equality exactly at `p=q=r=0`.

### 2.1. Separable base plus monotone interaction

Let

```text
Beta=n!(n+p+q)!/[(n+p)!(n+q)!]
    =product_(k=1)^q (n+p+k)/(n+k)
    =sum_(k=0)^min(p,q)
       binom(p,k)binom(q,k)/binom(n+k,k)>=1,           (14)

tau=(n+p+q+1)
       [1/(n+1)+1/(p+1)+1/(q+1)]-1.
```

Direct normalization gives `C_n=Beta tau`.  Split

```text
c_n(p,q)=(n+1)[1/(p+1)+1/(q+1)],
E_n(p,q)=C_n(p,q)-c_n(p,q).                           (15)
```

The boundary and boundary increment are

```text
E_n(0,q)
 =q(nq+2n+2q+3)/[(n+1)(q+1)]>=0,                     (16)

Delta_qE_n(0,q)
 =[nq^2+3nq+3n+2q^2+6q+5]/
   [(n+1)(q+1)(q+2)]>0.                               (17)
```

The mixed first difference is

```text
Delta_pDelta_qE_n(p,q)
 =Beta U_n(p,q)/[(n+1)(n+p+1)(n+q+1)]>0,             (18)

U_n(p,q)=
 n^2p+n^2q+2n^2
 +np^2+2npq+6np+nq^2+6nq+6n
 +p^2q+3p^2+pq^2+6pq+7p+3q^2+7q+4.
```

Symmetry and telescoping `(17)--(18)` show that `E_n` is nonnegative and
coordinatewise nondecreasing.

### 2.2. Effective-mean compensation

Define `F_n` by

```text
C_n(p,q)m_n(p,q)
 =c_n(p,q)(p+q+1)+F_n(p,q).                           (19)
```

For `s=p+q`, set

```text
Q=(n+1)^2(n+2)(s+2)^2,

S_n(p,q)=
 n^2p^2+2n^2pq+4n^2p+n^2q^2+4n^2q+4n^2
 +np^2q+4np^2+npq^2+6npq+11np+4nq^2+11nq+10n
 +2p^2+2pq+6p+2q^2+6q+6.
```

Every coefficient of `S_n` is positive, and exact simplification gives

```text
F_n-(s/2)E_n
 =[(Beta-1)Q+Beta sS_n]/
   [2(n+1)(n+2)(p+1)(q+1)]>=0.                        (20)
```

Assume `p<=q<=r`, write `x=q-p,y=r-q`, and abbreviate
`E_pq=E_n(p,q)`.  Monotonicity gives

```text
sum_cyc E_pq[(p+q)/2-r]
 =1/2{
   x[(E_qr-E_pq)+(E_qr-E_pr)]
  +y[(E_pr-E_pq)+(E_qr-E_pq)]
 }>=0.                                                (21)
```

The separable part cancels exactly:

```text
sum_cyc c_n(p,q)(p+q+1-r)=6(n+1).                     (22)
```

Equations `(19)--(22)` prove `(13)`.  The certificate gap `(20)` can
vanish on all three pairs only when `p=q=r=0`, proving the equality claim.

## 3. Strictness and a coercive margin

Dummy-index symmetrization gives

```text
B_(n+1)A_n-B_nA_(n+1)
 =1/3sum_(p,q,r)mu_pmu_qmu_rK_n(p,q,r)
 >=2A_n^3.                                             (23)
```

Hence

```text
R_(n+1)-R_n>=2A_n^2/A_(n+1)>0,                        (24)

R_m-R_i>=2sum_(t=i)^(m-1)A_t^2/A_(t+1),       m>i.   (25)
```

Equality in `(23)--(24)` occurs exactly for `V` proportional to `d_0`.
The increasing quotient is not generally convex:

```text
V=2d_0+d_10,
R_2-2R_1+R_0=-53721500/663<0.                         (26)
```

Thus `(24)`, not a hidden higher-order total-positivity assertion, is the
sharp structural output used below.

## 4. The separated-cone theorem

Return to `(2)`.  Since `d_j'=d_(j-1)` for `j>=1`, `V(0)=0`, and

```text
L(P)=integral_0^infinity P(s)e^(-s)ds,
```

integration by parts yields

```text
L(V^2)=2sum_jmu_jA_(j-1),
L(V^3)=3sum_jmu_jB_(j-1).                             (27)
```

For one lower atom `i<b`,

```text
D_i(V)
 :=2L(V^3)L(d_iV)-3L(d_iV^2)L(V^2)

 =6A_i sum_jmu_jA_(j-1)[R_(j-1)-R_i]>=0.              (28)
```

The quantitative version is

```text
D_i(V)>=12A_i sum_jmu_jA_(j-1)
                 sum_(n=i)^(j-2)A_n^2/A_(n+1),        (29)
```

with an empty inner sum at `j=i+1`.  Strictness shows that `D_i=0`
exactly when every occupied upper index is `j=i+1`.  Since `j>=b>i`,
this forces `i=b-1` and `V` proportional to `d_b`.  Finally,

```text
D(U,V)=sum_(i<b)lambda_iD_i(V),                        (30)
```

which proves `(3)--(4)`.

## 5. Stochastic transport and coefficient minors

For the general pair `(5)`, define probability laws

```text
alpha_i=lambda_iA_i/L(UV),
beta_i=mu_(i+1)A_i/[L(V^2)/2].                         (31)
```

Equations `(8),(27)` give

```text
D(U,V)
 =3L(UV)L(V^2)[E_betaR-E_alphaR].                     (32)
```

If `beta` first-order stochastically dominates `alpha`, then

```text
E_betaR-E_alphaR
 =sum_(t>=0)[R_(t+1)-R_t]
   [sum_(k>=t+1)beta_k-sum_(i>=t+1)alpha_i]>=0.        (33)
```

Every ratio increment is strict, so equality holds exactly when
`alpha=beta`.

There is also an exact coefficient-level Cauchy--Binet formula:

```text
D(U,V)
 =6sum_(i<j)A_iA_j(R_j-R_i)
   [mu_(j+1)lambda_i-mu_(i+1)lambda_j].               (34)
```

Condition `(6)` makes every summand nonnegative.  It also implies the
tail order directly, because

```text
sum_(j>=t)beta_j-sum_(j>=t)alpha_j
 =sum_(i<t<=j)(beta_jalpha_i-beta_ialpha_j)>=0.        (35)
```

Equality under `(6)` forces every occupied coefficient minor in `(34)` to
vanish, so `lambda=c nu`, equivalently `U=cV'`.

This extension is genuinely interlaced:

```text
U=2d_0+d_1,          V=d_1+2d_2,          D=17460>0,  (36)

U=d_0+d_2=V',        V=d_1+d_3,           D=0.        (37)
```

The order is load-bearing:

```text
U=d_2,                V=d_1+d_3,           D=-33540.  (38)
```

Formula `(34)` identifies the remaining arbitrary-positive-cone problem
precisely: shifted coefficient-ratio reversals contribute with opposing
signs and may cancel.

## 6. Derivative-flag holotopy

Suppose

```text
V=sum_(j>=r)mu_jd_j,                 U=V^(r)
  =sum_i mu_(i+r)d_i,                r>=1,             (39)
```

where `mu_j>=0` is a nonzero log-concave sequence with interval positive
support.  Its adjacent ratios are nonincreasing, hence

```text
mu_(j+1)mu_(i+r)>=mu_(i+1)mu_(j+r),       i<j.        (40)
```

This is `(6)`, so every independent derivative-flag plane
`span{V^(r),V}` satisfies `(3)`.  For `r=1`, `U=V'` is exactly the
equality family.  For every `r>1`, the finite supports of `V^(r)` and
`V'` have different extrema, so the orientation is strict.

The Gaussian corollary below additionally requires
`min supp(mu)>=r+1`, so that `V^(r)` has no `d_0` term.  Without this
extra leading zero, the factorial theorem still applies.

## 7. Moment-three and Gaussian consequences

For `H=alpha U+beta V`, let

```text
Q(alpha,beta)=L(H^2),               C(alpha,beta)=L(H^3).
```

If `U,V` are independent, `Q` is positive definite.  THM-2828 gives
`t111=L(U^3)>0`.  The THM-2824 division-free cubic remainder is

```text
I2
 =3t122g11g22-2t222g12g11-t111g22^2
 =-g11D(U,V)-t111g22^2<0.                             (41)
```

Thus `Q,C` have no common complex projective zero:

```text
L(H)=L(H^2)=L(H^3)=0                  implies H=0.     (42)
```

If `lambda_0=0`, both directions are divisible by `s`.  Put

```text
h=H/s,                 P=W+Zh(ZW),          W=conj(Z).
```

Charge balance gives

```text
E[P^(2k)]=binom(2k,k)L(H^k),          k=1,2,3,
E[P^(2k-1)]=0.                                        (43)
```

Hence every constant-free plane above gives a many-versus-many two-charge
Gaussian envelope detected by moment at most six.

## 8. Independent beta/Abel proof sidecar

The incoming concurrent proof of `(9)` is mathematically independent and
is preserved, rather than overwritten, in

```text
04-computation/gmc_disjoint_cone_matching_reduction_thm2830.py
05-knowledge/results/gmc_disjoint_cone_matching_reduction_thm2830.out
```

with LF hashes

```text
ce490561fbc888a94e58ba987912983280c286eba13b2638b5bd8c7a432560d4
5c54b52f8e04fbc363a0ce6f6026f03460bd9bf2053a5358ea116e767a5fc99c.
```

It starts from

```text
d_pd_q
 =binom(p+q,p)f_(p+q)
  +binom(p+q+2,p+1)d_(p+q+1),                         (44)
```

splits the adjacent-row cyclic determinant into high and low pieces
`H_*` and `L_*`, and proves the stronger pair

```text
H_*>0,                         H_*+4L_*>0.             (45)
```

After sorting `x<=y<=z`, a beta integral produces decreasing kernels
`phi_x>=phi_y>=phi_z`.  Abel summation reduces `(45)` to six prefix
polynomials in nonnegative gap variables.  Their exact nonzero-term counts
are

```text
high:       46, 51, 46,
quarter:   140,145,136,
```

and every coefficient is nonnegative with a positive constant.  The same
companion contains a second general-jump Newton certificate and the older
`19,800` matching-forward checks.  This independent proof gives strong
audit redundancy; the simpler floor `(13)` is the primary proof used by
the transport extension.

## 9. Stronger sidecars and exact scope

The individual four-upper-label matching inequality from the earlier
candidate remains open.  It would prove a stronger coefficientwise
statement before cyclic symmetrization, but neither proof of `(9)` needs
it.  Two other stronger shortcuts remain false: pairwise fourth
falling-factorial dominance fails by `-217452/901`, and raw
coefficientwise Laguerre quotient monotonicity fails for
`V=d_1+td_3` when `t>0` is small.

The primary exact companion verifies the normalized kernel identities,
`2,197` algebra cells, `6,591` residual monotonicity cells, `6,160` cyclic
cells, `60,000` random ratio cells, `1,452` separated orientation cells,
and `2,000` interlaced MLR cells.  It checks the sharp kernel equality,
coercive determinant, convexity hostile, Cauchy--Binet transport, and all
three interlaced controls.  Normal, optimized, and stored transcripts
agree exactly.

This theorem covers separated, stochastically ordered, and
MLR/derivative-ordered positive cones.  It does not cover arbitrary
incomparable positive cones, signed radial coefficients, all of HYP-8765,
SFC(3), NC2, or GMC2.

## 10. Independent hostile audits

Independent audits rederived the beta normalization, boundary and mixed
differences, effective-mean certificate, cyclic rearrangement, exact
`6(n+1)` cancellation, determinant symmetrization, integration by parts,
and equality classification.  Two further independent derivations checked
the stochastic-order and coefficient-minor extensions, derivative-flag
scope, negative unordered hostile, Gram/cubic-divisibility step, and
Gaussian `lambda_0=0` condition.  Both exact companions replay normally
and under `-O`; their LF-normalized hashes match.

**QED.**
