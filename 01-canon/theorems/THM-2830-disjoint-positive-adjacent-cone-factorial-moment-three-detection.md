---
id: THM-2830
title: "Disjoint positive adjacent-cone factorial moment-three detection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Two nonzero
  nonnegative combinations of adjacent factorial differences, separated
  by one support cut, have the required mixed cubic orientation.  Equality
  is exactly the two adjacent singleton rays at the cut.  The proof turns
  the quartic expression into strict monotonicity of a Pascal-kernel
  quotient.  Its adjacent-row determinant has positive cyclic cubic
  coefficients by a global beta/Abel prefix certificate.  Consequently
  factorial moments one through three detect every complex plane spanned
  by the two disjoint cones; when the constant slot is absent, the
  associated two-charge Gaussian envelope is detected by moment at most
  six.  A stronger individual four-label matching inequality remains open
  and is explicitly not used.
source: root/disjoint-adjacent-cone-factorial-orientation-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
related:
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_disjoint_cone_matching_reduction_thm2830.py
output: 05-knowledge/results/gmc_disjoint_cone_matching_reduction_thm2830.out
script_sha256: ce490561fbc888a94e58ba987912983280c286eba13b2638b5bd8c7a432560d4
output_sha256: 5c54b52f8e04fbc363a0ce6f6026f03460bd9bf2053a5358ea116e767a5fc99c
hash_basis: LF-normalized bytes
---

# THM-2830 -- two disjoint factorial prefix cones

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(s^n)=n!,                         f_n=s^n/n!,
d_i=f_(i+1)-f_i.                                          (1)
```

Fix `b>=1` and two finite nonzero coefficient families

```text
lambda_i>=0  (0<=i<b),             mu_j>=0  (j>=b),

U=sum_(i<b)lambda_i d_i,            V=sum_(j>=b)mu_j d_j. (2)
```

Then

```text
D(U,V):=
 2L(V^3)L(UV)-3L(UV^2)L(V^2) >=0.                       (3)
```

Equality holds exactly when, for positive constants `lambda,mu`,

```text
U=lambda d_(b-1),                   V=mu d_b.             (4)
```

This proves factorial-moment-three detection on every complex plane
spanned by the two disjoint cones.  If `lambda_0=0`, it also gives
two-charge Gaussian detection by moment at most six.

## 1. Positive tensors and the shifted derivative identity

Put

```text
H_n(k)=L(d_n d_k)=binom(n+k,n),
T_n(p,q)=L(d_n d_p d_q).                                (5)
```

Both tensors are strictly positive.  For the cubic tensor this is the
explicit THM-2828 identity

```text
T_n(p,q)
 =[(n+p+q)!/(n!p!q!)]
   [2(S+1)^2+S(np+nq+pq)-npq]/
   [(n+1)(p+1)(q+1)],            S=n+p+q.               (6)
```

For the fixed upper direction `V`, define

```text
A_n=L(d_nV)>0,                 B_n=L(d_nV^2)>0,
R_n=B_n/A_n.                                             (7)
```

The normalized differences satisfy

```text
d_j'=d_(j-1)                         (j>=1).              (8)
```

Also `V(0)=0`, and integration by parts against `e^(-s)` gives

```text
L(P)=L(P')                            when P(0)=0.         (9)
```

Therefore

```text
L(V^2)=2sum_j mu_j A_(j-1),
L(V^3)=3sum_j mu_j B_(j-1).                              (10)
```

For one lower atom write

```text
D_i(V)=2L(V^3)A_i-3B_iL(V^2).                            (11)
```

Substitution of `(10)` gives the load-bearing factorization

```text
D_i(V)
 =6A_i sum_j mu_j A_(j-1)[R_(j-1)-R_i].                 (12)
```

Thus the theorem follows once the quotient sequence is proved strictly
increasing:

```text
R_(n+1)>R_n                         for every n>=0.       (13)
```

Notice that `(13)` is global.  No assumption that the support of `V` lies
above the current row `n` will be used.

## 2. The adjacent-row cyclic determinant

The product identity

```text
d_p d_q
 =binom(p+q,p) f_(p+q)
  +binom(p+q+2,p+1)d_(p+q+1)                            (14)
```

and the hockey-stick sum

```text
S_n(N):=sum_(a<N)H_n(a)=binom(n+N,n+1)                  (15)
```

give, for `N=p+q`,

```text
T_n(p,q)
 =binom(N+2,p+1)H_n(N+1)+binom(N,p)S_n(N).              (16)
```

Define the labelled adjacent-row contribution

```text
Delta_n(p,q;r)
 =T_(n+1)(p,q)H_n(r)-T_n(p,q)H_(n+1)(r).               (17)
```

Its two terms from `(16)` are

```text
Delta_high
 =binom(N+2,p+1)H_n(N+1)H_n(r)(N+1-r)/(n+1),           (18)

Delta_low
 =binom(N,p)H_n(r)S_n(N)
   ( (n+1)(N-1)/(n+2)-r )/(n+1).                       (19)
```

For a labelled triple `(p,q,r)`, cyclically sum the three choices

```text
(p,q;r),                    (q,r;p),                    (r,p;q)
```

and call the results `E_high,E_low`.  Expanding `(7)` shows

```text
A_nB_(n+1)-A_(n+1)B_n
 =sum_(p,q,r)mu_pmu_qmu_r Delta_n(p,q;r).               (20)
```

For three distinct labels the corresponding monomial coefficient is
`2(E_high+E_low)`; for a multiset `(p,p,r)` it is the same cyclic sum
with its natural multiplicities; and for `(p,p,p)` it is one third of
the three equal cyclic terms.  It is consequently enough to prove

```text
E_high+E_low>0                                             (21)
```

for every positive integer triple `(p,q,r)`.

## 3. Exact beta/Abel normalization

Put

```text
x=p+1,              y=q+1,              z=r+1,
t=x+y+z,                                              (22)

C_0=1/[n!^2 x!y!z!(n+1)]>0.
```

Direct factorial cancellation in `(18)--(19)` gives

```text
E_high=C_0 H_*,
E_low =C_0 L_*,                                        (23)
```

where all sums over `{x,y,z}` count the labelled multiset, and

```text
H_*
 =sum_(u in {x,y,z})
    u(t-u)(t-2u)(n+u-1)!(n+t-u-1)!,                    (24)

L_*
 =xyz/[(n+1)(n+2)]
   sum_u (t-u-2)(n+u-1)!(n+t-u-2)!
    [(n+1)t-(2n+3)u-2n-1].                             (25)
```

We prove the stronger pair

```text
H_*>0,                         H_*+4L_*>0.              (26)
```

It then follows that

```text
H_*+L_*=[(H_*+4L_*)+3H_*]/4>0,                         (27)
```

which is `(21)`.

Sort `x<=y<=z` and write

```text
x=m+2,                y=x+a,                z=y+c,
m,a,c>=0.                                               (28)
```

For `u` among the three labels set

```text
W_u=(n+u-1)!(n+t-u-2)!,
phi_u(r)=r^(u-1)+r^(t-u-2).                            (29)
```

The beta integral, paired under `s<->1-s` and then changed by
`r=s/(1-s)`, gives the exact identity

```text
sum_u W_u Q(u)
 =(2n+t-2)! integral_0^1
   [r^n sum_u Q(u)phi_u(r)]/(1+r)^(2n+t-1) dr.          (30)
```

For `0<r<=1`,

```text
phi_u(r)
 =2r^((t-3)/2)
   cosh((u-(t-1)/2)log r).                             (31)
```

Since `x<=y<=z` and `x+y+z=t`,

```text
phi_x>=phi_y>=phi_z>0.                                 (32)
```

Indeed the successive squared-distance differences from `(t-1)/2` are
`(y-x)(z-1)` and `(z-y)(x-1)`.

For any three coefficients, Abel summation reads

```text
sum_u Q(u)phi_u
 =Q(x)(phi_x-phi_y)
  +[Q(x)+Q(y)](phi_y-phi_z)
  +[Q(x)+Q(y)+Q(z)]phi_z.                              (33)
```

It remains only to exhibit positive prefix polynomials.

### 3.1. The high part

Equation `(24)` is `(30)` with

```text
Q_n(u)=u(t-u)(t-2u)(n+t-u-1).                          (34)
```

After `(28)`, the three prefixes in `(33)` are polynomials in
`n,m,a,c` whose complete coefficient statistics are

| prefix | nonzero monomials | smallest coefficient | constant |
|---|---:|---:|---:|
| `Q_n(x)` | 46 | 1 | 48 |
| `Q_n(x)+Q_n(y)` | 51 | 1 | 96 |
| `Q_n(x)+Q_n(y)+Q_n(z)` | 46 | 1 | 144 |

Every coefficient is nonnegative.  Hence `(30)--(33)` prove `H_*>0`.

### 3.2. The quarter combination

Multiplying the second inequality in `(26)` by `(n+1)(n+2)` gives
`(30)` with

```text
P_n(u)
 =(n+1)(n+2)u(t-u)(t-2u)(n+t-u-1)
  +4xyz(t-u-2)
    [(n+1)t-(2n+3)u-2n-1].                            (35)
```

The three prefixes now have

| prefix | nonzero monomials | smallest coefficient | constant |
|---|---:|---:|---:|
| `P_n(x)` | 140 | 1 | 32 |
| `P_n(x)+P_n(y)` | 145 | 1 | 64 |
| `P_n(x)+P_n(y)+P_n(z)` | 136 | 1 | 96 |

Again every coefficient is nonnegative.  Equations `(30)--(33)` prove
`H_*+4L_*>0`.  The exact companion expands all six prefix polynomials
over `Z[n,m,a,c]`; this is a finite symbolic identity, not a bounded
sampling argument.

Combining `(20)--(27)` proves

```text
A_nB_(n+1)-A_(n+1)B_n>0.                              (36)
```

Since `A_n,A_(n+1)>0`, this is exactly `(13)`.

## 4. Equality and the cone theorem

Strict monotonicity in `(13)` turns `(12)` into a sum of nonnegative
terms.  Its `j=i+1` term is zero, and every `j>i+1` term is strict.
Therefore

```text
D_i(V)=0
 iff
 V is a positive multiple of d_(i+1).                 (37)
```

For `i<=b-2`, every occupied upper label satisfies `j>i+1`, so
`D_i(V)>0`.  For `i=b-1`, equality is possible exactly when
`V=mu d_b`.  Finally,

```text
D(U,V)=sum_(i<b)lambda_i D_i(V).                       (38)
```

Thus equality in `(3)` forces all occupied lower indices to equal
`b-1`, and `(4)` follows.  Conversely the adjacent pair in `(4)` gives
equality directly.

## 5. Factorial and Gaussian moment detection

For

```text
H=alpha U+beta V,                      alpha,beta in C, (39)
```

the real quadratic

```text
Q(alpha,beta)=L(H^2)                                   (40)
```

is positive definite: `U,V` are linearly independent real polynomials
and `L(F^2)=integral_0^infinity F(s)^2e^(-s)ds`.

Let `C(alpha,beta)=L(H^3)` and use the Gram/cubic notation

```text
g11=L(U^2),       g12=L(UV),       g22=L(V^2),
t111=L(U^3),      t122=L(UV^2),    t222=L(V^3).        (41)
```

THM-2824 proves that a common projective zero of `Q,C` exists exactly
when the real quadratic divides the real cubic.  One division-free
remainder invariant is

```text
I2
 =3t122 g11g22-2t222 g12g11-t111 g22^2
 =-g11 D(U,V)-t111 g22^2.                              (42)
```

The strictly positive cubic tensor `(6)` gives `t111=L(U^3)>0`.
Equations `(3)` and `(42)` therefore imply

```text
I2<0.                                                   (43)
```

Here `L(H)=0` identically by `(1)--(2)`.  Hence `Q,C` have no common
nonzero projective zero, and

```text
L(H)=L(H^2)=L(H^3)=0                  implies H=0.      (44)
```

If `lambda_0=0`, both `U,V` are divisible by `s`.  Put

```text
h=H/s,
P=W+Z h(ZW),                         W=conj(Z),         (45)
```

for a standard complex Gaussian `Z`.  Charge balance gives

```text
E[P^(2j)]=binom(2j,j)L(H^j),          E[P^(2j-1)]=0,
                                                       j=1,2,3. (46)
```

Thus every nonzero envelope `(45)` is detected by a Gaussian moment of
order at most six.

## 6. The stronger matching sidecar remains open

The earlier reduction polarized `(3)` into four upper labels.  For
`a,b,c,d>i`, put

```text
P_i(a,b,c,d)
 =H_i(a)T(b,c,d)+H_i(b)T(a,c,d)
  +H_i(c)T(a,b,d)+H_i(d)T(a,b,c),                     (47)

N_i(a,b;c,d)
 =T(i,a,b)H(c,d)+T(i,c,d)H(a,b).                      (48)
```

The individual matching conjecture

```text
P_i(a,b,c,d)>=3N_i(a,b;c,d)                            (49)
```

for each of the three complementary matchings would imply coefficientwise
positivity before the cyclic averaging used above.  After multiplying by
`i!a!b!c!d!`, `(49)` becomes

```text
sum_(x in {a,b,c,d})
 (i+x)!(a+b+c+d-x)! tau({a,b,c,d}\{x})

 >=3[
  (i+a+b)!(c+d)!tau(i,a,b)
 +(i+c+d)!(a+b)!tau(i,c,d)
 ],                                                     (50)

tau(a,b,c)
 =(a+b+c+1)
   [1/(a+1)+1/(b+1)+1/(c+1)]-1.
```

The companion verifies `(49)--(50)` and all coordinate forward
differences in the stated finite universe, but there is no universal
proof.  The theorem above does not require `(49)`: the beta/Abel proof
retains only a three-label cyclic coefficient and deliberately forgets
the individual matching.

Two still-tempting stronger shortcuts are false.  In the exact pair block

```text
i=0,                 (y,z)=(3,44),                 a=4,
```

the tilted fourth falling-factorial mean obeys

```text
Phi-[(y)_4+(z)_4]/2=-217452/901<0.                     (51)
```

And in the sign-Laguerre basis, for `V=d_1+t d_3`, the `k=3`
coefficientwise adjacent determinant is

```text
2t(5640t^2+371t-3),                                   (52)
```

which is negative for small positive `t`.  Pairwise convex dominance and
raw coefficientwise Laguerre monotonicity therefore remain invalid proof
routes.

## 7. Exact evidence and independent audits

The companion uses integers, rational numbers, and exact SymPy
polynomials with explicit exception gates; it has no truth-bearing Python
assertions.  It verifies:

1. the factorial tensors, pair product, polarization, and Pascal
   identities;
2. `1,080` exact high/low cyclic normalizations and quarter signs;
3. all six global beta/Abel prefix polynomials, including the complete
   `46/51/46` and `140/145/136` coefficient certificates;
4. `1,375` strict quotient comparisons and `3,016` cone orientations;
5. an independent direct general-jump Newton expansion in `2,000` source
   cells and `1,400` repeated-multiset aggregations;
6. the direct general-jump proof's exact `27`, `424`, `164`, and `474`
   term nonnegative-coefficient certificates;
7. all `19,800` forward matching checks, retained only as evidence for
   the stronger open sidecar; and
8. both false-route hostiles `(51)--(52)`.

One independent audit derived the global beta/Abel proof and equality
classification.  A second independently derived the general-jump
Newton/cyclic proof, including the repeated-label multiplicities and the
fact that only the first Newton order is needed for strictness.  Normal,
optimized, and stored transcripts agree, and the LF-normalized hashes
match.

This theorem enlarges THM-2828 from one arbitrary lower cone direction
against one interval to arbitrary nonzero positive directions on both
sides of a support cut.  It does not prove arbitrary signed radial
coefficients, the stronger matching sidecar, general HYP-8765, or any
new form of unrestricted GMC2 beyond the already closed THM-2022.

**QED.**
