---
id: THM-3136
title: "One-sided fixed-reference elementary-tail Hasse no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  On the certified support (a,b)=(1,2), the I1 prefix obtained by subtracting
  the top reduced pole 5 has a strictly positive scalar row response in every
  degree N>=5, but its fixed-reference normalized Young current violates the
  bottom-complement upset facet for every even N>=6.  In the complete
  THM-3120 active-prefix census through N=9, all 8,241 degree-five currents
  pass, while 43 later failures require genuinely nonprincipal upsets.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
depends_on:
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3128-active-pole-prefix-labelled-deletion-no-positive-preimage
  - THM-3129-bounded-poset-upset-facet-irredundancy
script: 04-computation/gmc_one_sided_fixed_reference_hasse_no_go_thm3136.py
output: 05-knowledge/results/gmc_one_sided_fixed_reference_hasse_no_go_thm3136.out
script_sha256: a2742adcf256a974be9947b45ccf9ddf2838ffee107ee667d143fceae96e4de0
output_sha256: b9895befd396a0878416cd41f8cef6cfc62c2978c3fe6092ff6ed6854449cfc4
semantic_sha256: 8a81f4b9cbd76aca6695deb6a685bf698336aa00d83a12eb63af48e2a1aeaa0a
hash_basis: LF-normalized bytes
---

# THM-3136 -- one-sided fixed-reference elementary-tail Hasse no-go

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-3128 rules out repairing the first active pole-prefix hostile while
keeping its labelled-deletion image fixed.  There is a natural transverse
move: shift only the signed bank functional by the pole prefix and hold the
distinguished physical reference alphabet fixed.  This does change the
deletion image, and in degree five it repairs every current in the finite
THM-3120 bank.  The repair does not continue.  A one-pole elementary tail
alternates forever and crosses a load-bearing upset facet in every even
degree.

## 1. The one-sided current

For one THM-3120 signed bank `Phi`, a reduced-pole prefix `R`, and its
distinguished residual alphabet `Q`, put

```text
Phi^R(f)=sum_X epsilon_X f[X-R].                              (1)
```

Only `Phi` is shifted; `Q` is not.  In degree `N`, define the raw Young-gap
coefficient vector

```text
H^R_(N,mu)
 =Phi^R(h_N)m_mu(Q)-Phi^R(m_mu)h_N(Q),      mu partition N.  (2)
```

Since `h_N=sum_(mu partition N)m_mu`, one has

```text
sum_mu H^R_(N,mu)=0.                                           (3)
```

Thus `(2)` lies in the zero-mass hyperplane on which the THM-3127 upward
Hasse-boundary cone lives.  It is Hasse-positive exactly when every
coarsening upset `U` satisfies

```text
sum_(mu in U) H^R_(N,mu)>=0.                                  (4)
```

The bottom-complement upset is

```text
U_bot=P_N\{(1^N)}.                                            (5)
```

By `(3)`, its inequality is the exact ratio test

```text
H^R_(N,1^N)
 =Phi^R(h_N)e_N(Q)-Phi^R(e_N)h_N(Q)<=0.                       (6)
```

## 2. A scalar-positive infinite hostile

Take support `(a,b)=(1,2)`, bank `I_1`, and subtract only its largest
remaining reduced pole `R=(5)`.  The exact reduced scalar row fraction is

```text
sum_(N>=0) Phi(h_N)t^N = t^5 P(t)/D(t),                       (7)

P(t)=36-108t-72t^2,
D(t)=(1-5t)(1-4t)(1-3t)^2(1-2t)^3(1-t)^4.                   (8)
```

Plethystic subtraction by `5` multiplies the complete-homogeneous series by
`1-5t`.  Re-flagging `P` against the remaining descending poles starts with
`4,3` and gives the exact positive identity

```text
P(t)
 =36(1-4t)(1-3t)+144t(1-4t)+72t^2.                           (9)
```

All later denominator poles are positive.  The same suffix-denominator
argument as THM-3120 therefore proves

```text
Phi^(5)(h_N)>0                    for every N>=5.             (10)
```

So the obstruction below is not caused by a scalar-negative prefix.

Now use the elementary generating series

```text
E_X(t)=sum_(N>=0)e_N[X]t^N.                                  (11)
```

Direct exact summation over the `I_1` signed bank gives

```text
Phi(E_X(t))=16t^5+40t^6.                                     (12)
```

Because `E_(X-5)(t)=E_X(t)/(1+5t)`, equations `(1)` and `(12)` yield

```text
sum_(N>=0)Phi^(5)(e_N)t^N=(16t^5+40t^6)/(1+5t),              (13)

Phi^(5)(e_5)=16,
Phi^(5)(e_N)=-40(-5)^(N-6),                 N>=6.            (14)
```

The fixed distinguished alphabet is

```text
Q={1,1,2,3,4,5},                                             (15)

E_Q(t)=1+16t+100t^2+310t^3+499t^4+394t^5+120t^6.            (16)
```

Hence `e_N(Q)=0` for `N>=7`, while `h_N(Q)>0` for every `N`.  For every
even `N>=8`, `(6)` and `(14)` give

```text
H^(5)_(N,1^N)=40*5^(N-6)h_N(Q)>0.                            (17)
```

Thus the upset mass in `(5)` is strictly negative.  Degree six is already
hostile before the finite-length vanishing takes over:

```text
Phi^(5)(h_6)=612,           h_6(Q)=297412,
e_6(Q)=120,                 Phi^(5)(e_6)=-40,

H^(5)_(6,1^6)=11969920>0.                                   (18)
```

Combining `(10)`, `(17)`, and `(18)` proves:

```text
the scalar prefix is positive for every N>=5,
but H^(5)_N is not Hasse-positive for every even N>=6.        (19)
```

This is an all-degree stopping theorem for the one-sided fixed-reference
proposal.

## 3. Complete active-prefix census through degree nine

The exact companion also scans both banks on all `115` THM-3120 supports,
every active reduced-pole prefix, and degrees `5<=N<=9`.  There are

```text
8241 active prefixes,                 41205 currents,
6798825 exact upset sums,             41205 independent max-flow tests. (20)
```

The pass and scalar-sign counts are

```text
N    currents   Hasse pass   scalar positive   positive/pass
5      8241        8241           8241             8241
6      8241        2600           5298             2600
7      8241        2820           7041             2820
8      8241        1731           5638             1731
9      8241        1466           6944             1466
----------------------------------------------------------------
all   41205       16858          33162            16858.       (21)
```

There are `8043` scalar-negative currents and none is Hasse-positive.  This
has a structural explanation.  The common-prefix row marginal still gives
`Phi^R(p_N)=0`, while `m_(N)=p_N`; hence

```text
H^R_(N,(N))=Phi^R(h_N)p_N(Q).                                (22)
```

Since `p_N(Q)>0`, every scalar-negative current fails the top upset
`{(N)}`.  More importantly, `(21)` contains `16304` scalar-positive failures,
so restricting attention to positive scalar flags does not repair the route.

The two extreme tests catch much, but not everything:

```text
all failures                                             24347
top-facet failures                                        8043
bottom-complement failures                               15384
union of the two extreme families                        17849
failures surviving both extremes                          6498
caught there by another principal upset                   6455
requiring a genuinely nonprincipal upset                    43. (23)
```

The `43` cases split as `35` in degree seven and `8` in degree nine.  They
are a live product-Gamma realization of THM-3129's facet irredundancy, not
only an abstract poset warning.

The first such case is again `(a,b)=(1,2),I_1,R=(5)`, now in degree seven.
Both extreme facets and every principal upset pass, but the nonprincipal
upset with minimal antichain

```text
{(3,1,1,1,1),(2,2,1,1,1)}                                  (24)
```

has exact mass

```text
-1151498720.                                                 (25)
```

Its full upset is reconstructed and checked in the companion.

## 4. Degree-five positive control and the deletion-image boundary

Every one of the `8241` degree-five currents in `(20)` is Hasse-positive.
In particular, at THM-3128's first co-shifted hostile

```text
(a,b)=(1,3), I_2, j=5, R=(8,7,6,5,5),                       (26)
```

holding `Q` fixed turns the hostile into a feasible Hasse current.  It does
so by changing the raw labelled-deletion target-`(4)` coordinate from

```text
-1779840  to  470835360.                                    (27)
```

Thus this control does not contradict THM-3128's no-positive-preimage
theorem, which fixes the deletion image.  It is the sharp degree boundary
that suggested the one-sided move, not a mechanism that survives in higher
degree.

## 5. Scope and reconstruction boundary

The current `(2)` is a derived one-sided wedge.  The pole prefix moves the
signed bank while the distinguished alphabet is held fixed.  No identity in
this theorem reconstructs the original nonrow normalized product-Gamma
response from these wedges.  Indeed `(27)` shows that the move changes a
raw labelled-deletion image.  Therefore even one of the `16858` finite positive
currents is not, by itself, a positive coefficient theorem for the original
response.

The theorem proves that scalar pole-flag positivity cannot be lifted by this
particular fixed-reference holotopy to an all-degree Hasse current.  It does
not rule out a different transverse selector, a degree-coupled combination
of prefixes, a larger positivity cone, or a direct arbitrary-support
product-Gamma argument.  It proves no Gaussian moment conjecture, NC2,
LRC(14), JC(2), or DC(2).

## 6. Exact companion

Run

```text
python 04-computation/gmc_one_sided_fixed_reference_hasse_no_go_thm3136.py
python -O 04-computation/gmc_one_sided_fixed_reference_hasse_no_go_thm3136.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_one_sided_fixed_reference_hasse_no_go_thm3136.out.
```

The companion reconstructs `(7)--(18)`, evaluates every current in `(20)`
by a bank-wide exact monomial recurrence, enumerates every upset, and checks
the equivalent max-flow certificate independently.  Only integers and exact
fractions are used.

**QED (candidate pending independent hostile audit).**
