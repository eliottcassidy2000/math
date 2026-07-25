---
id: THM-2216
title: "Residual-capacity hinge/Gram law"
status: >
  PROVED + VERIFIED-EXACT. For every finite labelled residual-cover problem,
  the top-p sum of the pointwise meet of two singleton capacity rows has an
  exact integer Ky--Fan hinge formula. Its tail is an exact positive
  semidefinite superlevel-incidence Gram kernel K_theta=R_theta R_theta^T.
  A lower-dimensional Hellinger Gram kernel gives a rigorous upper bound
  after integer upward rounding. On the scalar depth-(1,1,3) carrier, the
  threshold theta=2612 certifies all 514,605 unordered shallow-blocker
  pairs with minimum exact rounded-Gram certificate margin 78.6095718332.
  This independently recertifies that finite profile and provides a
  uniform truncated-correlation target for the corresponding all-depth
  residual-meet problem. No such all-depth bound is proved here, and this
  is not a proof of LRC(14).
source: klein-2026-07-24-residual-capacity-hinge-gram
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2205-scalar-depth-113-exact-lift-capacity-exclusion
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_depth113_hinge_gram_certificate_thm2216.py
output: 05-knowledge/results/lrc14_depth113_hinge_gram_certificate_thm2216.out
script_sha256: 6b0812be6ba3ad3f84cbc994cb3f8020d1972182b5056215414971abd60d4101
output_sha256: 675b03c1f9c329dbb8711c57445226a41ef1dbf40d7aa4fc499d25e24d3c390c
hash_basis: working-tree bytes (LF)
---

# THM-2216 -- residual-capacity hinge/Gram law

## 1. Abstract residual capacities

Let `X` be a finite phase set, let `(Y_x)_(x in X)` be pairwise disjoint
finite target fibres, and put `Y=disjoint_union_x Y_x`.  Let `Q` be a
finite set of terminal labels, with `M_q subset Y` the vertices covered by
terminal mask `q`.  Define

```text
h_q(x)=|M_q intersection Y_x| in {0,1,...,H}.       (1)
```

For a blocker-active phase set `A subset X`, put

```text
C_A(q)=sum_(x notin A) h_q(x),
F(q)=sum_(x in X) h_q(x),
X_A(q)=sum_(x in A) h_q(x).                          (2)
```

For two blockers `A,B`, write `C_AB=C_(A union B)`.  Then, label by
label,

```text
C_AB(q)<=min(C_A(q),C_B(q)),                         (3)

C_AB(q)
 =F(q)-X_A(q)-X_B(q)+X_(A intersection B)(q)
 <=F(q)-X_A(q)-X_B(q)+H|A intersection B|.          (4)
```

These are exact consequences of (1)--(2); no independence assumption is
made.  Put

```text
Y_AB=disjoint_union_(x notin A union B)Y_x,
W_AB=|Y_AB|.                                         (5)
```

If `p` terminal masks cover `Y_AB`, the ordinary union bound gives

```text
W_AB
 <=sum_(q among the chosen masks)C_AB(q)
 <=Top_p(C_AB),                                      (6)
```

where `Top_p(c)` is the sum of the `p` largest entries of `c`.

## 2. Ky--Fan hinge formula

For every nonnegative integer array `m:Q -> Z` and
`1<=p<=|Q|`,

```text
Top_p(m)
 =min_(theta in Z_{\ge 0})
    [p theta+sum_(q in Q)(m(q)-theta)_+].            (7)
```

Indeed, for every `p`-element set `T`,

```text
sum_(q in T)m(q)
 <=p theta+sum_q(m(q)-theta)_+.                      (8)
```

Taking the maximum over `T` gives one direction.  For the reverse
direction, take `theta` to be the `p`th decreasing order statistic of
`m`.  Every entry strictly above `theta` contributes its excess, and the
remaining places among the first `p` contribute exactly `theta`.
More precisely, if
`m_(1)>=...>=m_(|Q|)` and `m_(|Q|+1)=0`, the complete integer minimizer
set is

```text
m_(p+1)<=theta<=m_(p).                               (9)
```

For singleton capacity rows `C_u,C_v`, define their pointwise meet

```text
m_uv(q)=min(C_u(q),C_v(q))                          (10)
```

and, for an integer `theta>=0`,

```text
x_(u,theta)(q)=(C_u(q)-theta)_+,

K_theta(u,v)
 =sum_q min(x_(u,theta)(q),x_(v,theta)(q)).         (11)
```

Since

```text
(min(a,b)-theta)_+
 =min((a-theta)_+,(b-theta)_+).                     (12)
```

(7) becomes the exact meet-envelope identity

```text
Top_p(m_uv)
 =min_theta [p theta+K_theta(u,v)].                 (13)
```

The word "exact" in (13) concerns the meet envelope.  By (3), that
envelope can still exceed the actual two-blocker capacity row.

## 3. The exact meet-tail Gram kernel

Define the binary superlevel-incidence matrix

```text
R_theta[u,(q,s)]
 =1_{C_u(q)>=theta+s},             s=1,2,... .      (14)
```

Only finitely many columns are nonzero.  Integer layer cake gives

```text
(R_theta R_theta^T)[u,v]
 =sum_(q,s)
    1_{x_(u,theta)(q)>=s}1_{x_(v,theta)(q)>=s}
 =sum_q min(x_(u,theta)(q),x_(v,theta)(q))
 =K_theta(u,v).                                     (15)
```

Thus the entire pair table is the positive semidefinite Gram matrix

```text
K_theta=R_theta R_theta^T.                          (16)
```

Combining (3), (6), and (13), a `p`-mask cover would imply, for every
integer threshold,

```text
W_uv
 <=Top_p(C_uv)
 <=Top_p(m_uv)
 <=p theta+K_theta(u,v).                            (17)
```

Consequently

```text
W_uv>p theta+K_theta(u,v)                           (18)
```

is a rigorous sufficient one-threshold exclusion criterion, exact for the
meet envelope in (13).

## 4. A compact Hellinger upper kernel

The exact matrix (14) has one feature for every label and superlevel.
Define instead

```text
H_theta(u,v)
 =sum_q sqrt(x_(u,theta)(q)x_(v,theta)(q)).         (19)
```

Because `min(x,y)<=sqrt(xy)` for nonnegative `x,y`,

```text
K_theta(u,v)<=H_theta(u,v).                         (20)
```

If

```text
V_theta[u,q]=sqrt(x_(u,theta)(q)),
```

then

```text
H_theta=V_theta V_theta^T.                          (21)
```

This is a different, lower-dimensional Gram kernel, not the exact kernel
in (16).  It can be strictly weaker when the two tail coordinates are
imbalanced.  Equality in (20) holds exactly when, for every label `q`, the
two truncated coordinates are equal or at least one is zero.

It also admits a completely integral upper certificate.  Fix an integer
scale `S` and put

```text
r_S(n)=ceil(sqrt(S^2 n)),
V_int[u,q]=r_S(x_(u,theta)(q)).                     (22)
```

Then

```text
H_theta(u,v)
 <=S^(-2)(V_int V_int^T)[u,v].                     (23)
```

Every square root in (23) is bounded using integer square root, and
every later operation is an integer matrix product.

## 5. The depth-`(1,1,3)` certificate

Specialize to the scalar branch of THM-2192 and THM-2198 with

```text
N=13^4,       Q=13^3.                               (24)
```

The exact primitive carrier has

```text
2028 quotient phases,
1014 depth-one blocker sign classes,
13182 terminal unit sign classes.                   (25)
```

For a blocker class `u`, `C_u(q)` counts the guard-safe dangerous root
sheets of terminal label `q` over phases avoiding `u`.  For a pair
`u,v`, `W_uv` counts all guard-safe root sheets over phases avoiding both.
A terminal mask has at most two active sheets on each thirteen-root
fibre, so these are precisely the labelled capacities required by (1)--(6).

The companion script reconstructs all singleton rows directly from the
strict torsion inequalities and applies (23) with

```text
p=5,       theta=2612,       S=100000.              (26)
```

All `C(1014+1,2)=514605` unordered pairs satisfy

```text
(W_uv-5 theta)S^2-(V_int V_int^T)[u,v]>0.           (27)
```

The unique minimum is at shallow labels `(5,1098)`:

```text
W_uv=13580,
(V_int V_int^T)[u,v]=4413904281668,

minimum numerator=786095718332,
minimum margin=786095718332/10^10
              =78.6095718332>0.                     (28)
```

Thus five terminal masks cannot cover the residual for any shallow pair.
This is an independent fixed-threshold recertificate of the finite
depth-`(1,1,3)` exclusion already proved in THM-2205.

The script uses no floating point and no `assert` for load-bearing checks.
Normal and `python -O` runs are byte-identical.  The singleton-capacity
table digest agrees with the independently frozen table behind THM-2205.

## 6. All-depth target and paper connection

The exact kernel (16) identifies the missing all-depth statistic: not a
family average, but the intersection of the two truncated singleton
capacity stalks.  If a depth parameter `m` admits thresholds `theta_m`
such that, uniformly in every blocker pair,

```text
W_uv>=beta 13^m-o(13^m),
theta_m<=alpha 13^m+o(13^m),
K_(theta_m)(u,v)<=gamma 13^m+o(13^m),
beta>p alpha+gamma,                                 (29)
```

then (18) excludes every sufficiently large depth.  The stronger
replacement of `K` by the compact `H` in (29) also suffices.

This is the precise connection to the structured `XX^T` viewpoint routed
in
`05-knowledge/reference/CORE-PAPERS-COMPOSITIONAL-RELATIONS.md`:

```text
source object: truncated singleton-capacity stalks,
map:           superlevel incidence or square-root embedding,
target object: a structured Gram matrix,
preserved:     the pairwise tail correlation controlling Top_p,
lost by H:     the exact min-tail geometry,
needed input:  a uniform off-diagonal/correlation estimate,
decisive test: compare (29) with the residual coefficient beta.
```

Fast Gram multiplication can accelerate a finite census, but speed alone
does not prove (29).  Positive semidefiniteness alone also gives no useful
uniform off-diagonal bound.

## 7. Scope

The theorem proves the abstract hinge and Gram identities at every finite
depth and verifies one complete depth-four certificate.  It does **not**
prove that the fixed threshold `2612` scales to later depths, that
Hellinger loss stays bounded, or that every scalar valuation profile with
deepest valuation at least four is empty.  LRC(14) remains open.  QED.
