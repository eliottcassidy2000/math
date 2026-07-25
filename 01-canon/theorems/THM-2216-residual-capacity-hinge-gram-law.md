---
id: THM-2216
title: "Residual-capacity hinge/Gram law"
status: >
  PROVED + VERIFIED-EXACT. For every finite labelled residual-cover problem,
  the top-p sum of the pointwise meet of two singleton capacity rows has an
  exact integer Ky--Fan hinge formula. Its tail is an exact positive
  semidefinite superlevel-incidence Gram kernel K_theta=R_theta R_theta^T.
  The same embedding has an integral l1 cut metric D_theta, and exact
  polarization recovers K_theta from D_theta plus its diagonal energy.
  This yields a metric-separation form of the exclusion criterion and
  rigorous landmark lower certificates for D_theta.
  A lower-dimensional Hellinger Gram kernel gives a rigorous upper bound
  after integer upward rounding. On the scalar depth-(1,1,3) carrier, the
  threshold theta=2612 certifies all 514,605 unordered shallow-blocker
  pairs with minimum exact rounded-Gram certificate margin 78.6095718332.
  At the same threshold, seven farthest-first cut-metric landmarks certify
  every pair through triangle inequalities; the first landmark alone
  certifies 513,548 pairs.
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
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2238-terminal-cover-moment-hinge-gram-and-fibrewise-convex-hierarchy
script: 04-computation/lrc14_depth113_hinge_gram_certificate_thm2216.py
output: 05-knowledge/results/lrc14_depth113_hinge_gram_certificate_thm2216.out
script_sha256: 48957157479a220fec1015e12e085267cab36a735a8e5005f15df7aba3f9c1fe
output_sha256: 1d5c5cba196c62c769a32b1e132b3c2bcf7cc557f18d61cc610a8c48b08df2cc
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
made.  The gap between the singleton meet and the actual pair row is also
exactly directed:

```text
C_A(q)-C_AB(q)=X_(B\A)(q),
C_B(q)-C_AB(q)=X_(A\B)(q),

min(C_A(q),C_B(q))-C_AB(q)
 =min(X_(B\A)(q),X_(A\B)(q)).                       (4a)
```

Thus the meet envelope loses precisely the cheaper of the two directed
exclusive blocker losses, not an unspecified correlation error.  Put

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

### Cut-metric polarization and landmark certificates

Write the same superlevel row as a vector

```text
Phi_theta(u)_(q,s)
 =1_{theta<s<=C_u(q)}.                              (16a)
```

Its diagonal energy, inner product, and squared distance are

```text
E_theta(u)
 =||Phi_theta(u)||_2^2
 =sum_q x_(u,theta)(q)
 =K_theta(u,u),

K_theta(u,v)=<Phi_theta(u),Phi_theta(v)>,

D_theta(u,v)
 =||Phi_theta(u)-Phi_theta(v)||_2^2
 =sum_q |x_(u,theta)(q)-x_(v,theta)(q)|.            (16b)
```

Because the features are binary, `D_theta` is simultaneously a squared
Hilbert distance, an integral `l1` metric, and a nonnegative sum of
elementary cut semimetrics. Exact polarization gives

```text
K_theta(u,v)
 =[E_theta(u)+E_theta(v)-D_theta(u,v)]/2.           (16c)
```

This is precisely THM-2221's context-cut construction applied to the
superlevel features `(q,s)`, with owners as the core labels. The two
theorems are therefore two views of one embedding, not merely analogous:
THM-2216 reads its Gram overlap, while THM-2221 reads its Hamming
transport response. The diagonal energy is a necessary sidecar.
Complementing a feature column preserves its cut metric but can change
both `E_theta` and `K_theta`.

Substitution into (13) yields the exact metric form

```text
Top_p(m_uv)
 =min_theta [
    p theta
    +(E_theta(u)+E_theta(v)-D_theta(u,v))/2
   ].                                               (16d)
```

Consequently the one-threshold exclusion criterion is equivalently

```text
D_theta(u,v)
 >E_theta(u)+E_theta(v)-2(W_uv-p theta).            (16e)
```

Unlike a raw overlap ordering, `D_theta` admits metric acceleration. For
any landmark owner `a`, the triangle inequality gives

```text
D_theta(u,v)
 >=|D_theta(u,a)-D_theta(v,a)|.                     (16f)
```

Thus a landmark bank `A` already certifies exclusion whenever

```text
max_(a in A)|D_theta(u,a)-D_theta(v,a)|
 >E_theta(u)+E_theta(v)-2(W_uv-p theta).            (16g)
```

This is a rigorous route from a small number of owner--landmark distances
to many pair exclusions. It does not assert that a uniformly small
landmark bank exists at every depth.

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

### Finite threshold search and comparison with marginal order statistics

For fixed rows `C_u,C_v`, define

```text
B_H(theta)
 =p theta+sum_q
   sqrt((C_u(q)-theta)_+(C_v(q)-theta)_+).          (21a)
```

On an interval containing no coordinate value from either row, each
nonzero summand has the form

```text
sqrt((a-theta)(b-theta)).
```

Its second derivative is

```text
-(a-b)^2 /
 [4((a-theta)(b-theta))^(3/2)]<=0.                 (21b)
```

Hence `B_H` is concave on every such interval.  Its global minimum over
`theta>=0` is therefore attained at

```text
theta in {0} union {C_u(q):q in Q}
                   union {C_v(q):q in Q}.          (21c)
```

This makes optimized exact Hellinger certification a finite breakpoint
search; it does not require a continuous optimization.

There are now three rigorous meet envelopes:

```text
Top_p(min(C_u,C_v)),
min_theta B_H(theta),
sum_(k=1)^p min(C_u,(k),C_v,(k)).                  (21d)
```

The first is at most each of the latter two.  The two upper bounds are
incomparable.  For `p=2`, padding with zero coordinates if desired,

```text
C_u=(0,1), C_v=(1,0): Hellinger=0 < marginal=1;
C_u=(0,1), C_v=(0,2): marginal=1 < Hellinger=sqrt(2). (21e)
```

Thus one should take the cheapest certified envelope row by row rather
than treating either relaxation as universally stronger.

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

The exact cut metric supplies a much smaller second certificate. Start
with shallow label `1`, and repeatedly add an owner maximizing its minimum
`D_theta`-distance from the current bank, using the first index to break
ties. The first seven landmark labels are

```text
[1,183,799,244,1098,659,824].                       (28a)
```

After each addition, the numbers of unordered pairs certified by (16g)
are

```text
[513548,514564,514588,514597,514603,514603,514605]. (28b)
```

Thus seven owner--landmark distance columns certify all `514605` pairs.
This avoids an all-pairs computation of the `13182`-coordinate truncated
label overlap; the separate residual-size table `W_uv` is still retained.
At the binding Hellinger pair `(5,1098)`, direct polarization gives

```text
(E_5,E_1098,K,D)=(654,1458,302,1508),
2K=E_5+E_1098-D.                                    (28c)
```

The landmark result is a finite exact compression theorem for this
carrier, not yet an all-depth uniform landmark bound.

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
