---
id: THM-3886
title: "Cusp residual equality-seam second-layer trichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On
  the THM-3884 equality seam, the next filtration layers have three exact
  regimes.  For deg f>=3 they force a second full lift-gauge jet.  For
  deg f=2 the first collision leaves one Kummer-twisted symbol, but the next
  two homogeneous equations have incompatible x-adic orders, so the cell is
  empty.  For deg f=1 homogeneous square-root recursion forces the entire
  pair into a forbidden constant pure gauge, so that cell is empty.  The
  stable-range conclusion is associated-graded; no gauge-invariant descent
  or termination is claimed.
source: jc_quartic_c3_construct / post-THM-3884 second-layer filtration, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS, ONE WITH PROOF-SCOPE REPAIR.
  tournament-jc-breakthrough independently rederived the five-vertex
  Plucker gain packet and complete tie-depth envelope, including the n>=3
  collision exhaustion, second gauge jet, n=2 x-adic orders 8 and 3 for both
  Kummer signs, and the n=1 square-root recursion.  A separate sparse-
  convolution audit found that the original proof silently imposed the
  THM-3881 address in degrees one and two although the statement did not.
  Retaining the arbitrary constants repairs and strengthens both closures:
  the quadratic witnesses are unchanged, while the linear cell dies at a new
  degree-five obstruction.  Together the audits cover all seven stable degree
  envelopes, the associated-graded-only scope, and an explicit hostile showing
  that square survival is not invariant under gauge peeling.  Normal,
  optimized, and frozen streams byte-match their frozen transcripts.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
related:
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
script: 04-computation/jc2_cusp_residual_equality_seam_trichotomy_thm3886.py
output: 05-knowledge/results/jc2_cusp_residual_equality_seam_trichotomy_thm3886.out
script_sha256: 83898bc719a4fd2e520fcb002018347f1294a77d88db9d124b9a6c06fe901853
output_sha256: dc18c004abfdedf53e7910811b05c289c83c0b011dffb3208362a36515561a2a
semantic_sha256: f48033db03f9afa6de04998aa497a50260fa8271c720a10b0f5bc208a341da33
independent_audit_script: 04-computation/jc2_cusp_residual_equality_seam_trichotomy_independent_audit_thm3886.py
independent_audit_output: 05-knowledge/results/jc2_cusp_residual_equality_seam_trichotomy_independent_audit_thm3886.out
independent_audit_script_sha256: 5f4bc278a6a2d54167aa1337489409ac7d979e0f7474ef29736d323e9222c55b
independent_audit_output_sha256: ec5fb28b2e9869fc056250d269ed5499390dab9ecfe00199a89fee2f36230200
independent_audit_semantic_sha256: ab9e880e0d288913aea88525c79ac4202129e605fa69970590ab2b094c144543
hash_basis: raw LF bytes
---

# THM-3886 -- the equality seam closes in degrees one and two

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Use the
THM-3881 notation

```text
a=x+1,                       L=9x+4,
K=K_2+K_1+K_0,
K_2=y^2-15x^2,               K_1=-15x,             K_0=-4,  (1)
```

and its residual `S(T,f)`.  Suppose `f` is nonconstant, `S(T,f)` is a square,
and the THM-3884 equality case occurs:

```text
deg f=n>=1,                  deg T=n+1,
f_n=xq,                      T_{n+1}=-K_2q,                 (2)
```

where `q` is nonzero homogeneous of degree `n-1`.  Put

```text
s=f_{n-1},                   t=T_n,
R=K_2s+xt-y^2q.                                             (3)
```

Then the next layer obeys the following exact trichotomy.

1. If `n>=3`, then `R=0`.  Equivalently there is homogeneous
   `p` of degree `n-2` with

   ```text
   s=q+xp,                  t=-K_2p+15xq.                  (4)
   ```

   Thus the first two homogeneous layers of `(T,f)` are the negative of the
   first two layers of the full gauge vector `(K,-a)(q+p)`.

2. If `n=2`, the first collision necessarily has nonzero constants
   `u,epsilon in k` such that

   ```text
   q=xu^2,                  R=epsilon*x^3*u,
   epsilon^2=-216.                                          (5)
   ```

   but this Kummer-twisted symbol cannot lift through the next two homogeneous
   equations.  Hence no square survivor exists.

3. If `n=1`, no square survivor exists.

## 1. The intrinsic next symbol

Let

```text
r=aT+Kf,                    B=aL^2f^2-T^2.                 (6)
```

The leading degree-`n+2` part of `r` vanishes by `(2)`.  Its next part is

```text
r_{n+1}
 =xT_n+T_{n+1}+K_2f_{n-1}+K_1f_n
 =K_2s+xt-y^2q
 =R.                                                        (7)
```

Also

```text
B_{2n+3}=x(9x)^2(xq)^2=81q^2x^5.                          (8)
```

The quartic part `3r^2B` of the residual therefore contributes

```text
243q^2x^5R^2                                                  (9)
```

in degree `4n+5`.

The only possible competitor at this level is the cubic-coordinate term
`8AB`, where `A=KT+a^2L^2f` has degree `n+4`.  Its degree is `3n+7`.
Consequently:

```text
n>=3:       4n+5>3n+7,
n=2:        4n+5=3n+7,
n=1:        4n+5<3n+7.                                    (10)
```

All terms with fewer coordinate factors have still smaller degree.  This is
the source of the trichotomy, rather than an exceptional finite computation.

## 2. Stable range: a second gauge jet

For `n>=3`, degrees `4n+7` and `4n+6` vanish after `(2)`, and `(9)` is the
complete degree-`4n+5` form:

```text
S_{4n+5}=243q^2x^5R^2.                                    (11)
```

If `R` were nonzero, its valuation at the prime `x` in `(11)` would be
`5+2v_x(q)+2v_x(R)`, which is odd.  The highest form of a square must be a
square, so `R=0`.  Reducing `(3)` modulo `x` gives

```text
x divides s-q.                                             (12)
```

Write `s=q+xp`; substituting back gives `(4)`.  Direct expansion of the full
gauge vector confirms

```text
a(q+p) = xq+(q+xp)+lower,
-K(q+p)=-K_2q+(-K_2p+15xq)+lower.                         (13)
```

This is a two-jet alignment only.  Applying the gauge changes the residual
through its `Delta`-debt and is not asserted to carry a survivor to a
survivor.

## 3. The exceptional quadratic Kummer seam does not lift

For `n=2`, the leading part of `8AB` is

```text
52488q^3x^10=243q^2x^5(216qx^5).                          (14)
```

Hence the complete first nonzero form is

```text
S_13=243q^2x^5(R^2+216qx^5).                              (15)
```

Its degree is odd, so a square survivor requires

```text
R^2=-216qx^5.                                              (16)
```

Unique factorization now gives `(5)`: away from `x`, every valuation of `q`
is even, while `v_x(q)+5` is even.  Thus `q=xu^2` after absorbing a square
unit.  Since `q` has degree one, `u` is a nonzero constant, and `(16)` gives
`R=epsilon*x^3u`, `epsilon^2=-216`.  Both signs remain live at this filtration
level.

They die at the next level.  Write the remaining homogeneous pieces as

```text
f_1=A x,                 f_0=b,
T_2=(u^2-A)y^2+(15A+epsilon*u)x^2,
T_1=m x+rho y,           T_0=w,                          (17)
```

Indeed, substituting general linear `f_1` and quadratic `T_2` in
`R=epsilon*x^3u` and comparing the four degree-three monomials gives exactly
the first two formulas in `(17)`; the constant `w` is unrestricted.
After imposing `epsilon^2=-216`, direct homogeneous expansion gives

```text
S_13=0,
v_x(S_12)=8,             [x^8 y^4]S_12=-648u^6,
v_x(S_11)=3,             [x^3 y^8]S_11=8u^6.             (18)
```

The displayed coefficients are nonzero and are the unique terms at their
respective minimum `x`-orders, independently of `A,b,m,rho,w`.  If
`S=(g_6+g_5+...)^2` in homogeneous pieces, then `S_12=g_6^2`; hence
`v_x(g_6)=4`.  The next equation `S_11=2g_6g_5` would force
`v_x(S_11)>=4`, contradicting `(18)`.  Thus the apparent Kummer seam is a
first-collision artifact, not a liftable square symbol.

## 4. The linear cell closes by address-free square-root recursion

Let `n=1`; now `q in k^*`.  Keep the constant term unrestricted:

```text
f=qx+c,
T=-qK_2+ux+vy+w.                                          (19)
```

The residual has degree ten and

```text
S_10=52488q^3x^10.                                        (20)
```

Choose `mu in k^*` with `mu^2=52488q^3`.  If
`S=(g_5+g_4+...+g_0)^2` in homogeneous pieces, then

```text
g_5=mu*x^5.                                                (21)
```

The degree-nine equation uniquely determines `g_4`, because `S_9` is
divisible by `x^5`.  At degree eight, polynomiality of

```text
g_3=(S_8-g_4^2)/(2g_5)                                    (22)
```

requires divisibility by `x^5`.  After clearing the nonzero scalar
denominator, the coefficients at `y^8` and then at `x^4y^4` are respectively

```text
-9q(q-c)^4/32,                    -9qv^4/32.               (23)
```

Thus `c=q` and `v=0`.  Write `u=15q+d`; then

```text
f=aq,                  T=-Kq+dx+e,          e=w-4q.       (24)
```

The degree-seven equation uniquely determines `g_2` only if
`S_7-2g_4g_3` is divisible by `x^5`.  Its complete part below `x^5` is

```text
-(d^4q/288)x^3y^4.                                        (25)
```

Therefore `d=0`.  The unrestricted endpoint is

```text
(T,f)=(-Kq+e,aq),                         e in k.           (26)
```

For either choice of the leading root `g_5`, the equations in degrees nine
through six successively define polynomial `g_4,g_3,g_2,g_1`.  At degree
five, defining `g_0` would require divisibility by `x^5`, but the complete
part of its numerator below that order is

```text
(S_5-2g_4g_1-2g_3g_2)|_{v_x<5}
   =(243q/2) x y^2(y^2-30x^2).                            (27)
```

This is nonzero and independent of `e`, so the required division is
impossible.  This closes `n=1` without an address hypothesis.  At the
special value `e=0`, equation `(26)` is the constant pure gauge
`(KQ,-aQ)`, `Q=-q`; the THM-3881 divisibility gate supplies an independent
positive control for that one endpoint.

## 5. Gauge peeling really is non-invariant

The warning after `(13)` has a minimal explicit witness.  The base pair

```text
(T,f)=(0,0)
```

has square residual `L^4`.  Adding the unit gauge vector gives

```text
(T,f)=(K,-a),                    aT+Kf=0.
```

At `a=0`, equivalently `x=-1`, its residual specializes to

```text
625-8(y^2-4)^4.                                           (28)
```

The derivative is `-64y(y^2-4)^3`; its only possible common roots with
`(28)` would have `y=0` or `y^2=4`, where `(28)` equals `-1423` or `625`.
Thus `(28)` is squarefree and nonconstant, hence not a square.  A polynomial
square would remain a square after specialization.  Gauge peeling therefore
does not preserve the square equation even though it preserves the mixed
coordinate `aT+Kf`.

## 6. Scope

The theorem classifies the first layer after THM-3884 and, in the exceptional
small degrees, the additional layers needed for closure.  It closes
`deg f=1,2` and shows that all `deg f>=3` candidates begin with two gauge jets.
The stable `deg f>=3` implication uses only the seven degree envelopes and
the nonzero coefficient `243`; it remains valid over any integral domain of
characteristic different from three.  The exceptional degree-one and
degree-two closures retain the stated algebraically closed
characteristic-zero hypotheses.
It does not prove that gauge peeling preserves a square or that repeated
alignment terminates.  The `deg T>=deg f+2` region, alternate square/cube
lifts, a Keller atlas, and JC(2) remain **OPEN**.

Reproduce with

```bash
python3 04-computation/jc2_cusp_residual_equality_seam_trichotomy_thm3886.py
python3 -O 04-computation/jc2_cusp_residual_equality_seam_trichotomy_thm3886.py
python3 04-computation/jc2_cusp_residual_equality_seam_trichotomy_independent_audit_thm3886.py
python3 -O 04-computation/jc2_cusp_residual_equality_seam_trichotomy_independent_audit_thm3886.py
```

Each pair of streams must byte-match its corresponding frozen output in
`05-knowledge/results/`.  The independent companion uses its own sparse
convolution engine, retains the unrestricted constants in both exceptional
cells, exhausts all seven stable degree envelopes, and checks the gauge
hostile.  **QED.**
