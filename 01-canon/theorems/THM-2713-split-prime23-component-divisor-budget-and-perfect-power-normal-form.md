---
id: THM-2713
title: "Split prime-23 component divisor budget and perfect-power normal form"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every eta!=0
  and every b,c,d,e,w, the split even-Faber prime-23 weighted complete
  intersection is geometrically integral, not only generically integral.
  Its five old and three new branches exhaust the zero divisors of F1 and
  zeta: div_0(F1)=23O, div_0(zeta)=23N, and the rational chosen-sheet
  coordinate has div(q)=5N-3O.  A rational normalization would have a
  necessary perfect-power parametrization with binary-form degrees
  3,5,23,46.  Rationality is not excluded; the residual is exactly 89 extra
  delta units plus the displayed polynomial identities.  This is only the
  nonzero-lambda split even-Faber subchart, not the odd-seed bank, the full
  split branch, or JC(2).
source: thm2704-hostile-audit-2026-07-28-component-budget
audit: lrc-narrow-debt-queue-2026-07-28
extension_audit: coordinate-first-audit-2026-07-28-arbitrary-b-scale
depends_on:
  - THM-2704-split-even-prime23-generic-genus-eighty-nine
related:
  - THM-2636-degree-twenty-two-BCD-triple-spectral-square-Kummer-closure
  - THM-2692-degree-twenty-two-full-support-terminal-toric-order-twelve-closure
script: 04-computation/jc2_split_prime23_component_divisor_budget_20260728.py
output: 05-knowledge/results/jc2_split_prime23_component_divisor_budget_20260728.out
script_sha256: 81cb01e4dff454fc18417b3cfa2136fcab97901317d2b8d2df805a5a0c677350
output_sha256: 483754f20a492b3a77892cfea2f6c06612459d4366dc52647ce3d07cd4b9468b
extension_script: 04-computation/jc2_split_prime23_bscale_extension_20260728.py
extension_output: 05-knowledge/results/jc2_split_prime23_bscale_extension_20260728.out
extension_script_sha256: 9f8d8d6d14eff6d121f5ac45f070d4489d71a72f788ceb5f73db75c4fb63d216
extension_output_sha256: 4ef3d75c3bf6e82f05cefc7a133b1efe4d5b79c87f0e165e1c3bcbc9c2f2460f
hash_basis: working-tree bytes (LF)
---

# THM-2713 -- the split prime-23 curve cannot break its `3:5` divisor balance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The theorem
includes the arbitrary-`b` scale extension of THM-2704's exact nonzero-`eta`
curve, including `b=0`, and does not exclude its abstract rational
specializations.

## 1. Statement

Give `(h,t,v,zeta)` weights `(1,1,2,3)`, and let

```text
C=C_(b,c,d,e,w,eta):
F2_b=0,
zeta F1_b^4=eta t^23,                     eta!=0,      (1)
```

be the split even-Faber complete intersection of degrees `(6,23)` in
`P(1,1,2,3)`.  The parameters `b,c,d,e,w` are arbitrary complex numbers.
On the affine chart `h=1`, the arbitrary-`b` fluxes are defined from the
previously normalized `b=1` pair by

```text
f1_b-f1_1=-616t^2(b-1)(4840v-1331zeta-40),

f2_b-f2_1=49280t^2(b-1)
                    (29282v^2-1452v+1331zeta+2).      (1a)
```

Here `F1_b,F2_b` denote their weighted homogenizations.  Equivalently, choose
any `rho in C*` and put

```text
B_0=b rho^2, C_0=c rho^3, D_0=d rho^4,
E_0=e rho^5, W_0=w rho^6,
t=rho/y, v=u/y^2, zeta=Z/y^3.                         (1b)
```

The chosen-sheet equation then has

```text
eta=7496192^4 lambda^4/rho^23!=0.                     (1c)
```

Thus `rho` is a freely chosen nonzero scale; it need not be a square root of
`B_0`, and `b=0` is allowed.
Then:

1. `C` is a reduced, geometrically integral projective curve for every such
   parameter value.
2. Let `nu:X->C` be its normalization.  Write `O` for the reduced sum of the
   five old `L5` points over `t=0`, `N` for the reduced sum of the three new
   `G3` points, and `L=nu^*O_C(1)`.  Then

   ```text
   deg L=23,
   div_0(t)=4O+N,
   div_0(F1)=23O,
   div_0(zeta)=23N.                                  (2)
   ```

   After harmless nonzero normalization of the chosen split sheet,

   ```text
   q=-t^5/F1,                 div(q)=5N-3O.            (3)
   ```

3. Every member is an lci/Gorenstein curve of arithmetic genus `254`.  The
   five fixed `(4,23)` cusps contribute delta `165`, so, if `Delta_extra`
   denotes the total delta away from them,

   ```text
   g(X)=89-Delta_extra.                               (4)
   ```

   In particular a rational specialization must have
   `Delta_extra=89` exactly.
4. If `X=P1`, there are squarefree coprime binary forms `alpha,beta` of
   degrees `3,5`, a binary form `H` of degree `23` coprime to `alpha*beta`,
   a binary form `V` of degree `46` (possibly the zero section), and nonzero
   constants such that
   the normalized map has the necessary form

   ```text
   (h,t,v,zeta)=(H, tau*alpha*beta^4, V, sigma*alpha^23),
   F1_b=phi*beta^23,           sigma*phi^4=eta*tau^23. (5)
   ```

   On the affine chart this becomes, up to nonzero constants,

   ```text
   q=alpha^5/beta^3,
   t=alpha*beta^4/H,
   zeta=alpha^23/H^3,
   f1_b=beta^23/H^5.                                 (6)
   ```

The forms in (5) must still satisfy the exact identities `F2_b=0` and the
displayed `F1_b` equation.  No assertion here says that those identities have
no solution.  In this rational case, `t` and `q` have degrees `23` and `15`;
after the visible fibres, Riemann--Hurwitz leaves exactly `29` and `6`
additional ramification units, respectively.

## 2. The curve is a pure complete intersection on the smooth ambient locus

For readability, suppress the subscript `b` on `F1_b,F2_b` throughout the
geometric proof.  Formula (1a) shows that every new `b`-term contains `t^2`.
It therefore changes neither ambient coordinate-point value nor any datum on
the fixed fibre `t=0`.

The only singular points of `P(1,1,2,3)` are the `v`- and `zeta`-coordinate
points.  Uniformly in the parameters, `F2` takes the nonzero values

```text
-1190488992,                     15944049             (7)
```

there.  Hence `C` lies in the smooth ambient locus and `O_C(1)` is a line
bundle.

The exact corner calculation inherited from THM-2704 gives

```text
{F2=zeta F1^4=h=t=0}=empty.                           (8)
```

Thus the two weight-one sections `h,t` have no common zero on `C`.

The equations in (1) have no common hypersurface component.  Indeed, a
projective hypersurface component common to them would meet the positive
weight hyperplane `t=0` in a curve.  But the full intersection with `t=0`
is the finite affine `G3+L5` support below, and (8) excludes a component at
the corner.  Therefore `(F2,zeta F1^4-eta t^23)` is a regular sequence.
The curve is lci, Cohen--Macaulay, pure one-dimensional, and has no embedded
associated point.

The basepoint-free pencil defines

```text
pi:C -> P1,                         p |-> [h(p):t(p)], (9)
```

with `O_C(1)=pi^*O_P1(1)`.  Its restriction to every geometric irreducible
component is nonconstant: otherwise the pullback in (9) would be trivial on
that component, contradicting ampleness and positive degree of `O_C(1)`.
This excludes both affine vertical components and components contained in
`h=0` without a parameter-generic infinity assumption.

## 3. The eight fixed local domains

Every parameter-bearing term of `F1,F2`, including the arbitrary-`b`
deformation (1a), contains `t`.  Consequently the fibre `t=0` and all its
unit checks are parameter-independent.  Its reduced support is

```text
new:  zeta=0, G3(v)=0,                 three points,
old:  F1=0,    L5(v)=0,                five points,   (10)
```

where `G3,L5` are squarefree and coprime.

At a new point, `F1` and the relevant `v`-derivative of `F2` are units.  The
local curve is

```text
zeta=unit*t^23.                                      (11)
```

It is regular, with one branch and

```text
ord(t)=1,                  ord(zeta)=23, ord(F1)=0.   (12)
```

At an old point, `zeta` and
`det partial(F1,F2)/partial(v,zeta)` are units.  On the smooth local surface,
`t,F1` are coordinates and the curve is

```text
F1^4=unit*t^23.                                     (13)
```

Because `gcd(4,23)=1`, (13) is an analytically irreducible reduced plane
branch.  On its normalization,

```text
ord(t)=4,                  ord(F1)=23, ord(zeta)=0.   (14)
```

Thus every point in (10) belongs to exactly one geometric global component;
two components cannot share one of these local domains.

## 4. Opposing section degrees force every component to contain everything

Let `Y` be the normalization of one geometric irreducible component.  Since
`pi|Y` is nonconstant and projective, it is surjective.  Hence `Y` contains
some points from (10).  Let

```text
r=# old points on Y,              s=# new points on Y.       (15)
```

No other point lies over `t=0`, so (12)--(14) give

```text
d=deg(pi|Y)=deg(L|Y)=4r+s.                            (16)
```

Neither `F1` nor `zeta` vanishes identically on `Y`: by (1), either would
force the nonzero generic section `eta*t^23` to vanish.  They are sections
of `L^5` and `L^3`, respectively.  Their visible zeros in (12)--(14) imply

```text
23r <= deg div_0(F1|Y)=5d,
23s <= deg div_0(zeta|Y)=3d.                         (17)
```

Substituting (16) gives the opposing inequalities

```text
3r <= 5s,                         5s <= 3r.           (18)
```

Therefore `3r=5s`.  Since `0<=r<=5`, `0<=s<=3`, and `(r,s)!=(0,0)`, the
only possibility is

```text
(r,s,d)=(5,3,23).                                    (19)
```

Every geometric component would have to contain all eight fixed branches.
The local-domain uniqueness after (14) therefore permits only one component.

That component is generically reduced because it passes through the regular
new points (11).  The lci curve is Cohen--Macaulay, hence satisfies `S1` and
has no embedded associated points.  A generically reduced `S1` curve is
reduced.  This proves uniform geometric integrality.

This is the central mechanism: the naive branch grammar allows every degree
`1,...,23`, but the degree-three `zeta` section and degree-five `F1` section
push in opposite directions and leave only the full `3:5` packet.

## 5. Exact divisors and arithmetic genus

At (19), both inequalities in (17) are equalities:

```text
deg_visible div_0(F1)=23*5=5*23=deg div_0(F1),
deg_visible div_0(zeta)=23*3=3*23=deg div_0(zeta).    (20)
```

Thus the visible zeros exhaust the two sections, proving (2).  The rational
chosen-sheet reconstruction from the split first flux is, after rescaling by
a nonzero constant,

```text
q=-t^5/F1.                                           (21)
```

More precisely, before removing the scalar it is

```text
q_phys=-7496192 lambda t^5/(rho^5 f1_b).             (21a)
```

This formula is valid for every `b=B_0/rho^2`, including zero.

Equations (2) and (21) give

```text
div(q)=5(4O+N)-23O=5N-3O,                            (22)
```

so `q` has exactly three zeros of order five and five poles of order three,
with no additional zero or pole.

Because the equations are a regular sequence of degrees `(6,23)` wholly in
the smooth ambient locus, weighted adjunction applies for every parameter:

```text
omega_C=O_C(22),
deg O_C(1)=6*23/(1*1*2*3)=23,
p_a(C)=1+22*23/2=254.                                (23)
```

The five local branches (13) have

```text
delta=(4-1)(23-1)/2=33                               (24)
```

each.  Since `C` is integral, normalization gives (4).  In particular the
exceptional rational locus is no longer allowed to escape through
reducibility or nonflatness: it is exactly the locus carrying `89` additional
delta units.

## 6. The rational residual is a `3,5,23` perfect-power map

Assume only for this section that `X=P1`.  Pulling back `O_C(1)` gives a
degree-`23` line bundle, hence `L=O_P1(23)`.  Choose squarefree homogeneous
binary forms

```text
alpha, beta,                     deg(alpha,beta)=(3,5), (25)
```

whose zero divisors are `N,O`.  They are coprime because (10) is disjoint.
The zero-divisor identities (2) force, up to nonzero scalars,

```text
t=alpha*beta^4,
F1=beta^23,
zeta=alpha^23.                                      (26)
```

The remaining weight-one section is a degree-`23` binary form `H`.  The
basepoint-free statement (8) gives `gcd(H,alpha*beta)=1`.  The weight-two
coordinate pulls back to a binary form `V` of degree `46` (with the zero form
allowed).  Restoring
the scalars in (26) gives (5), and (1) gives their single displayed scalar
relation.  Dividing by the appropriate powers of `H` and using (21) yields
(6).

This form is necessary, not asserted sufficient.  The exact residual is the
simultaneous binary-form system

```text
F2(H,tau*alpha*beta^4,V,sigma*alpha^23)=0,
F1(H,tau*alpha*beta^4,V,sigma*alpha^23)=phi*beta^23. (27)
```

It retains the degrees `(3,5,23,46)`, all scalar units, and the omitted third
flux as a future sidecar.  It is a much smaller target than an arbitrary
genus-drop discriminant, but this theorem does not solve it.

There are also two sharp ramification invoices.  The degree-`23` map `t` has
visible contribution `5*(4-1)=15`; genus zero would leave `44-15=29` units.
By (22), `q` has degree `15`, and its zero and pole fibres contribute

```text
3*(5-1)+5*(3-1)=22.                                  (28)
```

Riemann--Hurwitz gives total ramification `28`, leaving exactly `6` other
units.  Thus the rational residual is simultaneously a `t`-cover with a
`29`-unit residual and a `q`-cover with only a `6`-unit residual.

## 6A. Independently audited arbitrary-`b` addendum

The extension companion derives (1a)--(1c) from the exact unscaled integer
fluxes, not by interpolating the normalized family.  Since both differences
in (1a) contain `t^2`, the following inputs to Sections 2--6 are literally
independent of `b`:

```text
the two ambient coordinate-point values;
the h=t=0 corner and basepoint-free pencil;
the complete t=0 support G3+L5 and every local unit;
the old/new orders (4,23,0) and (1,0,23);
the weights deg(F1,zeta)=(5,3).                       (28a)
```

The regular-sequence proof uses the finite fixed fibre plus the empty corner.
The component proof uses only (28a), ampleness of `O_C(1)`, and the opposing
section budgets.  The reducedness and adjunction arguments then use only the
resulting lci/integral curve and its unchanged degrees.  Therefore every
conclusion of the theorem transfers to arbitrary `b`, including `b=0`.

The full `h=0` face away from `t=0` is **not** asserted independent of `b`:
some homogenized terms from (1a) survive there when `t!=0`.  No description
of that face is needed.  Basepoint freeness of `[h:t]` and ampleness exclude
vertical or `h=0` components uniformly.

The independent audit replayed both execution modes, checked the exact
deltas and scale cancellation, and re-audited every dependency in (28a).
This is a theorem-level extension of the response geometry, not a claim about
the separately omitted `y=0` or odd-Faber physical branches.

## 7. The normalized `b=1` norm/trace sidecar

This section records the exact norm/trace sidecar on the original normalized
section `b=1`.  It is **not** used in the component, divisor, genus, or
perfect-power arguments above, and the arbitrary-`b` addendum makes no claim
that the displayed term counts remain unchanged away from `b=1`.

Write the affine fluxes as

```text
f1=a*zeta+b,
f2=A*zeta^2+B*zeta+C0,                  A=15944049.   (29)
```

Set

```text
D=A*b^2-B*a*b+C0*a^2=Res_zeta(f2,f1).                (30)
```

The universal quadratic norm identity is

```text
Res_zeta(f2,zeta*f1^4-X)
 =C0*D^4-X*T+A^5*X^2,                                (31)
```

where `T=A^5 Tr(zeta*f1^4)` in the quadratic algebra of `f2`.  Substituting
`X=eta*t^23` gives

```text
R(t,v)=C0*D^4-eta*t^23*T+A^5*eta^2*t^46.             (32)
```

The companion verifies exactly

```text
deg_(t,v) D=(10,5),          terms(D)=60,
deg_(t,v) T=(23,11),         terms(T)=1220,
C0(0,v)=-672 G3(v),
D(0,v)=-16071601392 L5(v),
gcd(T(0,v),G3 L5)=1.                                 (33)
```

Thus order `23` is literally the first coupling between the three new and
five old norm packets, and it is a unit on every fixed branch.  Formulas
(31)--(32) are the correct Hensel/subresultant instrument for (27); they are
not needed for the divisor-budget proof and are not themselves a rationality
obstruction.

## 8. Reproduction and exact boundary

Run

```bash
python3 04-computation/jc2_split_prime23_component_divisor_budget_20260728.py
python3 -O 04-computation/jc2_split_prime23_component_divisor_budget_20260728.py

python3 04-computation/jc2_split_prime23_bscale_extension_20260728.py
python3 -O 04-computation/jc2_split_prime23_bscale_extension_20260728.py
```

Both modes of both companions byte-match their declared outputs.  The first
script proves the universal norm identity over a polynomial ring, specializes
every exact `b=1` flux coefficient, verifies (31)--(33), exhausts all `23`
nonzero naive component degrees to the sole budget survivor (19), and checks
the exponent arithmetic in (6).  The extension companion derives (1a)
directly from the audited integer fluxes, verifies cancellation of arbitrary
nonzero `rho`, rechecks the fixed fibre and weighted corner, and includes
`b=0` in the same parameter-free budget.
An independent hostile audit separately reconstructed the weighted corner,
all fixed local units and completed branches, the nonvertical component
argument, the opposing divisor budgets, geometric reducedness, perfect-power
normal form, and both ramification invoices.  Normal and optimized replay
both emit the declared `19`-line transcript, and the companion contains no
Python `assert` nodes.

The theorem does **not**:

1. exclude solutions of the perfect-power system (27);
2. exclude `89` extra delta units on a special integral member;
3. restore or kill the eleven odd Faber seeds;
4. spend the third flux, the split `y=0` boundary, or every `A=0` physical
   trajectory consequence; or
5. prove the full split branch, `JC(2)`, or `DC(2)`.

The exact next target is (27), preferably through (32) at order `23`, with
the third flux retained before any coefficient quotient.

QED.
