---
id: THM-3621
title: "Russell-cylinder polynomial even-fold all-order vertical-cokernel rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED / ALTERNATIVE
  REDUNDANT CERTIFICATION OF THM-3619.
  Every polynomial even fold in the THM-3612 family is excluded: no pair of
  regular target functions pulls back to a nonzero constant source Jacobian.
  The all-order vertical-cokernel invoice forces the unique formal side germ
  Q_infinity=-3/4-9/(4x^2), which cannot be polynomial.  This closes only the
  stated even-fold family, not planar JC.
source: kps-s189 / THM-3619 all-order shifted-evaluation continuation, 2026-08-21
audit: >
  PASS -- an independent hostile reconstruction verified the shift
  filtration, branch signs, exact moving target and bivector identity, all
  target-displacement orders, normalized coefficient -16, recurrence,
  polynomial contradiction, and finite THM-3619 reconciliation.  Normal,
  optimized, and stored 752-gate transcripts agree after LF normalization.
  The theorem is a distinct proof, not a stronger closure result.
superseded_by: THM-3619-russell-cylinder-even-fold-higher-jet-staircase
depends_on:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
related:
  - THM-3619-russell-cylinder-even-fold-higher-jet-staircase
script: 04-computation/russell_even_fold_all_order_vertical_cokernel_thm3621.py
output: 05-knowledge/results/russell_even_fold_all_order_vertical_cokernel_thm3621.out
script_sha256: 456841200aeeb7985c8925dfa046edf6d12944a70ba54f62039624401c377656
output_sha256: 1500a76d4e7395fd2c2a7da1e71c0a70ed18d34eabbe13f9ddabf1c81311071b
semantic_sha256: b705c5c924d19e6042496420c9998291c77ef9f7c9d41bcaeea7bf0261f823c2
hash_basis: raw LF bytes
---

# THM-3621 -- Russell-cylinder polynomial even-fold all-order vertical-cokernel rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED / ALTERNATIVE
REDUNDANT CERTIFICATION OF THM-3619.**

All formal rings and differential forms below are over `C`.  This theorem
was developed concurrently with the promotion of THM-3619.  It gives an
independent shifted-evaluation/bivector proof of the same even-fold closure;
THM-3619 is the primary canonical statement.  The proof here is analytic
formal algebra, not an extrapolation from four computed rows.

## 1. Statement

Retain the stabilized Russell map and the polynomial fold of THM-3612,

```text
E_Q(x,t)=(x,Q(x)+t^2,t),
Q even,                 Q(0)=-3/4,          Q(+-1)=-3.   (1)
```

For every pair `F,G` of regular functions on the smooth Russell target,

```text
Jac_(x,t)(F o Psi o E_Q,G o Psi o E_Q) is not in C*.    (2)
```

Thus the doubly tuned locus left open by THM-3612 contains no polynomial
counterexample.  The result does not address non-even folds, folds with a
nonquadratic stable coordinate, implicit planes outside `(1)`, or arbitrary
polynomial maps of `A2`; planar JC remains open.

The exact all-order statement is the following.  Write

```text
q_k=Q^(k)(1),
c_n=2^(n+3)/(3^(n-1)(n-1)!).                            (3)
```

If a hypothetical pair has nonzero constant Jacobian and the side jets agree
with the unique survivor through order `n-1`, then for every `n>=3`,

```text
Delta_(2n-2) = -2^(n+3)/(3^(n-1)(n-1)!)
                * (q_n+(n+1)q_(n-1)).                    (4)
```

Here `Delta_N` is the unique vertical second-difference cokernel

```text
Delta_N(j)=[t^N]j_- -2[t^N]j_0+[t^N]j_+,               (5)
```

with `j_i dxi wedge dt` the three branch pullbacks of `dF wedge dG`.

## 2. Moving the side branches without changing the invoice

Suppose the constant Jacobian has been scaled to `12`, as in THM-3612.  Put

```text
n>=3,               m=n-1,               u=t^2,
h(u)=(1-4u/3)^(-1/2),                    rho(u)=1-h(u).  (6)
```

Use `xi=x+1,x,x-1` at the minus, middle, and plus preimages respectively.
If all source rows below total degree `2m` have been killed, then

```text
j_i(xi,t)=12+O((xi,t)^(2m)).                              (7)
```

The unshifted invoice `(5)` equals the shifted one

```text
Delta_(2m)
 =[u^m] Even_t(j_-(rho,t)-2j_0(0,t)+j_+(-rho,t)).        (8)
```

Indeed, a monomial `xi^a t^b` affected by the shift has `a>=1`.  If it
survived `(7)`, then `a+b>=2m`; after `xi=rho=O(t^2)` its order is

```text
2a+b=(a+b)+a >= 2m+1.                                   (9)
```

It cannot change `[t^(2m)]`.  Terms with `a=0` are literally unchanged.
This filtration argument is the all-order replacement for repeated matrix
elimination.

It also explains THM-3619's sparse quotient rows: their coefficients are

```text
d_a^(r)=-[u^r]rho(u)^a.                                  (10)
```

## 3. The rational germ has an exact moving triple collision

Introduce the comparison germ

```text
Q_infinity(x)=-3/4-9/(4x^2).                             (11)
```

At `x=+-h(u)` one has

```text
Q_infinity(h)+u=-3+4u,
D=1+h^2(-3+4u)=-2.                                      (12)
```

Both moving side points and the middle point `x=0,q=-3/4+u` therefore map,
in the smooth target coordinates `(C,Y,Z=S+3/4)`, to the same formal arc

```text
P(t)=(0,4t,7t^2/4).                                     (13)
```

Let `V_i` be the pushed source bivector at the corresponding point, written
in the basis

```text
(dC wedge dY, dC wedge dZ, dY wedge dZ).
```

Direct differentiation gives the identity in `C[[t]]^3`

```text
V_-^infinity+V_+^infinity
 =(24,21t,6t(11t^2-6))
 =2V_0.                                                  (14)
```

Since all three bivectors in `(14)` are based at the same target arc `(13)`,
pairing with **any** formal target two-form gives

```text
j_-^infinity(rho,t)-2j_0(0,t)+j_+^infinity(-rho,t)=0.    (15)
```

Closedness and decomposability are not used; this strictly enlarges the
actual forms `dF wedge dG` and makes the target-pair quantifier safe.

## 4. First departure and the universal invoice

Assume inductively that `Q` agrees with `(11)` through derivative `n-1` at
`x=1`, and put

```text
epsilon=Q-Q_infinity,
delta_n=q_n-q_n^infinity=q_n+(n+1)q_(n-1).              (16)
```

Because `h-1=(2/3)u+O(u^2)`,

```text
epsilon(h)=O(u^n),
epsilon'(h)=delta_n/(n-1)! * (2u/3)^(n-1)+O(u^n).        (17)
```

Evenness gives the same value and the opposite derivative at `-h`.  The
target-point displacement from `epsilon(h)` begins at `u^n`, too late for
`[u^(n-1)]`.  Only the derivative displacement remains.

Let `Phi=(C,Y,Z)` be the ambient Russell compiler and

```text
G_+-G_-=(Phi_q wedge Phi_w)|_(x=h,q=-3+4u,w=t)
        -(Phi_q wedge Phi_w)|_(x=-h,q=-3+4u,w=t).        (18)
```

At `t=0`, exact differentiation gives

```text
(G_+-G_-)(0)=(-16,0,0).                                  (19)
```

The middle constant row is `12`; because `V_0(0)=(12,0,0)`, the constant
coefficient of the target form on `dC wedge dY` is therefore `1`.  Every
positive-order term of the target form along `(13)` or of `(18)` lands above
`u^(n-1)` after multiplication by `(17)`.  Equations `(8),(15)--(19)` yield

```text
Delta_(2n-2)
 =-16*(2/3)^(n-1)/(n-1)! * delta_n,
```

which is exactly `(4)`.

The branch signs in `(18)` are forced: `epsilon'(-h)=-epsilon'(h)`, so the
sum of the two side perturbations is `epsilon'(h)(G_+-G_-)`, not a sum of the
`G` terms.

## 5. Induction and polynomial contradiction

THM-3612 already forces

```text
q_1=9/2,                         q_2=-27/2,               (20)
```

which are the first two derivatives of `(11)`.  A genuinely constant source
Jacobian makes every invoice `(4)` vanish.  Induction therefore forces

```text
q_n=-(n+1)q_(n-1),                  n>=3,                 (21)
q_n=(-1)^(n-1) 9(n+1)!/4.                                  (22)
```

Together with `Q(1)=-3`, all Taylor coefficients of `Q` at `1` equal those
of `(11)`.  Hence the polynomial

```text
x^2(Q(x)+3/4)+9/4                                      (23)
```

has zero Taylor germ at `1` and must vanish identically.  But its value at
`x=0` is `9/4`, a contradiction.  This proves `(2)`.

## 6. Boundary and reproduction

For `n=3,4,5,6`, `(4)` is respectively

```text
Delta_4 =-(32/9)(q_3-54),
Delta_6 =-(64/81)(q_4+270),
Delta_8 =-(32/243)(q_5-1620),
Delta_10=-(64/3645)(q_6+11340),                          (24)
```

exactly the provisional THM-3619 staircase.  The rational germ `(11)` is a
valid all-order formal side solution but is singular at the middle source
point; it is a comparison object, not a polynomial fold or a JC example.

The exact companion verifies `(10),(12)--(15),(18)--(19),(24)` and hostile
finite instances of the general coefficient formula with optimization-safe
gates.  The quantifier over all `n` rests on the formal proof in Sections
2--5, not on a bounded computation.

The companion's displayed ancestry `PIN` strings are historical provenance:
its executable gates check their syntax, not the current repository blobs.
The source/output hashes in this theorem and the 752 mathematical gates are
the verified package; no dependency claim rests on those provenance strings.

Reproduce with

```bash
python3 04-computation/russell_even_fold_all_order_vertical_cokernel_thm3621.py
python3 -O 04-computation/russell_even_fold_all_order_vertical_cokernel_thm3621.py
```
