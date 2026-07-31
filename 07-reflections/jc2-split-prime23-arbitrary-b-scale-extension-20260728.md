# The prime-23 `3:5` divisor packet does not require `B_0!=0`

**Status:** PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED as the
arbitrary-`b` addendum to
`THM-2713-split-prime23-component-divisor-budget-and-perfect-power-normal-form.md`.
The result extends the response-curve geometry only.  It does not by itself
promote the reserved unified physical closure, treat `y=0`, restore an odd
Faber seed, or prove `JC(2)`.

Companion:

```text
04-computation/jc2_split_prime23_bscale_extension_20260728.py
05-knowledge/results/jc2_split_prime23_bscale_extension_20260728.out
script SHA-256 9f8d8d6d14eff6d121f5ac45f070d4489d71a72f788ceb5f73db75c4fb63d216
output SHA-256 4ef3d75c3bf6e82f05cefc7a133b1efe4d5b79c87f0e165e1c3bcbc9c2f2460f
```

## 1. The normalization had hidden a harmless parameter

THM-2704/2713 chose `rho in C*` with `rho^2=B_0`, then normalized the
weight-two Faber coefficient to one.  That is convenient on `B_0!=0`, but
the geometric proof never used the value one away from positive powers of
the signed scale.

Choose instead any `rho in C*` and retain

```text
b=B_0/rho^2,
c=C_0/rho^3,        d=D_0/rho^4,
e=E_0/rho^5,        w=W_0/rho^6,
t=rho/y.                                                    (1)
```

This is valid whenever `y` is not identically zero, including `B_0=0`.
The old normalization is the section `b=1`; the omitted boundary is `b=0`.

## 2. Exact deformation from `b=1`

Rescaling the audited integer fluxes with (1) gives

```text
f1_b-f1_1
 =-616 t^2(b-1)(4840v-1331zeta-40),

f2_b-f2_1
 =49280 t^2(b-1)(29282v^2-1452v+1331zeta+2).          (2)
```

The arbitrary nonzero `rho` cancels.  The chosen-sheet relation remains

```text
zeta f1_b^4
 =(7496192^4 lambda^4/rho^23)t^23.                    (3)
```

Thus `eta!=0` whenever `lambda!=0`, with no condition on `b`.

Every new term in (2) contains `t^2`.  Consequently the complete fixed fibre
at `t=0` is literally unchanged:

```text
new: zeta=0, G3(v)=0,             three branches;
old: f1_b=0, L5(v)=0,             five branches;
fixed eliminant=G3 L5^4.                                  (4)
```

The imported audited scout rechecks all squarefree, coprime, Jacobian, and
unit controls before the extension assertions run.

## 3. Why the global component proof transfers

The geometry needed by THM-2713 consists of four parameter-insensitive
gates.

1. At the two quotient-singular coordinate points of `P(1,1,2,3)`, `t=0`;
   the values of `F2_b` remain `-1190488992` and `15944049`.
2. At the corner `h=t=0`, the top surviving forms remain

   ```text
   F1_b=-1449459 v zeta,
   F2_b=15944049 zeta^2-1190488992 v^3,
   ```

   so the degree-23 equation and `F2_b` have no common point there.
3. The old and new completed local equations remain

   ```text
   F1_b^4=unit*t^23,             zeta=unit*t^23,
   ```

   with order words `(4,23,0)` and `(1,0,23)`.
4. The line-bundle degrees of `F1_b` and `zeta` remain five and three.

The regular-sequence, Cohen--Macaulay, basepoint-free pencil, and reducedness
arguments therefore transfer verbatim.  On every geometric component the
same counts `(r,s)` satisfy

```text
d=4r+s,
23r<=5d,                    23s<=3d.                  (5)
```

Hence `3r=5s`, and the only nonempty bounded solution remains

```text
(r,s,d)=(5,3,23).                                    (6)
```

Every member is geometrically integral, and equality in (5) still exhausts
the two zero divisors.

## 4. The chosen sheet retains five poles

The physical split coordinate reconstructs as

```text
q_phys=-7496192 lambda t^5/(rho^5 f1_b).              (7)
```

After removing its nonzero scalar, THM-2713's normalized `q` therefore has

```text
div(q)=5N-3O                                         (8)
```

for every `b`, including zero.  The arithmetic genus, fixed cusp delta,
genus invoice, and rational perfect-power normal form use only the same
weights and divisors, so they also transfer unchanged.

This supplies the missing response-space input for extending the five-pole
physical no-go across `B_0=0`; that synthesis remains separately reserved
until its dependencies are integrated.

## 5. Sharp warning about infinity

The conclusion is **not** that the whole `h=0` face is independent of `b`.
After weighted homogenization, some `t^2(b-1)` terms survive when `h=0` and
`t!=0`.  What is unchanged is exactly what the proof consumes:

```text
t=0 fixed fibre;
h=t=0 corner;
ambient coordinate-point values;
old/new local units;
section weights and component budgets.                       (9)
```

The basepoint-free map `[h:t]` excludes vertical and `h=0` components by
ampleness; no generic or parameter-independent description of the remaining
infinity fibre is needed.

## 6. Audit and boundary

The independent audit replayed the companion in normal and optimized modes,
verified (2)--(3), and rechecked every gate in (9).  It specifically tested
that no hidden use of `B_0!=0` remains: `rho` is now a freely chosen nonzero
scale, while `b` records the actual weight-two coefficient.

The addendum covers

```text
b,c,d,e,w arbitrary,          eta!=0,
prime-23 even-Faber response curve.                           (10)
```

It does not cover `y=0` through the coordinate `t=rho/y`, `lambda=0`, the
odd bank, or the full split branch.  The separate exact `y=0` resultant is
the complementary boundary mechanism.
