---
id: THM-3627
title: "Russell-cylinder non-even hostile degree-six closure"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.
  The explicit non-even polynomial Q_h from THM-3624 survives the enlarged
  arbitrary-target-two-form collision system through total source degree five,
  but fails exactly at degree six.  A sparse exact left-cokernel certificate
  kills all 252 degree-six target columns and has nonzero value on the common
  constant column.  This closes Q_h only; the general non-even-fold family
  remains OPEN.
source: root / audit_thm3624 higher-order hostile continuation, 2026-08-21
audit: PENDING -- provisional theorem and exact companion require hostile audit
depends_on:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
related:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
  - THM-3619-russell-cylinder-even-fold-higher-jet-staircase
script: 04-computation/jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.py
output: 05-knowledge/results/jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.out
script_sha256: 0f79dda958cab53d59965427032610142d024e23028343b9797a26af94569020
output_sha256: c5e8785c03da19b19c22b930c88e908117a0964e82e683b3ace9c5c17096244f
hash_basis: raw LF bytes
---

# THM-3627 -- Russell-cylinder non-even hostile degree-six closure

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.**
This result continues the exact hostile polynomial from THM-3624 by two more
source degrees.  It closes that polynomial at degree six in a system strictly
larger than the decomposable system `dF wedge dG`.  It does not close the
general non-even family.

All rings, germs, derivatives, and differential forms are over `C`.

## 0. Setup and scope

Use the THM-3624 compiler coordinates

```text
D=1+x^2q,
c=xD(D+2),
b=(D-1)(D+2)^2,
e=q(D+3),                         c^2e=b(b+4),              (1)
```

and the quadratic stable-coordinate fold

```text
E_Q(x,t)=(x,Q(x)+t^2,t).                                  (2)
```

The hostile polynomial is

```text
Q_h(x)=1/5408 (
  44069x^11+112059x^10-154749x^9-406377x^8
 +188107x^7+513081x^6-82835x^5-230931x^4
 +5408x-4056).                                             (3)
```

It is non-even and has the THM-3624 collision jets

```text
x=-1:  (Q_h,Q_h',Q_h'',Q_h''')=(-3,-9/2,0,0),
x= 0:  (Q_h,Q_h',Q_h'',Q_h''')=(-3/4,1,0,0),
x= 1:  (Q_h,Q_h',Q_h'',Q_h''')
        =(-3,9/2,-243/13,10449/169).                       (4)
```

Thus the three source germs at `(-1,0),(0,0),(1,0)` map to the same smooth
target point.  Use its regular local coordinates

```text
(c,epsilon,w)=(c,e+3,w).                                   (5)
```

Write `Phi_Q` for the composite map from the source fold `(2)` through the
compiler `(1)` into these local target coordinates; explicitly,

```text
Phi_Q(x,t)=(c(x,Q(x)+t^2), e(x,Q(x)+t^2)+3, t).           (6)
```

The word *closure* in the title means closure of this finite hostile
candidate, not the differential condition `d Omega=0`.  The obstruction below
allows **every** target two-form, whether closed or decomposable or neither.

## 1. The enlarged arbitrary-two-form system

Let

```text
Omega=P(c,epsilon,w) dc wedge d epsilon
     +K(c,epsilon,w) dc wedge dw
     +R(c,epsilon,w) d epsilon wedge dw.                   (7)
```

For `N>=0`, put

```text
T_N=(C[c,epsilon,w]_(degree<=N))^3,
S_N=direct_sum_(i=-,0,+) C[xi,t]_(degree<=N).              (8)
```

Pull `(7)` back along the three branch germs

```text
x=x_i+xi,                 q=Q_h(x_i+xi)+t^2,              (9)
```

take the coefficient on `d xi wedge dt`, and truncate at total source degree
`N`.  This defines the exact linear map

```text
P_N:T_N -> S_N.                                           (10)
```

Let `tau_N` be `12` in each of the three constant rows and zero in every
positive-degree row.  The normalization `12` is harmless: any common nonzero
constant can be rescaled to it.

Every coordinate in `(5)` vanishes to source order at least one on every
branch.  Hence target coefficient terms of degree greater than `N` cannot
alter the pullback through source degree `N`.  The finite matrix `(10)` is
therefore exhaustive at that order.

## 2. Exact rank boundary

Exact rational arithmetic gives

| `N` | matrix shape | `rank P_N` | `rank[P_N|tau_N]` |
|---:|---:|---:|---:|
| 0 | `3 x 3` | 2 | 2 |
| 1 | `9 x 12` | 7 | 7 |
| 2 | `18 x 30` | 15 | 15 |
| 3 | `30 x 60` | 26 | 26 |
| 4 | `45 x 105` | 40 | 40 |
| 5 | `63 x 168` | 57 | 57 |
| 6 | `84 x 252` | 77 | 78 |

Thus the THM-3624 arbitrary-form survivor actually extends one step farther,
through total source degree five.  It fails for the first time at total source
degree six.

The equality at `N=5` asserts only an arbitrary target-two-form jet.  It is
not a closedness or decomposability claim.  The strict rank jump at `N=6`, by
contrast, obstructs the larger arbitrary-form system and therefore also every
closed or decomposable subsystem.

## 3. Sparse exact left-cokernel certificate

Here is a direct certificate independent of the rank comparison.  Order the
three branches as `(-,0,+)`.  Within each branch order source monomials first
by total degree and then as

```text
t^d, xi t^(d-1), ..., xi^d.                               (11)
```

For a pulled coefficient germ `j_i(xi,t)`, define

```text
L(j)=sum lambda_(i,a,b) [xi^a t^b] j_i.                   (12)
```

All omitted weights are zero.  The 23 nonzero weights are

| branch | source monomial | weight `lambda` |
|:---:|:---:|---:|
| `-` | `1` | `242002972/2313441` |
| `-` | `xi` | `-319114/59319` |
| `-` | `t^2` | `6189/2197` |
| `-` | `xi^2` | `1504/1521` |
| `-` | `xi t^2` | `-412/507` |
| `-` | `xi^3` | `-40/351` |
| `-` | `t^4` | `36/169` |
| `-` | `xi^2 t^2` | `20/117` |
| `-` | `xi t^4` | `-10/39` |
| `-` | `t^6` | `5/13` |
| `0` | `1` | `141280234/257049` |
| `0` | `xi` | `-16512/2197` |
| `0` | `t^2` | `-6189/2197` |
| `0` | `xi t^2` | `-216/169` |
| `0` | `t^4` | `-36/169` |
| `0` | `t^6` | `-18/13` |
| `+` | `xi` | `13016/4563` |
| `+` | `xi^2` | `128/117` |
| `+` | `xi t^2` | `32/39` |
| `+` | `xi^3` | `8/27` |
| `+` | `xi^2 t^2` | `4/9` |
| `+` | `xi t^4` | `2/3` |
| `+` | `t^6` | `1` |

Direct exact multiplication gives

```text
L(P_6 v)=0                         for every v in T_6,      (13)

L(tau_6)=465700024/59319 !=0.                              (14)
```

Thus `(12)` kills all 252 arbitrary target-two-form columns but not the common
constant column.  Deleting just the final `+`-branch `t^6` weight makes `(13)`
false; this is an active certificate rather than a vacuous zero-row check.

## 4. Consequence for the hostile polynomial

There is no formal target two-form `Omega` whose pullback along the three
`Q_h` branch germs is one common nonzero constant through all source orders:
such a form would already contradict `(13)--(14)` at order six.

In particular, there are no regular target functions `F,G` for which

```text
Phi_(Q_h)^*(dF wedge dG)=J dx wedge dt,             J in C*. (15)
```

Indeed `dF wedge dG` is a special case of `(7)`, and its local coefficient
series truncates to an element of `T_6`.  Via the Russell-cylinder
isomorphism, the same conclusion holds in the stabilized target coordinates
`(B,C,Y,S)` used in THM-3624.

This is a no-pair result for the **single explicit polynomial** `Q_h`.  It is
not a no-Keller theorem for every non-even fold.  Different side and higher
jets may alter the degree-six quotient, and the general non-even-fold problem
remains **OPEN**.

## 5. Exact reproduction and hostile controls

The deterministic companion verifies

- all derivatives of `Q_h` through order seven at the three collision points;
- every exact rational rank in the table, including the positive `N=5`
  survivor and the hostile `N=6` jump;
- all 252 annihilation identities in `(13)` and the nonzero debt `(14)`;
- a one-weight mutation that destroys annihilation; and
- an AST gate excluding inactive Python `assert` statements.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.py
python3 -O 04-computation/jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.out`.
