# Fixed-head affine paired paths: exact reduction and coherent all-word support

**Status: PROVED LEMMA CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED; not yet a canon theorem.**

This memo closes the finite mathematical gap behind the fixed-head
inverse-root experiment.  It assumes the scalar-row role thresholds already
proved in THM-2604 and the physical opposite-shift convention of THM-2599.
It proves a coherent statement that is stronger than the originally proposed
eight-stop census:

> For every canonical target role and its `13`-divisible blocker, at least
> nine invariant phases support the paired target-danger/blocker-safe event
> on **all thirteen physical inverse roots simultaneously**.  Consequently
> each such phase realizes every prescribed finite root word, with one affine
> shift law.  In particular it realizes any prescribed affine correction
> `q_j=q_0-j alpha` after the corresponding change of physical root.

The result is at the physical gate/cylinder level.  It does not identify the
physical root or shift with a THM-2542 chart coordinate, and it does not
produce a THM-2334 relation current or exclude a scalar row.

## 1. Inheritance pass and exact object

- **Closest proved mechanism:** THM-2599 proves that every physical root has
  some positive opposite-shift slice and that arbitrary root words have
  positive cylinder realizations once one may choose the shift root by root.
- **Canonical hostile:** THM-2602 shows that equal endpoint marginals do not
  determine an ordered transition; its identity and translation kernels have
  the same marginals and different sevenfold monodromy.
- **Corrected near miss:** MISTAKE-260 forbids reading an anchored nonzero
  mode from an unordered inverse-root count.  The invariant phase below is
  therefore retained explicitly rather than reconstructed after
  marginalization.
- **Least-used sidecar:** the lengths of the full danger teeth and safe gaps,
  including the two endpoint safe pieces when the blocker has one period,
  almost completely classify the possible null cells before computation.

Let

```text
d_L(x)=1_(||x||<L/14),                 u_1=1-d_1,
iota_r(z)=(z+r)/13,                    r in F_13,
a=13b,                                 b>=1.
```

For the THM-2599 paired gate with shift `q`, inverse-root substitution gives

```text
P_(r,q)(z)
 =d_L((kz+kr+q)/13) u_1(bz-q/13),      0<z<1.            (1)
```

Here `13` does not divide `k`; the canonical all-root thresholds are

```text
guard:      L=2, k>=10,
ordinary:   L=1, k>=12.                                     (2)
```

Put

```text
v=kr+q in F_13,
A_v={z:d_L((kz+v)/13)=1},
C_q={z:d_1(bz-q/13)=1}.                                  (3)
```

Then the cell `(r,q)` is null exactly when

```text
A_v subset C_q                         up to endpoints.    (4)
```

Call `(v,q)` **bad** when (4) holds.  The phase `v`, not `q` alone, is
the fixed-head coordinate.  The affine transport

```text
tau_Delta(r,q)=(r+Delta,q-kDelta)                         (5)
```

preserves it exactly:

```text
k(r+Delta)+(q-kDelta)=kr+q=v.                             (6)
```

This is the connection whose bad graph is classified below.

## 2. Interval geometry

The components of `A_v` are the intersections with `(0,1)` of

```text
((13m-v)/k-13L/(14k), (13m-v)/k+13L/(14k)),              (7)
```

and the components of `C_q` are the intersections with `(0,1)` of

```text
((m+q/13)/b-1/(14b), (m+q/13)/b+1/(14b)).                (8)
```

Thus the relevant full lengths are

```text
target danger tooth:       13L/(7k),
target safe gap:            13(7-L)/(7k),
blocker danger tooth:       1/(7b),
blocker safe gap:           6/(7b).                       (9)
```

If (4) holds, every connected target-danger component lies in one
blocker-danger component.  Taking complements, every full blocker-safe gap
lies in one target-safe component.  For `L=1`, every target-safe component,
including either component clipped by `0` or `1`, has length at most

```text
78/(7k).                                                   (10)
```

Indeed an endpoint component is only a subinterval of one gap in the
periodic target pattern.

## 3. Guard reduction

Suppose first that `L=2`.

### 3.1 Speeds `k>=13`

The variable `(kz+v)/13` traverses an interval of length `k/13`.  By
periodicity, its integral decomposes into `floor(k/13)` full-period
integrals, each with danger mass `2/7`, plus one nonnegative residual-arc
integral.  Hence

```text
mu(A_v)>=(13/k)(2/7) floor(k/13).                         (11)
```

Write `k=13n+s`, where `n>=1` and `1<=s<=12`.  Then

```text
(13/k)(2/7)n>1/7
iff 26n>13n+s
iff 13n>s.                                                (12)
```

On the other hand multiplication by the integer `b` preserves Haar measure,
so

```text
mu(C_q)=1/7.                                               (13)
```

Equations (11)--(13) rule out `A_v subset C_q` for every `b,v,q`.

### 3.2 Speeds `k=10,11,12`

For these three speeds, exact endpoint clipping over the thirteen values of
`v` gives

```text
min_v max{|J|:J a component of A_v}
 =3/35, 6/77, 13/84                 for k=10,11,12.       (14)
```

Every number in (14) is greater than `1/14`.  If `b>=2`, every component of
`C_q` has length at most `1/(7b)<=1/14`, so containment is impossible.
Only the `3*13*13=507` cases with `b=1` remain.  The two-route exact census
in Section 6 finds no bad pair.

Therefore the guard bad graph is empty for every canonical parameter.

## 4. Ordinary reduction for `b>=2`

Now let `L=1`.

### 4.1 Small all-root speeds

The only canonical speeds below `15` are `12` and `14`; `13` is excluded.
Their exact component floors are

```text
min_v max{|J|:J a component of A_v}
 =13/168, 13/98                       for k=12,14.         (15)
```

Both exceed `1/14`, so the same component-length argument rules out every
`b>=2`.

### 4.2 Speeds `k>=15`

Since

```text
k/13>=15/13>1+1/7,                                        (16)
```

the target parameter interval contains a complete ordinary danger tooth.
Thus `A_v` contains a component of length `13/(7k)`.  If (4) holds, it fits
inside a blocker component, giving

```text
13/(7k)<=1/(7b),                    hence k>=13b.          (17)
```

When `b>=2`, the interval `(0,1)` contains a complete blocker-safe gap of
length `6/(7b)`.  By complement containment and (10),

```text
6/(7b)<=78/(7k),                    hence k<=13b.          (18)
```

Together (17)--(18) force `k=13b`, contradicting `13` not dividing `k`.
So no ordinary bad cell has `b>=2`.

The opposing tooth/gap inequalities are the main infinite reduction.  The
large exact box in Section 6 is only a hostile control for them, not their
proof.

## 5. Ordinary reduction for `b=1`

For `q=0`, the blocker danger wraps around the endpoints and its complement
is one interior interval of length `6/7`.  If (4) holds, that interval lies
in one target-safe component, so (10) gives

```text
6/7<=78/(7k),                       hence k<=13.           (19)
```

Thus `q=0` is never bad for `k>=15`.

For `q!=0`, the blocker danger is one interior component.  Its complement
has two endpoint pieces whose total length is `6/7`.  Each piece lies in
the corresponding endpoint target-safe component, and each target component
has length at most (10).  Therefore

```text
6/7<=2(78/(7k)),                    hence k<=26.           (20)
```

The endpoint `k=26` is again excluded by `13|k`.  Combining (15), (19), and
(20), the complete ordinary finite universe is

```text
b=1,       k in {12,14,15,16,...,25}.                     (21)
```

It has `13*13*13=2,197` phase/shift cells.  With the `507` guard cells, the
entire unresolved universe has exactly `2,704` cells.

## 6. Exact bad graph

The exact census gives the following complete table.  An arrow `v->q` means
that `A_v subset C_q` up to endpoints.

| `k` | bad pairs `v->q` |
|---:|:---|
| 15 | `5->7, 6->6` |
| 16 | `2->9, 3->8, 7->5, 8->4` |
| 17 | `1->9, 4->7, 5->6, 8->4` |
| 18 | `2->8, 3->7, 5->6, 6->5` |
| 19 | `1->8, 3->7, 4->6, 6->5` |
| 20 | `1->8, 2->7, 4->6, 5->5` |
| 21 | `2->7, 3->6` |
| 22 | `1->7, 3->6` |
| 23 | `1->7, 2->6` |

There are no other pairs.  In particular, the graph is empty for the guard;
for ordinary `k=12,14,24,25`; for every `b>=2`; and for all larger canonical
speeds.  The table has `28` edges on nine speed rows, and every speed row has

```text
|Dom(G_k)|<=4.                                             (22)
```

The companion checks the reduced universe by two structurally independent
exact routes (both necessarily start from the same rational wall equations):

1. build the rational interval unions (7)--(8) and test containment up to
   their finitely many endpoints;
2. build all rational walls directly, evaluate both Boolean gates on every
   intervening midpoint, and sum the exact `Fraction` lengths on which
   `A_v=1` and `C_q=0`.

The routes agree on all `2,704` reduced cells.  As a non-proof hostile
control, direct interval containment was also checked on `1,774,500` cells
with `k<=200` and `b<=30`; it found exactly the same 28 pairs.

## 7. The all-word fixed-head theorem

Fix a canonical `(L,k,b)` and choose

```text
v in F_13 minus Dom(G_k).                                 (23)
```

There are at least nine choices by (22).  Define the coherent shift section

```text
q_v(r)=v-kr.                                               (24)
```

For every root `r`, the pair `(v,q_v(r))` is not bad.  Hence

```text
mu{z in (0,1):P_(r,q_v(r))(z)=1}>0                       (25)
```

simultaneously for all thirteen roots.  This is the key strengthening:
THM-2599 no longer needs thirteen unrelated choices `s_r`; one horizontal
phase (23) gives all of them through (24).

Choose one positive open chamber in (25) for each `r`, map it by `iota_r`,
and choose a common-depth base-13 cylinder `D_r` inside the resulting subset
of `I_r`.  Exactly as in THM-2599, if `T(x)=13x mod 1`, then for every finite
root word `(r_0,...,r_m)` and every block separation `N` at least the common
cylinder depth,

```text
mu(intersection_(j=0)^m T^(-jN)D_(r_j))
 =13^(-(m+1)ell)>0.                                       (26)
```

Every stop in (26) carries the shift `q_v(r_j)`.  Thus one good `v` embeds
the full one-sided thirteen-shift together with the affine label cocycle
(24).  The conclusion holds for arbitrary finite path length and arbitrary
repeated root word, not just for an affine eight-stop path.

For a prescribed initial root `r_0`, the map

```text
q_0 -> v=kr_0+q_0                                         (27)
```

is a bijection.  Therefore at least nine initial shifts `q_0` work for every
prescribed finite root path.

## 8. Sign and twisted-return audit

Write an affine physical root path as

```text
r_j=r_0-j delta.                                          (28)
```

Equations (24) and (27) give

```text
q_j=v-kr_j=q_0+kj delta.                                  (29)
```

This sign is forced by the invariant `v=kr+q`.  To obtain the correction
required in THM-2602's algebraic twisted-return test,

```text
q_j=q_0-j alpha,                                          (30)
```

choose

```text
delta=-k^(-1)alpha.                                       (31)
```

Equivalently, each physical root step is
`Delta=k^(-1)alpha` in (5), and (5) decreases `q` by `alpha`.  After seven
steps, (30) gives `q_7=q_0-7alpha`, the exact endpoint in THM-2602 equation
(14).

As a frozen diagnostic, the eight-stop census over the nine exceptional
speeds, all `13` starts, and all `12` nonzero root steps has good-start
histogram

```text
{9:86, 10:268, 11:584, 12:372, 13:94}.                    (32)
```

The lower bound nine is sharp in this diagnostic.  One control is

```text
(k,r_0,delta)=(16,0,8),
q_0 in {0,1,4,5,6,9,10,11,12}.                            (33)
```

The proof of the all-word result is (22)--(27), not the eight-stop census.

## 9. Connection audit and stopping boundary

The precise connection exposed here is:

```text
source object:       physical inverse-root/shift pairs (r,q),
target object:       horizontal phase fibres kr+q=v,
map:                 tau_Delta(r,q)=(r+Delta,q-kDelta),
preserved predicate: paired target-danger/blocker-safe positivity,
preserved quantity:  v,
forgotten data:      chamber location, owner/repair semantics, endpoint X,
needed sidecar:      lawful THM-2542 chart-clock identification and gluing,
cheapest next test:  realize tau_(k^-1 alpha) on each adjacent retained
                     chart edge and verify one positive sevenfold fibre
                     product before marginalization.                      (34)
```

Thus the memo supplies an exact physical analogue of THM-2602's desired
translation path.  It does **not** yet prove that the physical root `r` or
the dipole shift `q` is THM-2542's chart variable, nor that the cylinders
retain the selected source packet and its owner/repair fields through seven
adjacent charts.  Those are typing and gluing claims, not further interval
inequalities.

Positivity and ordering also do not imply nonzero `H`-drift.  THM-2367 shows
that an isolated lawful tensor is circulant when the blocker is strictly
`7`-deeper than the graft, and that extra Boolean masks can restore
circulancy even in an escaping valuation profile.  Therefore this affine
connection must not be promoted to a row exclusion without a noncancellation
sidecar.

The result does not recover a THM-2334 relation residue, a left endpoint
mode, a full-`X` current, or LRC(14).  What it removes is the narrower
obstruction that a single opposite-shift phase might fail somewhere along a
long prescribed physical root path.

## 10. Reproduction

Run

```bash
python3 04-computation/lrc14_fixed_head_affine_paired_path_referee.py
python3 -O 04-computation/lrc14_fixed_head_affine_paired_path_referee.py
```

The two executions byte-match

```text
05-knowledge/results/lrc14_fixed_head_affine_paired_path_referee.out
```

The checker is dependency-free and uses `Fraction` throughout.  Hashes on
the recorded bytes are

```text
script: 29be125d8174804cf688e7fd7d6be8ebb70e8ddcf29ea4faf2e89a1ab5df2a02
output: 3ae8a95514d021a719c0f52c9c7728ef4604061974bd6a094dff4ce7449a7cf1
```

An independent hostile audit replayed both execution modes and hashes,
rederived the full-tooth/full-gap inequalities, checked both `b=1` endpoint
geometries, audited the `tau_Delta` sign, and verified the common-depth
cylinder specification.  It found no mathematical or computational defect.
