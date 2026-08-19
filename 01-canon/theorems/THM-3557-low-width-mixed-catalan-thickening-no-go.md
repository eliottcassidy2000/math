---
id: THM-3557
title: "Low-width mixed Catalan thickening no-go and degree-six frontier"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT.  For polynomial mixed
  thickenings of the self-intersection (v^2,v^3-v), transverse widths one
  and two admit no constant-Jacobian solution in any coefficient degree.
  At width three, the affine coefficient variety is exactly empty when all
  coefficient degrees are at most 3, 4, or 5.  Thus the first unclosed
  width-three cell has coefficient degree six; width at least four is also
  open.  Every map in the ansatz retains the selected collision, so any
  solution would be a planar counterexample.  No solution is claimed.
source: kps-s187/kps-s188
depends_on:
  - THM-3545-catalan-self-intersection-keller-thickening-boundary
related:
  - THM-3548-planar-keller-conductance-shadow-gates
  - THM-3549-torus-quotient-correction-no-go
  - THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor
  - THM-3555-catalan-thickening-universal-cubic-root-cover
companion: 04-computation/jc2_catalan_mixed_thickening_recurrence_kps_s187.py
output: 05-knowledge/results/jc2_catalan_mixed_thickening_recurrence_kps_s187.out
proof_note: 07-reflections/jc2-catalan-mixed-thickening-recurrence-kps-S187.md
script_sha256: 0444ad61a0bb2cd165243db1a97f0cb0b299eb19263378c97d1dee9ff39a7e1e
output_sha256: cd23937963962bc4f43b83fc2f3ab6b477970c76f3151b7247c0025145380832
hash_basis: LF-normalized bytes
---

# THM-3557 -- low-width mixed Catalan thickening no-go

**PROVED + VERIFIED-EXACT + FINITE-EXACT.**  Mixed dependence on the boundary
parameter does not repair THM-3545 at shallow transverse width.  Width two
fails by a cap-free Wronskian/degree argument; width three first acquires a
genuinely new branch at coefficient degree six.

The field has characteristic zero.

## 1. Exact mixed recurrence

Put

```text
p(v)=v^2,                         q(v)=v^3-v,
P=p+sum_(j=1)^N a_j(v)w^j,
Q=q+sum_(j=1)^N b_j(v)w^j.                           (1)
```

Set `a_0=p`, `b_0=q`.  Collecting the coefficient of `w^k` in
`Jac_(v,w)(P,Q)` gives

```text
E_k=sum_(i+j=k+1) [j a_i' b_j-i a_i b_j'].           (2)
```

Therefore

```text
Jac(P,Q)=1
iff E_0=1 and E_k=0 for 1<=k<=2N-1.                  (3)
```

Every map in `(1)` restricts on `w=0` to

```text
v -> (v^2,v^3-v),                                     (4)
```

so `(v,w)=(1,0),(-1,0)` both map to `(1,0)`.  A solution of `(3)` would be a
noninjective planar Keller map.

## 2. Width one is impossible in all degrees

For `N=1`, write `(a_1,b_1)=(a,b)`.  The equations are

```text
p'b-aq'=1,                       a'b-ab'=0.            (5)
```

Neither `a` nor `b` can vanish.  The second equation gives `b=lambda a` in
`k(v)`, hence in `k[v]`, for a scalar `lambda`.  The first becomes

```text
a(1+2lambda v-3v^2)=1,                                (6)
```

which is impossible in `k[v]`.

## 3. Width two is impossible in all degrees

Write

```text
P=p+aw+cw^2,                    Q=q+bw+dw^2.           (7)
```

The four equations from `(2)` are

```text
p'b-aq'=1,
2(p'd-cq')+a'b-ab'=0,
2a'd-ad'+c'b-2cb'=0,
2(c'd-cd')=0.                                        (8)
```

There are three exhaustive top-row cases.

### 3.1 `c!=0`

The last equation gives `d=lambda c`.  Put

```text
s=b-lambda a,              H=lambda p'-q'=1+2lambda v-3v^2. (9)
```

Then `(8)` becomes

```text
p's+aH=1,
2cH+a's-as'=0,
c's-2cs'=0.                                           (10)
```

If `s=0`, the first equation asks `aH=1`.  Otherwise the last equation gives
`c=mu s^2`.  With `r=deg s`, the first row forces `r>=1` and
`deg a=r-1`.  In the second row, `2cH` has degree `2r+2`, while
`deg(a's-as')<=2r-2`, a contradiction.

### 3.2 `c=0,d!=0`

The third equation gives `d=mu a^2`.  If `r=deg a`, the first row forces
`deg b=r+1`.  The second row is

```text
4mu v a^2+a'b-ab'=0.                                  (11)
```

Its first term has degree `2r+1`, while its Wronskian has degree at most
`2r`; contradiction.

### 3.3 `c=d=0`

This is the width-one system.  Hence no width-two solution exists, without a
degree cap.

## 4. Width three through coefficient degree five

For `N=3`, write

```text
(a_1,a_2,a_3)=(a,c,e),       (b_1,b_2,b_3)=(b,d,f).   (12)
```

The six rows are

```text
E0=p'b-aq'=1,
E1=2(p'd-cq')+a'b-ab'=0,
E2=3(p'f-eq')+2a'd-ad'+c'b-2cb'=0,
E3=3a'f-af'+2(c'd-cd')+e'b-3eb'=0,
E4=3c'f-2cf'+2e'd-3ed'=0,
E5=3(e'f-ef')=0.                                     (13)
```

The bottom row forces `deg b=deg a+1`.  The top row splits into the four
exhaustive cases

```text
(e,f) both nonzero,       e=0/f!=0,
f=0/e!=0,                 e=f=0.                      (14)
```

For coefficient cap `D<=5`, the exact degree ledger and lower-row equations
close every branch:

| top branch | surviving common-power type | terminal obstruction |
|---|---|---|
| `e,f!=0` | `f=lambda e`; then `2 deg e=3 deg s` | only `(deg s,deg e)=(0,0),(2,3)` occur; `E3,E4` have incompatible top degrees |
| `e=0,f!=0,c=0` | one-sided rows | top coefficient `3` in `E3` |
| `e=0,f!=0,c,f` constant | constant common-power row | two forced constants differ by `9/4` |
| `e=0,f!=0,(deg c,deg f)=(2,3)` | square/cube branch | saturated rational ideal is `(1)`; terminal constant `-9/4` |
| `f=0,e!=0,d=0` | finite degree states | `E3` has nonzero leading coefficient `-7` |
| `f=0,e!=0,d,e!=0` | `2 deg e=3 deg d` | no simultaneous `E1,E2` degree survivor |
| `e=f=0` | width two | Section 3 |

The exceptional square/cube branch is transparent enough to record.  Its
top equations force

```text
c=h_1^2(v-rho)^2,             f=h_1^3(v-rho)^3,
a=1+2v(h_0+h_1v).                                     (15)
```

Divisibility after solving `E1` gives

```text
4h_0^2-4h_1^2rho^2-2h_1+3=0.                          (16)
```

The `v^2` coefficient of `E2` forces `h_0+h_1rho=0`, then `(16)` gives
`h_1=3/2`; the constant coefficient of `E2` is `-9/4`, contradiction.
Independently, adjoining the saturation inverse for `h_1` makes the complete
coefficient ideal have Groebner basis `[1]` over `Q`.

The companion implements every zero/nonzero branch and all degree states.
For each `D=3,4,5`, the affine coefficient variety is exactly empty.  This is
not a failure of a heuristic solver and makes no statement about projective
points at infinity, which do not represent affine polynomial coefficients.

## 5. Controls and the first open cell

The exact recurrence has two positive controls:

```text
(P,Q)=(v,w)                    -> (E0,E1)=(1,0),
three-row Catalan truncation   ->
 (1,0,0,-135/16,-405/64,-729/128).                    (17)
```

The second control preserves the genuine Catalan prefix and locates its first
leak at the correct row.  All 86 truth gates remain active under optimized
execution, and ordinary/optimized transcripts equal the stored output.

At `D=6`, new common-power cases appear:

```text
(deg s,deg e)=(4,6),
(deg c,deg f)=(4,6),
(deg d,deg e)=(4,6).                                  (18)
```

They are absent from the cap-five proof and cannot be discarded by repeating
its ledger.  Thus the first open width-three target is exactly coefficient
degree six.  Width four is independently open.

## 6. Cubic-cover reframe and scope

Under THM-3555's quadratic coordinate

```text
r^2=1-3w,                 w=(1-r^2)/3,                 (19)
```

every polynomial mixed thickening in `(1)` becomes even in `r`, and

```text
Jac_(v,r)(P,Q)=Jac_(v,w)(P,Q) dw/dr=-2r/3.             (20)
```

The actual Catalan polynomialization also has Jacobian `-2r/3`, but uses odd
terms in `r`; it is not an even pullback.  Hence the construction problem is
equivalently to find an even polynomial deformation with `(20)` that retains
the selected fiber.  Width is degree in `r^2`, so `(18)` is the first new
even `r^6` cell.

This theorem excludes only ansatz `(1)` at widths `N<=2`, and at `N=3` only
through coefficient degree five.  It does not exclude the degree-six cell,
width four, other boundary curves, or planar counterexamples.  It proves no
projective-closure statement and makes no inference from the interrupted raw
Groebner attempt recorded in the research note.
