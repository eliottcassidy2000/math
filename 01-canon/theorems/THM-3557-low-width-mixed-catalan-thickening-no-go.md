---
id: THM-3557
title: "Low-width mixed Catalan thickening no-go and degree-eight frontier"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT.  For polynomial mixed
  thickenings of the self-intersection (v^2,v^3-v), transverse widths one
  and two admit no constant-Jacobian solution in any coefficient degree.
  At width three, the affine coefficient variety is exactly empty when all
  coefficient degrees are at most seven.  The degree-six square/cube branch
  is forced to C=t^2,F=t^3 and ends in two exact lower-row contradictions;
  the five cells newly admitted at cap seven each have a unique nonzero top
  coefficient.  Thus the first unclosed width-three cell has coefficient
  degree eight; width at least four is also open.  Every map in the ansatz
  retains the selected collision, so any solution would be a planar
  counterexample.  Separately, the CITED 2022 sub-125 degree classification
  forces D+N>=108, so the internal degree-eight frontier remains globally
  excluded through D=104 at width three.  No solution is claimed.
source: kps-s187/kps-s188 + mixed_catalan_d7 / incoming-frontier extension, 2026-08-23
audit: >
  DEGREE-SEVEN EXTENSION HOSTILE AUDIT PASS (root + mixed_catalan_d7,
  2026-08-23).  The extension replays the complete 64-gate degree-six
  companion and its 86 inherited THM-3557 gates, enumerates the cap-six and
  cap-seven ledgers independently, isolates exactly five new cells, and
  extracts their unique first failed coefficients.  The only two terminal
  cells are literally the inherited exact-empty degree-six cells.  The
  38-gate companion has no Python assert; normal and optimized LF-normalized
  streams match the frozen transcript exactly.
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
degree_six_companion: 04-computation/jc2_catalan_mixed_thickening_degree6_exact.py
degree_six_output: 05-knowledge/results/jc2_catalan_mixed_thickening_degree6_exact.out
degree_six_proof_note: 07-reflections/jc2-catalan-mixed-thickening-degree-six-exact-closure-20260818.md
script_sha256: 0444ad61a0bb2cd165243db1a97f0cb0b299eb19263378c97d1dee9ff39a7e1e
output_sha256: cd23937963962bc4f43b83fc2f3ab6b477970c76f3151b7247c0025145380832
degree_six_script_sha256: a3f94ab9eb4157bec7effc15ac1f1a8ab0842814c3a4f4ae4f6705a97fe4346f
degree_six_output_sha256: 02d23d80e4c50edc1323105a9d6fcc2cf5c77d729731a4be354c415251c31616
degree_seven_companion: 04-computation/jc2_catalan_mixed_thickening_degree7_exact.py
degree_seven_output: 05-knowledge/results/jc2_catalan_mixed_thickening_degree7_exact.out
degree_seven_script_sha256: 6201e6b1cd5e199d5e65fff6118cdba586e638342f8a134a10c8d66fbb7de959
degree_seven_output_sha256: 2765d1d693129c282e26973b3c6ff5768a96ad2119639f23baf1041aa0d834b5
degree_seven_semantic_sha256: 8cab44bfd360cf0f80094eb76c130613110e9471bcd4b3ba82c8038713963b93
hash_basis: LF-normalized bytes
---

# THM-3557 -- low-width mixed Catalan thickening no-go

**PROVED + VERIFIED-EXACT + FINITE-EXACT.**  Mixed dependence on the boundary
parameter does not repair THM-3545 at shallow transverse width.  Width two
fails by a cap-free Wronskian/degree argument; width three is now empty
through coefficient degree seven and first remains open at degree eight.

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

## 4. Width three through coefficient degree six

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

## 5. Degree-six closure and the first open cell

The exact recurrence has two positive controls:

```text
(P,Q)=(v,w)                    -> (E0,E1)=(1,0),
three-row Catalan truncation   ->
 (1,0,0,-135/16,-405/64,-729/128).                    (17)
```

The second control preserves the genuine Catalan prefix and locates its first
leak at the correct row.  All 86 truth gates remain active under optimized
execution, and ordinary/optimized transcripts equal the stored output.

At `D=6`, the new common-power cases are

```text
(deg s,deg e)=(4,6),
(deg c,deg f)=(4,6),
(deg d,deg e)=(4,6).                                  (18)
```

They are absent from the cap-five proof and cannot be discarded by repeating
its ledger.  The degree-six companion performs the exhaustive refined ledger.
Every branch closes except

```text
e=0,                deg(a,b,c,d,f)=(3,4,4,5,6),       (19)
c=C h^2,            f=F h^3,             deg h=2.
```

An exact parametrization of `E0,E1`, saturated by every prescribed leading
coefficient, reduces the top coefficients of `E2,E3` to

```text
-3Ct+2F+t^3=0,              2Ft-C^2-Ct^2=0.           (20)
```

Their eliminant is `(C-t^2)^2`, so `C=t^2,F=t^3`.  The next rows force a
common coefficient `U=A` and then split as

```text
(B-R)(3(B-R)+2t)=0.                                   (21)
```

The first branch has terminal coefficient `-6t^2`.  In the second, putting
`X=48Bt+16t^2` forces simultaneously

```text
X=81,                         2X=135,                  (22)
```

which is impossible in characteristic zero.  Independently, both branch
ideals saturated by `t!=0` have Groebner basis `[1]` over `Q`.  Thus the
entire affine width-three coefficient space is empty through `D=6`.  The
next cap is closed below.  Width four is
independently open.

## 5A. Degree-seven closure and the first open cell

Cap seven introduces no new common-power type: the full list is still

```text
(0,0),                  (2,3),                  (4,6). (22a)
```

Relative to cap six, the branch `e,f!=0` admits exactly three new states,
written `(deg(s),deg(e);deg(t),deg(c))`:

```text
(2,3;7,1),             (4,6;2,7),             (4,6;7,3). (22b)
```

In the transformed row

```text
E3=2(c's-cs')+e't-3et',                                (22c)
```

their uniquely highest coefficients are respectively
`-18e_0t_0`, `6c_0s_0`, and `-15e_0t_0`, at degrees `9,10,12`.
Every displayed leading coefficient is nonzero by the degree cell, so none
survives.

In the branch `e=0,f!=0`, the only new states are

```text
(deg(c),deg(f);deg(a),deg(d))=(2,3;4,7),(4,6;4,7).     (22d)
```

In `E2`, the `a,d` Wronskian has unique degree ten and coefficient
`(2deg(a)-deg(d))a_0d_0=a_0d_0`; the other blocks have degrees at most eight.
All remaining cap-seven ledgers equal their cap-six counterparts.  The
terminal list is exactly the inherited exact-empty cells `(2,3;2,3)` and
`(4,6;3,5)` from Sections 4--5.  Hence the affine width-three coefficient
space is empty through `D=7`; the first internal open cap is `D=8`.

## 6. Cubic-cover reframe

Under THM-3555's quadratic coordinate

```text
r^2=1-3w,                 w=(1-r^2)/3,                 (23)
```

every polynomial mixed thickening in `(1)` becomes even in `r`, and

```text
Jac_(v,r)(P,Q)=Jac_(v,w)(P,Q) dw/dr=-2r/3.             (24)
```

The actual Catalan polynomialization also has Jacobian `-2r/3`, but uses odd
terms in `r`; it is not an even pullback.  Hence the construction problem is
equivalently to find an even polynomial deformation with `(24)` that retains
the selected fiber.  Width is degree in `r^2`; the first genuinely new even
`r^6` cell is precisely the degree-six branch closed by `(19)--(22)`.

## 7. Classical degree transfer -- CITED input

If all coefficient polynomials in `(1)` have degree at most `D`, then

```text
deg P,deg Q <= max(3,D+N).                             (25)
```

Any solution of `(3)` retains the collision `(4)`, so it would be a
nonautomorphic planar Keller pair.  The
[Guccione--Guccione--Horruitiner--Valqui](https://arxiv.org/abs/2204.14178)
classification of hypothetical counterexample degree pairs below
height `125` is a **CITED** input: it gives height at least `108`, and below
`125` leaves only the reduced target-pencil degree pair `(72,108)` and its
transpose.  Combining it with `(25)`
proves

```text
D+N>=108.                                              (26)
```

If `D+N<125`, a linear target reduction would furthermore have to expose
the pencil degree values `72` and `108`; the displayed components themselves
may both have degree `108`.  Thus the first internally open cell `(N,D)=(3,8)`
is far below the global degree floor: at width three one needs `D>=105`, and
at width four one needs `D>=104`.  This does not solve those high-degree
cells; it prevents a low-width search from confusing an internally new
recurrence branch with a globally admissible counterexample.

## 8. Scope

This theorem excludes only ansatz `(1)` at widths `N<=2`, and at `N=3` only
through coefficient degree seven.  It does not exclude degree eight or higher,
width four, other boundary curves, or planar counterexamples.  It proves no
projective-closure statement; affine coefficient-space emptiness is the
relevant quantified claim.
