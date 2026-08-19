# Mixed polynomial thickenings of the Catalan self-intersection

**Exact research note, 2026-08-18.**  This is deliberately not a theorem-ID
reservation and does not edit canon.  Within the displayed mixed-thickening
ansatz it proves a cap-free obstruction for transverse width `N<=2`, and an
exact affine coefficient-space obstruction for `N=3` through `v`-degree five.
It does not prove `JC(2)` and does not exclude width at least four or the first
new width-three degree-six branches.

Companion:
[`jc2_catalan_mixed_thickening_recurrence_kps_s187.py`](../04-computation/jc2_catalan_mixed_thickening_recurrence_kps_s187.py),
with stored output
[`jc2_catalan_mixed_thickening_recurrence_kps_s187.out`](../05-knowledge/results/jc2_catalan_mixed_thickening_recurrence_kps_s187.out).

```text
source_sha256 = 0444ad61a0bb2cd165243db1a97f0cb0b299eb19263378c97d1dee9ff39a7e1e
output_sha256 = cd23937963962bc4f43b83fc2f3ab6b477970c76f3151b7247c0025145380832
hash_basis    = LF-normalized bytes
```

## 1. Outcome and precise scope

Let `K` be a characteristic-zero field, put

```text
p(v)=v^2,                         q(v)=v^3-v,
P=p+sum_(j=1)^N a_j(v)w^j,
Q=q+sum_(j=1)^N b_j(v)w^j,       a_j,b_j in K[v].       (1)
```

The following statements are proved.

1. The exact coefficient recurrence for `Jac(P,Q)=1` is

   ```text
   E_k=sum_(i+j=k+1) (j a_i' b_j-i a_i b_j'),           (2)
   a_0=p, b_0=q,
   E_0=1,                     E_k=0  (1<=k<=2N-1).       (3)
   ```

2. There is no solution of `(1)--(3)` for `N=1`.
3. There is no solution for `N=2`, with no degree cap on the four correction
   polynomials.
4. For `N=3`, there is no solution satisfying

   ```text
   max_j(deg a_j,deg b_j)<=D                              (4)
   ```

   for `D=3,4,5`.  This is an exact emptiness statement for the affine
   coefficient variety, not a search that failed to find points.

At `w=0`, every map in `(1)` restricts to

```text
v |-> (v^2,v^3-v),
```

so the distinct points `(v,w)=(1,0),(-1,0)` retain the collision `(1,0)`.
Any polynomial solution would therefore be a noninjective planar Keller map.
The result says that this particular repair must either have `N>=4`, or, at
width three, contain a coefficient of `v`-degree at least six.

## 2. Inheritance pass and concept board

- **Closest proved mechanism:**
  [THM-3545](../01-canon/theorems/THM-3545-catalan-self-intersection-keller-thickening-boundary.md)
  gives the unique separated formal thickening and its nonterminating Catalan
  tail.
- **Canonical hostile:** the first three Catalan rows are a genuine partial
  solution: they satisfy `E_0=1,E_1=E_2=0` and first leak at
  `E_3=-135/16`.  A valid obstruction must preserve this prefix.
- **Corrected near miss:** polynomial truncation of the Catalan branch fails,
  but THM-3545 explicitly leaves mixed `v,w` corrections open.
- **Least-used sidecar:** the top coefficient equations are one-variable
  Wronskians.  Their zero set is not generic: unique factorization forces
  common powers before any Gröbner calculation.

The active board was:

1. the `w`-graded Jacobian recurrence;
2. Wronskian-zero/common-power rigidity;
3. the self-intersection collision at `w=0`;
4. affine coefficient varieties versus their projective closures;
5. the Catalan formal branch as a positive prefix control;
6. incoming coefficient-fibre/Segre and Kummer-puncture gates.

## 3. Derivation of the recurrence

Differentiate `(1)`:

```text
P_v=sum_i a_i'w^i,          P_w=sum_i i a_i w^(i-1),
Q_v=sum_j b_j'w^j,          Q_w=sum_j j b_j w^(j-1).
```

Collecting `w^k` in `P_vQ_w-P_wQ_v` gives `(2)` immediately.  For `N=3`,
write

```text
(a_1,a_2,a_3)=(a,c,e),       (b_1,b_2,b_3)=(b,d,f).
```

Then `(3)` is

```text
E0 = p'b-aq' = 1,
E1 = 2(p'd-cq')+a'b-ab' = 0,
E2 = 3(p'f-eq')+2a'd-ad'+c'b-2cb' = 0,
E3 = 3a'f-af'+2(c'd-cd')+e'b-3eb' = 0,
E4 = 3c'f-2cf'+2e'd-3ed' = 0,
E5 = 3(e'f-ef') = 0.                                  (5)
```

The companion derives these rows twice: once from `(2)`, and independently
by differentiating the full symbolic polynomials and collecting powers of
`w`.

## 4. Width one is impossible

For `N=1`, the equations are

```text
p'b-aq'=1,                         a'b-ab'=0.            (6)
```

Neither `a` nor `b` vanishes.  In `K(v)`, the second equation says

```text
(b/a)'=0,
```

so `b=lambda a` for `lambda in K`.  The first equation becomes

```text
a(lambda p'-q')=a(1+2lambda v-3v^2)=1,                 (7)
```

impossible because the displayed quadratic is not a unit.  This also closes
every later branch in which all rows above `(a,b)` vanish.

## 5. Width two is impossible without degree caps

Write

```text
P=p+aw+cw^2,                    Q=q+bw+dw^2.
```

Equation `(2)` gives exactly the four equations in the user's reduction:

```text
E0=p'b-aq'=1,
E1=2(p'd-cq')+a'b-ab'=0,
E2=2a'd-ad'+c'b-2cb'=0,
E3=2(c'd-cd')=0.                                     (8)
```

There are three exhaustive top-row cases.

### 5.1 `c!=0`

Equation `E3=0` gives `d=lambda c`, including `lambda=0`.  Put

```text
s=b-lambda a,                    H=lambda p'-q'
                                  =1+2lambda v-3v^2.   (9)
```

The first three rows become

```text
p's+aH=1,
2cH+a's-as'=0,
c's-2cs'=0.                                          (10)
```

If `s=0`, the first equation is the impossible factorization `aH=1`.
Otherwise

```text
(c/s^2)'=0,
```

so `c=mu s^2` for `mu!=0`.  If `r=deg s`, the first equation in `(10)`
forces `r>=1` and

```text
deg a=r-1.                                            (11)
```

In the second equation, `2cH` has degree `2r+2`, whereas

```text
deg(a's-as')<=2r-2.                                   (12)
```

The four-degree gap is a contradiction.

### 5.2 `c=0,d!=0`

Now `E2=0` gives `(d/a^2)'=0`, hence `d=mu a^2`; `a` cannot vanish by
`E0`.  If `r=deg a`, the equation `E0=1` forces

```text
deg b=r+1.                                             (13)
```

But `E1=0` reads

```text
4mu v a^2+a'b-ab'=0.                                  (14)
```

The first term has degree `2r+1`, while the Wronskian has degree at most
`2r`.  This is impossible, including `r=0`.

### 5.3 `c=d=0`

This is exactly the width-one system `(6)`.  Thus `N=2` is impossible in
all degrees.  This upgrades the requested bounded scans at `D=3,4` to a
cap-free theorem within the ansatz.

## 6. Width three through degree five

The proof of the bounded result is a finite degree decomposition of `(5)`.
Set

```text
A=deg a, B=deg b, C=deg c, Dd=deg d, E=deg e, F=deg f.
```

The bottom equation `E0=1` always forces

```text
B=A+1,                   0<=A<=D_cap-1.                (15)
```

Indeed the leading terms of `2vb` and `(3v^2-1)a` must cancel.  The top
Wronskian `E5=0` gives four exhaustive branches.

### 6.1 Both `e` and `f` are nonzero

Write `f=lambda e`, then put

```text
t=b-lambda a,                s=d-lambda c,
H=lambda p'-q'.                                           (16)
```

Equations `(5)` become the exact hierarchy

```text
p't+aH=1,
2p's+2cH+a't-at'=0,
3eH+2a's-as'+c't-2ct'=0,
2(c's-cs')+e't-3et'=0,
2e's-3es'=0.                                           (17)
```

Let `r=deg t`; then `deg a=r-1` by the first line.

If `s=0`, the fourth line gives `e=mu t^3`.  Caps at most five leave only
`r=1`; the second line then asks a nonzero constant Wronskian to cancel a
multiple of the quadratic `H`, which is impossible.

If `s!=0`, the last line gives

```text
2E=3 deg(s).                                            (18)
```

Through cap five the only possibilities are `(deg s,E)=(0,0),(2,3)`.

- In the constant case, the fourth line integrates to
  `c=(3e/(2s))t+constant`, so `deg c=r`.  The three degrees in the third
  line are `2`, at most `r-2`, and exactly `2r-1`; one is uniquely largest
  for every `r>=1`.
- In the `(2,3)` case, comparison in the third line leaves only the finite
  pairs printed in the transcript.  The fourth line kills every one: its
  `e,t` term has a strictly larger degree than its `c,s` Wronskian (or vice
  versa after the sole leading cancellation).  The survivor lists after
  that row are empty for caps three, four, and five.

### 6.2 `e=0,f!=0`

Now `E4=0` is

```text
3c'f-2cf'=0.                                           (19)
```

If `c=0`, rows `E1,E2` force

```text
Dd=2A-1,                    F=3A-3,                    (20)
```

but the top coefficient of `E3=3a'f-af'` is
`3A-F=3`, a contradiction.

If `c,f` are nonzero, `(19)` gives `3C=2F`.  Through cap five this leaves
only `(C,F)=(0,0),(2,3)`.

- For constants, `E3` gives `d=(3f/(2c))a+constant`.
  Degree comparison leaves only `A=1`.  Writing
  `a=1+2hv`, direct coefficients force simultaneously

  ```text
  c=h^2+3/4,                    c=h^2-3/2,              (21)
  ```

  whose gap is `9/4`.
- For `(C,F)=(2,3)`, the degree ledger leaves exactly

  ```text
  (A,C,Dd,F)=(2,2,3,3).                                (22)
  ```

  This is the exceptional square/cube branch treated in Section 7.

### 6.3 `f=0,e!=0`

The equation `E4=0` becomes

```text
2e'd-3ed'=0.                                           (23)
```

If `d=0`, rows `E1,E2` leave only

```text
(A,C,E)=(2,2,2)              at caps three and four,
(A,C,E)=(2,2,2),(3,4,5)      at cap five.              (24)
```

In every case the sole term `e'b-3eb'` in `E3` has nonzero leading
coefficient `E-3B=-7`.

If `d,e` are nonzero, `(23)` gives `2E=3Dd`.  Constants leave no state
through cap three; at caps four and five only `A=3` survives `E1`, after
which `E2` has the unique degree-seven term coming from
`c'b-2cb'`.  The alternative `(Dd,E)=(2,3)` has no simultaneous degree
survivor of `E1,E2`.  Thus this branch is empty.

### 6.4 `e=f=0`

This is the cap-free width-two obstruction in Section 5.

These cases exhaust every affine coefficient tuple.  No genericity,
nonvanishing of an unexamined coefficient, or rational-point assumption is
used.

## 7. The exceptional square/cube branch

Let `alpha,C_0,D_0,F_0` be the leading coefficients of `a,c,d,f` in `(22)`.
The degree-four parts of `E1,E2,E3` are

```text
4D_0-6C_0-(3/2)alpha^2=0,
6F_0+alpha D_0-6alpha C_0=0,
3alpha F_0-2C_0D_0=0.                                 (25)
```

Eliminating the first two linear unknowns from the third gives

```text
-(3/16)(alpha^2-4C_0)^2=0.                             (26)
```

Hence

```text
C_0=alpha^2/4,       D_0=3alpha^2/4,       F_0=alpha^3/8.  (27)
```

Equation `(19)` also gives `c^3/f^2=constant`.  Unique factorization, together
with degrees two and three, forces one common linear base.  Write

```text
a=1+2v(h_0+h_1v),
c=h_1^2(v-rho)^2,                  f=h_1^3(v-rho)^3,   (28)
```

where `h_1!=0`; equation `(27)` is exactly the leading normalization in
`(28)`.  Solve `E1=0` for `d`.  Divisibility by `v` requires

```text
(4h_0^2-4h_1^2rho^2-2h_1+3)/2=0.                      (29)
```

The coefficient of `v^2` in `E2` is

```text
9h_1(h_0+h_1rho)^2.                                   (30)
```

Thus `h_0+h_1rho=0`; equation `(29)` then forces `h_1=3/2`.  The constant
coefficient of `E2` becomes

```text
-9/4,                                                  (31)
```

the final contradiction.

As an independent algebraic replay, the companion adjoins
`zeta*h_1-1`, forms all coefficients of `(29),E2,E3`, and computes

```text
GroebnerBasis([1], zeta,rho,h_0,h_1; QQ, grevlex).     (32)
```

The same basis is `[1]` modulo `5` and `7`.  The rational basis `(32)` is
the proof; the modular runs are controls, not a modular-to-characteristic-zero
inference.

## 8. Solver status and controls

The exact status distinctions are:

- **Affine coefficient spaces:** empty exactly for `N=3,D<=5`, by the
  exhaustive branch proof above.
- **Exceptional ideal:** unit ideal over `QQ` after the necessary saturation
  `h_1!=0`.
- **Constant-top exceptional ideal:** also a unit ideal over `QQ` after
  saturating `h*c*f!=0`; this independently checks the `9/4` contradiction
  in `(21)`.
- **Projective closure:** not computed.  It may have points at infinity; they
  are not affine polynomial coefficients and are not needed for this claim.
- **Raw full Gröbner ideal:** an exploratory unreduced SymPy calculation did
  not complete on a useful timescale and was interrupted.  It contributes no
  verdict.  Triangularization reduced the only nonlinear survivor to `(32)`.
- **Positive controls:** `(P,Q)=(v,w)` returns rows `(1,0)`.  The cubic Catalan
  truncation returns

  ```text
  (E0,E1,E2,E3,E4,E5)=(1,0,0,-135/16,-405/64,-729/128). (33)
  ```

Ordinary and optimized Python runs are byte-identical after line-ending
normalization and match the stored output.  All 86 gates pass.

## 9. Connections to incoming work

### 9.1 Coefficient fibres and Segre leakage

[THM-3548](../01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md)
organizes the full Jacobian by equal exponent-sum fibres and warns that
fibrewise cancellation must retain global Segre cycle equations.  Recurrence
`(2)` is the `w`-graded one-dimensional slice of exactly that organization.
The source-to-target map is

```text
full exponent pair -> common w-sum k+1 -> row E_k.
```

It preserves the exact signed Jacobian contribution but forgets the separate
`v`-exponents until each row is expanded in `v`.  The `v`-degree ledger is the
required sidecar.  The top Wronskian is a two-channel fibre; its forced
proportionality is the analogue of a closed two-edge polygon.  The lower rows
then expose the cross leaks that a single cancelled top fibre cannot see.

### 9.2 One common-power mechanism in three places

The equations

```text
mA'B-nAB'=0                                             (34)
```

force exponent-ratio common powers in a UFD.  Here they create the
square/cube rows `(c,f)=(unit*u^2,unit*u^3)`.  In
[THM-3549](../01-canon/theorems/THM-3549-torus-quotient-correction-no-go.md),
the first open transverse cell has quartic/quintic leading rows with one
common base.  In
[THM-3550](../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md),
the top homogeneous Jacobian equation similarly forces two output leaders to
be powers of one form.  This is a genuine transfer:

```text
source: highest nonzero graded Jacobian row
target: one-variable logarithmic-derivative equation
preserved: vanishing of that row
lost: every lower Jacobian row
sidecar: the complete recurrence E_(top-1),E_(top-2),...
test: substitute the common powers and inspect the first leak.
```

In the present problem the first leak is decisive: it becomes the fixed
`-9/4` in `(31)`.

### 9.3 Kummer puncture versus mixed filling

[THM-3554](../01-canon/theorems/THM-3554-punctured-kummer-collision-surface-normal-form.md)
shows that the fixed three-variable collision restricts to the finite etale
cover `(s,b)->(b,4s^2)` only after puncturing `s=0`; filling the puncture
restores ramification.  Mixed corrections are one proposed way to change that
function field before filling.  The present obstruction says that the
smallest such correction is not a shallow perturbation: width two never
works, and width three must cross `v`-degree six.  This does not identify the
two constructions or prove a Kummer obstruction for `(1)`; it gives a concrete
complexity tariff on the suggested escape.

[THM-3555](../01-canon/theorems/THM-3555-catalan-thickening-universal-cubic-root-cover.md)
arrived during the final rebase and is now `PROVED + VERIFIED-EXACT`.  It
polynomializes THM-3545 by adjoining

```text
r^2=1-3w                                                   (35)
```

and identifies the resulting map with the universal marked-root cubic cover
`(t,p)->(p,-t^3-pt)`.  This gives a precise new reading of the present
calculation.  Pull any hypothetical polynomial solution `(P,Q)` of `(1)` back
along `(35)`.  Since `w=(1-r^2)/3`, both components become even polynomials in
`r`, while the chain rule gives

```text
Jac_(v,r)(P,Q)=Jac_(v,w)(P,Q) dw/dr=-2r/3.              (36)
```

The actual Catalan polynomialization has the same Jacobian divisor but uses
terms linear in `r`; it is not even.  Thus the finite mixed-thickening problem
can be reframed as follows:

```text
find an even-in-r polynomial deformation carrying Jacobian -2r/3,
then quotient the involution r->-r without losing the selected collision. (37)
```

The map in `(35)--(37)` preserves the exact Jacobian equation and
polynomiality after pullback.  It forgets the chosen local square-root sheet;
the sheet sign is the required sidecar.  The width is exactly the degree in
`r^2`.  Consequently the present results say that an even solution of `(37)`
cannot have `r^2`-degree at most two, and at `r^2`-degree three it needs a
coefficient of `v`-degree at least six.

THM-3555's first-jet theorem excludes corrections that fix the whole cubic
ramification line pointwise.  The current ansatz instead fixes the boundary
section `w=0` (equivalently `r^2=1`) and allows the ramification geometry to
move.  Neither theorem subsumes the other; together they isolate sheet parity
and ramification-line motion as the two sidecars a successful construction
must carry.

## 10. Sharp next targets

1. **Width three, degree six.**  This is the first cap where `(18),(19),(23)`
   acquire genuinely new multiples:

   ```text
   (deg s,deg e)=(4,6),
   (deg c,deg f)=(4,6),
   (deg d,deg e)=(4,6).                                (38)
   ```

   The next exact task is to parameterize the resulting fourth/sixth powers
   and run the lower-row leak map.  In THM-3555 coordinates this is the first
   new even `r^6` cell.  Repeating the cap-five ledger without these new
   branches would be invalid.
2. **Width four with low `v`-degree.**  The top ratio now has exponent pair
   `(3,4)` after the penultimate row.  Build the analogous common-power
   normal form before attempting a raw coefficient Gröbner basis.
3. **Automated UFD branch compiler.**  Given a terminal equation
   `mA'B-nAB'=0`, emit all zero/nonzero branches, common-power
   parametrizations, degree constraints, saturation variables, and the first
   two leak rows.  This is the reusable part of the calculation.
4. **Compare the first degree-six survivors with THM-3549's `(4,5)` cell and
   THM-3550's minimal `(6,8)` pencil skeleton.**  The comparison must retain
   total degree, transverse degree, common-base multiplicity, collision, and
   the lower-row leak; matching only the power exponents is insufficient.

The main methodological gain is triangular: solve the bottom Bézout row,
read the top Wronskian as a common-power law, and make the two fronts meet in
the middle.  A raw coefficient ideal hides both pieces and pays a much larger
Gröbner cost.
