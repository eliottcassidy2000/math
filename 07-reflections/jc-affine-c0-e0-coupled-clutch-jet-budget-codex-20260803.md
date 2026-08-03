# Coupled affine clutch: one more residual degree, no smaller obstruction

**Status: NON-CANONICAL SYNTHESIS AROUND PROVED, VERIFIED-EXACT,
INDEPENDENTLY HOSTILE-AUDITED THM-3289.**  This reflection is not a proof
source.  Every exact claim below must be checked against
[THM-3289](../01-canon/theorems/THM-3289-affine-transverse-c0-e0-coupled-clutch-critical-no-go.md)
and its frozen primary
[`companion`](../04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py)
and
[`output`](../05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out),
plus the independent
[`audit`](../04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.py)
and its
[`output`](../05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.out).

## Why this was the next coordinate

THM-3279 closes the affine `C_0` lane when `E_0` is constant.  Its mechanism
is a two-jet budget at the simple owner branch `S`: the value of `C_0` kills
one boundary coefficient, its slope kills the next, and a fifth-jet invariant
then stops further contact.  The least-used coordinate in the same local
gradient system was the slope of `E_0`.  It enters the reduced second gradient
equation as a genuinely new term,

```text
R_2=V^3 E_0'+V^2y+L(-V'y+2V^2C_0'),                    (1)
```

so it was not merely a reparametrization of the old `C_0` walls.

The tested enlargement is

```text
C_0=c_0+dx,               E_0=e_0+kx,                  (2)
```

with `dk!=0`.  The constant-`C_0` and constant-`E_0` faces are already owned
by THM-3212 and THM-3279.  Thus `(2)` isolates the first simultaneous
transverse deformation rather than re-running a saturated parameter family.

## One gate controls both the finite fibres and all T branches

At a root `alpha` of `g=ST`, the only zero of `P_z` is
`z=-C_0(alpha)`, and the remaining derivative is

```text
Delta(alpha)=k-A'(alpha)C_0(alpha).                     (3)
```

This is the exact finite-clutch gate.  What is structurally satisfying in the
proved calculation is that `(3)` reappears without translation in every
`T`-boundary leading row.  If `ord_alpha(V)=m` and
`A'(alpha)=2/(2-m)`, then

```text
[t^(3m-1)]K
 =16m(m-1)/(m-2) v^3 Delta(alpha).                      (4)
```

So there are not two unrelated arguments, one local and one at the
saturated resultant.  The source is the residue of the finite gradient; the
target is the first unremoved `T` coefficient of the resultant; the map is
local elimination in `y=Vz`; and the preserved predicate is nonvanishing of
the same scalar `Delta`.  What elimination destroys is the actual critical
`z`-coordinate, but the finite formula supplies it whenever the gate fails.

This is a useful procedural pattern: before launching a separate boundary
PRS, ask whether its leading coefficient is already the finite clutch
sidecar.  Here that observation disposes of all four `T` branches at once.

## The correct S wall and why its apparent exception is hostile

At the simple `S` branch, write

```text
V=v_1t+v_2t^2+...,
C_0=c+dt.                                                (5)
```

The first nontrivial coefficient factors as

```text
q_3=(8/3)v_1^2(2c-k)(6cv_1^2+3kv_1^2+4v_2).            (6)
```

The first factor is not a new wall: `A'(S)=2`, so `k=2c` is exactly failure
of `(3)`.  The genuine post-clutch wall is

```text
k=-2c-4v_2/(3v_1^2).                                   (7)
```

Keeping both factors visible matters.  A premature mental simplification of
`(6)` suggested `k=-2c`; recomputing the literal factor caught the missing
`v_2` shift before any artifact was frozen.  The reliable procedure was to
name the finite gate first and only then solve the other factor.

On `(7)`, `q_4` is affine in `d`, with slope proportional to

```text
3cv_1^2+v_2.                                             (8)
```

At first `(8)=0` looks like an exceptional branch on which the slope control
is lost.  Substituting it back into `(7)` gives `k=2c`, however, so the entire
exceptional slope lies on the already excluded finite-clutch wall.  This is a
small but reusable hostile test:

```text
when a wall-solving denominator vanishes, substitute it into every earlier
gate before treating it as a new case.                                  (9)
```

In the live lane, `q_4=0` therefore determines a unique `d`.

## The fifth/sixth-jet sidecar

After the live `k` and `d` walls, the remaining two coefficients have the
form

```text
q_5 = unit * (3cv_1^2+v_2) Q_5(c),
q_6 = unit * Q_6(c),                                    (10)
deg((3cv_1^2+v_2)Q_5)=3,       deg Q_6=4.
```

Exact Euclidean algorithms over each of the two cubic accessory fields give

```text
(4,3) -> (3,2) -> (2,1) -> (1,0).                       (11)
```

After saturating by the excluded factor `(8)`, the live profile is

```text
(4,2) -> (2,1) -> (1,0).                                (12)
```

Thus the proved mechanism is not that one of `q_5,q_6` is universally
nonzero.  Each can vanish on a genuine algebraic locus.  The obstruction is
compatibility: their zero loci have empty intersection, witnessed by a unit
PRS terminal.  This is precisely the kind of two-polynomial sidecar that a
single discriminant or a numerical sample would erase.

The hostile checks are layered:

1. retain the excluded slope factor and prove the full degree-`3/4` gcd is
   one;
2. remove that factor and independently prove the live degree-`2/4` gcd is
   one;
3. perform both calculations in both cubic accessory fields; and
4. exhibit `C_0=1+x,E_0'=1` controls whose degree-53 residuals are squarefree
   and disjoint from `ST`.

## What the E slope buys, exactly

For affine `C_0` with constant `E_0`, THM-3279 has

```text
deg K=96,       deg H=52,       ord_S H<=2.              (13)
```

Turning on `k` introduces the unique infinity monomial

```text
128 C_0 V^6 d k^2,                                      (14)
```

whose affine leading term is `128d^2k^2V_lead^6`.  The coupled invoice is
therefore

```text
deg K=97,       deg H=53,       ord_S H<=3.              (15)
```

The new slope buys exactly one additional possible unit of `S` contact, but
it also creates exactly one additional residual degree at infinity.  The
off-owner lower bound does not move:

```text
old lane: 52-2=50,
new lane: 53-3=50.                                      (16)
```

This explains the no-go more sharply than “another parameter failed.”  The
parameter has real geometric effect—it changes both the resultant degree and
the contact budget—but those effects balance.  Low-order transverse freedom
is being charged one-for-one by the infinity divisor.

The mechanism suggests a candidate meta-pattern, not yet ready for promotion:

```text
an added jet that raises owner contact and infinity degree by the same amount
cannot improve the residual critical-point invoice.                    (17)
```

Evidence currently comes only from the constant-`E_0` and affine-`E_0`
clutch lanes, so `(17)` should remain a local research heuristic rather than
a universal principle.

## What information is still missing

The resultant calculation preserves the response identity, the finite owner
gate, and existence of common gradient zeros.  It destroys or never records:

- a cofactor or Bezout certificate for the two gradient equations;
- the branchwise data needed to assemble a second polynomial coordinate;
- the Jacobian equation for an unknown mate `Q`;
- an inverse cover or global chart; and
- deformations of `V,A` themselves that preserve the response geometry.

Consequently the next useful coordinate is unlikely to be another scalar
coefficient of `C_0` or `E_0`.  A quadratic `E_0` might buy another local wall,
but `(14)--(17)` predict another infinity charge and give a cheap reason not
to privilege that route.

## Follow-up: the rational controls are global critical graphs

The finite-exact
[linear-subresultant scout](jc-affine-c0-e0-linear-subresultant-critical-section-codex-20260803.md)
retains the first lost cofactor for `C_0=1+x,E_0'=1` in both accessory fields.
Its chosen standard sequence has degree profile `(3,3,2,1,0)`, and the linear
row has the exact form

```text
S_1=unit * boundary_i * (a(x)y+b(x)),
deg(a)=36,       deg(b)=44.
```

On the squarefree degree-`53` residual `H`, the boundary, `V`, and `a` are all
units.  Exact substitution gives

```text
y=-b/a,       R_1=R_2=0 in K_i[x]/(H).
```

Thus each reduced residual critical scheme is one global graph, with a
degree-`52` representative for `y` and hence also `z=y/V`.  No rootwise
selector split is needed.  This upgrades resultant existence to a coherent
critical section for two controls only; it is not an inverse cover, mate, or
statement about every affine parameter.

## Next problem: generalize the section or retain a second coordinate

The next underexplored object should live before scalar elimination.  Two
concrete versions are:

1. **Parameter-uniform subresultant sidecar.**  Retain the affine-parameter
   version of

   ```text
   U R_1+W R_2=boundary * (a y+b),                      (18)
   ```

   and compute the degeneracy divisor `Res_x(a,H)` in parameter space.  Where
   `a` vanishes, keep the homogeneous pair `(a,b)` rather than divide it away.
   The rational controls are positive controls; the finite-clutch and
   exceptional-slope walls are hostiles.

2. **Second-coordinate jet.**  Introduce a bounded `z`-degree ansatz for `Q`
   and impose

   ```text
   P_xQ_z-P_zQ_x=1.                                     (19)
   ```

   Reduce `(19)` first modulo the owner branches and then through the same
   `y=Vz` localization.  The cheapest decisive probe is whether the finite
   residues of the unknown `Q` coefficients can satisfy all five branches
   while respecting one infinity normalization.  A failure would identify a
   missing branchwise coordinate; a positive signal would finally carry more
   information than critical-point nonexistence.

The source of either map is the coupled gradient pair `(R_1,R_2)`.  Its target
is, respectively, a syzygy module or a truncated mate module.  The preserved
predicate is the exact response and owner geometry.  The essential retained
information is the cofactor/branch label that the scalar resultant destroys.

The stopping reason for the present lane is therefore structural, not merely
computational: simultaneous affine `C_0,E_0` freedom spends its new jet on a
matching infinity degree.  Retaining one earlier object recovers a global
critical graph for the controls; progress now requires its parameter-uniform
degeneracy locus or a genuine second-coordinate ansatz, not one more scalar
coefficient.
