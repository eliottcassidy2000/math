# The divergence obstruction is one step above the unit class in an exact torsion ladder

**Status:** PROVED + GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE BANK; PROMOTED AS
[THM-3318](../01-canon/theorems/THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z.md).
Over any characteristic-zero field, the family

```text
P_(r,lambda)=x+lambda x^r z,       r>=2, lambda!=0,
```

has unimodular gradient but nonzero canonical divergence class.  More
precisely, its divergence class generates a length-`r` `K[P]`-torsion block,
the ordinary unit-response obstruction is its first `P`-multiple and has
length `r-1`, and both become exact after deleting the special fibre `P=0`.
This classifies the requested `lambda=1` family and its nonzero scalar
deformation.  It proves no new case of `JC(2)`.

The exact companion is
[`jacobian_divergence_family_agent_exact.py`](../04-computation/jacobian_divergence_family_agent_exact.py),
with frozen
[`output`](../05-knowledge/results/jacobian_divergence_family_agent_exact.out).

## Inheritance pass

The closest mechanism is the canonical class introduced in
[`jc-critical-inverse-cover-cofactor-jacobian-probe-agent-20260803.md`](jc-critical-inverse-cover-cofactor-jacobian-probe-agent-20260803.md):
after a Bezout row `AP_x+BP_z=1` is found,

```text
mu(P)=[A_x+B_z] in coker(D_P),
D_P=P_x partial_z-P_z partial_x,
```

is independent of the row, and `mu(P)=0` is equivalent to the existence of a
polynomial mate.  That reflection handled only `P=x+x^2z` as a hostile.

The historical provenance is the current active lens
[`HYP-8950`](../05-knowledge/hypotheses/INDEX.md), which asks how Hamiltonian
cokernels connect local face obstructions to global polynomial exactness.
Its source reflection is
[`S109`](jc2-progress-the-hamiltonian-cokernel-realizes-the-manufactured-valuation-and-the-dvdk-face-deathstar-S109.md).
The DvdK-face wording there is explicitly corrected by
[`MISTAKE-244`](../01-canon/MISTAKES.md); only the Hamiltonian-cokernel and
weight viewpoints are inherited here.

Three proved neighbors fix the scope.

- [`THM-2045`](../01-canon/theorems/THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate.md)
  already excludes this family by a weighted-sector argument: take its
  parameters `(a,b,u,s)=(1,-lambda,r-1,1)`.  The new result is not a broader
  no-mate theorem; it computes the exact cokernel classes and their torsion.
- [`THM-2230`](../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md)
  identifies polynomial response fibres once a polynomial mate exists.  Here
  there is no such mate.  A Laurent response fibre will instead be computed
  directly, without extending THM-2230 beyond its hypotheses.
- [`THM-3306`](../01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md)
  and its
  [exceptional-quadratic blow-up note](jc-affine-c-exceptional-quadratic-c2-blowup-codex-20260803.md)
  carry a nonsplit `C_2` root deck at a critical base locus.  The localized
  family below has a Kummer deck, but only after the typed comparison near the
  end; no common geometric carrier is assumed.

The live board was:

| object | operation | retained fact | cheapest hostile |
|---|---|---|---|
| gradient Bezout row | divergence | canonical class `mu` | reverse the bracket sign |
| `coker(D_P)` | `K[P]` action | exact annihilator | inspect the last uncancelled pole |
| localization at `x` | integrate at fixed `P` | Laurent mate and correction | demand polynomial extension |
| local inverse map | forget the `x` root | Kummer deck of degree `r-1` | `r=2`, where the deck is trivial |

## 1. Gradient unit and the canonical representative

Let `K` be a characteristic-zero field, `R=K[x,z]`, `r>=2`, and
`lambda in K*`.  Put

```text
P=x+lambda x^r z,
P_x=1+lambda r x^(r-1)z,       P_z=lambda x^r.          (1)
```

The following polynomial row is an exact Bezout identity:

```text
A=1-lambda r x^(r-1)z,
B=lambda r^2 x^(r-2)z^2,
AP_x+BP_z=1.                                             (2)
```

Its divergence is

```text
A_x+B_z=lambda r(r+1)x^(r-2)z.                          (3)
```

Consequently the gradient is unimodular and

```text
mu=[lambda r(r+1)x^(r-2)z] in C_P:=R/D_P(R).            (4)
```

The quotient `C_P` is a `K[P]`-module, not a quotient ring, because

```text
D_P(F(P)h)=F(P)D_P(h).                                  (5)
```

For sign control, another Bezout row has the form

```text
(A+hP_z, B-hP_x),
```

whose divergence is `(A_x+B_z)-D_P(h)`.  Thus `(4)` uses
`D_P(f)=P_xf_z-P_zf_x=Jac(P,f)`.  Reversing the bracket changes the response
sign; the companion tests both conventions.

## 2. The localization computes every possible correction

Localize at `x`.  There is an exact coordinate presentation

```text
S=R[x^(-1)]=K[P,x,x^(-1)],
z=(P-x)/(lambda x^r).                                   (6)
```

Holding `P` fixed gives

```text
D_P=-lambda x^r partial_x,       ker(D_P:S->S)=K[P].    (7)
```

The kernel statement follows by differentiating a Laurent expansion
`sum_n a_n(P)x^n`; characteristic zero forces every coefficient with `n!=0`
to vanish.

Two exact Laurent primitives are

```text
h_0=r z/x-lambda^(-1)x^(-r),
Q_0=x^(1-r)/(lambda(r-1)),                              (8)

D_P(h_0)=lambda r(r+1)x^(r-2)z,
D_P(Q_0)=1.                                             (9)
```

The correction is not merely existential.  It turns the polynomial Bezout
row into the closed Laurent row of the mate:

```text
A+h_0P_z=0=Q_(0,z),
B-h_0P_x=lambda^(-1)x^(-r)=-Q_(0,x).                   (10)
```

Equations `(7)--(9)` classify the entire localized response fibres:

```text
{h in S:D_P(h)=A_x+B_z}=h_0+K[P],
{Q in S:D_P(Q)=1}=Q_0+K[P].                            (11)
```

Neither fibre meets `R`.  In the first line, every element has the uncancelled
ambient pole `-lambda^(-1)x^(-r)`; in the second, every element has the pole
`x^(1-r)/(lambda(r-1))`.  Adding `H(P)` is polynomial in the original
coordinates and cannot cancel either pole.  Hence `mu!=0`, and no polynomial
mate exists, uniformly for every `r>=2`.

## 3. Exact annihilators and the one-step bridge

Let

```text
theta=[1] in C_P.                                       (12)
```

The two localized primitives give exact polynomial torsion witnesses.  With

```text
v=lambda x^(r-1)z,             P=x(1+v),                (13)
```

one has

```text
P^(r-1)Q_0 = (1/(lambda(r-1)))(1+v)^(r-1) in R,
P^r h_0    = (1/lambda)(1+v)^r(rv-1) in R.              (14)
```

Therefore `P^(r-1)theta=0` and `P^r mu=0`.  These powers are sharp.  If a
nonzero `F(T) in K[T]` has `T`-adic order `m`, write
`F(T)=T^mG(T)` with `G(0)!=0`.  By `(7)`, any polynomial primitive of
`F(P)` or `F(P)(A_x+B_z)` would have to be respectively

```text
F(P)Q_0+H(P),          F(P)h_0+H(P),       H in K[T].   (15)
```

Their boundary orders are `m+1-r` and `m-r`.  If these are negative, the
leading coefficient is a nonzero multiple of `G(0)` and cannot be cancelled
by `H(P)`.  Conversely `(14)` handles every `m` at or above the displayed
threshold.  Thus

```text
Ann_(K[P])(theta)=(P^(r-1)),
Ann_(K[P])(mu)=(P^r).                                  (16)
```

There is also an exact bridge between the two classes.  Put

```text
g_0=(r-1)z+lambda r x^(r-1)z^2.                        (17)
```

Directly,

```text
P h_0+(r-1)Q_0=g_0,
D_P(g_0)=P(A_x+B_z)+(r-1).                             (18)
```

Therefore

```text
P mu=-(r-1)theta.                                      (19)
```

The complete cyclic picture is

```text
K[P]mu  ~= K[T]/(T^r),
K[P]theta=P K[P]mu ~= K[T]/(T^(r-1)).                  (20)
```

This is the main structural gain.  The canonical divergence class is not an
unrelated second obstruction: inside its explicitly classified cyclic block,
it is a one-step torsion lift of the original unit-response class.

Both classes vanish after inverting `P`.  Indeed `P=x(1+v)`, so `x` is a unit
on `P!=0`, and `(8)` applies.  Hence the obstruction is an exact special-fibre
extension failure, not a class that survives on the generic fibre.  The
generic nonzero fibre is `G_m`, but its nonzero `H^1` should not be conflated
with these particular torsion classes.  This sharpens the local/global
language around HYP-8950: generic fibre cohomology and the location of the
specific class `[1]` are different data.

## 4. The Laurent mate has a Kummer deck

Although `Q_0` is not polynomial, the localized pair has

```text
Jac(P,Q_0)=1.                                           (21)
```

Let `p=P`, `q=Q_0`, and `w=x^(-1)`.  Then

```text
w^(r-1)=lambda(r-1)q.                                  (22)
```

Over `K(p,q)`, the polynomial in `(22)` is irreducible by Eisenstein at the
prime `q`.  After inverting `q`, its derivative `(r-1)w^(r-2)` is a unit.
Thus the localized map `(x,z)->(p,q)` is a finite etale Kummer cover of degree
`r-1`, with geometric deck group `mu_(r-1)`.

- At `r=2` the deck is trivial and the localized map is an isomorphism, yet
  the polynomial extension still fails.  This is the hostile showing that a
  deck is not the cause of `mu!=0`.
- At `r=3` the local deck is `C_2`.  This is the only literal group-level
  resemblance to the exceptional quadratic under THM-3306.
- At `r>3` the forgotten inverse label is an `(r-1)`-st root, while the pole
  obstruction remains the same typed polynomial-extension failure.

The typed comparison with the THM-3306 quadratic is:

| field | this family | THM-3306 exceptional quadratic |
|---|---|---|
| source | noncritical open `Spec R[x^(-1)]` | blow-up of a critical coefficient base ideal |
| target map | `(x,z)->(P,Q_0)` | strict transform of the quadratic subresultant |
| forgotten label | the Kummer root `w=x^(-1)` | one of two common critical `y` roots |
| preserved predicate | local Jacobian response `1` | common-gradient-root incidence |
| destroyed information | extension across `x=0` | individual root over the residue field |
| sidecar | an `(r-1)`-st root from `(22)` | `sqrt(delta)` |
| decisive hostile | `r=2`: obstruction with no deck | nonsquare norm and first-normal motion |

There is no morphism between these carriers and no transfer of THM-3306's
degree-72 result.  The comparison only identifies the common operation:
forgetting a root label creates a finite cover that must be retained as a
sidecar.  Here that cover lives after localization and is separate from the
failure to extend polynomially.

## 5. Exact audit

The companion performs two independent kinds of checks.

1. It verifies the formulas `(2)--(3)`, `(8)--(10)`, `(14)`, `(18)`, and
   `(21)--(22)` symbolically with formal `r` and `lambda`.
2. It checks `r=2,...,9` at `lambda in {1,-1,2,-3}`, including the reversed
   bracket sign, sharp boundary orders, polynomial killing powers, and exact
   Kummer factor degrees.  For `r=2,...,8` it also builds the full matrix of
   `D_P` on the rectangular universe

   ```text
   0<=deg_x<=2r+2,             0<=deg_z<=4,
   ```

   verifies that adjoining the `mu` representative raises the exact rank by
   one, and verifies a known image target as a positive control.  This bounded
   matrix bank is diagnostic only; the localization-and-pole proof supplies
   the all-degree result.

The source has zero assertion nodes and zero floating literals.  Normal and
optimized runs match the frozen transcript line for line.  The LF-normalized
digests are

```text
script  fc95ebad308c9d953c5a331a6c545c273907665cbc24bb0b32639775492f3731
output  004b82f379ae9726cb569f561349f89f23b41737945d1a581bbaa5e35675a22a
```

## Scope and nonconsequences

This result classifies one explicit non-Keller family.  The Jacobian
conjecture starts with a **polynomial** pair already satisfying constant
Jacobian.  Here every polynomial first coordinate has been proved to admit no
polynomial mate; the displayed mate uses negative powers of `x`, and for
`r>2` its localized map is not even birational.  A unit Jacobian on the
non-affine chart `x!=0` is neither a Keller map of `A^2` nor a counterexample.

The useful export is instead a testable search object.  For a future
gradient-unimodular candidate, one should determine whether `theta=[1]` and
`mu` lie in a finite `K[P]`-torsion block, whether one is a `P`-multiple of
the other, and which deleted divisor trivializes them.  That refines a blind
mate ansatz into a torsion-order and extension problem, while `JC(2)` and
`DC(2)` remain open.
