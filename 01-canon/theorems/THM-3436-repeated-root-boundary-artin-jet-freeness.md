---
id: THM-3436
title: "Repeated-root boundary Artin-jet freeness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For
  P=ax+b+g(x)z^d in characteristic zero, base-change to a finite normal
  splitting field and fix a repeated root alpha_i.  In target character
  sigma-1, every boundary quotient modulo (P-beta_i)^q is free over the
  Artin jet ring K'[lambda]/lambda^q: its rank is N=deg(rad(g)) exactly when
  alpha_i is THM-3433-selected, and N-1 otherwise.  The compatible formal
  completion is free of the same rank.  This is a local completion packet,
  not a direct-sum description of the integral sector or a new JC(2) case.
source: root boundary-jet session, 2026-08-15
audit: two-pass independent immutable-package and promotion-delta audit CLEAN at 7e44c5141b after the explicit C tensor K' repair; CRT, Neumann contraction, horizontal unit conjugacy, local-system rank, Artin freeness, formal-completion compatibility, wrap/nonsplit type, derived-evidence relabel, normal/-O/stored replay, pinned hashes, AST/security, documentation, and routing all checked
depends_on:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3419-generic-kummer-response-regular-sector-rank
  - THM-3433-all-sector-multiroot-primary-torsion-classification
related:
  - THM-3430-nonlinear-wrap-linearization-and-canonical-vertical-torsion
  - THM-3422-one-root-nonlinear-integral-hamiltonian-response
script: 04-computation/jc_multiroot_boundary_jet_packet_probe_20260815.py
output: 05-knowledge/results/jc_multiroot_boundary_jet_packet_probe_20260815.out
script_sha256: 3dcfc4f3f36ed658a64f3d45b7e055eb6bab10fa536522368a1aa3e32eab6332
output_sha256: 0d9908f5695a64895fd1d78cace03340b35bf9831e3ebd1eb53747eb563fc01c
hash_basis: LF-normalized bytes
---

# THM-3436 -- repeated-root boundary Artin-jet freeness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Let `K` be a field of characteristic zero, let `d>=2`, and put

```text
P=ax+b+g(x)z^d in K[x,z],              a in K*,           (1)
D=Jac(P,-),
C=K[x,z]/D(K[x,z])=direct_sum_(sigma=1)^d C_(sigma-1).    (2)
```

Assume that `g` is nonconstant.  Choose a finite normal extension `K'/K`
which splits `g` and contains the `d`th roots of unity, and write

```text
g=gamma product_(j=1)^N (x-alpha_j)^e_j,
beta_j=a alpha_j+b,                 e_j>=1,
N=deg(rad(g)).                                             (3)
```

Fix a repeated root `e_i>1`, a target character `1<=sigma<=d` (`sigma=d`
is wrap), and `q>=1`.  Base-change before using the geometric boundary
coordinate:

```text
Cbar_(sigma-1)=C_(sigma-1) tensor_K K',
lambda=P-beta_i,
R_q=K'[lambda]/(lambda^q).                                (4)
```

Define

```text
selected(i,sigma) iff
  d divides sigma(e_i-1), and
  d divides sigma e_j for every j!=i.                     (5)
```

Then there is an `R_q`-module isomorphism

```text
Cbar_(sigma-1)/lambda^q Cbar_(sigma-1) ~= R_q^c,          (6)

c=N     if selected(i,sigma),
c=N-1   otherwise.                                        (7)
```

The isomorphisms may be chosen compatibly in `q`, using one formal binomial
gauge.  Consequently

```text
lim_q Cbar_(sigma-1)/lambda^q Cbar_(sigma-1)
 ~= K'[[lambda]]^c.                                       (8)
```

For wrap, every repeated root is selected and `c=N`.  For a nonsplit root,
the statement is geometric after the declared faithful base change; Galois
conjugation transports the packet to its conjugates.  No canonical split
coordinates over `K` are asserted.

The theorem is deliberately local.  A Pruefer arm from THM-3433 is
`lambda`-divisible and disappears in every quotient `(6)`, so boundary jets
alone do not recover primary torsion or an integral direct sum.

## 2. Connection and loss ledger

| item | exact content |
|---|---|
| source | the integral Hamiltonian target sector after `K'/K` base change |
| target | a free packet over the principal Artin jet ring `R_q` |
| map | `lambda=yH` CRT, a vertical finite-Neumann contraction, and a horizontal formal unit gauge |
| preserved | character, boundary root, jet order, selected congruence, and free rank |
| destroyed | every `lambda`-divisible channel, the Pruefer arm, nonsplit primary coordinates, and the global torsion-free complement |
| needed sidecar | THM-3433 DeathBars/primary torsion, not fibre dimension alone |
| cheapest hostile | a simple root `g=x`, where `(x,a+z^d)` is not the unit ideal and the CRT mechanism fails |

The closest proved mechanism is THM-3433's selected-root law.  The canonical
hostile is its locally resonant but globally blocked profile
`(d,sigma;e_i,e_j)=(4,2;3,1)`.  The corrected near miss is inference of
freeness from residue dimension and total length: here freeness is instead a
consequence of an operator-level conjugacy.  The least-used sidecar is the
vertical nilpotent contraction, which prevents a thickened component from
contributing hidden extension classes.

## 3. Exact cokernel reduction and CRT

Work in `A=K'[x,z]`.  Put

```text
y=x-alpha_i,                 e=e_i,
g=y^e u,
lambda=a y+g z^d=yH,
H=a+y^(e-1)u z^d.                                      (9)
```

Because `e>1`, `H=a mod y`; hence `(y,H)=A`.  Therefore

```text
A/(lambda^q) ~= A/(y^q) direct-sum A/(H^q).              (10)
```

Moreover `D(lambda)=0` and

```text
D(y)=-d g z^(d-1),
D(H)=-H D(y)/y.                                          (11)
```

The quotient in the second line is polynomial because `y|g`.  Thus both CRT
ideals are `D`-stable.  Since `D(lambda^q A)=lambda^qD(A)`, right exactness
gives, sector by sector,

```text
Cbar_(sigma-1)/lambda^q Cbar_(sigma-1)
 =coker(D_sigma:A_sigma/lambda^q -> A_(sigma-1)/lambda^q).
                                                               (12)
```

There is no finite-versus-infinite-dimensional interchange in `(12)`; it is
the literal quotient of `A_(sigma-1)` by `D(A_sigma)+lambda^qA_(sigma-1)`.
Equation `(10)` splits this two-term complex into vertical and horizontal
factors.

## 4. The vertical factor is contractible

Modulo `y^q`, write

```text
D=a partial_z+R,                                          (13)
```

where every term of `R` raises `y`-adic order by at least `e-1`.  Let `J` be
coefficientwise integration in `z`, divided by `a`.  Then

```text
D J=1+R J.                                                (14)
```

The operator `RJ` is nilpotent on `A/(y^q)`.  Hence the finite sum

```text
J(1-RJ+(RJ)^2-...)                                        (15)
```

is a right inverse to `D` and respects the character grading.  The vertical
target cokernel is zero for every `q`, including the first nonreduced case
`q=2`.

## 5. The horizontal complex is constant over the jet ring

On `A/(H^q)`, the elements `y`, `u`, and `z` are units: their product
`y^(e-1)u z^d=-a+H` is a unit.  Eliminating `z^d` identifies the coefficient
ring in a fixed character with

```text
B_q=R_q[x,rad(g)^(-1)].                                  (16)
```

For an input `p z^sigma`, direct differentiation gives the target
coefficient operator

```text
A_lambda(p)
 =a sigma p+(lambda-a y)[sigma(g'/g)p-dp'].              (17)
```

Set

```text
s=1-lambda/(a y),
r=s^(-sigma/d) in B_q^*.                                 (18)
```

The binomial expansion of `r` truncates modulo `lambda^q`.  Differentiating
it with `lambda` held constant proves the exact identity

```text
A_lambda(rp)=r s A_0(p).                                 (19)
```

Multiplication by `r` on the domain and by `rs` on the target are units.
Thus the horizontal two-term complex is isomorphic to its special-fibre
complex tensored from `K'` to `R_q`:

```text
coker(A_lambda) ~= coker(A_0) tensor_(K') R_q.           (20)
```

This is the step that proves freeness.  Neither a residue-field dimension nor
an invariant-factor count is being substituted for `(19)`.  The formal
series in `(18)` has compatible truncations, proving `(8)` as well.

## 6. The special-fibre rank

At `lambda=0`, the horizontal Kummer component is

```text
z^d=-a/[y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j].        (21)
```

Its base is the affine line with the `N` roots of `g` deleted.  In character
`sigma`, the rank-one local system has monodromy exponents

```text
sigma(e_i-1),             {sigma e_j:j!=i} mod d.        (22)
```

It is trivial exactly under `(5)`.  A one-vertex, `N`-edge cochain model of
the punctured line has

```text
(dim H^0,dim H^1)=(1,N)       in the trivial case,
(dim H^0,dim H^1)=(0,N-1)     otherwise.                 (23)
```

Combining `(12)`, `(15)`, `(20)`, and `(23)` proves `(6)--(7)`.  Faithfully
adjoining roots of unity or enlarging the splitting field does not change
the rank; conjugate gauges are carried to one another by Galois action.

## 7. Equality and failure boundaries

- `sigma=d` makes every exponent in `(22)` divisible by `d`; the gauge is
  simply `s^(-1)`, so wrap has rank `N` at every repeated root.
- At `q=2`, `r=1+[sigma/(ad y)]lambda mod lambda^2`; the linear correction
  cancels the actual first extension term.  Nilpotent thickness is not being
  ignored.
- If `e_i=1`, then `H mod y=a+u(alpha_i)z^d`, so `(y,H)` is generally a
  proper ideal.  For `g=x` the components meet at the roots of `a+z^d`; the
  CRT and Neumann proof stops exactly there.  No simple-root analogue is
  claimed.
- Distinct boundary values cannot collide:
  `beta_i-beta_j=a(alpha_i-alpha_j)` and `a!=0`.
- The theorem does not identify the integral completion with a canonical
  direct summand, recover THM-3433 torsion, construct a polynomial mate, or
  settle any open part of `JC(2)`.

## 8. Exact companion

The companion verifies the factorization and `D`-stability, `28` explicit
finite-Neumann inverses through jet order four, `35` horizontal gauges through
order five, `7` explicit `q=2` coefficients, `10` wrap gauges, and `126,576`
special-fibre graph ranks.  It enumerates `506,304` **derived** Artin-packet
instances across `21,686` selected and `104,890` unselected geometric
root-character profiles.  That last loop is bookkeeping evidence; the
universal identity `(19)` is the proof of freeness.  Split, wrap, blocked, and
nonsplit controls are included.

Reproduce with

```bash
python3 04-computation/jc_multiroot_boundary_jet_packet_probe_20260815.py
python3 -O 04-computation/jc_multiroot_boundary_jet_packet_probe_20260815.py
```

Both modes reproduce the stored transcript byte-for-byte.

**QED.**
