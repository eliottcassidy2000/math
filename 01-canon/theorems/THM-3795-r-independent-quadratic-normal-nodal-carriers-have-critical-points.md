---
id: THM-3795
title: "R-independent quadratic-normal nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  THM-3785 cubic pseudo-plane, every carrier with nodal arm value e^2,
  nonzero first-normal constant, and r-independent canonical normal form
  A=e^2+z f(e)+z^2 h(e) has a critical point.  On the torus curve
  c^3+2re=0 its derivative is a nonmonomial Laurent polynomial whose three
  sources occupy disjoint exponent classes modulo 3.  Hence a surviving
  carrier must use genuine r-coupled data in its canonical normal form.
  This is a carrier obstruction, not a Darboux nonexistence theorem.
source: jc_zero_debt_lift / nodal quadratic-normal Laurent obstruction, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The monic z-normal
  form, exact torus parametrization and converse, all Hamiltonian signs,
  chain-rule proportionality, residue-class noncollision, Laurent-root
  reconstruction, algebraic-closure step, sharp z^3e^3 boundary, and the
  constant-r resultant (including finite-root and nonzero-coordinate seams)
  were rederived.  The surface is smooth and its Poisson packet is
  symplectic.  The deterministic companion uses active gates only; normal
  and optimized runs byte-match the frozen transcript and all recorded
  hashes match.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
related:
  - THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry
script: 04-computation/jc2_cubic_pseudoplane_quadratic_normal_carrier_thm3795.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_quadratic_normal_carrier_thm3795.out
script_sha256: 5e19cff995177659cdad8ba726e88636657ae2c50b119e22bff7b6d8714ee8ed
output_sha256: 472474bdf09236746f7c02549ed204436472220904032eb73a0a58f64da4f58f
semantic_sha256: 7c86c5fc372ca92cfa7b000b333dc38155131cd1c0166f3394e2b3144262a494
hash_basis: raw LF bytes
---

# THM-3795 -- an r-independent quadratic-normal carrier is critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
an algebraically closed field `k` of characteristic zero, fix `c in k*`, and
retain the smooth symplectic cubic pseudo-plane

```text
Y=Spec k[r,z,e]/(r^2e-z^3+c^3r)                       (1)
```

with Poisson packet

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re.  (2)
```

For arbitrary `f,h in k[e]` with `f(0)!=0`, put

```text
A=e^2+z f(e)+z^2 h(e).                               (3)
```

Then `A:Y->A1` has a critical point.

The phrase *r-independent canonical normal form* is intrinsic.  Since the
relation in `(1)` is monic in `z`, division gives a unique expression

```text
B=k[Y]=k[r,e] direct_sum z k[r,e] direct_sum z^2 k[r,e].             (4)
```

A regular function with arm value `e^2` whose three coefficients in `(4)`
are independent of `r` is exactly `(3)`.  Thus the theorem does not merely
test one quadratic ansatz: it closes the whole `r`-independent cell of the
canonical normal form with the nodal arm value.  In the Darboux problem the
first-normal Bezout law of THM-3790 forces

```text
f(0)=-1/(3c^3),                                       (5)
```

so every such cell is covered.

## 1. The dangerous curve is a one-dimensional torus

Put

```text
K=c^3+2re.                                            (6)
```

The curve `V(K)` on `Y` is parametrized by `z in k*` as

```text
e=-c^6/(4z^3),             r=2z^3/c^3.               (7)
```

Both `(1)` and `(6)` vanish identically under `(7)`, and conversely `K=0`
forces `r,e,z` to be nonzero and gives these formulas.  Write

```text
a=-c^6/4,
Phi(z)=A|_(K=0)=a^2z^(-6)+z f(az^(-3))+z^2h(az^(-3)).                (8)
```

Because `A` is independent of `r`, its Hamiltonian components are

```text
{A,r}=-3r^2 A_z-9z^2 A_e,
{A,z}=-3K A_e,
{A,e}= 3K A_z.                                       (9)
```

The last two vanish on `(7)`.  Moreover

```text
de/dz=3c^6/(4z^4),             r^2(de/dz)=3z^2,       (10)
```

so the remaining component is exactly

```text
{A,r}|_(K=0)=-3r^2 Phi'(z).                           (11)
```

Consequently every nonzero root of `Phi'` reconstructs an actual critical
point by `(7)`; no elimination converse is being assumed.

## 2. Three exponent classes cannot cancel

Expand

```text
f(e)=sum_(j>=0) f_j e^j,       h(e)=sum_(j>=0) h_j e^j.
```

Differentiating `(8)` gives the finite Laurent polynomial

```text
Phi'(z)=
 -3c^12/8 z^(-7)
 +sum_(j>=0) (1-3j)f_j a^j z^(-3j)
 +sum_(j>=0) (2-3j)h_j a^j z^(1-3j).                 (12)
```

The three displayed sources occupy respectively the exponent classes

```text
2 mod 3,                    0 mod 3,                    1 mod 3.       (13)
```

They therefore cannot cancel across sources.  Within either sum each
exponent occurs only once, and characteristic zero makes every factor
`1-3j` and `2-3j` nonzero.  In particular `(12)` contains both of the
nonzero terms

```text
-3c^12/8 z^(-7)                         and             f(0).         (14)
```

Thus `Phi'` is not a Laurent monomial.  Multiply it by the inverse of its
lowest power of `z`.  The result is a nonconstant polynomial with nonzero
constant term.  Algebraic closedness supplies a root `z_0`, and the
constant term makes `z_0!=0`.  Equations `(7),(9)--(11)` then give an actual
critical point of `A` on the smooth symplectic surface `Y`.

## 3. Why the next coordinate really is r-coupled

The modulus three in `(13)` is the ramification index of the torus law
`e=az^(-3)`.  A monomial `z^q e^j` contributes derivative exponent
`q-1-3j`.  The first higher `z`-power able to collide with the mandatory
`-7` exponent is

```text
q=3, j=3:             d/dz(z^3e^3)|_(K=0) has exponent -7.           (15)
```

Indeed the hostile boundary

```text
A_sharp=e^2-z/(3c^3)+(4/c^6)z^3e^3                  (16)
```

restricts on `(7)` to `-z/(3c^3)` and has no critical point on this torus.
But `(1)` rewrites its new term as

```text
(4/c^6)z^3e^3=(4/c^6)r^2e^4+(4/c^3)re^3,            (17)
```

which has genuine `r`-dependence in the unique normal form `(4)`.  Thus
the theorem's boundary is exact: a smooth nodal carrier, if one exists,
must use `r`-coupled data (or change the arm profile).  Formula `(16)` only
defeats this particular torus witness; no smoothness or Darboux claim is
made for it.

The cheapest `r`-coupled correction is not yet enough.  As a hostile
control set `c=1` and

```text
A_b=e^2-z/3+b r,                         b!=0.         (18)
```

For a critical point with `z!=0`, put `u=r/z` and recover
`e=(z^2-u)/(u^2z)` from the surface equation.  Two Hamiltonian equations
are then equivalent to

```text
p_1=18z^2-u^4z-18u=0,
p_2=81bz^2-u^3z-9=0.                                  (19)
```

Their exact `z`-resultant is

```text
81(26244b^2u^2+9bu^8-5832bu-2u^7+324).               (20)
```

It has degree eight and nonzero constant term for every `b!=0`, hence a
nonzero root `u`.  The two leading `z`-coefficients are `18` and `81b`, so
there is no degree-drop root at infinity: resultant zero gives a common
finite root, and that root has `z!=0` because `p_2(u,0)=-9`.  The recovered
point lies on `Y`, and the third Hamiltonian follows from the surface
Casimir relation.  Thus even constant `r`-repair remains critical.
This control does **not** close nonconstant `r g(e)`, which is the live
carrier lane.

The deterministic companion named in the metadata checks `(7)--(20)`,
3,535 active residue-class, coefficient, reconstruction, and hostile gates,
the THM-3792 `h=0` boundary, and the sharp cubic collision.  Normal and
optimized executions byte-match the frozen transcript, and the source has
no inactive Python assertions.
**QED.**
