---
id: THM-3552
title: "Two-channel cyclic-fiber holomorphic-exactness obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let p<q be coprime
  integers at least two, T=x^p y^q, phi a nonconstant squarefree polynomial
  with phi(0) nonzero, and Psi'=phi H nonzero.  Then
  P=x phi(T)+Psi(T) is a polynomial submersion but has no rational, hence no
  polynomial, constant-Jacobian mate.  On the normalized generic Kummer
  fibre every mate would make y dT/(q phi T) exact, while an exhaustive
  valuation ledger makes that differential nonzero and holomorphic.  The
  explicit p=2,q=5 example has a four-vertex Newton polygon and generic
  genus ten.  This is a first-coordinate no-mate theorem, not JC(2).
source: root-2026-08-18-planar-jacobian-counterexample-hostiles
depends_on: []
related:
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
  - THM-3551-one-ray-planar-jacobian-mate-no-go
script: 04-computation/jc_two_channel_kummer_obstruction_thm3552.py
output: 05-knowledge/results/jc_two_channel_kummer_obstruction_thm3552.out
script_sha256: cef53b6f0d12c0a9a4c9a74547a33f2db68ec82f84f602d75903147bd12ccbe2
output_sha256: ed088e1d47fb46328d2b34429a293db5d5b0059ce351f2ce3213f5e82043ac22
semantic_sha256: 12cf6002c32d9569d68ea0457ea0b28677a6141854b2470adf1171f4c3c1e4e0
hash_basis: LF-normalized bytes
---

# THM-3552 -- two-channel cyclic-fiber holomorphic-exactness obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Let `k` be a field of characteristic zero.  Choose coprime integers

```text
2 <= p < q,
```

put `T=x^p y^q`, and take polynomials `phi,Psi,H in k[T]` satisfying

```text
deg(phi)>=1,       phi(0)!=0,       phi squarefree,
Psi'=phi H!=0.                                             (1)
```

Define

```text
P=x phi(T)+Psi(T).                                        (2)
```

Then:

1. `(P_x,P_y)=(1)`, so `P` has no affine critical point after any algebraic
   extension of `k`;
2. there is no `Q in k(x,y)` and no `kappa in k*` for which
   `Jac(P,Q)=kappa`.

In particular, `(2)` has no polynomial Jacobian mate.  The second additive
channel `Psi(T)` can make the Newton polygon genuinely two-dimensional and
the generic fibre arbitrarily complicated, but it does not provide an exact
Hamiltonian response.

## 1. The gradient gate really passes

Write

```text
A=x phi'(T)+Psi'(T).
```

Direct differentiation gives

```text
P_x=phi(T)+(pT/x)A,             P_y=(qT/y)A.              (3)
```

Because `p,q>=2`, the polynomial interpretations of the quotients in `(3)`
vanish appropriately on either coordinate axis, where `P_x=phi(0)!=0`.
In the torus, a common zero of `P_x,P_y` would give

```text
A=0,             phi(T)=0.                                (4)
```

At a root of `phi`, condition `(1)` makes `Psi'=0`, so `(4)` reduces to
`x phi'(T)=0`.  Both factors are nonzero: `x` is a torus coordinate and
`phi` is squarefree.  This contradiction proves the unit-gradient claim.

The divisibility in `(1)` is load-bearing.  Merely adding an unrelated
`Psi` can restore critical points; the exact companion includes
`Psi'=1` as a hostile control.

## 2. The generic Kummer fibre

The cancellation in the tangential minor is exact:

```text
Jac(P,T)=q phi(T)T/y.                                     (5)
```

Let `z` be transcendental over `k`.  In the function field of `P=z`,

```text
x=(z-Psi(T))/phi(T),
y^q=R(T):=T phi(T)^p/(z-Psi(T))^p.                       (6)
```

The order of `R` at `T=0` is one, so `(6)` is not a proper-power Kummer
equation.  It therefore defines the generic irreducible fibre after passing
to `kbar(z)` and normalization.

If a rational `Q` satisfied `Jac(P,Q)=kappa`, equations `(5)` and `(6)`
would force on that normalized fibre

```text
dQ=kappa omega,             omega:=y dT/(q phi(T)T).      (7)
```

Thus the mate problem has become an exactness problem for one explicit
meromorphic differential.

## 3. Exhaustive valuation ledger

Set

```text
d=deg(phi),       e=deg(Psi),
m=p(e-d)-1,       g_inf=gcd(q,m).                         (8)
```

Since `Psi'=phi H!=0`, one has `e>=d+1` and hence `m>=p-1>0`.
For generic `z`, the roots of `z-Psi` are simple and disjoint from the roots
of `phi`.  The orders of `omega` at every possible nonordinary place are:

| place below the Kummer cover | order of `omega` above it |
|---|---:|
| `T=0` | `0` |
| a root of `phi` | `p-1` |
| a root of `z-Psi` | `q-1-p` |
| each point over infinity | `(m+dq)/g_inf-1` |

Here is the calculation.  At `T=0`, full ramification gives
`ord(y),ord(dT),ord(T)=(1,q-1,q)`.  At a root of `phi`, those three relevant
orders are `p,q-1,q`, giving `p-1`.  At a pole of `R` coming from
`z-Psi=0`, they are `-p,q-1,0`, giving `q-1-p`.  At infinity,

```text
ord(y)=m/g_inf,       ord(T)=-q/g_inf,
ord(phi)=-dq/g_inf,   ord(dT)=-q/g_inf-1,                (9)
```

which gives the last row.

Every displayed order is nonnegative: `p>=2`, `p<q`, and
`m+dq>q>=g_inf`.  At all remaining places the Kummer map is unramified and
`omega` is regular.  Therefore `omega` is a nonzero holomorphic differential
on the smooth projective normalization.

In characteristic zero the differential of a nonconstant rational function
has a pole at every pole of that function.  A rational function with no pole
on a projective curve is constant.  Consequently a nonzero holomorphic
differential cannot equal `dQ`.  Equation `(7)` is impossible, proving the
theorem.

## 4. A four-vertex, genus-ten hostile

Take

```text
p=2,                    q=5,
phi=1+T+T^2,
Psi=T+T^2/2+T^3/3.                                      (10)
```

Then

```text
P=x+x^2y^5+x^3y^5+(1/2)x^4y^10+x^5y^10+(1/3)x^6y^15.  (11)
```

Its Newton hull is

```text
(1,0), (5,10), (6,15), (2,5),                            (12)
```

with normalized area `20`; its gradient Groebner basis is `[1]`.  The
generic fibre is

```text
y^5=T(1+T+T^2)^2/(z-T-T^2/2-T^3/3)^2.                  (13)
```

The seven branch orders on the `T`-line are

```text
1;  2,2;  -2,-2,-2;  1 at infinity.                     (14)
```

All are coprime to five, so Riemann--Hurwitz gives genus `10`.  The divisor
of `omega` has orders

```text
0;  1,1;  2,2,2;  10,                                   (15)
```

whose degree `18` is exactly `2g-2`.  This is not a bounded-degree accident:
the valuation proof closes every rational mate.  The exact computation also
finds rank-one augmented inconsistency for every polynomial response space
of total degree at most `15`.

The example is outside the literal families of THM-2045 and THM-3418, and
its `x`- and `y`-fibre degrees exceed the inherited degree-three walls.  Its
mechanism is closest to the holomorphic-exactness obstruction in THM-2784,
but the present `C_5` Kummer family and two-channel formula are different.

## 5. What this does and does not challenge

This theorem refutes three common construction heuristics:

1. a unit gradient ideal is only the entrance gate;
2. positive Newton area is not the same as independent response channels;
3. high generic-fibre genus can make the mate equation harder, not easier,
   because it creates holomorphic cohomology classes that are not exact.

The internal cyclic cover in `(6)` is not the generic degree or monodromy of
a polynomial map `(P,Q)`--no such `Q` exists.  The theorem is therefore not
a counterexample and supplies no Keller pair.  A live construction must
alter the divisor of the forced response form, for example through two
genuinely independent invariants or mixed terms that create cancellable
poles rather than another holomorphic class.

## 6. Exact verification

Reproduce with

```bash
python3 04-computation/jc_two_channel_kummer_obstruction_thm3552.py
python3 -O 04-computation/jc_two_channel_kummer_obstruction_thm3552.py
```

The ordinary and optimized transcripts agree byte-for-byte.  The companion
checks the gradient ideal and identity `(5)`, Newton support and area,
squarefreeness and generic branch separation, the genus and complete divisor
ledger, bounded response systems through degree `15`, four general valuation
packets, the broken-divisibility hostile, and a positive affine/shear mate.
An independent audit rederived `(3)`--`(9)` and the explicit genus-ten
divisor without importing the companion.

**QED.**
