---
id: THM-3272
title: "Tame pure-C3 centered-norm packet integral fixed projector"
status: >
  PROVED + VERIFIED-EXACT.  Let O be a henselian DVR with 2 and 3 units,
  and let a monic separable quartic factor integrally as
  f=(X-a)g, with g an irreducible cubic, g mod m=(X-m)^3, and a-m a unit.
  For the universal centered-norm packet nu of THM-3271, the fixed component
  reduces to zero while every moving component reduces to
  2(a-m)^3/27.  Hence its characteristic polynomial reduces to
  Z(Z-2(a-m)^3/27)^3, and the singleton spectral denominator has residue
  -8(a-m)^9/19683, a unit.  The unique simple residue root Hensel-lifts the
  fixed centered norm without a sheet mark, and the spectral projector lies
  in the monogenic order O[nu]; it is exactly the fixed-factor idempotent.
  Thus the centered-norm packet automatically closes its own local integral
  marking gate in the tame pure-C3/unit-cross-resultant scope.  Residue
  characteristic two, centering at three, and nonunit fixed/moving
  cross-resultant are sharp boundaries.  The theorem supplies no affine-open
  incidence, chain-rule cofactor, inverse-different unit, cross-place base
  section, C3/S4 exclusion, or JC(2) theorem.
source: jc-tame-c3-integral-packet-2026-08-02
depends_on:
  - THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector
related:
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3230-marked-c3-trace-centered-norm-and-terminal-prefactor-recovery
script: 04-computation/jc_tame_c3_centered_norm_integral_projector_thm3272.py
output: 05-knowledge/results/jc_tame_c3_centered_norm_integral_projector_thm3272.out
script_sha256: bfe07eb604e1689e660210ac5256e1182ae85ef64e4bdf8aa888139f25af848c
output_sha256: eef61796ad1dbac37836d87262dfa1daf1121a281c5d92347b974902d70baba3
hash_basis: LF-normalized bytes
---

# THM-3272 -- the tame centered-norm packet generates its fixed idempotent integrally

**PROVED + VERIFIED-EXACT.**

## 1. Exact local hypotheses

Let `O` be a henselian discrete valuation ring with fraction field `K`,
maximal ideal `m_O`, and residue field `k`.  Assume

```text
2,3 in O*.                                                  (1)
```

Let

```text
f(X)=(X-a)g(X) in O[X],
a in O,
g monic cubic, separable and irreducible over K.            (2)
```

Assume there is `m_0 in O` such that

```text
bar g(X)=(X-bar m_0)^3 in k[X],
bar a != bar m_0.                                          (3)
```

The last condition is equivalent to

```text
g(a) in O*,                                                (4)
```

because reduction of `g(a)` is `(bar a-bar m_0)^3`.
Thus the fixed and cubic factors are comaximal in `O[X]`.

These are exactly the order-theoretic hypotheses used here.  In the tame
pure-`C3` application, `(3)` says the three moving roots form one totally
ramified cubic orbit and coalesce in the special fibre, while `(4)` is the
unit fixed/moving cross-resultant.  No assertion is made for a pure-`C3`
field before such an integral split normalization has been supplied.

Let

```text
A_O=O[X]/(f),                     xi=X mod f,              (5)
```

and let `nu=N_f(xi)` be the universal centered-norm packet from
[THM-3271](THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector.md).
Because `3` is a unit, its explicit denominator `27` is harmless and

```text
nu in A_O.                                                 (6)
```

## 2. The special-fibre packet has one separated value

Pass to an integral closure in a splitting field of `g`.  Every root
`beta_i` is integral because `g` is monic.  Equation `(3)` has only the root
`bar m_0` over an algebraic closure of `k`, so

```text
bar beta_1=bar beta_2=bar beta_3=bar m_0.                 (7)
```

This splitting base change is used only to calculate symmetric residues;
the resulting characteristic polynomial and denominator lie in `O`, so
unitness descends faithfully.

Let `n_a` be the packet value belonging to the fixed root `a`.  Its
complementary roots are `beta_1,beta_2,beta_3`.  Their centered differences
all reduce to zero, hence

```text
bar n_a=0.                                                (8)
```

Now choose a moving root `beta_i`.  Its three complementary roots reduce to

```text
{bar a,bar m_0,bar m_0}.                                  (9)
```

Their mean is `(bar a+2bar m_0)/3`.  Writing

```text
d=bar a-bar m_0 in k*,                                    (10)
```

the three centered factors are

```text
2d/3, -d/3, -d/3.
```

Therefore every moving packet value has the same nonzero residue

```text
bar n_(beta_i)=2d^3/27 !=0.                               (11)
```

The factor `2` and the denominator `3` in `(11)` are the reason both units
in `(1)` are load-bearing.

## 3. Automatic unit spectral denominator

Let

```text
P_nu(Z)=Norm_(A_K/K)(Z-nu) in O[Z]                        (12)
```

be the packet characteristic polynomial.  Equations `(8)--(11)` give the
complete special fibre

```text
bar P_nu(Z)=Z (Z-2d^3/27)^3.                              (13)
```

The root zero is simple.  Hensel's lemma therefore identifies `n_a` without
first naming the root `a`: it is the unique root of `P_nu` in `O` reducing
to zero.  Put

```text
H(Z)=P_nu(Z)/(Z-n_a),
Delta_nu=H(n_a)=P_nu'(n_a).                               (14)
```

Taking the product of the three fixed/moving gaps in `(11)` gives the exact
residue

```text
bar Delta_nu=(-2d^3/27)^3
            =-8d^9/19683 in k*.                           (15)
```

Thus `Delta_nu` is automatically a unit.  THM-3271's spectral projector

```text
e_fix=H(nu)/Delta_nu                                      (16)
```

therefore belongs to the **monogenic suborder** `O[nu]`, not merely to the
generic quartic algebra.  It satisfies

```text
e_fix=(1,0) under
A_O ~= O x O[X]/(g),                                      (17)
```

so it is exactly the fixed-factor idempotent.  In particular,

```text
O[nu] ~= O x (1-e_fix)O[nu].                              (18)

```

Comaximality `(4)` already put a fixed idempotent in the full graph order
`A_O`.  The new point is stronger and more specific: the single centered-
norm packet itself generates that idempotent integrally, and its scalar
characteristic polynomial selects the fixed centered norm by a simple
residue root.

## 4. Relation to the pure-`C3` cubeclass recovery

At a THM-3230 pure-`C3` place satisfying `(1)--(4)`, the scalar `n_a` is
exactly the trace-centered cubic norm used there.  Equations `(13)--(16)`
recover it, and its factor idempotent, from the unmarked packet.  Its first
nonzero Laurent coefficient can therefore recover the tame cubeclass and
the intrinsic terminal class `Lambda` without choosing a sheet label.

This is a **local order** statement.  Across different places the simple
root of `P_nu` may belong to a different quartic sheet.  The global object is
still the four-component `nu` on the quartic spectral cover, as in
THM-3271; the local idempotents do not become one idempotent over a connected
base.

The theorem also does not construct the primitive-element chain-rule
cofactor from THM-3064.  The derivative star and its inverse-different pole
are separate from the centered norm, and Keller constancy compares their
**product** with the unknown cofactor.  Integral fixed-factor selection is
only the first gate.

## 5. Sharp failure boundaries

### 5.1 Residue characteristic two

Take `O` to be the ring of integers of the unramified quadratic extension
`Q_2(zeta_3)`, put `a=1`, and take `g=X^3-2`.  The cubic is Eisenstein and,
because `zeta_3` is present, its generic extension is cyclic of degree three.
Its roots reduce to zero, and `g(1)=-1` is a unit.  But `(11)` becomes zero
in the residue field because `2=0`.  Hence

```text
bar P_nu(Z)=Z^4,
bar Delta_nu=0.                                           (19)
```

The packet no longer separates the fixed factor integrally.  Thus the
assumption `2 in O*` is sharp for this mechanism.

### 5.2 Centering at residue characteristic three

If `3` is not a unit, the mean in `(5)` and the universal formula with
denominator `27` need not lie in the order.  This theorem makes no wild
characteristic-three replacement claim.

### 5.3 Nonunit cross-resultant

Let `O=C[[t]]`, `a=m_0=0`, and `g=X^3-t`.  The generic cubic is irreducible
and cyclic after adjoining the already available cube roots of unity, but

```text
g(a)=-t notin O*.                                         (20)
```

All four roots reduce to zero.  Both sides of `(11)` collapse, and again
`bar P_nu=Z^4`.  Thus pure-`C3` inertia without fixed/moving separation does
not imply the integral projector theorem.

These hostiles isolate three different failures: coefficient `2`, division
by `3`, and the cross-resultant distance `a-m_0`.

## 6. Exact Kummer control

The family

```text
f=(X-a)((X-m_0)^3-t)                                     (21)
```

has fixed packet value exactly

```text
n_a=t.                                                    (22)
```

At `t=0`, its complete packet polynomial is

```text
P_nu(Z)=Z (Z-2(a-m_0)^3/27)^3,                            (23)
```

and its spectral denominator specializes to `(15)`.  For the concrete
choice `(a,m_0)=(1,0)`, the packet projector reduces in the quartic algebra
to

```text
e_fix=(X^3-t)/(1-t),                                      (24)
```

the ordinary factor idempotent.  Formula `(24)` directly checks that the
packet-generated and factor-generated projectors agree.

## 7. Exact companion and stopping point

Run

```bash
python3 04-computation/jc_tame_c3_centered_norm_integral_projector_thm3272.py
python3 -O 04-computation/jc_tame_c3_centered_norm_integral_projector_thm3272.py
```

Both modes byte-match the stored transcript.  The companion checks the
symbolic means and signs in `(8)--(11)`, the exact residue `(15)`, the full
Kummer identities `(21)--(24)`, the spectral/factor idempotent equality,
`172` separated controls over `F_5,F_7,F_11`, both characteristic-two
collisions, and `23` nonseparated controls.

No conclusion is made outside `(1)--(4)`.  Even inside that scope, the
theorem proves neither that the affine Zariski-main open contains the fixed
boundary trace, nor that the chain-rule cofactor has the required normalized
inverse-different unit.  It gives no cross-place base section, polynomial
Keller realization, pure-`C3` exclusion, `A4/S4` exclusion, `JC(2)`, or
`DC(2)`.

**QED.**
