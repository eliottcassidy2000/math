---
id: THM-2686
title: "Coprime Galois invariant Kummer vanishing and the A4 standard-plane gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let E/K be a
  finite Galois extension of degree prime to a prime ell over the function
  field of complex affine n-space, let R be the affine-target normalization
  in E, and put U=R_reg.  Then H^1_et(U,mu_ell)^Gal(E/K)=0.  An invariant
  nontrivial ell-cover descends to a direct-product extension over the target;
  its cyclic ell-subfield is unramified at every target divisor and is
  therefore trivial by purity and affine-space simple connectedness.  For
  ell=2 and Gal(E/K)=C3, finite mod-2 Kummer cohomology is consequently a sum
  of standard two-dimensional F2[C3]-planes.  For a quartic A4 Keller candidate,
  THM-2655's module test therefore reduces to H^1_et(U,mu_2)!=0, equivalently
  a nonzero unit squareclass or Cl(R)[2] class.  A dense regular chart A^n or
  G_m x A^(n-1) excludes the candidate because it cannot carry the restricted
  connected V4 torsor.  No general A4, quartic, or Jacobian exclusion follows.
source: root-2026-07-28-a4-cyclic-resolvent-mod2-gate
depends_on:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
script: 04-computation/jacobian_a4_cyclic_resolvent_mod2_gate.py
output: 05-knowledge/results/jacobian_a4_cyclic_resolvent_mod2_gate.out
script_sha256: ca8aa0cf99ac2f4fd85606c78c8083e7c076c8bc744be2a3b4b3ffb4d4fc21ac
output_sha256: 4d7ef2c104f6dcd069bb474bb9e0ebcbaa355844cf046e324dbfd84bf408c8a0
hash_basis: LF-normalized bytes
---

# THM-2686 -- coprime Galois normalizations have no invariant Kummer character

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The general mechanism is coprime invariant Kummer vanishing: a prime-order
character on the regular locus of a Galois normalization cannot be invariant
when the Galois degree is prime to that character order.  The descended cyclic
cover would be unramified over every divisor of affine space, contradicting
purity and simple connectedness.

THM-2655 forces a quartic `A4` Keller map to carry a connected `V4` torsor on
the regular locus of its cyclic cubic resolvent normalization.  Here the
general theorem applies with `(ell,Q)=(2,C3)`.  Since `C3` acts on the three
nonzero `V4` characters, the absence of invariant quadratic characters makes
every nonzero mod-two character generate the required standard plane.

This does not exclude `A4` by itself.  It turns the equivariant module test
into a scalar nonvanishing test and makes simple regular charts decisive.

## 1. Setup

Let

```text
X=A^n_C with n>=1,          K=C(X),
ell a prime,                E/K finite Galois,
Q=Gal(E/K),                 gcd(ell,|Q|)=1.                (1)
```

Let `R` be the normalization of `X` in `E`, and put

```text
U=R_reg.                                                   (2)
```

The `Q`-action preserves `U`.  We claim

```text
H^1_et(U,mu_ell)^Q=0.                                    (3)
```

All cohomology and fundamental groups below are geometric over `C`.

## 2. An invariant Kummer cover descends as a direct product

Suppose that a nonzero class

```text
alpha in H^1_et(U,mu_ell)^Q                              (4)
```

exists.  It gives a connected etale `C_ell`-cover `T->U`: because `ell` is
prime, a nonzero torsor over connected `U` cannot split into smaller orbits.
Let `L/E` be its function-field extension.  Equivalently, restriction to the
generic point is injective,

```text
H^1_et(U,mu_ell) -> H^1(E,mu_ell)=E^*/E^(*)ell,          (5)
```

Thus `alpha` gives a nonzero `Q`-invariant field Kummer class.  Inflation--
restriction for `E/K`, together with

```text
H^1(Q,mu_ell)=H^2(Q,mu_ell)=0,
```

shows that this class descends uniquely to a nonzero class in
`K^*/K^(*)ell`; the group cohomology vanishes because multiplication by
`|Q|` is invertible on the `ell`-torsion module.  Let `D/K` be the corresponding
cyclic degree-`ell` extension.  Then `E intersect D=K`, `L=ED`, and

```text
Gal(L/K)=C_ell x Q.                                      (6)
```

Equivalently, invariance identifies all conjugates of `T`, and the vanishing of
`H^2(Q,mu_ell)` linearizes the cover.  The generic Kummer route makes explicit
that no descent across the possibly branched geometric quotient `U/Q` is
assumed.  Coprimality in (1) is load-bearing.

## 3. The cyclic quotient is everywhere unramified

Every height-one point of the normal variety `R` is regular and therefore
lies in `U`.  Since `T->U` is etale,

```text
L/E is unramified at every height-one prime of R.          (8)
```

Let `v` be a target divisor of `X`, choose a prime of `L` above it, and let
`I<=C_ell x Q` be its inertia group.  In the tower `L/E/K`, (8) says

```text
I intersect C_ell=1.                                     (9)
```

Projection `I->Q` is therefore injective, so `|I|` divides `|Q|`.  Let

```text
D=L^Q                                                    (10)
```

be the cyclic degree-`ell` subfield.  The inertia of `D/K` is the image of `I`
in `(C_ell x Q)/Q=C_ell`; its order divides both `|I|` and `ell`, hence is
trivial by (1).  Therefore the
normalization of `X` in `D` is unramified in codimension one.

Affine space is excellent, so its normalization in `D` is finite.
Zariski--Nagata purity for a normal source over the regular scheme `X`
promotes codimension-one unramifiedness to a finite etale cover of all affine
space.  But

```text
pi_1^et(A^n_C)=1.                                        (11)
```

The field `D` must therefore equal `K`, contradicting its degree `ell`.
This proves (3).

## 4. Every nonzero character generates the standard plane

Now specialize to `(ell,Q)=(2,C3)`.
The finite-coefficient etale cohomology of `U` is finite dimensional.  Since
`2` does not divide `|C3|`, Maschke semisimplicity applies to

```text
H=H^1_et(U,mu_2)                                         (12)
```

as an `F2[C3]`-module.  Over `F2`, the only simple `C3`-modules are the
trivial line and the standard irreducible plane on which `C3` cycles the
three nonzero vectors.  Equation (3) eliminates every trivial summand.  Thus

```text
H is a direct sum of standard planes;                    (13)
dim_F2 H is even;
H!=0 iff H contains the standard character plane W.      (14)
```

The singular locus of normal `Spec R` has codimension at least two (and is
empty in dimension one), so `Gamma(U,O_U)=R` and `Pic(U)=Cl(R)`.  The Kummer
sequence therefore gives the equivariant exact sequence

```text
0 -> R^*/R^(*)2 -> H -> Cl(R)[2] -> 0.                   (15)
```

Semisimplicity splits (15) as a `C3`-module.  Hence (14) is equivalent to

```text
R^*/R^(*)2 !=0              or              Cl(R)[2]!=0. (16)
```

Any nonzero class in either term has a three-element nonzero orbit and spans
a standard plane.  No separate action-matrix computation is needed in the
cyclic branch.

## 5. Consequence for quartic `A4`

For a quartic `A4` Keller map, THM-2655 constructs a connected `V4` torsor on
`U`; dualizing injects its standard character plane into (12).  Equations
(14)--(16) therefore give the exact necessary test

```text
H^1_et(U,mu_2)!=0.                                       (17)
```

Conversely, nonvanishing in (17) supplies the **abstract equivariant Kummer
carrier** required by THM-2655, because every nonzero summand is already the
standard plane.  It does not construct a Keller map or show that the carrier
is the actual quartic normalization.

There is a useful chart exclusion.  Let `V` be a nonempty open subset of `U`.
The actual `V4` torsor supplied by a Keller map has integral total space, so
its restriction over `V` remains connected.  If

```text
V is isomorphic to A^n
or V is isomorphic to G_m x A^(n-1),                     (18)
```

this is impossible: the first chart has trivial etale fundamental group and
the second has procyclic fundamental group, neither of which has a `V4`
quotient.  Thus an actual dense chart of either form excludes the `A4`
resolvent branch.

## 6. Controls and scope

Run

```bash
python3 04-computation/jacobian_a4_cyclic_resolvent_mod2_gate.py
python3 -O 04-computation/jacobian_a4_cyclic_resolvent_mod2_gate.py
```

Both modes byte-match

```text
05-knowledge/results/jacobian_a4_cyclic_resolvent_mod2_gate.out.
```

The companion checks `19` exact finite group/module identities.  It realizes
the standard plane as `F2^2` with an order-three matrix, reconstructs
`V4 semidirect C3=A4`, checks the relevant `C6` subgroup lattice, and uses two
geometric controls.  The smooth cyclic quotient `(t,u,v)->(t^3,u,v)` has
normalization `A^3` and zero Kummer rank but is ramified, hence not Keller.
The toric `d^2=abc` hostile has `Cl[2]=(F2)^2` with `C3` cycling its three
nonzero classes, showing the standard-plane carrier can occur on a singular
normalization.  Purity, (11), and the field descent in Sections 2--3 are
mathematical inputs rather than executable claims.

An independent hostile audit checked the generic-point injection and
inflation--restriction descent, every inertia subgroup and quotient, the
excellence/finite-normalization input to purity, the `Pic(U)=Cl(R)` passage,
the semisimple module criterion, and connectedness after chart restriction.
It also byte-matched the normal, optimized, and stored executable outputs and
their declared hashes.  The same coprime-order audit yields the general
`(ell,Q)` statement in Sections 1--3; the executable controls specialize to
`(2,C3)`.

The theorem is a necessary resolvent obstruction and an abstract-carrier
criterion.  It does not exclude every cyclic cubic normalization, the general
`A4` branch, degree four, `JC(2)`, general JC, or DC(2).

QED.
