# The determinant gate, without the Heegner transfer

> **CORRECTED 2026-07-21 by MISTAKE-229 / THM-2055.** The determinant maximum
> is a polygonal support norm `h_conv(+-c_i)(Rd)`, not a single binary
> quadratic form of discriminant `-7`. The Euclidean comparison form has
> discriminant `-4`; gate failure is only uncertified. The Heegner and
> isotropic/rank-versus-anisotropic/Euler claims below are retracted. The
> surviving route is the signed-column normal fan, owner tangent disks, and
> Farey/Klein-sail enumeration with phase-height sidecars.

*boxeph-2026-07-21-S216, corrected after audit. The original reflection made
several unsupported identifications. This file preserves the useful geometric
observation and records exactly where the proposed synthesis fails. The
companion script is
`04-computation/anisotropic_determinant_gate_heegner_residual_boxeph_S216.py`.*

## Status correction

The original S216 conclusion is **REFUTED**. In particular, THM-2053 does not
attach discriminant `-7` to an LRC(14) relation plane, the Paley form supplies
no such attachment, and class number one supplies no compression of the
THM-2052 atlas. The words *anisotropic*, *rank*, *Euler*, and *lonely* referred
to different predicates that had not been connected by a map or theorem.

What survives is narrower and useful:

- `a z_i-b u_i` is the determinant of the parameter direction `(a,b)` and the
  plane column `(u_i,z_i)`;
- the exact THM-2053 sufficient certificate has a quadratic right side and a
  linear left side, so every chosen plane has a finite uncertified tail;
- failure of that certificate is exactly a union of `2n=26` open tangent disks
  in the chosen parameter lattice; and
- reduced-form enumeration independently verifies the stated classical class
  numbers, including `h(-7)=1`.

None of those four statements proves that an uncertified direction is unsafe,
or that a class group organizes the uncertified directions.

## What THM-2053 actually says

Choose a saturated basis `u,z` for a rational two-plane and write a primitive
parameter as `d=(a,b)`. THM-2053 proves the sufficient implication

```text
max_i |a z_i-b u_i| <= (a^2+b^2)/91
    ==> M(a u+b z) >= 1/14.                            (1)
```

The left side of (1) is indeed a maximum of `2x2` determinants. The right side
is the ordinary Euclidean norm in these **chosen basis coordinates**. As a
binary quadratic form it is `[1,0,1]`, with discriminant `-4`, not `-7`.
Changing the saturated basis by `GL(2,Z)` changes its coefficient matrix in
the new coordinates. Thus (1) does not define a plane-intrinsic discriminant
that could be read from the number `14=2*7`.

For `c_i=(u_i,z_i)` and `Jc_i=(z_i,-u_i)`, failure of (1) is exactly

```text
d in union over i and sigma in {+1,-1} of
  {x : ||x-(91 sigma/2)Jc_i||^2 < (91^2/4)||c_i||^2}. (2)
```

There are 26 disks when `n=13`; each is translated away from the origin and
has the origin on its boundary. The union in (2), not one quadratic form or
one centered ball, is the exact carrier. The round bound
`||d|| < 91 max_i||c_i||` is only an envelope for it.

The certificate is one-way. A point of (2) is merely **not certified by this
estimate**. It may still be lonely for a different reason.

## Why the Heegner and Paley transfers fail

Classical reduction theory does give

```text
h(-3)=h(-4)=h(-7)=h(-8)=h(-11)=h(-19)=h(-43)=h(-67)=h(-163)=1.
```

But `h(-7)=1` classifies proper equivalence classes of primitive positive
definite binary quadratic forms of discriminant `-7`. The coordinate norm in
(1) has discriminant `-4`, while (2) is a union of translated disks determined
by all thirteen columns. Neither the column data nor the translations are
encoded by a form class. Consequently class number one does not identify,
rigidify, or reduce the actual THM-2053 residual atlas.

For odd primes `p=3 mod 4`, the nonprincipal Paley eigenvalues satisfy a
quadratic polynomial whose homogeneous form has discriminant `-p`. That is a
classical spectral identity in a different object. S216 supplied no map from
Paley eigenvectors, characters, or form classes to a saturated LRC relation
plane preserving the loneliness predicate or the disk union (2). Equality of
a displayed discriminant would not by itself provide such a map; here even
the THM-2053 coordinate discriminant is `-4`.

The local-symbol argument also failed at its first calculation. If
`D=-p`, then

```text
(D/p)=(-p/p)=0,
```

because `p` divides `D`. The separate value `(-1/p)` is `-1` for
`p=3 mod 4` and `+1` for `p=1 mod 4`, but it is not `(D/p)`. It cannot be used
as written to classify the discriminant-`-p` form over `Q_p`.

## The rank/Euler identification has no established implication

The earlier reflection asserted

```text
isotropic direction => resonance => rank jump,
anisotropic direction => no resonance => lonely => Euler survivor.
```

No definitions or maps made those arrows valid. A null vector of a quadratic
form is not automatically a new bounded integer relation on thirteen speeds.
Conversely, anisotropy of a form only forbids the relevant form from
representing zero; it says nothing by itself about the maximin
`M(v)=max_t min_i ||t v_i||`. In particular, **anisotropic does not imply
lonely**. The rank and Euler alternatives remain distinct proposed proof
mechanisms, not synonyms for local isotropy and anisotropy.

## Repaired research question

The determinant view is still a productive coordinate description. For each
actual THM-2052 two-plane, one can ask whether the 26 tangent disks have
overlap, symmetry, boundary, or residue patterns that permit a smaller exact
discharge. Any proposed arithmetic organizer now has to provide:

1. a map from its objects to the chosen plane and its thirteen columns;
2. invariance or controlled transformation under a saturated basis change;
3. preservation of the exact disk predicate (2), not just a discriminant; and
4. an implication from the organized cases to an actual `1/14` loneliness
   witness.

Until such a map is supplied, class-number computations and Paley spectra are
comparative prompts only. The current concrete LRC(14) task remains the one in
THM-2053: construct or compress the finite plane atlas and exactly discharge
the primitive lattice points in its tangent-disk unions.

Links: HYP-8865 (refuted synthesis), THM-2053, THM-2052, HYP-8841,
[[each-prime-is-its-paley-tournament-a-periodic-table-of-2-3-5-7-11-for-lrc14-boxeph-S215]].
