---
id: HYP-8846
title: "LRC14 pointed transport on the rank-eleven relation planes"
status: >
  OPEN finite completion program after a proved large-direction theorem.
  THM-2052 reduces every hypothetical counterexample to a finite atlas of
  rational planes. THM-2053 now proves that every two-dimensional atlas cell
  has torus margin at least 1/13. Its exact safe gate is
  max_i|a z_i-b u_i|<=(a^2+b^2)/91; the round corollary is ||(a,b)||>=91L.
  The missing work is controlled basis/atlas compression and exact discharge
  of the finite anisotropic residuals, not an infinitary transport theorem.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-2050
  - THM-2052
  - THM-2053
  - THM-2054
  - THM-2055
  - HYP-2108
  - HYP-1842
  - HYP-1977
  - HYP-2647
  - HYP-2896
  - HYP-2986
  - HYP-3267
  - HYP-4336
  - HYP-4346
  - HYP-8841
  - HYP-8845
  - HYP-8850
  - HYP-8860
  - HYP-8871
---

# HYP-8846 -- pointed transport, not an unpointed plane escape

The original proposal in this file had a real quantifier gap:

```text
known:  target v lies in a rational plane ker(W);
known:  rank rigidity supplies some escaping direction;
missing: "some direction is safe" does not imply the specified v is safe.
```

THM-2053 repairs the quantifier on all sufficiently long directions. If
`u,z` is a saturated integer basis of the plane, put

```text
L=max_i sqrt(u_i^2+z_i^2),
v(a,b)_i=a u_i+b z_i.
```

For every primitive `(a,b)`, the corresponding subcircle is
`1/(2sqrt(a^2+b^2))`-dense in the parameter two-torus. The phase-height
objective is `L`-Lipschitz. Keeping the actual normal displacement gives the
sharper error

```text
E(a,b)=max_i |a z_i-b u_i|/(2(a^2+b^2)).
```

A planar-direction argument produces a
full-support projection with a repeated absolute speed, so settled LRC for at
most twelve distinct speeds gives two-torus margin at least `1/13`. Therefore

```text
max_i |a z_i-b u_i| <= (a^2+b^2)/91
                         ==> M(v(a,b))>=1/14.           (1)
```

In particular, `sqrt(a^2+b^2)>=91L` implies (1). This is pointed: it applies
to the specified `(a,b)`, not merely to a different direction in the same
plane.

## The revised terminal

The rank-eleven program now has two finite terminals rather than one:

```text
rank 12       -> maximal-minor box (THM-2052),
rank 11/plane -> primitive parameter disk ||(a,b)||<91L (THM-2053).
```

The intended per-target routing remains useful inside that disk:

1. some peel has the HYP-2108 endpoint functional `P_w>=0`;
2. active owner data supplies a bounded relation outside `W`; or
3. an exact resolved phase/topological certificate discharges the row.

There is no longer a need to postulate an infinite wall-crossing transport to
close large directions. Wall transport is a finite-disk accelerator.

## Concrete completion algorithm

For each rank-eleven bounded triple code from THM-2052:

1. compute a saturated kernel basis by Smith/Hermite normal form;
2. lattice-reduce the basis to minimize the determinant gate and its round
   envelope `91L`;
3. enumerate only primitive parameter pairs violating (1): equivalently,
   lattice points in the `26` open disks centered at
   `+-(91/2)(z_i,-u_i)` with radii `(91/2)||(u_i,z_i)||`; use the round disk
   only as a coarse outer bound;
4. reject zero, repeated, nonpositive, non-distinct, or imprimitive speed rows;
5. run the exact pair-sum phase-height test;
6. retain `P_w`, endpoint-owner, peel-tax, and relation-rank labels for any
   survivors, then quotient only by transformations preserving those fields.

Each residual circle has the origin on its boundary. This tangent-disk picture
is the preferred parameter carrier: intersect it with the star's positivity
cone and primitive lattice before any phase-height work.

This is structural compression, not the first proof of finiteness: THM-763
already supplies the global ceiling `sum v_i<=91^12`. The potential gain is
that the plane gate discards whole projective tails before raw speed
enumeration and retains the geometry needed for resonance fans.

The finite atlas is astronomically large if generated as all possible bounded
triple rows. The theorem-facing next step is therefore an atlas-compression
lemma: enumerate rank-two **column configurations/matroids**, not raw relation
matrices, and bound a reduced-basis `L` directly from the sparse triple code.

## What the unrelated repo work contributes

- HYP-4342 supplied the `(1,N)` Lipschitz/net mechanism; THM-2053 extends it to
  every primitive direction.
- HYP-4346 supplied the rank-two algebra but also exposed the wrong-quantifier
  trap. It is now optional acceleration, not the bridge.
- HYP-2896 is the exact one-tail model of the finite-disk fan. On
  `span((1,...,13),e_12)`, target rows split by the resonance walls `12|w` and
  `14|w`; the cells carry `q=12`, `q=14`, or the affine binding phase
  `(35m+2)/(84m+5)`. In the basis `u=(1,...,13),z=e_12`, writing `w=12+b`,
  THM-2053's determinant is exactly `13|b|` for `b!=0`, so it certifies every
  integer `|b|>=1183`. The HYP-2896 fan discharges the finite residue
  symbolically. That geometry-then-arithmetic split is the desired output
  format for a general plane.
- HYP-2986 supplies the faithful three-state terminal: open tope, boundary
  cocircuit, or forbidden wall packet.
- HYP-2647 supplies the addressed wall-transport matrix for moving between
  neighboring parameter cells. Its scalar signed total is not enough.
- HYP-1842 suggests private endpoint pivots as relation-rank suppliers, but a
  binary private pivot does not yet imply a rational relation outside `W`.
- HYP-1977 proves that even observer-pointed tournament classes can mix safe
  and unsafe cells; endpoint lengths/owners must remain in the quotient.
- HYP-3267's contact holonomy is a useful connection coordinate, but its
  empty/full collision prevents it from replacing the endpoint cell.
- HYP-8845/HYP-8850 parity duplicates a first survivor into its mirror. It
  cannot create that first survivor and therefore belongs after the pointed
  gate.
- THM-2054 is the incoming analytic complement: relative Fejer decorrelation
  along a character line could certify whole off-resonance cells before exact
  lattice-point enumeration. It is currently claimed/in progress and is not
  used in THM-2053.
- HYP-8860's Paley-prime table usefully assigns roles to moduli `3,7,11`
  (resonance atom, period-14 apex, rank scale). It is a modulus-selection lens,
  not a carrier: Paley orientation discards the signed coefficients and
  endpoint owners needed to decide a tangent-disk point.
- THM-2055 replaces the raw tangent-disk union by the signed column polygon's
  normal fan. Only hull vertices own the determinant maximum; each owner cone
  has one disk and an owner-local radius. HYP-8871 proposes ordinary
  Stern--Brocot/Klein-sail traversal of those rational sectors. MISTAKE-225
  records why Heegner form classes cannot replace this polygonal carrier.

## Tournament analysis

Vertices are proof carriers, not runners or raw relation rows:

```text
torus_geodesic_disk
rank12_maximal_minor
resolved_phase_height
endpoint_owner_tope_cocircuit
addressed_wall_transport
private_pivot_rank_supplier
contact_holonomy
unpointed_rank_rigidity
raw_relation_count
```

Pairwise observable: `(pointed target retention, finite terminal, exact LRC
predicate, owner/magnitude retention, cost)`. The gauge prefers a carrier only
when it preserves the specified row and supplies either a finite domain or an
exact exit. The tie path begins

```text
torus_geodesic_disk
> rank12_maximal_minor
> resolved_phase_height
> endpoint_owner_tope_cocircuit
> addressed_wall_transport
```

and places unpointed rank rigidity and raw relation count last. This quotient
preserves the target-row LRC predicate and destroys raw relation-matrix identity;
the destroyed identity must be recoverable through the saturated kernel basis.

## Honest remaining statement

LRC(14) is still open. THM-2053 proves that the two-dimensional branch is
finite in parameter space, not that the finite disks are empty of
counterexamples. The next decisive target is:

> Bound and classify reduced saturated bases of THM-2052's bounded
> support-three rank-eleven codes strongly enough that every residual to (1)
> is either covered by the existing exact finite window or splits into
> finitely many HYP-2896-style resonance fans with explicit phase-height
> certificates.
