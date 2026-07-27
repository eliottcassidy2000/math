---
id: THM-1375
title: "Reduced-Jacobian lattice correction, and the degree-three discriminant/reflection law"
status: >
  (0) REFUTED-BY-WITNESS, exact: sigma.F verified Keller (det JG = 2 constant),
  non-injective (the three collision points still share an image), and untwisted
  (target weight vector = source weight vector = (1,-1,-2)).  My own hypothesis,
  killed before it was claimed.
  (I) CORRECTED / PARTLY SUPERSEDED by THM-2598 and MISTAKE-297.  The
  geometric-to-Galois chain and Campbell's Galois exclusion survive, but the
  claimed unconditional Smith self-normalizing rule does not: THM-1365 assumes
  that generic deck transformations extend polynomially.  D4 therefore remains
  a live degree-four monodromy type alongside A4,S4.
  (II) PROVED (modulo the cited Campbell 1973), and it strictly upgrades the
  salvage catalog's grade for the same statement from VERIFIED (0/8316) to PROVED.
  (III) PROVED -- elementary group theory (A_3 regular; A_4 on 4 points transitive
  but not regular).
  (IV) The D_3 split and the Campbell reading are PROVED.  The explicit lambda = -1
  transposition is VERIFIED-EXACT in rational arithmetic.  NOT CLAIMED: that the
  torus reflection coincides with a monodromy reflection -- one is a symmetry of
  the special fibre over the collision point, the other is a loop in the base, and
  their relation is OPEN.
  Nothing here decides JC_2, DC_1 or DC_2, which remain the live frontier.
source: kind-pasteur-2026-07-20-S128c105 (owner: creatively construct a reduced Jacobian Conjecture that does hold; consider tournaments and dihedral groups)
depends_on:
  - THM-1300    # the counterexample, its det, its collision, and its torus weights
  - THM-1310    # fibre geometry: S_3 resolvent, cyclic 3-tournament, Jelonek quartic
  - THM-1330    # Keller monoid; the cusp selection rule
  - THM-127     # Paley D_2p: rotations are automorphisms, reflections anti-automorphisms
related: [THM-1305, THM-1345, THM-2598, MISTAKE-196, MISTAKE-197, MISTAKE-297]
script: 04-computation/reflection_driven_keller_kps_S128c105.py (+ .out), equivariant_reduced_jc_kps_S128c105.py
---

# THM-1375 — the reduced-Jacobian lattice, and the reflection that drives every counterexample

> **CORRECTION (THM-2598 / MISTAKE-297).**  The original §1 treated the
> generic deck group `N_G(H)/H` as a polynomial automorphism group of affine
> space.  That extension across the Jelonek set is not automatic and is an
> explicit hypothesis in THM-1365.  The unconditional self-normalizing rule,
> its “maximal element” claim, and its deletion of `D4` are retracted.  The
> corrected current statements appear below; §§2--4 are unaffected.

## 0. Self-refutation, stated first

THM-1300 records that F is ℂ*-equivariant with weights `(1,-1,-2)` on the source
and `(-2,-1,1)` on the target — the same multiset, related by the transposition
`(1 3)`. I proposed:

> **(RJC-untwisted, FALSE).** F Keller and ℂ*-equivariant with the *same weight
> vector* on source and target ⟹ F invertible.

It is false. Let `σ(a,b,c) = (c,b,a)` and `G := σ ∘ F`. Then `det JG = -det JF = 2`
(still Keller), `G` has the same fibres as `F` (still non-injective), and the
target weight vector of `G` is `(1,-1,-2)` — **untwisted**. Verified exactly.

The reason is structural and worth keeping: *Keller* and *non-injective* are both
invariant under composing with a linear automorphism, and the weight twist is not.
No coordinate-dependent quantity can ever be the operative hypothesis.

This is the same lesson as my LRC `(1,2,3)` episode — an apparent symmetry breaking
that was really an artifact of an ordering convention. I have now made that error
in two unrelated threads, which suggests it is a habit rather than an accident:
**when a structure looks distinguished, first ask what group is allowed to act.**

## 1. The proved lattice, and the missing extension arrow

The geometric chain survives:

```text
injective ==> automorphism ==> proper ==> nodal/empty Jelonek
           ==> abelian complement monodromy ==> Galois ==> automorphism.
```

The cited inputs remain Ax--Grothendieck, Zariski--Lefschetz plus
Deligne--Fulton in the nodal step, and Campbell's Keller-Galois criterion at
the last step.  What does **not** survive is the proposed stronger group
branch.  For monodromy `G` and a point stabilizer `H`, group theory gives

```text
H not self-normalizing ==> N_G(H)/H is a nontrivial generic deck group.
```

THM-1365 can continue from there only under the extra statement

```text
every element of N_G(H)/H extends to a polynomial automorphism of affine space.
```

For a nonproper map that extension across the Jelonek set is not automatic.
Accordingly “every Keller stabilizer is self-normalizing” remains the
Deck-Poverty **conjecture**, not a proved reduced Jacobian conjecture.

THM-2598 now identifies the missing arrow geometrically.  In the canonical
Zariski-main factorization `A^n -> Xbar -> A^n`, a generic deck element always
extends to the finite normalization `Xbar`; it is polynomial exactly when it
preserves the open affine source.  The ramification support in `Xbar` is
intrinsic and invariant.  Thus the `D4` involution can evade the Smith
fixed-point obstruction only if an unramified missing divisor is exchanged
with an included divisor over the same target component.  Normalizer data
alone cannot see this sidecar.

The corrected small-degree table is:

| d | monodromy | `|G|` | `|H|` | `N_G(H)/H` | Galois | proved verdict |
|---|---|---:|---:|---|---|---|
| 2 | `S2` | 2 | 1 | `S2` | yes | excluded by Campbell |
| 3 | `A3` | 3 | 1 | `A3` | yes | excluded by Campbell |
| 3 | `S3=D3` | 6 | 2 | 1 | no | not excluded; §2 forces this type for a counterexample |
| 4 | `C4,V4` | 4 | 1 | `C4,V4` | yes | excluded by Campbell |
| 4 | `D4` | 8 | 2 | `C2` | no | **OPEN: generic involution need not extend polynomially** |
| 4 | `A4,S4` | 12,24 | 3,6 | 1 | no | not excluded by this test |

THM-2598 independently enumerates this table and proves that the `D4`
matching resolvent has type `1+2`.  Thus the honest degree-four live list is
`D4,A4,S4`, while at degree three the only non-Galois transitive type remains
`S3`.

## 2. A status upgrade: the discriminant law is PROVED, not verified

The salvage catalog grades "the fibre cubic's quadratic resolvent is a non-square"
as VERIFIED for F (`0/8316` sample points). It is forced:

> **Proposition.** Every degree-3 Keller counterexample has fibre discriminant a
> **non-square** in `ℂ(F)`; equivalently its monodromy is exactly `S₃`.
>
> *Proof.* Non-injective ⟹ not an automorphism ⟹ (Campbell 1973, contrapositive)
> the cover is not Galois. A degree-3 extension is Galois iff its Galois group is
> `ℤ/3`, iff the discriminant is a square. ∎

So canon's `8316`-point computation is a **corollary**, not evidence for a pattern.

## 3. Why degree 3, and only degree 3, is character-detectable

`A₃` has order 3 = the degree, so a transitive subgroup of `A₃` is *regular*, hence
the cover is Galois. For `d ≥ 4` this fails — verified: at `d = 4` the witness is
`A₄`, at `d = 5` there are two (`D₅` and `A₅`), all transitive, all even, none
regular. Therefore

> **the sign (= discriminant) character detects Galois-ness exactly at `d = 3`.**

With the degree-two Galois impossibility, `d = 3` is simultaneously the minimal
counterexample degree *and* the unique degree at which the reduced JC is decidable
by a single quadratic character. That is why the resolvent `√(-L)` is the right
invariant for F and why no analogous scalar exists one rung up — a concrete reason
the `2-jet d = 4` programme (PROBLEM-ATLAS target 2) will need a different tool.

## 4. The tournament reading — the two flagship threads are one

THM-1310 shows the generic fibre is a **cyclic 3-tournament**. Its symmetry group
is
```
   D_3 = S_3  =  {3 rotations}        u  {3 transpositions}
               =  Aut(3-cycle)         u  ANTI-automorphisms (arc-reversing)
```
which is exactly the split THM-127 proves for the Paley tournament `T_p`: rotations
`v -> v+1` are automorphisms, the reflection `v -> -v` is an anti-automorphism onto
`T_p^op`. Under this dictionary Campbell's theorem reads:

> **Rotation-only Keller covers are trivial.**

and its contrapositive is the structural statement:

> **Every degree-3 Keller counterexample is REFLECTION-DRIVEN: its monodromy must
> contain an arc-reversing element of the fibre tournament.**

The discriminant character is precisely the rotation/reflection indicator — the
same parity bookkeeping as Rédei's sign on the tournament side.

**An explicit reflection.** `λ = -1` in F's own torus fixes the collision image
`(-1/4,0,0)` (target weights `(-2,-1,1)`, so `λ^{-2},λ^{-1},λ` act as `+,-,-` and
the last two coordinates are `0`), and permutes its three preimages as
```
   (0,0,-1/4) -> itself ;   (1,-3/2,13/2) <-> (-1,3/2,13/2)
```
— a transposition, fixing one point and swapping two. Verified in exact rationals.

**Not claimed.** This torus reflection is a symmetry of the *special* fibre over
the collision point; the reflection Campbell forces lives in the *monodromy*, a
loop in the base. Whether they coincide is **open**, and asserting it would be
precisely the kind of coincidence-into-mechanism slide that MISTAKE-197 warns
against.

## 5. Instrument gate (MISTAKE-196)

Any search over untwisted equivariant maps must be able to rediscover `G = σ∘F`.
In the exact parametrisation `G₁ = xA(s,t)`, `G₂ = yB + xzC`, `G₃ = zD + y²E`
(`s = xy`, `t = x²z`, and this is a complete description of the weight pieces, not
an ansatz), `G` needs `deg D = 3`. **My first run used ansatz degree 1 and was
therefore vacuous**; its "no non-linear solutions" line is not reported as a
finding here.

## Named next

- **Do the two reflections coincide?** Compute the monodromy of a loop around the
  Jelonek quartic and compare with the `λ = -1` action. If they agree, the torus
  *is* the monodromy generator and the counterexample is "reflection = torus".
- **Repair the missing extension arrow.** Determine when a generic element of
  `N_G(H)/H` extends across the Jelonek divisor to a polynomial automorphism.
  The first sharp test is the `D4` involution in THM-2598.
- `d = 4`: since no quadratic character detects Galois there, the `2-jet d = 4`
  target must distinguish the `D4` `1+2` matching algebra, the `A4/C3` cubic,
  and the `S4/S3` cubic; a discriminant alone cannot do so.
