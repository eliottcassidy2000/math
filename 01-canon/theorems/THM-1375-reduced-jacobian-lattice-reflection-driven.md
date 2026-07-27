---
id: THM-1375
title: "The corrected reduced-Jacobian lattice and the degree-three reflection law"
status: >
  (0) REFUTED-BY-WITNESS, exact: sigma.F verified Keller (det JG = 2 constant),
  non-injective (the three collision points still share an image), and untwisted
  (target weight vector = source weight vector = (1,-1,-2)).  My own hypothesis,
  killed before it was claimed.
  (I) CORRECTED BY MISTAKE-297.  The geometric chain through Galois is a
  synthesis of proved/classical implications, but the former terminal Smith
  self-normalizing-stabilizer rung is RETRACTED: deck transformations live on
  the punctured source and were never proved to extend freely across the
  Jelonek divisor.  D4 is therefore open at degree four, alongside A4/S4.
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
related: [THM-1305, THM-1345, MISTAKE-196, MISTAKE-197, MISTAKE-297]
script: 04-computation/reflection_driven_keller_kps_S128c105.py (+ .out), equivariant_reduced_jc_kps_S128c105.py
---

# THM-1375 — corrected reduced-Jacobian lattice and the degree-three reflection law

> **MISTAKE-297 CORRECTION.** The original Section 1 applied Smith theory to
> deck transformations on the punctured source as if they extended freely to
> affine space. That implication is retracted. The geometric chain through
> Galois and every degree-three reflection statement survive; the alleged
> self-normalizing-stabilizer rung and its exclusion of `D4` do not.

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

## 1. The corrected lattice and the retracted Smith rung

The proved chain stops at Campbell's Galois criterion:

```text
injective -> proper -> nodal Jelonek -> abelian monodromy -> Galois -> d=1.
```

The last implication is Campbell. A transitive abelian monodromy action is
regular, so the original extension is Galois. What does **not** follow is
that every non-self-normalizing point stabilizer is excluded. The historical
discussion and table below are retained to show the failed route; wherever
they say that Smith excludes `D4` or proves `(RJC-max)`, read that assertion as
**RETRACTED BY MISTAKE-297**.

The repo holds two catalogs of true restricted JCs (the salvage catalog,
kind-pasteur-S128c99 §2; the true-Jacobian ledger, boxeph-S144) but never orders
them. Ordered by *weakness of hypothesis* (weaker hypothesis = stronger theorem):

```
injective ==> proper ==> nodal Jelonek ==> abelian monodromy ==> GALOIS ==> d=1
(Ax-Groth.)  (classical)  (Zariski-Lefschetz + Deligne-Fulton)  (Campbell)
```

The proved branch is:

- **Geometric branch.** Empty Jelonek set (= proper) ⟹ vacuously nodal ⟹ by
  Zariski–Lefschetz the generic plane section controls `π₁(ℂⁿ∖A)`, and by
  Deligne–Fulton a nodal plane curve has abelian complement group ⟹ the monodromy
  is abelian ⟹ transitive abelian ⟹ regular ⟹ the cover is Galois.
- **Galois terminal.** Transitive abelian monodromy is regular, so the original
  function-field extension is Galois; Campbell then forces invertibility.

The former extra `(RJC-max)` rung is retracted.  Nontrivial deck
transformations act on the open source above the target complement, not on
contractible affine space.  Without a regular extension across the Jelonek
divisor, Smith theory supplies no self-normalizing-stabilizer criterion.
Degree two is still impossible because its extension is automatically Galois.

**F violates every rung, maximally.** Its Jelonek set is Zariski's three-cuspidal
quartic (THM-1330) — not merely non-nodal but the worst rational quartic; its
monodromy is `S₃`, the largest possible at degree 3; its point stabiliser
`⟨transposition⟩` is self-normalising in `S₃`.

### The ordering is strict, and it reproduces canon's table

Computed over the transitive subgroups of `S_d`, `d ≤ 5` (normalisers by direct
enumeration):

| d | monodromy | \|G\| | \|stab\| | self-norm | Galois | verdict |
|---|---|---|---|---|---|---|
| 2 | `S₂` | 2 | 1 | no | yes | excluded (Campbell) |
| 3 | `A₃` | 3 | 1 | no | yes | excluded (Campbell) |
| 3 | **`S₃ = D₃`** | 6 | 2 | yes | no | **ALLOWED** |
| 4 | `Z₄`, `V₄` | 4 | 1 | no | yes | excluded (Campbell) |
| 4 | **`D₄`** | 8 | 2 | **no** | **no** | **OPEN** |
| 4 | `A₄`, `S₄` | 12, 24 | 3, 6 | yes | no | ALLOWED |

The group census itself is exact, but its old verdict column was not: `D₄` is
not killed by Campbell and the Smith extension step was unavailable.  The
honest wild degree-four list is therefore `D4,A4,S4`.

Note also that at `d = 3` the only allowed monodromy is `S₃` — which is what §2
proves must occur.

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

With Campbell's degree-two Galois impossibility, `d = 3` is simultaneously the minimal
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
- **D4 pole tax.** A `D4` cover has a deck involution on the punctured source.
  Prove and exploit the fact that it cannot extend regularly across the Jelonek
  divisor; its unavoidable rational pole is the corrected Smith-side invariant.
- `d = 4`: since no quadratic character detects Galois there, the `2-jet d = 4`
  target needs a genuine `D₄/A₄/S₄` lane invariant, not a discriminant alone.
