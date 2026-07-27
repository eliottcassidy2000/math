---
id: HYP-9031
title: "The D5 dictionary: LRC word-current and JC flux as explicit H^1 classes -- anti-parallel quantifiers, twin orientation sidecars"
status: >
  OPEN (a typed cross-domain GRAMMAR with explicit classes on both
  sides, three falsifiable predictions, and named hostile controls --
  NOT a reduction; per the bridge discipline this file names map,
  preserved predicate, loss, sidecar, and cheapest test for each entry).
source: opus-2026-07-26 (successor to the shape-confirmation in
  07-reflections/spectral-vs-geometric-rank-why-LRC-stays-open-and-JC-fell-opus-20260726.md)
related:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2337 / THM-2356 / THM-2334 (word-current, Bockstein, 169 twists)
  - THM-2389 / THM-2406 / THM-2463 / THM-2468 (flux, pole-descent, square lifts)
  - HYP-9030-keller-degree-semigroup
---

# HYP-9031 -- the D5 dictionary, written out

Agent D5 (opus reflection, 2026-07-26) predicted that the two hardest
live frontiers are the *same* `H^1`-shape with different coefficient
arithmetic. This file writes the two classes explicitly and turns the
shape into a typed dictionary. The headline sharpening:

> **JC(2) and LRC(14) are the two DIRECTIONS of one cohomological
> question.** JC(2) asserts a universal VANISHING (every planar Keller
> map has trivial monodromy class: injectivity); LRC(14) asserts a
> universal NONVANISHING (every scalar row has nonzero word-current).
> The sporadic core (THM-1300/2473) is precisely a *computed nontrivial
> class* on the geometric side -- which is why that side FELL while the
> spectral side stays rigid: exhibiting one nontrivial class is
> geometry; proving all classes nonzero is Diophantine.

## The geometric class (JC side), now fully explicit

For the sporadic Keller map `F` (THM-2473): let `L = 27a^2c^2 - 18abc +
16a + b^3c - b^2` and `U = C^3 \ Z(L)`. Over `U`, `F` restricts to a
degree-3 finite etale cover, classified by

```text
rho_JC in H^1(U; S3)   (nonabelian; = Mon(F) = S3, surjective, PROVED)
```

- **Distinguished cycle:** a loop around the fold locus; its image is
  the transposition swapping the folded pair (`4ax^2+1` on the x-axis
  line, disc `-16a`).
- **Abelianization:** `sgn(rho_JC) = [-L] in H^1(U; Z/2) = C(a,b,c)^*/squares`
  (from `disc_x E = -4(27ac^2-9bc+8)^2 L`): the Kummer class of `-L`.
- **The JC(2) program in this language:** every degree-18/22 plane
  closure (THM-2463/2468/2469/2470/2472 ...) restores a discarded
  physical relation as a CYCLIC KUMMER LIFT (`Y^2 = 1/p`, `y^3 = C/p`)
  and derives genus >= 2: i.e. proves the would-be planar class cannot
  exist. Pole-descent = computing the class's local restriction;
  square/cube lifts = un-abelianizing it.

## The spectral class (LRC side)

On the oriented `C_13` root chart of a canonical blocker word
(THM-2305/2334/2337): the word-current assigns to the deepest
primitive-91 `c_3`-edge cycle the coefficient sum `C(q) = sum_z A(q,z)`
valued in the mod-91 character group; CRT splits every class

```text
rho_LRC in H^1(chart; Z/91) = H^1(.; Z/7) (+) H^1(.; Z/13),
```

and the crux (kps-S131b; the `7 (x) 13` tensor) is that BOTH components
of a mixed unit character must be certified jointly nonzero.

- **Distinguished cycle:** the primitive-91 deepest `c_3` edge loop
  (THM-2293-type edges; the trichotomy of MSG-2153, under audit this
  session, would make the 91-unit *marked* edge exist unconditionally).
- **Sidecar:** the first Bockstein `H^1 -> H^2` (THM-2337/2356): the
  refined-to-coarse gluing lives one degree up; the 169-twist variance
  (THM-2334 eq 42) is the concrete nonvanishing test.

## The dictionary (typed, with losses)

| slot | JC / geometric | LRC / spectral | loss if identified |
|---|---|---|---|
| base | `U = C^3 \ Z(L)` | oriented `C_13` chart minus danger combs | no common site; GRAMMAR only |
| coefficients | `S3` monodromy (nonabelian) | `mu_7 (x) mu_13` characters (abelian) | rank-type: collision vs no-cancellation |
| class | `rho_JC` (computed, nontrivial) | word-current `C(q)` (nonvanishing OPEN) | anti-parallel quantifiers |
| abelianized invariant | Kummer class `[-L]` | mixed unit character mod 91 | `Z/2` vs `Z/91`: no map |
| distinguished cycle | loop at the fold (transposition) | primitive-91 `c_3` edge loop | none (both honest cycles) |
| degeneration locus | `Z(L)`: even pair-escape; empty-fiber curve (F not surjective) | comb zero sets (`7`-line invisibility, MSG-2152 boundary) | escape vs cancellation |
| orientation sidecar | `F o sigma = tau o F`, `sigma != tau`: NO depth-2 involution (THM-2473(6)) | tournament directed only up to global converse; orientation representative consumed (THM-2504) | **twin phenomena, proved independently 2026-07-26/27 on both sides** |
| one degree up | disc-square surface (projection collision, order-2) | first Bockstein mod 91 (THM-2337/2356) | different connecting maps |

## Three falsifiable predictions

1. **(spectral)** If the MSG-2153 trichotomy survives audit, the
   H^1-existence half of the LRC crux is DONE at the marked-edge level,
   and everything remaining is the `H^2`/incidence half (Bockstein,
   address participation, THM-2321 (51)). Test: the audit (this
   session) + the 169-twist variance computation: a nonconstant twist
   bank = first nonzero-target survival.
2. **(geometric)** Each of the three remaining degree-22 support-two
   planes -- `(B,C)`, `(B,E)`, `(C,E)` -- closes by restoring one
   discarded ROOT coordinate and forcing genus on the resulting curve:
   either a cyclic Kummer lift (`Y^2 = 1/p` as in THM-2463/2468, the
   un-abelianized class) or a root-ratio normalization (`r = E/(Dy)`
   style as in THM-2469/2470/2472, the class's local restriction).
   Seven of ten planes already followed one of these two flavors; the
   prediction is that no third mechanism will be needed, because both
   flavors are the same act -- putting the discarded `H^1` datum back.
3. **(orientation)** Any attempt to canonicalize the LRC endpoint
   current without carrying an explicit orientation representative
   (chart + word + converse-gauge) will fail the way THM-2504 shows;
   dually, any Keller-tower argument that assumes the collision
   involution acts at depth >= 2 is FALSE (THM-2473: no sign-diagonal
   `rho` with `F o rho = sigma o F`). Both sides force the same sidecar
   type: a chosen representative of an orientation 2-torsor per level.

## Hostile controls

- MSG-2152's sharp signed-step boundary (spectral): disjointness +
  endpoint Prony CANNOT force `gcd(m,91)=1`; the 13-line can hide.
  Any claimed "class nonvanishing" proof must exhibit the extra
  indicator/root-sheet structure it consumes.
- THM-2473's empty-fiber curve (geometric): over the omitted rational
  curve the cover does not exist at all; "the class" is only defined on
  `U`, and arguments must not quantify over all of `C^3`.
- MISTAKE-226/235 (shared kernels are not bridges): this file is a
  grammar with typed losses, NOT an LRC-to-JC reduction; no theorem may
  cite it as a dependency.

## Cheapest decisive next steps

(i) finish the MSG-2153 audit and the 169-twist variance run (both in
flight this session); (ii) on the JC side, check prediction 2 against
the next plane closure (THM-2469/2472 just landed: verify the lift type
used); (iii) formalize the orientation-sidecar parallel: state the
"orientation 2-torsor per level" object once, in a form both THM-2504's
chart gauge and THM-2473's sigma/tau bookkeeping instantiate.
