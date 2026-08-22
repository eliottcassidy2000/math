---
id: HYP-9031
title: "The D5 dictionary: LRC word-current and JC flux as explicit H^1 classes -- anti-parallel quantifiers, twin orientation sidecars"
status: >
  OPEN only as a typed cross-domain GRAMMAR. The former direct-map reading is
  REFUTED: THM-2542 constructs an F13-valued graph-Cech class, while C91 is
  its mapping-torus carrier length, not a mixed Z/91 coefficient class.
  Berggren ancestry has zero graph H1, and the odd LRC coefficients admit no
  nonzero homomorphism to or from the JC S3, V4, or characteristic-zero
  response carriers. THM-3354 proves the direct-map no-go and records the
  corrected comparison cospan as DEFINITIONAL bookkeeping. THM-3431 now
  proves a secondary refinement: a deck-group H1 transgression on the LRC
  side and a vertical local-H1 injection on the selected one-root JC side,
  with zero additive cross-maps. THM-3496 adds one sharply marked
  coefficient-changing exception: the graph seam maps isomorphically to the
  mu_13 Kummer line of an algebraically closed punctured normal slice and
  commutes with degree pullback, but it stops before additive flux and every
  physical/Keller predicate. This file is NOT a reduction.
source: opus-2026-07-26 (successor to the shape-confirmation in
  07-reflections/spectral-vs-geometric-rank-why-LRC-stays-open-and-JC-fell-opus-20260726.md)
related:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2337 / THM-2356 / THM-2334 (word-current, Bockstein, 169 twists)
  - THM-2389 / THM-2406 / THM-2463 / THM-2468 (flux, pole-descent, square lifts)
  - THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction
  - HYP-9030-keller-degree-semigroup
---

# HYP-9031 -- the D5 dictionary, written out

## TYPE CORRECTION -- 2026-08-12

The historical proposal below overidentified objects which merely share an
obstruction/realization grammar. The current typed facts are:

- THM-2542 proves [g]=7a in H^1(C_7^graph;F_13). Its skew product is one
  cycle of length 91, but that orbit length does not supply a nonzero C7
  coefficient or a primitive mixed Z/91 class.
- THM-3336's Gaussian charges live on the multiplication group. Their pullback
  to Berggren endpoint paths is a coboundary because the ancestry graph is a
  tree (THM-3345).
- The sporadic S3 Galois-closure torsor, quartic V4/mu2 Kummer plane, and
  Hamiltonian response/de Rham module are three different JC carriers. None
  is a common "JC flux" object.
- Every coefficient homomorphism between C13 or C91 and S3 or V4 is trivial
  in both directions; the same holds between those finite odd groups and the
  additive characteristic-zero response module.

Therefore there is no unmarked map between the original geometric sites and
no additive LRC-class-to-Hamiltonian-flux map in the present canon.  THM-3496
does construct a marked coefficient-changing map to the selected divisor's
`mu_13` normal Kummer line; its algebraically closed unit-root hypothesis,
oriented uniformizer, deck generator, and exponent-one normalization are all
load-bearing.  The broader lawful survivor is a **typed comparison cospan** recording site,
coefficient object, distinguished class/observer, target predicate, lost
information, missing realization sidecar, and quantifier. The historical
same-H1 wording below is retained only for provenance and must not be cited as
current truth. THM-3354 gives the proved direct-map no-go, the integral versus
generic response split, and the exact hostile controls. THM-3431 proves the
strongest current secondary replacement: chart holonomy transgresses to the
degree-13 cover primitive's deck defect, a selected one-root JC observer
embeds as `[lambda^(-q)]` in vertical local cohomology, and the two classes
still admit no nonzero additive map in either direction. Their shared
DeathBar record is explicitly lossy and has no universal-property claim.

## 2026-08-16 refinement -- the marked normal Kummer line is reachable

[THM-3496](../../01-canon/theorems/THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction.md)
proves the normalized map

```text
[g] |-> (sum_i g_i)[y^13=lambda]
```

from the oriented seven-chart graph line to the `mu_13` Kummer line of one
oriented punctured formal normal slice.  Degree-`k` graph pullback matches
`lambda=t^k`; degree thirteen kills and degree fourteen restores both sides.
This is not a reversal of THM-3354: changing coefficients to `mu_13` and
marking the normal parameter creates a new correspondence object.  There are
twelve nonzero unmarked scalar choices, and over `Q((lambda))` unit classes
such as `[2]` make the full Kummer group larger than the valuation line.

The same theorem sharpens the target obstruction.  For `P=x+x^2z`, the
nonzero characteristic-zero unit response becomes exact over `Z/13^r` for
every finite `r`.  This is finite-coefficient extinction, not a universal
Bockstein or derived-completion no-go.  The physical word-current-to-chart
arrow and the Kummer-to-additive-flux arrow both remain open/blocked in their
previous precise senses.

## Historical proposal -- superseded as a direct map

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

**Scorecard (updated 2026-07-27):** prediction 1's existence half is
fulfilled (THM-2479, countersign pending) and then STRENGTHENED beyond
the prediction: the THM-2334 (42) variance object is positive on the
canonical typed row and (klein + opus exact inverse DFT; two independent
reruns) ALL 168 nonzero target aggregates `A(q)` are nonzero -- full
target-plane support. Prediction 2 is CONFIRMED 10/10: every degree-22
support-two plane closed by restoring a discarded root coordinate
(`Y^2=1/p` Kummer lifts: THM-2463/2468; root-ratio normalizations
`r = E/(Dy)`-style: THM-2429/2437/2469/2470/2472/2475/2476/2480; the
'Hensel' steps are irreducibility subroutines, not new closure
mechanisms). Prediction 3 (orientation 2-torsor per level) remains open
and is now the live prediction.

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

## 2026-07-27 convergence: the LRC cocycle is now CONSTRUCTED

Independently of this file, the overnight fleet built exactly the
dictionary's LRC column: THM-2542 (codex) computes a **nonzero Cech
1-class on the seven transported selector/scheduler charts** (the C91
mapping torus of the root-deck local system; audit in flight), and
[root-frontier-cont16]'s C91 ledger proves the representation-theoretic
no-go `Hom_G(Q(zeta_91), Q(zeta_13)) = 0` -- no equivariant map
produces a pure 13-root dipole from the primitive conductor-91 defect.
Both slot directly into the dictionary: the constructed Cech class is
the spectral `rho_LRC` made explicit (chart cover replacing the
abstract cycle), and the Hom-vanishing is the coefficient-mismatch row
sharpened to a theorem (the `7 (x) 13` tensor admits no equivariant
collapse to its 13-line -- THM-2506/2508's cut/F_7-character necessity
is the same fact in transform coordinates). The named residual on the
LRC side is now the **semantic 2-cell** (identify the abstract mapping
torus with a physical endpoint/current) -- i.e. the `H^2`/incidence
half, exactly as prediction 1 framed it. The dictionary's "one degree
up" slot is where both sides' open problems now live: JC's higher mixed
strata vs LRC's semantic-arrival identification.

**Hostile correction (THM-2547 / MISTAKE-281).**  The later `432/432`
cut/root pairing does **not** discharge this semantic 2-cell.  It multiplies
a synthetic THM-2471 companion control by THM-2506's high-branch defect only
after an untyped numerical character-label alignment; no common ancestry
base or physical endpoint/current is present.  The Boolean `C=delta_2`
hostile has full support in both factors but exact cross-colour cancellation.
Thus the convergence statement above remains deliberately one degree short:
the missing incidence map must be built atomwise on each live row before
marginalization and integration.

## 2026-07-27 refinement: the channels separate in JC but MERGE in LRC

mac-mini S143 (rank-12 box anatomy, refereed): pair-relation saturation
makes the rank-12 box the all-small height box, and total mod-p
invisibility (the support-3 zero-trick restores full sparse rank at
p = 7, 13) means NO CRT/mod-p instrument sees it -- so there is no
independent finite discharge of LRC's "geometric channel": its depth
RE-MERGES with the spectral 91-line. This sharpens the dictionary's
divergence axis: in JC the geometric channel genuinely fell separately
(a collision realized the class in dim 3); in LRC the geometric and
spectral channels are ONE at depth, and the only forcing primes are 7
and 13. Same session: the mod-2 (palindrome/mirror) forcing of the
91-line is IMPOSSIBLE (Smith theory: Z/2-equivariance with odd-order
7 (x) 13 coefficients vanishes identically) -- the repaired forcing is
Z/7 x Z/13 ROTATION equivariance, whose localization lemma (S144:
covariant counts localize onto duty carriers, N(v) = 10 z_13(v) mod 13)
is the correctly-primed Redei mechanism. Dictionary consequence: the
"abelianized invariant" row's LRC entry should be read at the rotation
fixed loci (degenerate-word/duty contributions), not at any 2-adic
fixed locus.

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
