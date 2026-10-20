# Inheritance, concept board, and connection contracts

Status: current session routing, 2026-10-20. The proofs are in README;
this file records why these particular routes were pursued.

## Inheritance pass

* Closest proved mechanism: [THM-4512, coefficient-descent classes](../../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
  and [the Oct4 boundary compiler, BC1--BC2](../collatz_boundary_compiler_20261004.md).
  The old actual-carry threshold already leaves at most one representative
  per coarse cylinder at every length. We do not present that as new.
  The new threshold bounds the stronger positive-only charge and has a
  different sharp constant, `161/256`.
* Canonical hostile: the [adaptive-budget packet](../adaptive_budgets_smooth_interfaces_20261020/README.md)
  retains source63, five necessary first-exit sorting blocks, and the
  common-triangular-lattice obstruction. We retain its actual word.
* Corrected near miss: exact word classes use `2^(K+1)`, not `2^K`.
  The positive unit emission uses the post-step denominator. Both are
  checked here by a coarse/exact split and a separate suffix formula.
* Least-used sidecar: [THM-2646, central braid height](../../../01-canon/theorems/THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre.md)
  (July28) explains why a permutation or modular graph forgets a knot's
  integer winding information. [The K4 defect web](../defect_boundary_20261020/web.md)
  already distinguishes graph complement, tetrahedral opposition and
  line-graph loss; the new game consumes this distinction.
* Geometric guardrail recovered before iterating: [THM-4558](../../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md)
  bounds the Moser field by four colours. [Boundary palette compression](../boundary_palette_compression_20261011.md)
  already gives the general port-relation formalism and classical spindle
  controls. Our contribution is the exact two-wheel embedding and its
  six-colour missing-palette kernel, not a new spindle or a new general
  boundary calculus.

## Portfolio and live board

Anchor: the user-assigned first-exit positive-charge question. Niche:
non-lattice boundary composition. Wildcard: complementary spanning games
and a genuine line-graph/knot comparison.

| Object / representation | Predicate | Operation | Retained obstruction / next test |
|---|---|---|---|
| Actual odd orbit / affine word | positive-only payment | exact cylinder refinement | source, K, m, first-exit history; least representatives beyond logarithmic depth |
| Contracting coarse cylinder | every source above its head pays | isolate one representative | extra final divisions; do not identify coarse and exact endpoints |
| Two hinged wheels /17 ports | extension in six colours | rotation plus one new hub edge | actual missing colour names; seek a further field escape |
| Complementary graphs | both nontraceable | edge deletion / separator | component count after deleting named vertices; connected examples from n7 |
| Torus braid / step permutation | one component and knot type | reduce q modulo p | retain integer q and strand order; counts alone fail at T(2,1),T(2,3) |

Existing META-PATTERNS cards used: **Search the statement before the
method**, **Find the hidden second coordinate in a nearly true theorem**,
and **Re-evaluate a certificate after every fibre-changing operation**.
No duplicate method card is added.

## Incoming work read as mathematical signal

Startup fast-forwarded seven commits from `ffc5216306` to `7f8df95a96`.
A later fetch during the session found no further incoming commits at
that checkpoint.
The checkpoint push subsequently integrated `c402d288ac`, `9852737522`,
`411c1fd35a` and `5e93d8adf6`. Their game and geometry proofs were then
read, and the stronger results below supersede the former open game target.

1. [Neighboring disjoint two-swap exclusion](../collatz_two_swap_neighbor_disjoint_exclusion_20261020/NEIGHBOR_DISJOINT_EXCLUSION.md)
   closes a precisely specified all-length cycle-word family, using a
   content bound, logarithmic separation and a finite residual certificate.
   This motivated reusing the existing effective power-gap estimates for
   a charge threshold. The route does not exclude arbitrary cycles or
   provide first exit for our marked source.
2. [Golden inventories and outer detours](../spectrum_basins_20261020e/README.md)
   now connect all active `{1,C3,R4}` basins and handle nonstrong mixed
   blue detours. Strong-to-strong detours remain open. These are actual
   legal move sequences, unlike sorting a Collatz word for an upper bound.
3. [The local matching interaction](../spectrum_basins_20261020e/arithmetic_interaction.md)
   retains E and C, with `H=3^(E-2C)5^C`; matching size alone loses the
   adjacent-row interaction. The geometric analogue has an exact map:
   separate extension counts `|M|,|N|` lose `|M intersection N|`, while
   the named incidence vectors recover (8). The two C's are unrelated
   quantities; their shared feature is an explicitly recovered interaction.
4. [The reviewed complement-game packet](../complement_path_games_20261020/AUDIT_AND_PROVENANCE.md)
   independently agrees on the four existence thresholds and connected-both
   family, proves the remaining sharp W/H costs, and supplies prescribed-
   start and exact join/path-cover criteria. Our open general W/H budget
   question is therefore closed by this incoming work.
5. [Fixed-deficiency kernels and forcing](../kernels_games_detours_20261020f/README.md)
   now give a polynomial-times-double-factorial formula at every fixed
   deficiency D, with degree floor(5D/2), while leaving the zero-heavy
   boundary unrestricted. That quantifier is analogous to our bounded-
   depth certificate: every fixed regime is different from a uniform
   entry theorem. Its sixteen-point palette gate is a concrete new input;
   [PALETTE-GLUE.md](PALETTE-GLUE.md) composes two actual copies at one
   root, preserving all first-copy ports and producing the exact
   intersection kernel. The all-order rotation-avoidance proof and the
  31-point certificate were developed after reading that incoming result.

## Connection contracts

**Carry to positive charge.** Source: the ordered exact word. Target:
`(A,d,C)` and the threshold R. Map: retain only exponent-one emissions
in shifted coordinates. Preserves a sufficient upper bound and actual
integer clocks. Destroys the negative reset dividend. Sidecar: D and
the first-exit predicate. Decisive hostile: source55's later descent
to53 fails, whereas the short first exit pays. Transfer: the inherited
all-length B bound proves the new all-length C bound, and the effective
power gap proves the logarithmic-depth region.

**Geometric composition.** Source: two actual point sets with named ports.
Target: two missing-colour masks and the hub edge. Map: forget interiors
only after recording their neighbour sets and every cross contact.
Preserves extension count exactly. Destroys the individual rim assignment
and coordinate lift. Sidecar: full points and port assignments. Test:
equal singleton masks give zero; distinct singleton masks give one.
The rank-six kernel is indefinite; a PSD assumption would erase the
negative overlap term. This is not a proof about the ROOT Gram operator.

**Hamilton to Euler.** Source: directed edge-words. Target: vertex-words
of the directed line graph. Map: each labelled edge becomes a vertex.
Preserves ordered traversal and completeness. Destroys the old vertex
identity unless endpoints are retained. Test: claw edges give an
undirected-line-graph false positive; directed endpoint compatibility
restores the correspondence. The two-player game itself has symmetric
edge colours, so no artificial tournament orientation is introduced.

**Torus braid to finite graph.** Source: p ordered strands and integer q.
Target: the permutation `i->i+q mod p`. Preserves component count and
cyclic orbit length. Destroys full twists. Sidecar: integer q, over/under
crossings, strand order. Test: q1 versus q3 with p2. No claim identifies
Hamilton-path counts or Euler lengths with topological winding numbers.

## Literature pin and scope

* Ellison, *On a Theorem of S. Sivasankaranarayana Pillai*, 1970--1971,
  [primary PDF, Theorem3, printed p12-05](https://www.numdam.org/article/STNB_1970-1971____A10_0.pdf).
  Read the theorem, exceptional exponent set and subsequent continued-
  fraction reduction. Imported only the stated exponential separation
  of `2^x` and `3^y`. This input already appears in THM-4512's Oct7
  alternative proof. It does not assert anything about source residues
  or Collatz entry. The proof of the source/depth cone is local here.
* Rozier--Terracol, [*Paradoxical behavior in Collatz sequences*,
  arXiv:2502.00948v2](https://arxiv.org/html/2502.00948v2), Definitions1.1--1.2,
  Lemma2.1 and Theorem3.2 were read. The coefficient/actual stopping-time
  equality is a separate open conjecture; half-step and odd-step clocks
  must be translated. Their adjacent parity swap changes a carry;
  it does not authorize rearranging a source's actual orbit. No numerical
  verification range from that paper is imported into our charge result.
* [Madore, arXiv:1509.07023](https://arxiv.org/abs/1509.07023) was located
  as the classical field-colouring source. The specific four-colour
  carrier bound is inherited with its attribution from THM-4558, not
  claimed to follow from the abstract alone.

## Remaining questions ranked by consequence

Final incoming pass through b9cc570666 read the strong-blue successor's
README, pinned-stress proof and boundary-transport identities. Its
two-block repairs and role-exchanging boundary hubs strengthen the
tournament lane without establishing a Collatz map. The source spectrum
packet now routes to this successor; its earlier statements are unchanged.
The sharper pinned-stress bound applies exactly to our two-hub gate:
h=2, a=12, t=10 and sigma_B=9. The derivation in PALETTE-GLUE.md retains
all17 port positions and counts repeated forbidden colours. This adds
a quantitative geometric screen, not global six-colour forcing.

1. At first coefficient crossing, can the least source r satisfy
   `r<=R(w)`? The symbolic certificate rules it out through18; the
   analytic cone rules it out for `m<=4 log r`, r>=65536. The unrestricted
   first-exit version is still open, and existence of a crossing is an
   additional obligation.
2. Can the reciprocal reset budget in (6) be controlled during a long
   first excursion without assuming its length? Proper-prefix values
   matter: an arbitrary later descent already fails the target at55.
3. Can an authenticated context in a field beyond the known four-colour
   carrier force the same missing singleton at the adjacent hubs? A
   boundary relation can be restrictive while the full graph remains
   four-colourable; the current gadget demonstrates both facts.
4. The disconnected/hamiltonian minimum cost is now closed by incoming
   work: `2n-5` except n6,8. The minimum under the additional restriction
   that **both** graphs remain connected is still a separate target;
   our family gives the upper bound2n-4. No sharpness claim is made.
