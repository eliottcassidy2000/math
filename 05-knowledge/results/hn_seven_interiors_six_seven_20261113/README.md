# Six and seven: seven private interiors cannot force, and where 6|7 thresholds come from

mac-mini-2026-11-13-cuts (second prompt) · **Hadwiger–Nelson: χ(R²) ∈ {5,6,7} unchanged; Chan–Pak OPEN.**
Owner prompt: use "30000-like numbers relating to 6- and 7-vertex strong components" to improve the repository's
Hadwiger–Nelson advances, noting that 30000 is half of 60000, and relate this to Fermat blocks from past sessions.

Status: **PROVED, computer-assisted** (§1, using finite enumeration plus exact planar necessary conditions; an
independent replay is welcome but was not done here). **FINITE-EXACT** (§2, §3). The remaining items are ideas and
conjectures (§4). Scripts are in [`scripts/`](scripts/). They need nauty `geng`/`gentourng` and networkx;
`horizon.py` needs the `hl*.txt` data of
[busch_staircase_fermat_block_cuts_20261113](../busch_staircase_fermat_block_cuts_20261113/scripts/README_SCRIPTS.md).

## 1. Seven-interior extension theorem (HN; extends the five-interior theorem)

Setting (exactly that of [FIVE_INTERIORS](../paired_budget_selectors_matrix_20261028/FIVE_INTERIORS.md)): eight
distinct points form four disjoint named pairs. A proper six-colouring of their unit-distance graph is prescribed,
constant on each pair. Then add private interior points.

**Theorem.** With at most **seven** private interiors, every such prescription extends to a proper six-colouring of
the whole unit-distance graph. Hence a forcing selector in this architecture needs at least **eight** private
interiors. The previous bound was six.

**Reduction.** Merge equal-coloured pairs into k ≤ 4 colour groups and form the critical quotient: the groups form a
clique, every interior unit edge is kept, and a group is joined to an interior if one of its points touches it. A
failed extension contains an induced vertex-7-critical subgraph Q, of minimum degree ≥ 6. Its complement therefore
has maximum degree ≤ |Q| − 7. Enumerating complements with `geng`:

| \|Q\| | 7-critical graphs | shape |
|---|---:|---|
| 7 | 1 | K7 |
| 8 | 0 | – |
| 9 | 1 | K4 ∨ C5 |
| 10 | **7** | all are **K3 ∨ G7**, with G7 one of the seven 4-critical 7-vertex graphs whose complement has max degree ≤ 3. One of them is **K3 ∨ (Moser spindle)**, another K3 ∨ C̄7 |
| 11 | **237** | – |

**Planar necessary conditions**, all exact:
- no unit K4;
- no K_{2,3}, since two points have at most two common unit neighbours; this is applied to interior pairs (counting
  port points), to interior–port pairs, and to port–port pairs;
- each port point's interior contacts lie on its unit circle, so their induced unit graph is a linear forest or a
  regular-hexagon 6-cycle;
- a colour group has two points, or four when two pairs share a colour. WLOG each contact is made by one point,
  since every filter is monotone in the contacts.

**Cases** (scripts `six_interiors.py`, `seven_interiors.py`, `seven_k3_fix.py`):

| interiors | quotient | groups | survivors |
|---|---|---|---:|
| ≤ 5 | K7, K4 ∨ C5 | – | 0 (five-interior theorem) |
| 6 | 10-vertex K3 ∨ G7 | 4 pairs | **0** |
| 6 | K4 ∨ C5 | 3 groups | 0: the interiors would contain a unit K4 |
| 7 | 10-vertex K3 ∨ G7 | 3 groups, sizes (4,2,2) | **0** |
| 7 | 11-vertex (237 graphs) | 4 pairs | **0** |
| 7 | K7 or K4 ∨ C5 with ≤ 2 groups | – | 0: unit K4 or K5 among the interiors |

**Why it fails.** For K3 ∨ G7, the three universal vertices of K3 cannot all be port groups unless one G7 vertex
joins them. If one universal vertex is an interior, it is a hub: every other interior in Q lies on its unit circle,
so G7 minus the ports must be a 60°-linear forest. That is incompatible with 4-criticality once K_{2,3} and the
shared-centre restrictions are imposed. Two universal interiors give K_{2,3}, and three give a unit K4. The 11-vertex
graphs die to the same filters at the contact-split stage.

**Next step.** Eight interiors require the 12-vertex 7-critical graphs, whose complements have max degree ≤ 5. The
pipeline is unchanged, but `geng 12 -D5` is large, so it should be generated with a criticality prefilter. The
10-vertex step shows what a selector must eventually contain: a 7-vertex 4-critical core joined to universal
contacts. Of the seven cores, three pass the unit-distance necessary conditions (no K4, no K_{2,3}): the Moser
spindle (11 edges) and two 12-edge cores. The other four contain K_{2,3}. So a genuine selector is likely to be
*spindle-like inside* and *hub-free outside*.

## 2. The tournament 6|7 thresholds (FINITE-EXACT; the first item is classical)

- **Dichromatic number.** Every tournament on ≤ 6 vertices splits into two transitive halves. Exactly 4 of the 456
  tournaments on 7 vertices do not (`dichro.py`), and 258 of 6880 on 8 vertices do not. This matches Neumann-Lara:
  the least 3-dichromatic order is 7, with four tournaments, and the least 4-dichromatic order is 11. Their path
  counts are 99, 151, 171 and 189, the last being the Paley tournament P_7.
- **Multiples of 7.** No strong 6-tournament has H divisible by 7. Its counts are 15, 17, 19, 23, 25, 27, 29, 31,
  33, 37, 41, 43, 45. Order 7 is where the 7-lane is born: 35 (2 classes), 49 (7), 77 (10), 91 (15), 105 (4),
  133 (6), 147 (1), 175 (1), 189 (P_7). The values 63, 119 = 7·F2 and 161 = 7·23 first need order 8 or more.
- **21** is the gap of the Busch step [15, 25) at order 6
  ([staircase](../busch_staircase_fermat_block_cuts_20261113/README.md)).
- **P_7 = Cay(Z_7, ⟨2⟩).** It exists because ord_7(2) = 3 is odd, so {1, 2, 4} ∩ −{1, 2, 4} = ∅. This is the same
  fact that makes 7 divide no Fermat-block count 2^j + 1. Also 7 = Φ3(2) = N(2 − ω) and 21 = Φ3(4) = N(4 − ω).

**HN dictionary (PROVED, elementary).** The hexagonal 7-colouring of the plane is the Eisenstein residue map
Z[ω] → Z[ω]/(2 − ω) ≅ F_7, with ω ↦ 2. The six unit directions go to ±⟨2⟩ = F_7^×. Orienting the three directions
1, ω, ω² forward, the colouring is a homomorphism of the oriented triangular lattice onto P_7. Every
digraph that maps homomorphically to a tournament on ≤ 6 vertices is 2-dicolourable: pull back the two transitive
halves. So, for any *geometrically prescribed* orientation, an oriented six-colouring forces a 2-dicolouring. An
oriented-colouring obstruction at six colours therefore reduces to a two-colour (no monochromatic directed cycle)
problem. This is a variant of HN, not HN itself, which allows the orientation to depend on the colouring.

## 3. The 30000 / 60000 trail and Fermat blocks

- **THM-5070** (port chains) glues blocks after a port, with a single back arc w → u. These are Fermat blocks
  generalized, with a nonnegative 3-dimensional state (I, E, O) in any base x. Its gap sets of I(Ω(T), b) in
  **bases 6 and 7**, up to α-sum **30,000**, are {3, 10, 14, 21, 46} and {3, 10, 14, 29}. The gap m = 3 (exactly
  three pairwise-meeting odd cycles) is universal across bases.
- **Halving.** Deleting the ready vertex of a Fermat block sends H to (H+1)/2
  ([mutation_collatz](../mutation_collatz_chanpak_20261112/arithmetic/PROOF.md)). So 2^16 + 1 = 65537 = F4 halves
  to 2^15 + 1 = 32769, "60000 → 30000". Both exceed H_max(10) = 15745, so both first occur at order ≥ 11, which is
  the least 4-dichromatic order. The Fermat count packet records least orders 3, 4, **6**, 8, 11 for F0–F4, so F2 = 17
  is born at order 6, from R_4.
- **Cut horizons.** The least odd N > 67 that is not a single Fermat-block cut h + 2^j L of a source on ≤ m vertices
  is X(m) = 71, 275, 1439, 22229, 398131 for m = 5, …, 9 (`horizon.py`). These are the numbers of the right scale
  produced by order-6 to order-9 sources.

## 4. Ideas toward HN

1. **A critical-quotient staircase.** Five-, six- and seven-interior exclusions are the steps K4 ∨ C5, K3 ∨ G7 and
   the 237 eleven-vertex graphs. This parallels Busch's staircase, where a count fixes how large a strong realizer
   can be. Pushing to 8 or 9 interiors is mechanical. The first step with a survivor is where to build a selector,
   and survivors come with an explicit quotient to realize.
2. **Spindle cores.** Since K3 ∨ Moser spindle is a 10-vertex critical quotient, try selectors whose interior graph
   contains a spindle. Port pairs should touch the spindle without creating hubs; this is what the filters forbid
   at 7 interiors. With more interiors, contacts can be spread across distinct circles.
3. **Positivity barriers.** Positive ear factors preserve boundary support (two-port unit-distance packet), and
   nonnegative port-chain transfers (THM-5070) cannot create zeros. In both problems, forcing needs a
   sign-changing or support-killing gadget. On the tournament side the only zero-creating events are Busch-floor
   rigidity and the 7, 21 gaps.
4. **"At most one more" model.** The idoneal planes (THM-4566, checked for |D| ≤ 30000) have at most one further
   exception by Weinberger's theorem, via Tatuzawa's at-most-one-Siegel-zero result. That is the precise template
   for the owner's dichotomy wish ("one more, or a family"). An analogous Chan–Pak statement would need realizer
   counts with a uniform lower bound failing for at most one target.
