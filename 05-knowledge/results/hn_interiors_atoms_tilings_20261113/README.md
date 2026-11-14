# Nine private interiors, the four-vertex atom, and 3|5 from tilings to Chan–Pak

mac-mini-2026-11-13-cuts (third prompt) · **Hadwiger–Nelson and Chan–Pak OPEN.** Owner prompt: make progress on the
eight-interior (12-vertex) case. Relate 12 and 8 to the 3-and-1 structure of the four-vertex tournament, using the
atom made of four arcs with the two remaining arcs as flip-able diagonals. Merge prior Collatz and hexagonal-tiling
work. Read the 3 regular tilings against the 5 Platonic solids (one object stays, the two subordinates each split)
and relate that to the 3s and 5s of Chan–Pak.

Status: **PROVED, computer-assisted** (§1; positive and negative controls included; an independent replay is
welcome). **CITED / ELEMENTARY** (§2, §4). **ANALOGY with exact arithmetic content** (§3). Script:
[`scripts/interiors_sat.py`](scripts/interiors_sat.py) `M [--nogeo]` (pysat + nauty `geng` + networkx).

## 1. Nine-interior extension theorem (HN)

Setting of [FIVE_INTERIORS](../paired_budget_selectors_matrix_20261028/FIVE_INTERIORS.md): four disjoint port pairs
carry a prescribed proper six-colouring, constant on pairs, and private interior points are added.

**Theorem.** With at most **nine** private interiors, every prescription extends. So a forcing selector in this
architecture needs at least **ten** private interiors. This improves the seven-interior theorem of
[hn_seven_interiors_six_seven](../hn_seven_interiors_six_seven_20261113/README.md), which itself improved on six.

**Method: interior-first CEGAR.** Rather than enumerating 12-vertex 7-critical graphs, fix the interior unit graph H
and let SAT choose the contacts t(p, x) of the port points.
- **Minimality.** Run m = 5, 6, … in order. Configurations with ≤ m − 1 interiors are then already excluded, so a
  failing m-interior configuration has a critical quotient using all m interiors.
- **Criticality bounds.** This gives deg_H(x) + |c(x)| ≥ 6 for every interior, where c(x) is the set of colour groups
  x touches. For each group, |S_g| ≥ 7 − k with k the number of groups used. Consequently H is connected, has
  min degree ≥ 6 − k, and is K4-free and K_{2,3}-free.
- **Planar clauses.**
  - Each port point's contacts induce a forest of max degree 2 in H, or exactly a hexagon.
  - Two port points share at most 2 contacts.
  - An interior and a port point share at most 2 interior neighbours.
  - Two interiors have at most 2 common neighbours, counting interiors and port points.
- **Counterexample loop.** List-colour H with lists [6] minus the colours of touched groups, preferring the two free
  colours. Each successful colouring φ yields the clause "some interior touches the group whose colour φ gives it".
- **Group sizes.** These are (2,2,2,2) for 4 colours; (2,2,2) or (4,2,2) for 3; and (2,2), (4,2), (4,4) or (6,2)
  for 2.

The clauses are relaxations of the true geometry (for example, the 300° arc budget on a circle is not imposed).
Relaxing only weakens the filters, so it can produce false survivors but never false exclusions.

| interiors | interior graphs (k = 4 / 3 / 2) | cases | survivors, geometry on | survivors, geometry off (control) |
|---|---|---:|---:|---:|
| 5 | – | 4 | **0** (reproves five-interior theorem) | 18 of 21 |
| 6 | – | 24 | **0** (agrees with §1 of the seven-interior note) | – |
| 7 | – | 99 | **0** (agrees) | – |
| 8 | 607 / 41 / 1 | 693 | **0** | 12,001 of 14,316 |
| 9 | 5,639 / 352 / 4 | 6,359 | **0** | – |

The search takes 8 seconds at eight interiors (at most 402 counterexample rounds) and 3.5 minutes at nine (at most 852). Without the planar clauses,
forcing configurations are abundant, so every exclusion is genuinely geometric.

**Twelve and eight.** The eight-interior quotient has 12 = 4 + 8 vertices: a K4 of colour groups plus eight
interiors. At ten vertices every critical quotient had the shape K3 ∨ G7, three universal vertices over a
seven-vertex 4-critical core. In the language of joins, the "3 and 1 structures" are the critical joins K3 ∨
(4-critical) and K1 ∨ (6-critical). The four-vertex atom (§2) is exactly K1 ∨ K3 as a graph. The interior-first
search never needs the 12-vertex join list, because the geometry already decides at the level of contacts.

## 2. The four-vertex atom (CITED: [trees5_tournaments4](../trees5_tournaments4_erdos_sos_20261005.md))

Fix four arcs of a 4-tournament: the path 1→2→3→4 plus the long arc 1→4, an undirected square. The two remaining arcs
are the diagonals {1,3} and {2,4}. The Klein group (Z_2)² of diagonal flips produces all four isomorphism classes:

| flipped diagonals u | class | H = 2u + 1 |
|---|---|---|
| none | TT4 | 1 |
| {1,3} | C3 → sink | 3 |
| {2,4} | source → C3 | 3 |
| both | strong R_2 | 5 |

Its counts 1, 3, 5 are the whole Chan–Pak story in miniature. 3 = F0 and 5 = F1 are the only strong counts on ≤ 4
vertices, and they are the Busch bricks: every Busch floor value is 3^e·5^k. The missing next value 7 = 2u + 1 at
u = 3 would need a third diagonal, which a 4-vertex atom does not have. That is the "3 and 1" decomposition
4 = 1 + 3 (one apex over a cyclic triangle) seen arithmetically.

## 3. Three tilings, five solids, and the 3|5 of Chan–Pak (ANALOGY with exact arithmetic)

| family | self-dual or self-converse | dual or converse pairs | classes after merging |
|---|---|---|---|
| regular tilings {3,6}, {4,4}, {6,3} | {4,4} | ({3,6}, {6,3}) | 2 |
| Platonic solids | {3,3} | ({3,4}, {4,3}), ({3,5}, {5,3}) | 3 |
| 4-tournaments | TT4, R_2 | (C3+sink, source+C3) | 3 |

The Schläfli product (p−2)(q−2) is 4 for every tiling, and 1, 2, 2, 3, 3 for the solids. The self-dual tetrahedron
gives 1, and the two dual pairs give 2 and 3. The one fixed object stays put while the dual pair of tilings
({3,6}, {6,3}) splits into two dual pairs of solids, with the 6 dropping to 4 or 5.

**The exact Chan–Pak counterpart** comes from the 2-adic lines 1 + c·2^s (THM-5126, HYP-9658):
- c = 1 is the line of the Fermat blocks themselves, 1 + 2^s = H(R_s). It is complete, the object that stays.
- c = 3 = F0 and c = 5 = F1 are the two subordinate lines. Each loses exactly one point:
  **7 = 1 + 3·2^1** and **21 = 1 + 5·2^2**.
- Both lost points sit at s = 2^r: 7 = 1 + F0·2^(2^0) = Φ3(2) and 21 = 1 + F1·2^(2^1) = Φ3(4).

So 3 and 5 enter Chan–Pak exactly as F0 and F1 enter the Busch floor, and the two exceptions are the first two values
Φ3(2^(2^r)) = F_{r+1} + F_r − 1 of the Fermat "join" sequence. The next value, 273 = Φ3(16), is attained.

## 4. Hexagons, the torus, and Collatz (CITED / ELEMENTARY dictionary)

- **Heawood.** The torus needs 7 colours. The Heawood map has seven mutually adjacent hexagons (cf.
  [natural_number_roles](../natural_number_roles_20261004.md)), and its dual is K7 on the torus. Lifting it to the
  plane gives exactly the hexagonal 7-colouring Z[ω]/(2 − ω) ≅ F_7 with ω ↦ 2. Orienting K7 by differences in
  {1, 2, 4} = ⟨2⟩ gives the Paley tournament P_7. The set {1, 2, 4} is the Fano perfect difference set (its
  differences cover 1–6 once).
- **Seven and twenty-one.** 7 = Φ3(2) is the number of Heawood hexagons, Fano points and P_7 vertices, and it is a
  Chan–Pak exception. 21 = Φ3(4) is the number of K7 edges, P_7 arcs and Fano flags, and it is also a Chan–Pak
  exception.
- **Collatz on the hexagonal lattice.** In Z[ω], 3 ramifies: 3 = −ω²(1 − ω)², the √3-rotated subdivision of index
  3. And 2 is inert: the ×2 subdivision of index 4. The two Collatz moves, 3x+1 and x/2, are therefore the two
  primitive self-similarities of the hexagonal tiling. The colouring sublattice (2 − ω) of index 7 is the first split
  prime, and its norm is the first Chan–Pak exception.
- **Prior hexagon threads.** THM-4663 ties unit-distance hexagons to level 6, with the four-bar second moment the eta
  product (η1η2η3η6)²; the torus-clock note's cellular link is a union of 60 hexagons.
- **Suggested probe.** The hexagonal "Collatz tower" (index 3 versus index 4 sublattices) composed with the
  index-7 colouring gives a natural family of plane colourings indexed by Löschian norms 3^a·4^b·7. Their
  conflict-free radii are an HN-flavoured test of how 3, 4 and 7 interact.

## 5. Next

1. Push the interior bound. At ten interiors, `geng 10 -c -d2` gives about 10^7 graphs before the K4 and K_{2,3}
   filters, so generate with an edge cap (K_{2,3}-free graphs on 10 vertices are sparse) and filter in C. The SAT
   cost per graph is milliseconds. Use symmetry breaking, and add the circle-arc budget if survivors appear.
2. If survivors ever appear, they come with explicit contacts. Realize them in exact coordinates with the repo's HN
   toolkit; a survivor that realizes is a genuine six-colour gate.
3. Chan–Pak: test whether "Φ3(2^(2^r)) is attained for r ≥ 2" has a uniform construction (273, 65793, …). That
   would show the Fermat-join sequence is exceptional only at its two bricks.
