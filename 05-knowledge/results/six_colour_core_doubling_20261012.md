# The six-colour core at the eight ports: nine double-port vertices, the doubling extension theorem, and why no eleventh-root context with a three-colourable exterior can force

2026-10-12, opus-2026-10-12-S24. **PROVED:** the census, the doubling lemma
and the extension theorem ([THM-4614](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)).
**FINITE-EXACT:** the 570-pattern confirmation and six three-colourable probe
contexts. **NUMERICAL with margin:** the Hoffman bound at the residue field
of 109. **OPEN:**
three-colourability of every finite unit graph of `Q(zeta55)`
([HYP-9247](../hypotheses/HYP-9247-every-finite-unit-distance-graph-in-q-zeta55-is-three-colourable.md)),
the field's chromatic number, six-colour forcing, and the plane. Independently
audited 2026-10-12 (CORRECT WITH FIXES: the census, witnesses, apexes, the
570 patterns, the probes and the Hoffman bound were reproduced with
independent code; a port quantifier, the residue-colouring scope and a
transfer condition were corrected; MISTAKE-593). No literature priority is
claimed for the lemma (it is the standard product-of-palettes trick, applied
to a fixed port set).

## Inheritance and the changed question

`L = Q(zeta55)`, `alpha = zeta5`, `gamma = zeta11`, `eta = (-3+sqrt(-55))/8`,
ports `B` = the unit pentagon `P_j = sum_(i<j) alpha^i` and the unit rhombus
`(0, -eta, 3/4, -bar eta)` sharing `0`. [The thirds packet](eleven_thirds_boundary_20261012.md)
proved that the old 8,436-vertex mixed-generator context, the 2,278-vertex
thirds context and their 10,598-vertex union (88 seam edges) all retain the
*entire* bare relation: every one of the 570 proper six-colourings of `G[B]`
(up to relabelling; 327,600 labelled) extends, and asked for an eleventh-root
context that removes at least one pattern while keeping source and ports.
[P2 of the port-layer packet](eleven_port_layers_20261011.md) proved this for
one-layer contexts on at most five K-ports, and [B1 of the pentagon packet](pentagon_context_number_fields_20261011.md)
for five-colourable contexts of the twenty-port gadget.

The question changed here is not *which generator* but *which property of
the exterior* decides. Closest mechanism: complete boundary equality-partition
gluing with a witness lift. Hostile: a 6-regular lattice patch, whose
extreme vertex has at most half its unit directions behind it and so peels.
Corrected near miss: reading "the context is three-colourable" as a weak
positive fact about the context rather than as the reason the relation is
full. Least-used sidecar: *the double-port vertices of the port set itself*,
a finite census independent of any context.

Board: double-port census / doubling lemma / 6-core / three-colourability of
`L` / residue colourings / the two number identities. Anchor here is the
niche of the [Collatz companion](collatz_band_budget_hostile_20261012.md),
which records the shared board.

## 1. The census: exactly nine points of L see two ports

For ports `p != q` with `D = |q-p|^2` (an element of the real subfield),
the two apexes `|x-p| = |x-q| = 1` are `x = p + (q-p)(1/2 +- s/(2D))` with
`s^2 = D(D-4)`; they lie in `L` iff `D(D-4)` is a square in `L`, and since
the ports lie in `K = Q(alpha, sqrt(-11))` and `[L:K] = 5` is odd, iff it is
a square in `K`. Squareness is decided exactly: a root proposed by an
eight-term integer relation is verified by squaring in the 40-dimensional
power basis, a non-square is certified by a quadratic non-residue at a split
prime `ell = 1 mod 55`.

| port pairs | `D` | verdict |
|---|---|---|
| the nine unit edges | 1 | no apex (`D(D-4) = -3`; `sqrt(-3)` not in `L`) |
| `(P_j, P_(j+2))`, five pairs | `phi^2` | apexes `P_(j+1)` (port) and `P_j + alpha^(j+1)` (new) |
| `(P_1,-eta)`, `(P_1,-bar eta)`, `(P_4,-eta)`, `(P_4,-bar eta)` | `5/4`, `5/4`, `3.995`, `0.468` | apexes `0` (port) and `1-eta`, `1-bar eta`, `P_4-eta`, `P_4-bar eta` (new) |
| `(0, 3/4)`, `(-eta, -bar eta)` | `9/16`, `55/16` | apexes are ports |
| `(P_2,-eta)`, `(P_3,-eta)` | `4.40`, `6.10` | distance above two |
| `(P_1,3/4)`, `(P_2,3/4)`, `(P_3,3/4)`, `(P_4,3/4)`, `(P_2,-bar eta)`, `(P_3,-bar eta)` | `1/16`, `1.217`, `2.431`, `2.026`, `0.873`, `0.390` | non-residue at 661, 661, 661, 2531, 991, 331 |

So the points of `L` at unit distance from two ports are the ports and nine
further points, every one the fourth vertex of a unit rhombus spanned by two
port edges at a common port (the five pentagram apexes and four completions
at `0`); each sees exactly two ports, no point outside the port set sees
three (the port `0` itself sees four), and no two of the nine are at unit
distance (complete unit graph on the seventeen points,
exact). The closed forms are checked against the census.

## 2. The doubling lemma, and the count fifteen

Let `H` be any finite graph, `B` marked, `c` a proper 3-colouring of `H-B`,
and let every vertex of `H-B` have at most two neighbours in `B`. Pair the
six colours into three disjoint pairs (15 matchings) and assign the pairs to
the three classes of `c` (6 bijections): 90 options. Given a proper
six-colouring `t` of `H[B]`, colour `x` in `H-B` by any colour of its class's
pair not used by its port neighbours. The only vertices that can fail are
those with two port neighbours of distinct colours whose pair *equals* the
class's pair, and such a vertex blocks exactly six options (three matchings
contain its pair, two bijections send its class there). Adjacent exterior
vertices lie in different classes, hence in disjoint pairs; port edges are
proper by construction. **Fewer than fifteen blocking vertices never block
every option.**

## 3. The theorem, and its exact confirmation

**THM-4614.** For every finite unit-distance graph `G` in `L` containing `B`
with `G-B` three-colourable, `R_6(G;B) = R_6(G[B];B)`: every proper
six-colouring of the ports extends to `G`.

*Proof.* By section 1, every exterior vertex has at most two port neighbours
and at most nine have two; by section 2, at most 54 of the 90 options are
blocked. ∎

Confirmation (FINITE-EXACT): for each of the 570 patterns a SAT instance over
class assignments to the nine apexes (with the 90 option clauses) is
unsatisfiable, so no pattern is blockable even by a context containing all
nine apexes at once; abstract controls show that fifteen star vertices, or
fifteen distinct colour pairs in one class, do block.

This explains every failure to date without reference to generators,
denominators, layer counts or seam edges: all those contexts were
three-colourable, and a three-colourable exterior's only lever at palette six
is a double-port vertex, of which the whole field has nine.

## 4. Richer eleventh-root contexts stay three-colourable; residue colourings cannot help

The apex context (ports, the nine apexes, and one layer of 52 directions:
the fifth and eleventh roots of unity, the ten Mersenne directions
`u_j = (2-gamma^j)/(2-gamma^(-j))` with denominators on `2047 = 23*89`, the
thirds unit `(-5+sqrt(-11))/6`, and their negatives) has 846 vertices and
2,187 complete unit edges, is three-colourable, peels to the eight ports at palette six, and has
counts 30/180/247/113, the full relation, as the theorem requires. Five balls
over subsets of these directions are three-colourable too:

| context | vertices | complete unit edges |
|---|---:|---:|
| radius 2, roots 10+11 and Mersenne | 1,304 | 2,534 |
| radius 2, roots 10+11, Mersenne, thirds unit | 1,408 | 2,738 |
| radius 3, roots 10 and Mersenne | 4,974 | 14,224 |
| radius 3, roots 10+11 and thirds unit | 6,000 | 17,170 |
| radius 4, roots 10 and Mersenne (declared `--large` run) | 38,514 | 143,614 |

Edges are found by a bucketed version of the engine's two-prime filter
(331, 661), every survivor receiving the exact norm; the bucketed finder is
checked against the engine's pairwise finder on the marked points.

Two no-go facts bound what a proof of three- (or six-) colourability of `L`
can look like. Additive colourings of the unit-path component need a modulus
supported on conjugation-fixed primes (the component contains `1/p` for every
`p` with `-1` outside `<p>` modulo 55, hence for 2, 3, 5, 7, 11, 13, 17), the
least being 19. The least fixed prime of residue degree two is 109
([C3 of the palette-escape packet](cyclotomic_palette_escape_20261011.md)),
and the residue colouring there is a proper colouring of the Cayley graph of
`F_(109^2)` on its 110 norm-one residues; that graph is triangle-free, its
least eigenvalue is `-19.64` (floating-point character sums; the bound 6.60
clears 6 with margin), and Hoffman's bound gives chromatic number at least 7.
So the residue colouring with the smallest residue field does not reach six
colours (the residue fields `F_(ell^10)` of 19, 29, 79 are unexamined), and
any colouring of `L` with at most six colours, in particular any
three-colouring, is a non-additive rule like the periodic thirds rule. This
is the content of HYP-9247.

## 5. The reroute: what a forcing context must contain, and what the identities do not do

A context that changes the relation must have a four-chromatic exterior.
Since `L` is triangle-free, that means a triangle-free four-chromatic
unit-distance graph with coordinates in `Q(zeta55)`; none is known, and the
repo's four-critical diamond chains use `zeta60` and triangles, so they do not
qualify. The cheapest decisive tests are listed with HYP-9247 (placing the
23-vertex fish or O'Donnell's graph in `L` by solving their rigidity
equations exactly; Mycielski-type configurations on unit five-cycles).

The two number identities of the prompt are consistent with this picture.
`2047 = N(2-gamma) = 23*89` makes the Mersenne directions `u_j` and their
denominator receipt (`1799 - 248 sum u_j = 1/2047`) available; the apex
context and the balls above use all ten of them, and they are absorbed by
the theorem like every other direction. `127 = 10^2 + 3^3` gives the
Eisenstein norm-one rotation `(73+60 sqrt(-3))/127`, which leaves `L`; a
context using it is a context in `Q(zeta165)`, which contains `zeta3`, unit
triangles and the Moser spindle, and in which the theorem's census does not
apply (the pentagon edges then have apexes). That field, not `L`, is where
four-chromatic exteriors are cheap, at the cost of the triangle-free port
geometry.

The general form of the lemma applies to any port set: with palette `m`,
exterior chromatic number `k` and at most `r` port neighbours per exterior
vertex, extension is unconditional when `m >= k(r+1)`; at `m = 2k` and
`r <= 2` the double-port census decides. The repo's six-port conditional
mechanism ([marked return thickness, section 4](marked_return_thickness_20261012a.md))
therefore needs, for each of its 151 separating partitions, an eliminating
fragment whose exterior is at least four-chromatic, whenever its six ports
have fewer than fifteen double-port vertices (automatic for any six of the
eight ports here; a different port set needs its own census).

## 6. Connection ledger

| Source object | Map and preserved predicate | Lost data / retained sidecar |
|---|---|---|
| Port pair `(p,q)` | `D(D-4)` square in `K` iff apexes in `L` | the sign of the root; both apexes constructed |
| Exterior three-colouring | classes to disjoint colour pairs | the blocked options; six per double-port vertex |
| Previous contexts | instances of THM-4614 | none; their lifts remain valid witnesses |
| Residue field at 109 | Cayley graph on 110 norm-one residues | Hoffman bound only; exact chromatic number not computed |
| Mersenne directions | absorbed by the theorem | no forcing; denominator receipt retained |
| 127 rotation | lives in `Q(zeta165)` | triangle-freeness and the census |

## 7. Reproduction

    python3 -B 04-computation/experiments/six_colour_core_doubling_20261012.py
    python3 -B -O 04-computation/experiments/six_colour_core_doubling_20261012.py

[Script](../../04-computation/experiments/six_colour_core_doubling_20261012.py),
[output](six_colour_core_doubling_20261012.json), and the declared large run
`python3 -B <script> --large` with [its output](six_colour_core_doubling_20261012_large.json)
(adds the radius-four ball). Ordinary and `-O` runs agree byte for byte. All coordinates are exact 40-vectors over `Q(zeta55)` in the
inherited engine; numerics (mpmath, 120 digits, PSLQ on the eight-term basis
of `K`) only propose roots, which are then squared exactly. The output keeps
the pair census with witnesses, the nine apex coordinates, their port pairs,
the empty apex adjacency, the 570-pattern result, the apex context's boundary
packet, the probe records and the eigenvalue control.
