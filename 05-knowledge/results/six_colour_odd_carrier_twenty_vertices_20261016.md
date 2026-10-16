# The cubic carrier and the six-colour boundary: what a coupled exterior must be, why twenty vertices, and what an equality tester would prove

2026-10-16, session opus-2026-10-16-S38 (second thread of the session). Continues clock_duality_oct15's
[elliptic bridge](hn_elliptic_bridge_20261016c.md) and [cubic carrier G7](elliptic_crossedges_20261016c.md).
Canonical statements: [THM-4661](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md),
[THM-4662](../../01-canon/theorems/THM-4662-the-fish-curve-has-no-point-over-the-real-quartic-field-so-the-fish-and-ei21-are-not-in-q-zeta55-or-any-odd-carrier.md).
Audit: independent adversarial referee, [report](../../04-computation/experiments/s38_audit_six_colour_20261016/REPORT.md):
no mathematical error in THM-4661 or THM-4662 (both reproduced independently, including the two-descent); errors in
remarks and in this note fixed below (MISTAKE-622).

**Owner's prompt (verbatim).** "continue as you left off and keep pushing the cubic-extension direction toward a
coupled unit-distance configuration with at least 18 exterior vertices, preserving an equal-endpoint six-colouring
while eliminating all separating colourings"

**Status.**
- **PROVED (THM-4661):** an exterior of at most **nineteen** points of an odd-degree extension preserves every old
  six-colouring (the inherited bound was seventeen; "at least eighteen" in the prompt becomes "at least twenty",
  and an exhaustive SAT search over the necessary conditions makes it **twenty-two**);
  the structure of a vertex-minimal failed exterior (lists of size five or six, low vertices span a forest of edges
  and odd cycles of length at least five, not three-colourable, not five-choosable, minimum degree five, a vertex of
  degree six); an equality tester is exactly six-chromatic and, in any field, decides the plane at seven colours;
  EI19 is excluded from every odd-degree extension of `Q(zeta55)` and a fish placement in any conjugation-stable
  odd-degree extension of `K = Q(zeta5, sqrt(-11))` needs a point of the genus-one curve `C3` over the real quartic field,
  as `L` already does; **EI21 carries the same genus-one curve**; and `C3` has no such point (**THM-4662**, a two-descent
  over `Q(sqrt5)`), so the fish and EI21 are excluded from `Q(zeta55)` and from every conjugation-stable odd-degree
  carrier over `K`.
- **FINITE-EXACT:** the arithmetic table of the parity step for every even palette from 4 to 40; non-residue
  certificates for `sqrt(-2)`, `sqrt(-35)` at `p = 661` in the cubic carrier; three-colourability of four probe
  contexts of `L' = Q(zeta55, theta)` built from the width-three unit and its conjugates (up to 13,287 vertices),
  which are, however, Cartesian products in disguise (section 5); clock_duality_oct15's circuit-based cubic
  assemblies (received this session) are three-chromatic as well.
- **OPEN:** a four-chromatic triangle-free unit-distance graph in `L'` (or in any odd-degree carrier over `K`); the
  realization fields of EI17 and EI21; `C3(F_4)`; the plane's chromatic number.

## 1. What the prompt asks for, in the repo's language

"Preserving an equal-endpoint six-colouring while eliminating all separating colourings" is the **equality tester**
of [existence_count_amplifiers](existence_count_amplifiers_20261011e.md): two non-adjacent ports `x, y` with a
positive equal-endpoint extension count and zero separating count. Colour symmetry makes the two-port relation
all-or-nothing: if one separating assignment dies, all do.

**THM-4661 statement 3.** An equality tester is exactly six-chromatic: if it were five-colourable, recolouring `y`
with the sixth colour would separate the ports. Two copies sharing `x`, rotated until the two images of `y` are at
unit distance (a chain of copies if `|x - y| < 1/2`), form a finite unit-distance graph that is not six-colourable.
This is the Moser spindle principle with six colours in place of three. Hence:

- an equality tester of any size, in any field, would prove `chi(R^2) = 7` and would in particular be a
  six-chromatic unit-distance graph, which nobody has exhibited (the record lower bound is de Grey's five);
- conversely `chi(R^2) = 7` gives one (an edge-critical seven-chromatic unit graph minus an edge);
- no equality tester lies in a five-colourable field plane, in particular not in `K` (residue colouring at 19).

So the prompt's target is not a size question: it is the six/seven question of
[the bridges note](constructive_interfaces_and_bridges_20261011.md) restricted to one carrier. What can be pushed is
the list of conditions any such configuration must meet, and how far the cubic carrier gets through them.

## 2. Nineteen exterior points never suffice

Setting (inherited from clock_duality_oct15): old points in a conjugation-stable field `F` without `sqrt(-3)`,
exterior points in an odd-degree extension `M`. Then (I1) an exterior point has at most one old neighbour, (I2) the
unit graph is triangle-free, (I3) two points have at most two common neighbours. For the cubic carrier take
`F = K' = K(theta)`, `M = L' = L(theta)`, `[L':K'] = 5`.

A vertex-minimal exterior `H` that kills an old colouring is `L`-critical for the lists `L(b) = [6]` minus the colour
of `b`'s old neighbour. The new ingredients compared with the seventeen-vertex theorem:

- **Degree-choosability (Borodin; Erdős–Rubin–Taylor).** The low vertices (`deg = |L|`) span components that are
  Gallai trees; triangle-free, so their blocks are edges and odd cycles of length at least five, and they have at
  most `(5/4)(ℓ - 1)` edges. In particular some vertex has degree above its list size, so degree at least six. The
  counting bound `1 + d + d(q-2)/2` at that vertex gives nineteen at once.
- **Nineteen is impossible.** At nineteen every degree-six vertex is tight: its neighbours have degree five, every
  second neighbour sees exactly two of them, and nothing lies further out. The low-degree count forces at least six
  high vertices. Rooting at one high vertex `z`, the other high vertices sit among the twelve second neighbours,
  and counting the edges inside that layer gives `2e(W_5) = 36 - 7|W_6|` with `|W_6| >= 5`: impossible
  (THM-4661, step (v)).

The same parity step excludes the extremal order for every even palette from 4 to 40
([table](six_colour_odd_carrier_checks_20261016.out)). The audit adds a stronger fact (E2 of the report): **no graph on
18 or 19 vertices is triangle-free, has at most two common neighbours per pair, and has minimum degree five**
(exhaustive SAT and backtracking, with the Clebsch graph recovered as a positive control at 16). Every set of at most
19 exterior points therefore spans a 4-degenerate unit graph, which gives the nineteen-vertex bound even without the
Gallai argument. Five-cores on 20 vertices do exist; none of 3,000 sampled ones passes the Gallai filter.

**Twenty-two, by exhaustive search (THM-4661 statement 5, FINITE-EXACT).**
[six_colour_exterior_twenty_sat_20261016.py](../../04-computation/experiments/six_colour_exterior_twenty_sat_20261016.py)
roots a SAT search at a high vertex of maximum degree. It enumerates the attachment patterns of its second
neighbourhood up to isomorphism and imposes every necessary condition above (triangle-free, two common neighbours at
most, degrees and lists, no low 4-cycle, the Gallai edge count, a high vertex), with non-three-colourability added
lazily. Results ([log](six_colour_exterior_twenty_sat_20261016.out)):

| `n` | patterns | verdict |
|---:|---:|---|
| 19 | 1 | UNSAT (re-derives the hand proof, without needing three-colourability) |
| 20 | 7 | UNSAT |
| 21 | 38 | UNSAT |
| 22 | 137 | 36 patterns (all with root degree six) have models; every model found is non-three-colourable |

So **every failed exterior has at least 22 vertices**: 21 field-escaping points never delete an old six-colouring.
At 22 the screen stops excluding, but the screen is weaker than the target. Each of the 36 models is
**five-choosable**, certified by an exact Alon–Tarsi coefficient
([alon_tarsi_five_choosability_20261016.py](../../04-computation/experiments/alon_tarsi_five_choosability_20261016.py)):
an orientation with outdegrees at most four whose monomial in `prod (x_u - x_v)` has nonzero coefficient. A failed
exterior is not five-choosable, so none of them is one. The obstruction that remains is choosability, not size: every
22-vertex model the search finds satisfies the necessary conditions N1-N5 and has choice number at most five.

## 3. The three-colour gate, now for arbitrary old sets

THM-4614's doubling lemma with (I1) in place of the nine-point census: if the exterior is three-colourable, give the
classes three disjoint colour pairs; each exterior point loses at most one colour of its pair. So a failed exterior
is **not three-colourable**, and also not five-choosable. In the cubic carrier, the very first object the direction
needs is therefore a four-chromatic triangle-free unit-distance graph inside `L' \ K'` with minimum degree five. In
`L` itself none is known and HYP-9247 conjectures none exists.

## 4. What odd degree cannot buy

Odd-degree extensions of `L` keep its quadratic subfields (`Q(sqrt5)`, `Q(sqrt(-11))`, `Q(sqrt(-55))`) and keep
every genus-one obstruction with an even-degree divisor:

| known triangle-free four-chromatic unit graph | status in `L` | status in odd-degree carriers (THM-4661 statement 4) |
|---|---|---|
| EI19 (19 vertices, rigid) | excluded: forced `sqrt(-2)`, `sqrt(-35)` | excluded from every odd-degree extension of `L`, including `L'` (also certified at `p = 661`) |
| fish (23 vertices, one-parameter) | needed `C3(F_4) != empty` (fish note (5)); **now EXCLUDED (THM-4662)** | needs the same condition in any conjugation-stable odd-degree `M ⊇ K`, hence **EXCLUDED** |
| EI17 (17 vertices, rigid) | excluded to coordinate height `1e63`; field OPEN | OPEN |
| EI21 (21 vertices, one-parameter) | OPEN in S24 (117-unit scan); carries the fish's `C3` (THM-4661 4(c)); **now EXCLUDED (THM-4662)** | **EXCLUDED** from every conjugation-stable odd-degree carrier over `K` |

The fish argument is the fish note's own: a point of `C3` over a field of odd degree over `F_4`, together with the
rational degree-eight fibre of `C3 -> P^1`, gives an `F_4`-rational divisor of degree one, hence a point over `F_4`.
The cubic `theta` is irrelevant to it.

**EI21 carries the same genus-one curve (new).** At the vertex `v_15` of EI21 (HoG 51378) four unit 4-cycles meet,
`(10,13,15,14)`, `(13,15,19,17)`, `(14,15,20,16)`, `(15,19,18,20)`. With `v_15 = 0` and `a, b, c, d = v_13, v_14, v_19,
v_20`, the parallelograms put `v_10 = a + b`, `v_17 = a + c`, `v_16 = b + d`, `v_18 = c + d`. The edges `10-18` and
`16-17` then say `|P + Q| = |P - Q| = 1` for `P = a - d`, `Q = b - c`. So `|P|^2 + |Q|^2 = 1` and `P ⊥ Q`, and with
`X = |P|^2` the three skew products give exactly the fish's square classes `X(1-X)`, `X(4-X)`, `(1-X)(X+3)` times
`11`. Every faithful numerical realization satisfies the two relations to `1e-15`, with `X` sweeping `(0.006, 0.998)`
along the family ([checks](six_colour_odd_carrier_checks_20261016.out)). S24 left EI21 open after a 117-unit scan; it
is now tied to the same single question `C3(F_4)` as the fish. It is excluded from `Q(zeta11, sqrt5)` and from every
carrier whose real field has odd degree over `Q(sqrt5)`. A census of the hub pattern (a hub, four neighbours, four
rhombus completions, two closing edges) finds it at vertex 15 of EI21 and vertex 12 of EI19, and nowhere in EI17 or
the fish (whose `C3` comes from its antipodal central rhombus instead). Of the four smallest known triangle-free
four-chromatic unit graphs, only EI17 is untouched.

**And `C3(F_4)` is empty (THM-4662, new).** [two_descent_q_sqrt5_fish_curve_20261016.py](../../04-computation/experiments/two_descent_q_sqrt5_fish_curve_20261016.py)
settles the question the fish note left open. `C3` is the two-covering of `E: y^2 = x(x-9)(x+7)` with class `(33, 33)`.
Over `F_4 = Q(sqrt5)(sqrt delta)`, `delta = (55 + 11 sqrt5)/2`, the rank of `E` is `rank E(Q(sqrt5)) + rank E^delta(Q(sqrt5))`.
A full two-descent over `Q(sqrt5)` gives two-Selmer groups of orders 4 and 16. Explicit points such as
`(484, 1573 - 3751 sqrt5)` and `((-65 + 19 sqrt5)/2, 160 + 400 sqrt5)` on the twist fill the whole Selmer group. So
`rank E(F_4) = 2`, and the descent image over `F_4` is the span of four explicit classes, which does not contain
`(33, 33)`. The controls reproduce the known ranks of six congruent-number curves over `Q` and of two curves over
`Q(sqrt5)`, and the Selmer group is reproduced with a second sampling seed. Consequence: **the fish and EI21 have no
placement in `Q(zeta55)` and none in any conjugation-stable odd-degree carrier over `K`**, the cubic carriers included.
Among the four smallest known triangle-free four-chromatic unit-distance graphs only EI17, which is rigid, survives,
and it is excluded from `L` only numerically. More generally, THM-4662 statement 5 forbids the nine-vertex hub
configuration `Gamma_9` itself (four rhombi at one vertex plus two closing edges) in `Q(zeta55)` and in every
conjugation-stable odd-degree carrier over `K`. Any four-chromatic exterior these fields could carry must avoid it.

A general remark explains part of the pattern. Suppose a graph can be built from one edge by repeatedly adding a vertex
with two placed neighbours, and it is placed inside an odd-degree extension `M` of a conjugation-stable field `F`
containing that edge. Then by (I1), step by step, it lies in `F`. So odd-degree carriers help only graphs whose
placements use genuinely new directions. The square-root conditions that come with an even-degree divisor argument
(as for `C3`) descend. Conditions without such an argument (EI17's rigid field, or the fish's genus-nine cover `C4`)
are not touched by this remark.

## 5. Probes of the cubic carrier, and why they are uninformative

[cubic_carrier_three_colour_probes_20261016.py](../../04-computation/experiments/cubic_carrier_three_colour_probes_20261016.py)
builds the width-three unit `v` of G7 exactly in the basis `zeta55^i theta^k` (norm one checked exactly), its ten
conjugates `sigma_m(v)`, `m ∈ (Z/11)^x`, and Cayley balls over them together with roots of unity. All unit pairs are
found by a three-embedding filter and confirmed by the exact norm in `L'`.

| probe | directions | radius | vertices | exact unit edges | triangles | three-colourable (verified) |
|---|---:|---:|---:|---:|---:|---|
| `mu10` + `v`-conjugates | 30 | 3 | 4,971 | 14,220 | 0 | yes |
| `mu22` + `v`-conjugates | 42 | 3 | 13,287 | 37,086 | 0 | yes |
| `mu110` + `v`-conjugates | 130 | 2 | 8,581 | 17,230 | 0 | yes |
| `mu10` + `v`-conjugates + `v`-conjugates times fifth roots | 110 | 2 | 6,161 | 12,430 | 0 | yes |

**These probes test nothing new, and this is the finding.** Two exact checks in the same script:

- every unit edge in all four probes is a generator step (zero non-generator unit pairs);
- the 55th roots of unity span a 40-dimensional `Q`-space, and adding the ten conjugates of `v` raises the rank to
  exactly 50; the fifth roots of unity and the ten rotated copies `v_m mu5` span exactly `4 + 10 * 4 = 44` dimensions.

So the conjugates of `v` satisfy no rational linear relation with the roots of unity: the probes are Cartesian because
their direction sets are linearly independent. The first three probes are induced subgraphs of the Cartesian product
of a cyclotomic Cayley graph with a bipartite lattice graph, the fourth of a product of eleven rotated pentagonal
lattices. Their three-colourings come from the factors: probes 1 and 2 map to `C5` and `C11` (residue maps plus the sum
of the `v`-coordinates), probe 4 to its pentagonal factors; probe 3's cyclotomic projection is three-coloured by SAT.
Coupling itself is short. **Rotated** first-power conjugates `zeta^k v_r` satisfy integer relations of length 20, 70,
80 and 88 for the triples `(1,3,4)`, `(2,9,3)`, `(1,10,2)`, `(1,2,3)` (kernel-LLL on the saturated relation lattice;
each relation is verified exactly). The length-20 one is two 10-step unit walks with a common endpoint (audit E1; an
earlier version of this section reduced an unsaturated lattice and reported lengths about 1000 times too long,
MISTAKE-622). clock_duality_oct15's C2 ([cubic cross-edge packets](cubic_crossedge_packets_20261016d.md), received
during this session) gives another circuit, of mass 64, in the powers `v^i xi^j`. Their five assemblies are
three-chromatic (rows 4 and 5 are built from that circuit, rows 1 to 3 are power, mixed and conjugate banks), and so
are the audit's walk banks built from the 20-step circuit (up to 1,279 vertices). The probes here are uninformative by
construction, while the circuit-based assemblies are three-colourable too: the three-colour gate of section 3 is
where the direction stops.

**First-moment calibration (HEURISTIC).** If an exterior behaved like a random graph with `n` vertices and `m` edges,
each with one forbidden colour, the expected number of extensions would be about `6^n (5/6)^(m+n)`. This falls below
one only at average degree about 17.7. The known triangle-free four-chromatic unit graphs have average degree 3.6 to
3.7 (EI17: 31 edges on 17 vertices; fish: 42 on 23), and the probes above stay below 6. In the plane the number of
unit distances among `n` points is at most of order `n^(4/3)`, and the dense finite unit graphs available in these
fields (cyclotomic balls) keep coming out three-colourable (HYP-9247). The proved bound of twenty is therefore a
floor, not an estimate. A coupled exterior, if one exists, is either large or
very special.

## 6. The board after this session

| Requirement on a coupled exterior in `L' \ K'` | Status |
|---|---|
| at least 20 vertices; at least 22 by exhaustive search over the necessary conditions | PROVED (20); FINITE-EXACT, solver verdicts (22) |
| minimum degree 5, degree 6 where there is no old contact, some degree-6 vertex | PROVED necessary |
| low vertices span a forest of edges and odd cycles (length >= 5) | PROVED necessary |
| not three-colourable (four-chromatic triangle-free unit graph in `L'`) | PROVED necessary; none known; EI19, fish, EI21 excluded (THM-4661, THM-4662); EI17 open |
| not five-choosable | PROVED necessary; it kills every 22-vertex model of the other conditions found (Alon–Tarsi) |
| all separating assignments killed (equality tester) | would prove `chi(R^2) = 7` (THM-4661 statement 3; the converse only gives a tester in some field) |

The canonical hostile for this board is a 22-vertex model of the search: it meets every necessary condition N1-N5,
yet is five-choosable and so cannot kill a colouring. (Five-cores of 18 or 19 vertices do not exist at all.) The
least-used sidecar is choosability: no earlier screen in the repo used it, and it is the condition that separates "has
a five-core" from "can kill a colouring".

## 7. Reproduction

```text
python -B 04-computation/experiments/six_colour_odd_carrier_checks_20261016.py
python -B 04-computation/experiments/cubic_carrier_three_colour_probes_20261016.py
python -B 04-computation/experiments/six_colour_exterior_twenty_sat_20261016.py 19 20 21 22
python -B 04-computation/experiments/alon_tarsi_five_choosability_20261016.py 05-knowledge/results/six_colour_exterior_models_n22_20261016.json
python -B 04-computation/experiments/two_descent_q_sqrt5_fish_curve_20261016.py
```

Outputs: [checks](six_colour_odd_carrier_checks_20261016.json) ([log](six_colour_odd_carrier_checks_20261016.out)),
[probes](cubic_carrier_three_colour_probes_20261016.json).
