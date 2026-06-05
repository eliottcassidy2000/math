# Concept Map — Tournament Parity Research

**Purpose:** Complete, structured database of every mathematical concept, object, technique, and connection in this project. Organized for rapid lookup by future Claude instances. Created by kind-pasteur-2026-03-07-S34.

**Last updated:** kind-pasteur-2026-03-07-S34 (continued session — 10+ commits, 30+ connections cataloged)

---

## I. CORE MATHEMATICAL OBJECTS

### Tournaments
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Tournament** T | Complete directed graph on n vertices | Every pair has exactly one arc | Everywhere |
| **Tournament Analysis** | Pipeline from pairwise/geometric data to binary comparators/gauges/switchboards, tournaments, invariants, and wall-crossing paths | Rankers summarize objects; switchboards summarize relations; edge-local/switch lifts preserve cyclic structure; THM-372/373/374 give the finite switchboard, runner-clock, and semicircle base cases; HYP-1970 corrects H as an exact half-turn detector but not a scalar max-gap meter; LRC hard rows now use two-clock corridor overlays, endpoint-labelled pressure DAGs, H as phase-loneliness entropy, marked H-loneliness masks, an arc-criteria metric vector, a ranker/entropy/pressure gauge trichotomy, endpoint-aware gauge fingerprints, S506c's open-arc/local-crowding/pressure/observer gauge bundle, a second-order pair-cell operation-grid vector, rooted/marked A000568-style quotient models, S510's sparse clock-image refinement, S511's operation-weighted danger fibers, observer-source reachability, threshold-decorated A000568 class fibers, an operation arc zoo splitting ledgers from phase perturbations and marked bundles, S513's add/multiply gauge stack over the odd-core/dyadic denominator grid, S514's pressure-core refinement of that stack, THM-384/HYP-1986 compactified source-gap forcing, HYP-1987's arc-confined source menu, THM-385/HYP-1988 observer-score blocker stratification, HYP-1989's n=145 CRT aperture route, THM-386/HYP-1990's two-gap race flow, and S520/HYP-1992's observer-source n=18 gate/debt audit | HYP-1931, HYP-1932, HYP-1940, HYP-1941, HYP-1942, HYP-1951, HYP-1952, HYP-1964, HYP-1965, HYP-1966, HYP-1967, HYP-1968, HYP-1969, HYP-1970, HYP-1971, HYP-1972, HYP-1973, HYP-1974, HYP-1975, HYP-1976, HYP-1977, HYP-1978, HYP-1979, HYP-1980, HYP-1981, HYP-1982, HYP-1983, HYP-1984, HYP-1985, HYP-1986, HYP-1987, HYP-1988, HYP-1989, HYP-1990, HYP-1991, HYP-1992, THM-372, THM-373, THM-374, THM-381, THM-382, THM-383, THM-384, THM-385, THM-386, S23, S454, S455, S480, S471, S493, S495, S502b, S504, S505, S506, S507, S509, S510, S511, S512, S513, S514, S516, S517, S518, S519, S520, S26 |
| **H loneliness meter** | Reading `H(T)` along runner-clock tournaments as a scalar shadow of circular geometry | Low `H` detects unanchored open-semicircle bunching; high `H` often tracks LRC-style local separation, but LRC also needs marked stationary vertex, safe-gap mask, and pressure/deletion data | HYP-1971, HYP-1970, HYP-1951, THM-370, THM-374, S506, `07-reflections/h-as-loneliness-meter-s506.md` |
| **LRC arc-criteria metric vector** | Panel of Tournament Analysis gauges for LRC loneliness | Phase-half `H` tracks global spread; close-sector tie rate tracks safe-gap count; local-moat and safe-deficit marked origin scores track origin clearance; pressure gauges are read by strict SCCs and directed triangles rather than completed `H` | HYP-1972, HYP-1971, HYP-1967, S506, `04-computation/lrc_arc_criteria_loneliness_s506.py` |
| **LRC arc-gauge bundle** | Parallel panel of completed tournament gauges for LRC selected rows | Phase/open-arc gauges record global spread, close/danger/kinetic switches record local `1/n` crowding, relief gauges record pressure peeling, observer outdegree records the marked bracket, and endpoint debt must stay labelled before projection | HYP-1975, HYP-1972, HYP-1973, HYP-1974, S506c, `04-computation/lrc_arc_gauge_zoo_s506.py` |
| **LRC pair-cell operation grid** | Second-order Tournament Analysis whose vertices are unordered runner pairs | Pair-cell danger-deficit tie rate tracks close-pair count; dyadic-row, odd-core, same-chain, and product-sum-defect criteria label the `x+2`/`x*2` arithmetic branch structure rather than serving as time-loneliness meters | HYP-1976, HYP-1972, HYP-1963, HYP-1965, HYP-1966, S509, `04-computation/lrc_pair_cell_operation_grid_s509.py` |
| **LRC over A000568 quotient** | Projection tower from runner time movies to observer-pointed tournament classes and then unmarked tournament isomorphism classes counted by A000568 | Exact small-cell audits show LRC safe/unsafe status is not a function of the unmarked class, and often not of the pointed class; the plausible analogy is a section/projection-defect problem over `G_N`, with endpoint/gap labels and HYP-1976 pair-cell data as fiber coordinates | HYP-1977, HYP-1976, HYP-1975, HYP-1967, HYP-1970, S509, `04-computation/lrc_a000568_iso_analogy_s509.py` |
| **LRC observer-coupled perspective gap** | Reinterpretation of rooted tournament perspective counts as the first observer lift above unmarked A000568 classes | S586 corrects the old `P(5)=48` versus `U(6)=56` reading: all six-classes are reachable under the full add-a-vertex extension, but the root is unused, so the gap is not orphan classes but missing incident-word/threshold payload. HYP-2120 identifies the exact source-perspective slice; HYP-2121 explains why the full rooted slice still needs incident threshold fibers. | HYP-2121, HYP-2120, HYP-1977, HYP-1978, HYP-1981, HYP-2113, HYP-2109, S586, S582, `04-computation/tournament_perspective_observer_coupling_s586.py` |
| **VT trienerment polygon rigidity** | Vertex-transitive local geometry split into cyclic polygon, dihedral bracelet, and abstract Cayley cascade | S589 repairs the slogan "VT trienerment iff regular polygon": cyclic sharp-transitivity forces a regular polygon, but dihedral VT allows imprimitive bracelets and abstract VT tournaments allow noncyclic Cayley meshes. In `Z/N` point-set audits through `N<=18`, `83` dihedral-VT sets split into `31` regular polygons and `52` bracelets, first `P=(0,1,3,4) in Z/6` with gaps `(1,2,1,2)`. Cayley `F21` is vertex-transitive and regular but has no element of order `21`, so rigidity cascades through a nonabelian relator mesh. HYP-2124 is the exact cyclic/AP unit-clock witness; HYP-2125 is the broader taxonomy. LRC should keep local rooted view plus cascade law plus source/endpoint/denominator labels. | HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, THM-052, THM-389, THM-381, S589, `04-computation/vt_trienerment_polygon_rigidity_s589.py` |
| **LRC orbit-sheaf monodromy rigidity** | Rigidity as a labelled section over orbit quotients, with defects measured by boundary maps and mixed fibers | S590 fuses HYP-2126's orbit-rigidity and fixed-point-propagation axes: local/global are only two types among static, dynamic, boundary, quotient, gluing, monodromy, spectral, isostatic, pincer, and automaton rigidity. AP unit clocks under doubling give the key seam: odd `n` keeps all units, even `n` sends every unit witness to gcd-2 boundary residues; at `n=14`, each doubled residue has one unit lift and one gcd-2 nonunit lift. HYP-2128 supplies the complementary additive `+2`/triangular face where `2n-1` is the odd-square root of `8*C(n,2)+1`. Dihedral point-set monodromy is trivial polygon or order-2 bracelet; rooted tournament quotients show fixed/parent labels mix when cross fibers are forgotten. | HYP-2127, HYP-2128, HYP-2126, HYP-2125, HYP-2124, HYP-2101, HYP-2020, HYP-1977, S590, `04-computation/lrc_orbit_sheaf_monodromy_s590.py` |
| **LRC observer-source marked tournament** | Mark the stationary observer and orient observer-to-runner arcs by the LRC threshold `||v_i t|| >= 1/n` | THM-381: observer loneliness is exactly the observer being a source; THM-385: observer indegree is exactly the blocker count, so source, almost-source, score layers, and side-defect/2-king repair predicates form an exact distance-to-source fiber around HYP-1981 | THM-381, THM-385, HYP-1988, HYP-1987, HYP-1981, HYP-1977, HYP-1979, S517, S511, `04-computation/lrc_observer_predicate_zoo_s517.py`, `04-computation/lrc_observer_source_tournament_s511.py` |
| **LRC sparse marked chamber walk** | Sparse-clock refinement of the A000568 quotient frame | Initial and hard-row runner clocks visit tiny arithmetic subsets of the A000568 chamber space; observer-marked safe witness samples identify target hits, while HYP-1977 explains why endpoint/gap fibers cannot be discarded | HYP-1979, HYP-1978, HYP-1977, HYP-1951, HYP-1967, HYP-1975, S510, `04-computation/lrc_a000568_marked_chamber_s510.py` |
| **LRC operation-weighted danger fiber** | Runner-level arc criteria induced from pair-cell operation-grid labels and current `1/N` danger mass | Static `x+2`/`x*2` labels are fixed coordinates, but additive, dyadic, multiplicative, and product-sum danger-weighted gauges become moving loneliness metrics over the marked A000568 projection fiber | HYP-1980, HYP-1979, HYP-1978, HYP-1977, HYP-1976, HYP-1972, S511, `04-computation/lrc_operation_grid_arc_criteria_s511.py` |
| **LRC threshold-decorated class fiber** | A000568 tournament class plus circular/tie-wall stratum, observer mark, and `1/n` threshold colors on gaps or pair-cells | THM-382: bounded `n=3,4` audits separate only after threshold decoration; THM-383: equality walls expand the bounded half-turn class menu; THM-384: the observer is lonely exactly when both observer-adjacent gaps are threshold-long, reducing the remaining proof to HYP-1986 compactified source-gap forcing | THM-382, THM-383, THM-384, HYP-1986, HYP-1982, HYP-1976, HYP-1951, S516, S512, `04-computation/lrc_gap_threshold_proof_route_s516.py`, `04-computation/lrc_iso_class_constraint_s512.py` |
| **LRC restricted quotient stack** | Stack of LRC-to-tournament maps: raw circular body, source-deleted target, observer-source mark, gap/apex threshold fiber, kinetic gap-flow, blocker-deficit shadow, and endpoint-pressure/channel labels | S535 shows raw phase classes remain mixed in small exact audits, while source-deleted and threshold-decorated maps certify all sampled systems; new metrics include source codimension, label tax, purity curvature, kinetic torsion, blocker viscosity, compactification index, pressure survival rank, monodromy defect, and quotient-lens resolvent | HYP-2020, HYP-2018, HYP-2019, HYP-1981, HYP-1982, HYP-1986, HYP-1988, HYP-1999, HYP-2001, HYP-2008, HYP-2009, HYP-2015, S535, `04-computation/lrc_restricted_tournament_mappings_s535.py` |
| **LRC section-boundary functors** | Non-runner vertex choices for LRC tournaments: fixed circle sections, fixed boundaries, void sectors, events, cover arcs, residue channels, Fourier modes, and proof obligations | S539 extends HYP-2022's sector-duality result and HYP-2023's event-media result: danger/occupancy-colored section tournaments are pure in bounded `n=4..6` open cells; boundary-flux and void-pressure become pure at `q=2n`; AP wall misses indicate compactification debt. The associated protocol asks future agents to justify the tournament vertex set before defaulting to runners. | HYP-2024, HYP-2023, HYP-2022, HYP-2021, HYP-2020, HYP-2018, THM-382, THM-383, THM-384, THM-387, S539, `04-computation/lrc_section_boundary_functors_s539.py` |
| **LRC support-flow / cut-flow duality** | Directed-graph reframing of LRC safe measure and danger cover | HYP-2025 identifies full-support resonance as a nowhere-zero flow on the speed dipole. S540 extends this to all Fourier supports as NZ flows on sub-dipoles, while the dual danger-cover cell cycle exposes lonely intervals as zero-flow cuts. AP walls have nowhere-zero open-cell coverage but empty strict endpoint cores; true counterexamples would need both support-flow cancellation and a nonpeeling cover-flow core. | HYP-2026, HYP-2025, HYP-2022, HYP-2024, HYP-2001, THM-379, THM-380, S540, `04-computation/lrc_flow_cut_support_s540.py` |
| **LRC operation arc zoo** | Runner-level Tournament Analysis switches using phase, endpoint, two-neighbor, dyadic/odd-core, product-sum, Goldbach/Lemoine, danger-grid, and bundle votes | Scalar endpoint and two-neighbor switches are ledgers that often collapse to transitive H; phase-overridden operation gates perturb the half-turn clock; majority bundles produce separate marked-fiber coordinates to read by score histograms, observer score, strict SCCs, and directed triangles | HYP-1983, HYP-1982, HYP-1980, HYP-1981, HYP-1976, HYP-1979, S511, `04-computation/lrc_operation_arc_zoo_s511.py` |
| **Double-round-robin vertex doubling** | Replace each vertex by a two-vertex fiber and each old edge by a `2 x 2` block whose winner and loser take complementary perfect matchings | THM-378: these are signed SC/voltage lifts; all quotient pair records are `2-2`, clone scores are universally `(n-1)^n,n^n`, and sheet-flip gauge classes are triangle parity vectors of dimension `binom(n-1,2)` | THM-378, T269, T508, S501c, `07-reflections/sc-blowup-and-twin-gaining.md` |
| **LRC first-even bridge seam** | First doubling `m -> 2m` in Lonely Runner denominators, with inherited odd-parent units and new nonunit bridge room | `n=14` has six even bridges; `n=18` has two triadic bridges `(6,12)`; gap-debt product and pressure peeling are key ledgers; hard positive gaps are short corridors through half-turn tournament-clock cells with boundary alignment ratios `3/7` and `3/11`; H drops at the first gate then plateaus under further dyadic scaling; S520 verifies the n=18 no-gate branch is killed by unit observer-source targets while the `9->18->36` ladder halves gap and doubles endpoint debt under the THM-386 gap-flow background | HYP-1905, HYP-1910, HYP-1920, HYP-1942, HYP-1950, HYP-1952, HYP-1967, HYP-1969, HYP-1981, HYP-1986, HYP-1987, HYP-1988, HYP-1989, HYP-1990, HYP-1991, HYP-1992, S453, S455, S490, S492, S493, S502b, S505, S520 |
| **LRC n=18 square-core packet** | Normalized endpoint-depth packet for the `18=2*3^2` first-even row | Row-parent/gate/double-gate ladders are dyadic translates of `96*(0,2)+16*(0,3)+64*(4,2)`; bridges `6` and `12` are locally identical but globally signed | HYP-1952, S493 |
| **LRC endpoint-pressure certificate** | Two-layer finite certificate combining endpoint protection cores with owner-pressure graphs | THM-379/THM-380: a nonempty pressure-realized endpoint core forces an owner cycle; the sieve/core/DAG trilemma narrows counterexamples to labelled pressure SCCs; S504 refines the induction: source blockers charge sink/private endpoint layers | THM-357, THM-359, THM-366, THM-369, THM-379, THM-380, HYP-1960, HYP-1961, HYP-1968, S503, S504 |
| **LRC tournament gauge trichotomy** | Menu of arc-assignment criteria for LRC snapshots | Ranker gauges identify lonely coordinates but H collapses; phase/entropy gauges make H a global loneliness metric; pressure gauges expose blocker-debt hierarchy via score/SCC/source-sink shape | HYP-1973, S506, `07-reflections/lrc-tournament-gauge-lab-s506.md` |
| **LRC A000568 marked quotient** | Reading LRC as a rooted/marked tournament-isomorphism problem over the ordinary A000568 quotient | Ordinary phase tournaments give the unrooted base; the stationary observer roots the class; endpoint masks and pressure labels form marked fibers; hard ladders can keep the same coarse rooted shadow while endpoint debt doubles | HYP-1978, HYP-1977, HYP-1973, HYP-1900, HYP-1901, S507, `07-reflections/lrc-a000568-isoclass-analogy-s507.md` |
| **LRC add/multiply gauge stack** | Denominator-level Tournament Analysis over the `x+2` odd-core and `x*2` dyadic grid | Addition contributes Goldbach/Lemoine representation fibers; multiplication contributes divisor/dyadic branches; product-sum equations label operation collisions; scalar gauges are transitive ledgers, while majority conflicts produce nontrivial H, c3, SCC, and edge-flip fingerprints; S514 shows the full hard-row stack has shape but the predicted bad conjunction fails at coarse pressure, so the pressure coordinate should be THM-380 owner-compatible endpoint cores | HYP-1985, HYP-1984, HYP-1983, HYP-1982, HYP-1981, HYP-1980, HYP-1979, HYP-1978, HYP-1977, HYP-1976, HYP-1975, S514, S513, `04-computation/lrc_add_mult_gauge_stack_s513.py`, `04-computation/lrc_three_layer_stack_audit_s514.py` |
| **LRC parity-ladder proof program** | HYP-2091 clean/wall simplex-polygon parity split overlaid with `2n-1` unit/nonunit shell strata and D/U/N obligation labels | S576 routes rows into clean unit, clean composite, wall unit, and wall composite proof lanes. The n=14 row is clean geometry but composite clock burden (`C=27=3^3`, no tie pairs, unit/nonunit shells `9/4`, D/U/N total `34`, `190` converse-merged round nodes), so the immediate target is a labelled fibre table over the round/converse seam with private pivots, gcd descent, and THM-397 endpoint owners. | HYP-2092, HYP-2091, HYP-2090, HYP-2089, HYP-2088, HYP-2087, HYP-2086, HYP-2083, THM-397, S576, `04-computation/lrc_parity_ladder_proof_program_s576.py` |
| **LRC n=14 unit-spine slack fibre** | Dihedral quotient of the strict sub-edge n=14 proof ledger after the nine unit shells modulo `27` are forced | S576 fixes the canonical unit spine `(1,2,4,5,7,8,10,11,13)` and scans four nonunit/zero-residue slack runners. S578 through slack `<=42` finds `531` full D/U/N covers, only AP and `V*` at the floor, and one-unit lifts through `81` all cheap. HYP-2162/THM-407 folds the raw `C=27` shell layer to gcd strata `{1,3,9}`; HYP-2163 supplies the efficient proof-pipeline route; HYP-2164 supplies the pinch certificate. S609 extends the canonical slack fibre itself through `81`: `9506` full covers, `9504/9506` cheap-pair exits, `2/2` no-cheap rows positive measure, only AP and `V*` measure-zero, and zero residual. S610 composes these into an 11-atom quotient tower whose remaining seam is lift/CRT conservativity. | HYP-2166, HYP-2165, HYP-2164, HYP-2163, HYP-2162, HYP-2100, HYP-2099, HYP-2098, HYP-2097, HYP-2096, HYP-2095, HYP-2094, HYP-2092, HYP-2088, HYP-2087, HYP-2084, THM-397, THM-407, S610, S609, S578, S576, `04-computation/lrc_n14_res27_quotient_tower_s610.py`, `04-computation/lrc_n14_res27_fixed_bridge_s609.py`, `04-computation/lrc_n14_dihedral_spine_pivots_s576.py`, `04-computation/lrc_n14_unit_spine_exchange_s578.py` |
| **LRC n=14 carry-fiber glue** | Lift/CRT descent datum over HYP-2166's `Res_27` quotient tower, with `v=r+27k` and `v == r-k (mod 14)` | S611 shows why HYP-2164 is the base quotient, not the whole lifted theorem: all `36` AP/`V*` unit scalar lifts remain floor, but their least-positive shadows are only `3` floor and `33` strict.  The three floor shadows are AP, `V*`, and nonprimitive `2*AP`, matching HYP-2164.  AP/`V*` local `+27` perturbations and S612 bounded carry fibers are strict away from zero carry.  S654 adds that the same carry controls parity (`v mod2=r+k`) and apex divisibility (`14|v iff k=r mod14`), so any new floor lift must use a globally coherent carry pattern or route to HYP-2165 owner/Cprime certificates. | HYP-2230, HYP-2167, HYP-2166, HYP-2165, HYP-2164, HYP-2163, HYP-2162, HYP-2161, THM-401, THM-407, S654, S612, S611, S610, S609, S608, `04-computation/lrc14_parity_carry_bridge_s654.py`, `04-computation/lrc_n14_carry_conservativity_s611.py`, `05-knowledge/results/lrc14_parity_carry_bridge_s654.out`, `05-knowledge/results/lrc_n14_carry_conservativity_s611.out` |
| **LRC n=14 parity-carry bridge** | Even/odd number-theory shadows and the `n=14` apex obstruction share the carry coordinate in the `C=27` lift | For `v=r+27k`, S654 records exact identities `v mod2=r+k` and `v mod14=r-k`, plus `14|v_i+v_j iff k_i+k_j=r_i+r_j`.  This imports HYP-2218's even/odd same-pair diagonal `p=7 -> (14,21)` into the LRC carry theorem: hitting the multiple-of-14 branch is not separate from parity toggling.  Minimal apex carries over AP and `V*` are all strict (closest `28/365` and `2/25`), and one/two-toggle apex neighborhoods have no near-floor or below-floor rows. | HYP-2230, HYP-2167, HYP-2218, HYP-2229, HYP-2220, THM-401, S654, S642, S653, `04-computation/lrc14_parity_carry_bridge_s654.py`, `05-knowledge/results/lrc14_parity_carry_bridge_s654.out`, `07-reflections/lrc14-parity-carry-bridge-s654.md` |
| **LRC n=14 sin/cos even-odd wall carrier** | Typed parity bridge from odd active walls to even derivative/slack layers | S655 reads `sin` as odd boundary data, `cos` as adjacent derivative/slack, and `cot` as a scalar pole ledger that forgets composite side channels, downstream of HYP-2230's carry coordinate tying parity to the `mod 14` obstruction.  The Goldbach/Lemoine duplicate branch `7+7=14`, `7+2*7=21` is a scalar shadow, while AP and `Vstar` hit the LRC `n=14` floor through off-diagonal odd pairs `(1,13),(5,9),(3,11)` and mirrors.  Both rows have exact `M=1/14`, the same six best times, and `C=27` gcd-shell counts `{1:9,3:3,9:1}`; the single-swap scan has `195` rows, `0` below, and only `12->24` tight.  The next target is an odd-wall/gcd-shell/carry-owner no-leak lemma. | HYP-2231, HYP-2230, HYP-2229, HYP-2222, HYP-2219, HYP-2218, HYP-2167, HYP-2164, HYP-2116, T724, S655, `04-computation/lrc_sin_cos_parity_carrier_s655.py`, `05-knowledge/results/lrc_sin_cos_parity_carrier_s655.out`, `07-reflections/lrc14-sin-cos-even-odd-wall-carrier-s655.md` |
| **LRC n=14 information bottleneck** | Sufficient-forgetting reading of HYP-2164..2170: compress raw rows only along channels independent of the floor predicate, then reattach relevant side information | S612 separates address entropy from relevance entropy. Raw `Res_27` row address entropy drops `23.3102 -> 3.4594 -> 1.0000` bits across the quotient tower and floor atoms, but rare proof predicates carry the action: pinch status has `0.0016` bits, owner route has `0.9981`, scalar carry route has `0.4138`, and the scalar carry indicator perfectly determines floor-vs-strict shadow route in S611 probes. The proposed sufficient statistic is `(Res_27 proof atom, owner route, carry cocycle, Cprime window)`, now explicitly including incoming HYP-2168's C2b `p_0>0` / `2n-1` tick window and HYP-2170's AP/`V*` bounded tight-family evidence. | HYP-2171, HYP-2170, HYP-2169, HYP-2168, HYP-2167, HYP-2166, HYP-2165, HYP-2164, HYP-2161, HYP-2157, HYP-2156, THM-401, THM-406, THM-407, S612, S611, S610, `04-computation/lrc_n14_information_bottleneck_s612.py`, `05-knowledge/results/lrc_n14_information_bottleneck_s612.out` |
| **LRC dimension-descent residue depth** | Cross-n proof-burden ledger for reducing n=14 toward n=13, n=12, n=11 by `C=2n-1` arithmetic rather than runner count | S613 audits the unit-shell quotient, D/N gates, and pair-sum pinch route for `3<=n<=14`.  It recovers n=14's `340928` unit rows and `27733` survivors exactly, then shows the real descent is `27=3^3 -> 25=5^2 -> 23 prime`, while `23 -> 21=3*7` grows the burden (`420 -> 3415` survivors) despite smaller n.  Floor atoms recover AP/`2*AP`, n=14 `V*`, and the p0-collapse sporadics `(1,3,4,7)`, `(1,3,4,5,9)`, `AP[6->12]`, and `(1,4,5,6,7,11,13)`; HYP-2169/HYP-2174 explain the doubled mod-3 subfamily. | HYP-2175, HYP-2174, HYP-2172, HYP-2169, HYP-2168, HYP-2167, HYP-2166, HYP-2156, HYP-2153, THM-401, THM-407, S613, S608, `04-computation/lrc_dimension_descent_salience_s613.py`, `05-knowledge/results/lrc_dimension_descent_salience_s613.out` |
| **Unit distance n=22 coimage frontier** | Small planar unit-distance analogue of the LRC sufficient-statistic problem | S614 records the current published frontier `60 <= u(22) <= 61`: the `F`-free graph coimage reaches `62`, but totally-unfaithful/geometry side channels kill both `62` candidates. Any `61`-edge solution has min degree `4` or `5` and must extend a `56/57`-edge `21`-core. The route tournament ranks Moser search, dense-core extension, and unfaithful obstruction mining above raw graph-only/grid routes. It also corrects the incoming S599/HYP-2170 `49` reading as a triangular/Harborth lattice baseline, not the planar optimum under the published `60`-edge lower bound. | HYP-2176, HYP-2171, HYP-2170, S614, S612, `04-computation/unit_distance_n22_tournament_lrc_s614.py`, `05-knowledge/results/unit_distance_n22_tournament_lrc_s614.out`, `07-reflections/unit-distance-n22-tournament-lrc-grid-disproof-s614.md` |
| **Unit distance Moser-carrier counting speedups** | Frontier-gain counting and dense-core extension ledgers for the `n=22` unit-distance search | S617 keeps the S614 rank-4 Moser unit shell (`18` directed vectors) and updates child counts by `E(S union {q})=E(S)+gain(q)` instead of recounting all pairwise distances. A width-1200 beam recovers `60` edges at `n=22` with `211.4x` edge-count-check speedup. The retained `57`-edge `21`-cores only gain `3`, and retained `56`-edge cores only gain `4`, so the Moser beam explains why it reaches the known `60` lane but not `61`. The useful quotient preserves unit-edge increments and destroys embedding/unfaithful side channels, so the next proof route is a `21`-core extension solver plus obstruction library. | HYP-2188, HYP-2176, HYP-2184, S617, `04-computation/unit_distance_moser_beam_speedups_s617.py`, `05-knowledge/results/unit_distance_moser_beam_speedups_s617.out`, `07-reflections/unit-distance-moser-beam-speedups-s617.md` |
| **Unit distance impairment spectroscopy** | Deliberately damage construction/proof channels to expose load-bearing side information before scaling | S622 runs controlled impairments on the Moser carrier: direction-pair dropout, observed direction-support masks, and gain-cap throttling. Exact values through `n=14` are recovered by the full width-260 beam; the best `n=14` witness uses seven direction pairs, dropping pairs `0,1,4` is still exact while dropping any of `2,3,5,6,7,8` costs one edge, and gain cap `4` is the first cap that recovers exactness through `n=14`. The resulting techniques are impaired-carrier spectroscopy, direction Helly certificates, shadow-price ledgers, gain-order sieves, and obstruction-first beams. | HYP-2194, HYP-2188, HYP-2176, HYP-2152, T691, T686, S622, `04-computation/unit_distance_impairment_lab_s622.py`, `05-knowledge/results/unit_distance_impairment_lab_s622.out`, `07-reflections/unit-distance-impairment-spectroscopy-s622.md` |
| **Unit distance impairment atlas** | Small controlled failures that identify load-bearing construction/proof side channels | S623 extends HYP-2194's direction-pair spectroscopy by deliberately impairing triangular and Moser carriers through width, policy, gain, and canonicalization before scaling. Moser target `14` jumps from `28/29` edges at widths `1/3/10` to `33` at width `30`, exposing future-gain shape as a retained invariant. Moser target `12` loses under sprawl bias (`27 -> 24`); direction-drop jackknife loses one edge for six of nine Moser antipodal directions at target `10`; gain ceilings reveal high-gain moves as load-bearing (`12,21,25,27` for ceilings `1,2,3,4`). The technique program is damage-response ledger, direction jackknife certificates, gain-threshold extension solver, deletion-core resilience, and orbit-budget accounting. | HYP-2196, HYP-2194, HYP-2188, HYP-2176, T692, T691, S623, S622, `04-computation/unit_distance_impairment_atlas_s623.py`, `05-knowledge/results/unit_distance_impairment_atlas_s623.out`, `07-reflections/unit-distance-impairment-atlas-s623.md` |
| **LRC / unit-distance channel ledger** | Cross-problem bridge that treats state-local frontiers as lossy unless load-bearing side channels are retained | S624 pairs an LRC witness-orbit jackknife over survival masks on `Z/(2n-1)` with a lightweight Moser direction/gain impairment. On LRC n=14, AP/Vstar/doubled-AP release `18/6/2` residues across gcd `1/3/9` channels when an orbit is omitted; redundancy prices are `18/6/8` because the two gcd-9 residues are killed four times. Prime `C=23` has one raw shell orbit, pushing side-channel scarcity into lift/carry data. On Moser target `9`, direction-drop loss histogram is `{0:3,2:6}` and critical directions have usage `3`, loss `2`, price `6`; gain ceilings `1,2,3,4` give `9,15,16,18` edges. Channel Tournament Analysis over proof obligations is transitive under proof and scaling gauges but flips `19/36` edges, so proof relevance and scaling value are different quotients. | HYP-2197, HYP-2196, HYP-2194, HYP-2188, HYP-2187, HYP-2175, T693, T692, S624, S623, S599w-x, `04-computation/lrc_unit_distance_channel_ledger_s624.py`, `05-knowledge/results/lrc_unit_distance_channel_ledger_s624.out`, `07-reflections/lrc-unit-distance-channel-ledger-s624.md` |
| **Unit-distance tournament Hamiltonian path flop** | Points-to-tournament mapping that asks whether Rédei's guaranteed Hamiltonian path is carried by unit pairs, nonunit pairs, or a base-order artifact | S625 tests unit-flip and nonunit-flip tiling rules on triangular rows through `n=22` and Moser rows through `n=14`, complementing HYP-2201's triangular-lattice traceability theorem. The graph-level unit Hamiltonian path persists in every tested row. The nonunit complement first has a Hamiltonian path at `n=6`, fails at compact hexagon `n=7` because the center is complement-isolated, then reappears from `n=8`. The canonical unit-flip directed all-unit path flops at `n=7`, so the first flop is order-level rather than geometric. Snake-base tiling restores an all-unit directed path whenever the unit graph has one; HP unit-step histograms are the richer invariant. | HYP-2202, HYP-2201, HYP-2197, HYP-2196, HYP-2188, T694, S625, `04-computation/unit_distance_tournament_hampath_s625.py`, `05-knowledge/results/unit_distance_tournament_hampath_s625.out`, `07-reflections/unit-distance-tournament-hampath-flop-s625.md` |
| **Unit distance Moser unit-spine traceability** | Separating non-lattice Moser traceability from point-order flip-gauge artifacts | HYP-2201 gives the lattice traceability theorem, and HYP-2202 separates graph-level unit HPs from canonical tiling-order flops; S626 checks the Moser-carrier lane where a true-optimum flop could first matter. In the width-1200 Moser beam, every exact witness through `n=14`, the exact `n=21` `57`-edge row, and the `n=22` `60`-edge lane has a unit spine; `n=14` has `55160` unit HPs. But fixed lexicographic point-order flip tournaments already fail to carry an all-unit directed HP at `n=7` (`5/6` unit arcs) even though the unit graph has `60` unit HPs. Thus the first observed flop is gauge/tie-order noise. | HYP-2203, HYP-2202, HYP-2201, HYP-2197, HYP-2196, HYP-2194, HYP-2188, HYP-2176, T695, S626, `04-computation/unit_distance_unit_spine_tournament_s626.py`, `05-knowledge/results/unit_distance_unit_spine_tournament_s626.out`, `07-reflections/unit-distance-moser-unit-spine-gauge-separation-s626.md` |
| **Unit-distance H-gap carrier reframing** | Treats `H=7`/`H=21` as forbidden scalar-collapse tests for unit-distance edge carriers rather than literal unit-distance tournament counts | S627 shows exact `u(5)=7` is legal because it decomposes as `4` unit-spine edges plus `3` tile/bulk edges, while the S599z literal unit-tournament count is `15`, not `7`. The triangular/Harborth `n=11` lattice echo `21=10+11` is likewise an edge-carrier scalar; exact planar `u(11)=23` already jumps past it. Carrier-packet Tournament Analysis uses proof-obligation vertices such as spine, tile/bulk, frontier, OCF H-gap guardrail, round-LRC channel, Moser/non-lattice carrier, literal H, raw edge quotient, and equidecomposability ledger. Geometry/H-gap/scaling gauges flip `23/36`, `17/36`, `22/36` pair edges; the majority tournament has a large SCC, `5` directed 3-cycles, and `25` Hamiltonian paths. | HYP-2204, HYP-2203, HYP-2202, HYP-2201, HYP-2200, HYP-2178, HYP-2181, T696, S627, `04-computation/unit_distance_hgap_carrier_s627.py`, `05-knowledge/results/unit_distance_hgap_carrier_s627.out`, `07-reflections/unit-distance-hgap-carrier-reframing-s627.md` |
| **Unit-distance n=21 centered-hex carrier split** | Moser `n=21` as a spine/bulk carrier, with SC coset checks used as a guardrail rather than a new H-value route | S630 treats the incoming S629/THM-409/HYP-2205 SC perspective theorem as primary, then extends the finite SC audit through `n=7`: `|Anti(T)|=|Aut(T)|`, scores reverse by `d -> n-1-d`, vertex anti-image sets are dual automorphism orbits, and `H=7,21` remain absent. The unit-distance claim is HYP-2206: exact/Moser `n=21` reads `57 = 20` unit-spine edges `+ 37` tile/bulk edges, where `37=C_hex(3)`. HYP-2207 records `Phi3(2)=7` and `3*Phi3(2)=21` as root-of-unity carrier guardrails, not raw equalities between tournament `H`, unit-distance edges, and LRC shells. | HYP-2207, HYP-2206, HYP-2205, HYP-2204, HYP-2203, HYP-2200, THM-409, T700, T699, S630, S631, `04-computation/sc_complement_perspective_s630.py`, `05-knowledge/results/sc_complement_perspective_s630.out`, `07-reflections/sc-complement-cosets-and-unit-distance-n21-s630.md` |
| **Applied perspective-carrier compression** | Use converse/perspective retention as a state-space reducer for LRC, rooted tournament counts, and unit-distance cores | S634 applies THM-409/HYP-2205 and HYP-2206 to concrete work queues.  LRC source targets and sink targets are converse-paired copies of `U(n-1)`, so the merged source/sink target list has size `U(n-1)`.  Rooted perspectives modulo all-edge reversal obey Burnside: `Q_n=(P_n+F_n)/2`, with `Q_1..Q_6=1,1,3,6,28,148`.  The exact-core audit verifies all five stored `57`-edge `n=21` graph6 cores are traceable, making `57=20+C_hex(3)` graph-real and pointing the `n=22` route toward endpoint-compatible ears. | HYP-2210, HYP-2206, HYP-2205, THM-381, THM-409, T703, S634, `04-computation/applied_perspective_carrier_s634.py`, `05-knowledge/results/applied_perspective_carrier_s634.out`, `07-reflections/applied-perspective-carriers-s634.md` |
| **Rational shadow carrier obstruction** | Treat rational/irrational as field-of-definition descent, parallel to parity and two-shadow carrier retention | S635 corrects the pi/e claim: it is known that `pi+e` and `pi*e` cannot both be algebraic, not that exactly one is rational.  The algebraic mechanism is Vieta reconstruction: if `S=x+y` and `P=xy` both descend to a tame field, then `{x,y}` descends via `T^2-S*T+P`; hence for `{e,pi}` at least one shadow must retain transcendence.  Finite-field unordered pairs show sum-only and product-only fibers while joint `(sum,product)` has max fiber `1`.  LRC analogue: classify by scale-free ratios, reset period/rational rank, relation lattice, and short circuits before raw time search. | HYP-2211, HYP-2210, HYP-2209, HYP-2154, HYP-2155, T704, S635, `04-computation/rational_shadow_carrier_s635.py`, `05-knowledge/results/rational_shadow_carrier_s635.out`, `07-reflections/rational-shadow-carriers-s635.md` |
| **Quadratic pi/e discriminant carrier** | Adds the branch sheet `D=e-pi` to the trace/norm carrier `S=e+pi`, `P=e*pi` | S636 proves a sharper unconditional statement: any two of `{S,P,D}` reconstruct `e` and `pi` algebraically, so no two can both be algebraic; therefore at least two of `e+pi`, `e*pi`, and `e-pi` are transcendental.  The discriminant identity is `S^2-4P=D^2`, or `S^2-D^2-4P=0`.  PSLQ finds no small scalar/pair relations in the tested range and detects the structural sheet relation.  Schanuel would make all nonconstant polynomial shadows of `(e,pi)` transcendental. | HYP-2212, HYP-2211, HYP-2210, HYP-2154, HYP-2155, T705, S636, `04-computation/quadratic_pi_e_carrier_s636.py`, `05-knowledge/results/quadratic_pi_e_carrier_s636.out`, `07-reflections/quadratic-pi-e-carrier-s636.md` |
| **Dehn/scissors side-channel addendum** | Extends HYP-2211/HYP-2212 through cuboid/simplex scissors congruence and unit-distance impairments | S637 records the Dehn addendum: cuboid/simplex comparisons mirror the rational-shadow lesson because volume is the scalar and Dehn data is the retained channel.  A regular tetrahedron has nonzero Dehn side data (`cos(theta)=1/3`, so `theta/pi` is irrational by Niven), while cuboids have zero Dehn invariant.  The repo dictionary maps this to tournament SCC/beta/c3 packets, LRC `C=2n-1` unit/nonunit/lift/carry/owner labels, and unit-distance spine/bulk/traceability/direction-support packets. | HYP-2213, HYP-2212, HYP-2211, HYP-2210, HYP-2209, HYP-2187, HYP-2186, T706, S637, `04-computation/dehn_scissors_shadow_s637.py`, `05-knowledge/results/dehn_scissors_shadow_s637.out`, `07-reflections/dehn-scissors-side-channel-addendum-s637.md` |
| **Pi/e lonely-shadow fallout** | Algebraic descent of one pi/e shadow forces transverse shadows to stay transcendental | S638 adds two fallout routes to S635/S636.  If `log(pi)` were algebraic, Lindemann-Weierstrass/Hermite-Lindemann would force `e+pi`, `e*pi`, and `e-pi` all transcendental, so algebraicity of any one proves `log(pi)` transcendental.  Newton power sums show that algebraic `S=e+pi` forces every `e^k+pi^k` for `k>=2` transcendental, while algebraic `P=e*pi` forces every `e^k+pi^k` for `k>=1` transcendental.  The repo analogue is the H=21 side-condition closure: scalar quotients become usable only with retained transverse packets. | HYP-2214, HYP-2213, HYP-2212, HYP-2211, HYP-2200, T707, S638, `04-computation/pi_e_lonely_shadow_s638.py`, `05-knowledge/results/pi_e_lonely_shadow_s638.out`, `07-reflections/pi-e-lonely-shadow-fallout-s638.md` |
| **Basel pi period carrier** | Even-zeta pi identities split into disjoint elementary packets, power-sum moments, and Bernoulli local-prime data | S653 treats Basel identities as one period-carrier family.  The sine product gives `zeta({2}^m)=pi^(2m)/(2m+1)!`; Newton identities/log-derivative recover `zeta(2k)=c_k*pi^(2k)`; Euler's Bernoulli formula gives the same coefficients through `k=12`; and von Staudt-Clausen records denominator primes, including the repo chain `6 -> 42 -> 1806 -> 1806`.  The bridge to tournaments is OCF: disjoint reciprocal-square packets parallel disjoint odd-cycle collections before scalar evaluation. | HYP-2229, HYP-2228, HYP-2227, HYP-2214, HYP-2211, HYP-1627, T722, T252, S653, `04-computation/basel_pi_family_carrier_s653.py`, `05-knowledge/results/basel_pi_family_carrier_s653.out`, `07-reflections/basel-pi-family-carrier-s653.md` |
| **Obstruction carrier atlas** | Shared side-channel-jackknife method across H=21, LRC, unit distance, A000568, Schanuel, twin primes, and Goldbach | S639 names a common carrier schema: raw scalar quotient, hidden witness, retained side channels.  It computes LRC `C=27` shell strata, twin-prime and Goldbach primorial survivor ledgers, carrier tag similarities, and a transfer-method tournament.  The tournament is transitive with one Hamiltonian path and ranks `side_channel_jackknife` above spine/bulk decomposition, finite-window extinction, local-prime ledgers, transverse-shadow fallout, Burnside transporters, quotient-collapse tests, Schanuel completions, and raw scalar numerology. | HYP-2215, HYP-2214, HYP-2213, HYP-2212, HYP-2211, HYP-2209, HYP-2206, HYP-2200, T708, S639, `04-computation/obstruction_carrier_atlas_s639.py`, `05-knowledge/results/obstruction_carrier_atlas_s639.out`, `07-reflections/obstruction-carrier-atlas-s639.md` |
| **Resonance product-ledger atlas** | Local obstruction-product ledger refining the side-channel jackknife method | S640 takes HYP-2215's side-channel jackknife as the parent move and turns it into a product-ledger target: for each domain, record local prime/modulus, forbidden residues or channels, surviving residues, local weight, and the lost side channel.  The atlas connects `H=21` strong/OCF packets, LRC `C=27` gcd/lift/carry/owner data, unit-distance `57=20+37` spine/bulk traceability, A000568 Burnside shadows, pi-e discriminant/transverse sheets, and twin-prime/Goldbach admissibility and singular-series data.  Its top-ranked next build is `local_obstruction_product_ledger`. | HYP-2216, HYP-2215, HYP-2214, HYP-2213, HYP-2212, HYP-2211, HYP-2210, HYP-2200, HYP-2164, T709, S640, `04-computation/resonance_product_ledger_s640.py`, `05-knowledge/results/resonance_product_ledger_s640.out`, `07-reflections/resonance-product-ledger-s640.md` |
| **LRC14 / UD21 difficulty bridge** | Shared `27`-quantum section/bulk obstruction connecting LRC `n=14` to unit-distance `n=21` | S641 sharpens the bridge: LRC `n=14` has shell modulus `C=2n-1=27=3^3`, with `<2,-1>` strata `1,3,9`, AP/V*/`2*AP` floor atoms, and open lift/CRT carry conservativity; unit-distance `n=21` is the Moser slab `P_2^-` with `E=2*27+3=57`, `20` unit-spine edges, and `37=C_hex(3)` bulk.  Monad S6's exhaustive `n=9` H-spectrum keeps low gaps `[7,21]`, so `21` is a durable tournament obstruction, but the shared LRC/UD difficulty is a retained section plus a small load-bearing side channel: LRC unit-shell section vs unit Hamiltonian spine, and LRC carry/owner lift vs Moser bulk/ear gluing. | HYP-2217, HYP-2216, HYP-2215, HYP-2206, HYP-2197, HYP-2188, HYP-2167, HYP-2166, HYP-2164, HYP-2200, THM-408, THM-407, THM-401, THM-115, T710, S641, `04-computation/lrc14_ud21_difficulty_bridge_s641.py`, `05-knowledge/results/lrc14_ud21_difficulty_bridge_s641.out`, `07-reflections/lrc14-unit-distance21-difficulty-bridge-s641.md` |
| **Goldbach/Lemoine same-pair projection** | Invertible pair-space carrier linking even `E=p+q` and odd `O=p+2q` targets | S642 turns the user's unordered-even/ordered-odd idea into a two-shadow projection: for odd primes `p,q`, the pair `(E,O)` reconstructs `q=O-E`, `p=2E-O`; swapping the pair fixes `E` and reflects `O -> 3E-O`; the duplicate diagonal maps `p -> (2p,3p)`, with `p=7` giving the visible `(14,21)` shadow.  This complements the LRC14/UD21 bridge: `14/21` is a prime-7 additive diagonal, while the hard LRC/UD proof channel remains the retained `27` shell/spine-bulk ledger. | HYP-2218, HYP-2217, HYP-2216, HYP-2215, HYP-2211, HYP-2208, HYP-2051, HYP-2049, HYP-2044, T711, S642, `04-computation/goldbach_lemoine_pair_projection_s642.py`, `05-knowledge/results/goldbach_lemoine_pair_projection_s642.out`, `07-reflections/goldbach-lemoine-same-pair-projection-s642.md` |
| **Goldbach/Lemoine same-pair bridge graph** | Finite even/odd companion graph refining the S642 pair projection | S643 builds on HYP-2218 by treating each unordered odd-prime pair `{p,q}` as an even node `E=p+q` with one or two odd companion nodes `E+p` and `E+q`.  The bridge graph exposes parent-rich odd nodes, the duplicate branch locus `(E,O)=(2p,3p)`, and the prime-2 apex channel: Lemoine reps with doubled prime `2` have no even same-pair companion; up to odd target `121`, only `7=3+2*2` is apex-only. | HYP-2219, HYP-2218, HYP-2215, HYP-2049, HYP-2044, HYP-1984, T712, S643, `04-computation/goldbach_lemoine_pair_bridge_s643.py`, `05-knowledge/results/goldbach_lemoine_pair_bridge_s643.out`, `07-reflections/goldbach-lemoine-same-pair-bridge-s643.md` |
| **Vieta/perfect aliquot carrier** | Triangular pair-count carrier `A=C(n,2)` with discriminant root `sqrt(8A+1)=2n-1` and proper-divisor shadow `s(A)` | S644 makes `C=2n-1` a Vieta square-root observer on pair counts.  Perfect numbers are triangular fixed controls (`6=C(4,2)`, `28=C(8,2)`, `496=C(32,2)`), while `C(14,2)=91=7*13` has exact aliquot shadow `s(91)=21`.  The semiprime family `n=2p`, `2p-1` prime gives `s(C(2p,2))=3p`, matching HYP-2219/S643's duplicate branch at `p=7`: `(2p,3p)=(14,21)`. | HYP-2220, HYP-2219, HYP-2218, HYP-2217, THM-401, THM-361, T713, S644, `04-computation/vieta_perfect_aliquot_carriers_s644.py`, `05-knowledge/results/vieta_perfect_aliquot_carriers_s644.out`, `07-reflections/vieta-perfect-aliquot-carriers-s644.md` |
| **Perfect-number aliquot fixed-point ecology** | Perfect numbers as loops of the divisor-sum dynamical map `s(n)=sigma(n)-n` | S645 builds on S644/HYP-2220's Vieta/perfect aliquot carrier by scanning the aliquot map itself: `D(n) -> sigma(n) -> s(n)`.  Perfect numbers are fixed points `s(n)=n`, equivalently `sigma(n)=2n` or abundancy `A(n)=2`; the proof side channel is the prime-power product for `sigma(n)/n`.  The finite scout through `100000` finds fixed points `[6,28,496,8128]`, `13` amicable pairs, one length-5 sociable cycle, no quasi-perfect defect `+1` rows, and the closest odd near-fixed row `32445` with `sigma(n)-2n=6`. | HYP-2221, HYP-2220, HYP-2216, HYP-2215, HYP-2211, HYP-2208, T714, S645, `04-computation/aliquot_fixed_point_carrier_s645.py`, `05-knowledge/results/aliquot_fixed_point_carrier_s645.out`, `07-reflections/perfect-aliquot-fixed-point-carrier-s645.md` |
| **LRC Pillai fixed-clock carrier** | The `C=2n-1` antipodal shell clock has a gcd-shell divisor mass `A(C)=sum gcd(a,C)`; `n=14` is a fixed clock | S646 transfers the perfect-number fixed-point lens to LRC.  For odd `C`, `A(C)=(P(C)-C)/2`, where `P(C)=sum_{r=1}^C gcd(r,C)` is Pillai's gcd-sum.  The fixed equation `A(C)=C` is equivalent to `P(C)=3C`, and the only odd fixed clocks are `C=15,27`, i.e. LRC `n=8,14`.  At `n=14`, AP is shell-fixed while `Vstar` is gcd-carrier-fixed: it moves shell `12` to shell `3` but preserves gcd counts `{1:9,3:3,9:1}`.  The next proof target is a no-leak lemma: mass-changing rows are loose; mass-fixed rows are AP/Vstar or killed by pair-pinch/carry-owner labels. | HYP-2222, HYP-2221, HYP-2220, HYP-2217, HYP-2177, HYP-2167, THM-401, THM-369, T715, S646, `04-computation/lrc_pillai_fixed_clock_s646.py`, `05-knowledge/results/lrc_pillai_fixed_clock_s646.out`, `07-reflections/lrc-pillai-fixed-clock-carrier-s646.md` |
| **Triangular fixed tournament controls** | Perfect triangular arc counts split into Hamiltonian spine plus off-path deformation fiber | S647 translates HYP-2220 into tournament language: every `m`-vertex tournament has `C(m,2)` arcs, and fixing one Hamiltonian path leaves `C(m-1,2)` off-path free arcs.  The second perfect control is `C(8,2)=28=7+21`, with `7` a section/spine length and `21` a fiber dimension.  Double-counting the full labelled `n=8` H-spectrum gives the exact fixed-base-path fiber of size `2^21`, where `H=7` and `H=21` have count `0`; this is a role-mismatch guardrail for LRC/unit-distance scalar transfers. | HYP-2223, HYP-2220, HYP-2200, HYP-2202, HYP-2203, HYP-2217, THM-115, T716, S647, `04-computation/triangular_fixed_tournament_controls_s647.py`, `05-knowledge/results/triangular_fixed_tournament_controls_s647.out`, `07-reflections/triangular-fixed-tournament-controls-s647.md` |
| **Unit-distance Moser fixed quantum** | THM-408's Moser add-one-slab operation has a stable `27` edge quantum split into spine and bulk carriers | S648 turns the fixed-carrier lens onto unit distance.  In both `P_m^+` and `P_m^-`, after the cap transient, adding one slab has direction increment `(0,1,8,4,0,4,4,2,4)`, sum `27`, split as `8` unit-spine edges plus `19` hidden bulk edges.  At the frontier `m=2` rows, `P_2^-` has `n=21,E=57,spine=20,bulk=37,pure_bulk=20`, while `P_2^+` has `n=22,E=60,spine=21,bulk=39,pure_bulk=21`.  The `n=22` row has one degree-3 cap vertex, so a hypothetical `61`-edge row inside this carrier must repair the cap endpoint or leave the fixed Moser quantum. | HYP-2224, HYP-2222, HYP-2217, HYP-2204, HYP-2203, HYP-2188, HYP-2176, THM-408, T717, S648, `04-computation/unit_distance_moser_fixed_quantum_s648.py`, `05-knowledge/results/unit_distance_moser_fixed_quantum_s648.out`, `07-reflections/unit-distance-moser-fixed-quantum-s648.md` |
| **Euler-Heegner prime-window boundary** | Prime-generating quadratics split into source, interior prime window, and forced square sink | S649 records the exact boundary for `f_p(x)=x^2+x+p`: the lucky primes `{2,3,5,11,17,41}` are precisely recovered in a search through `p<=500` by prime values `x=0..p-2` and forced failure `f_p(p-1)=p^2`.  The Heegner side channel is `d=4p-1`, sending `{7,11,19,43,67,163}` to the lucky primes.  The tournament dictionary makes the indexing useful: `p-1` is fixed Hamiltonian spine length, `p-2` is the interior/Moon `c3` floor, and `C(p-1,2)` is the off-path fiber. | HYP-2225, HYP-2224, HYP-2223, HYP-2222, HYP-2217, HYP-2215, HYP-2200, THM-115, T718, S649, `04-computation/prime_heegner_tournament_boundary_s649.py`, `05-knowledge/results/prime_heegner_tournament_boundary_s649.out`, `07-reflections/prime-heegner-tournament-boundary-s649.md` |
| **Heegner prime-polynomial horizon** | Euler/Rabinowitsch prime quadratics have a boundary/interior `p-2` witness horizon matching THM-410 interval reversals | S650 refines HYP-2225's Euler-Heegner boundary packet by adding the exact THM-410 pair-witness model.  The Euler primes `{2,3,5,11,17,41}` map by `d=4p-1` to the Heegner tail `{7,11,19,43,67,163}`.  For `Q_p(x)=x^2+x+p`, the exact horizon is `Q_p(0)=p`, `Q_p(1..p-2)` prime, and `Q_p(p-1)=p^2`; among primes `p<=500`, only the Euler set has first positive composite exactly at the boundary `x=p-1`.  Reversing one long edge in a `p`-vertex transitive tournament creates exactly `p-2` cyclic-triangle witnesses with split `sigma=0, lambda=p-2, delta=0`. | HYP-2226, HYP-2225, HYP-2224, HYP-2223, HYP-2222, HYP-1220, HYP-1202, THM-410, T719, S650, `04-computation/heegner_prime_horizon_tournament_s650.py`, `05-knowledge/results/heegner_prime_horizon_tournament_s650.out`, `07-reflections/heegner-prime-horizon-tournament-s650.md` |
| **Number-theory tournament carrier atlas** | Hard number-theory problems become finite local-witness tournaments when the right side channels are retained | S651 generalizes S650: tournaments fit number theory as local witness/obstruction carrier quotients, not raw integer graphs.  Paley tournaments are literal finite-field character orientations (`p=7` has `c3=14,H=189`; `p=11` has `c3=55,H=95095`).  Sieve obstruction tournaments show side-channel dependence: twin gap `2` vs Goldbach `N=210` has `3` local-prime edge flips, vs `N=2110` has `1`, Goldbach `210` vs `2110` has `4`, and twin vs Euler `p=41` has `62`.  The carrier supplement adds dynamic prime races: `q=11`, `X=100000` has `15` flips versus final count, `1` directed 3-cycle, SCC sizes `[3,1,1,1,1,1,1,1]`, and `H=3`. | HYP-2227, HYP-2226, HYP-2225, HYP-2224, HYP-2223, HYP-2217, HYP-2216, HYP-2215, THM-410, THM-115, T720, S651, `04-computation/number_theory_tournament_atlas_s651.py`, `05-knowledge/results/number_theory_tournament_atlas_s651.out`, `04-computation/number_theory_tournament_carriers_s651.py`, `05-knowledge/results/number_theory_tournament_carriers_s651.out`, `07-reflections/number-theory-tournament-carriers-s651.md` |
| **THM-410 square-blowup enumeration carriers** | Near-transitive and imprimitive tournament enumeration can be split into exact interval weights, upset bitsets, and block substitution templates | S652 combines THM-410 with modular decomposition.  Matching reversals from a transitive order have exact `c3=sum(b-a-1)`, giving a theorem-grade low-`c3` branch key (`188273/2390480` matchings at `n=14` have weight `<=10`); the supplemental scout adds the arbitrary-upset rule `r_xy=r_yz!=r_xz`.  Square-node substitution gives a Burnside companion to A000568 with values `1,1,12,704,2127984` through `n=5`, exact square formulas `score((i,a))=n*score_T(i)+score_T(a)` and `c3(Sq(T))=(n^3+n)c3(T)`, and H on lexicographic products is computed by block path-cover polynomials plus macro-word DP rather than naive multiplication (`cycle3[chain3,chain3,chain3]` has `H=2721`, not `3`; `C3[C3]` has `H=3159`, not `81`). | HYP-2228, HYP-2227, HYP-2226, HYP-2209, HYP-2200, THM-410, THM-115, T721, S652, `04-computation/thm410_square_blowup_enumeration_s652.py`, `05-knowledge/results/thm410_square_blowup_enumeration_s652.out`, `04-computation/thm410_square_blowup_speedups_s652.py`, `05-knowledge/results/thm410_square_blowup_speedups_s652.out`, `07-reflections/thm410-square-blowup-enumeration-speedups-s652.md` |
| **Unit distance Moser spine-ladder theorem** | Explicit traceable-section theorem for layered rank-4 Moser slabs | THM-408 proves two infinite unit-spine families behind the S626 witnesses. The plus family `P_m^+` has `8m+6` vertices and includes `n=14`/`n=22`; the shifted minus family `P_m^-` has `8m+5` vertices and includes the exact `n=21` row. Each full fiber is an 8-vertex Gray-code slab, caps truncate the final fiber, and the bridge vector `(-1,1,0,0)` glues adjacent layers. Abstractly this is a traceable section over a carrier quotient: point-order tournaments can forget the section, but slab/ear proof obligations retain it. | THM-408, HYP-2204, HYP-2203, HYP-2202, HYP-2188, T697, S628, `04-computation/unit_distance_spine_ladder_s628.py`, `05-knowledge/results/unit_distance_spine_ladder_s628.out`, `07-reflections/unit-spines-as-traceable-sections-s628.md` |
| **SC perspective-flip cyclotomic carrier** | Self-converse edge reversal acts canonically on rooted perspectives, not necessarily on vertices | THM-409 proves `Anti(T)` is a coset over `Aut(T)` and induces a canonical involution on `Aut(T)`-vertex orbits. S629 verifies this through `n=6`, including anti-automorphisms with vertex cycle type `(6,)`, so the rooted-perspective lift is the correct object. The same packet keeps `Phi_3(2)=7`/`Phi_3(4)=21` as cyclotomic carrier echoes: unit-distance edge counts can visibly hit `7` or `21` only with spine/bulk side channels retained, while `P_2^-` is a 21-vertex, 57-edge Moser slab rather than `H=21`. | THM-409, HYP-2205, HYP-2204, THM-408, HYP-2121, HYP-2181, T698, S629, `04-computation/sc_perspective_flip_cyclotomic_s629.py`, `05-knowledge/results/sc_perspective_flip_cyclotomic_s629.out`, `07-reflections/sc-perspective-flip-cyclotomic-carriers-s629.md` |
| **Anti-coset transporter atlas** | Common `Anti(X)=Iso(X,JX)` schema for complements, converses, reflections, and conjugations | S632 generalizes THM-409's SC anti-coset into a transporter atlas.  SC tournaments through `n=7` have no coset failures and still avoid `H=7/21`; LRC folded shells under `<2,-1>` for even `n=6..20` have doubled-shell transporters the same size as stabilizers, with `n=14`, `C=27` yielding gcd strata `1,3,9`; Eisenstein prefixes keep the full `D6` anti-coset at centered shells `n=1,7,19` but only a small or empty shard on most partial shells.  The operational rule is to replace raw complement merges by `(X,J,Aut,Anti)` records and name quotient-loss side channels. | HYP-2208, HYP-2207, HYP-2206, HYP-2205, THM-409, THM-407, T701, S632, `04-computation/anti_coset_everywhere_s632.py`, `05-knowledge/results/anti_coset_everywhere_s632.out`, `07-reflections/anti-cosets-everywhere-s632.md` |
| **Sequence-shadow recursion lab** | Extending hard counts by fixed, merged, bisection, q-shadow, skinny, and transporter companion sequences | S633 treats `SC=1,1,2,2,8,12,88,...` as one layer in a family around A000568 rather than a lone next-term target.  The exact scout records `T=A000568`, `SC`, `(T+SC)/2`, `(T-SC)/2`, q-deformed shadows `A(m,q)`, SC bisection ratios, round/LRC self-converse shadows, and folded shell transporter-orbit counts.  The useful recursion is that `SC(2m)=A(m,4)` while `SC(2m+1)` is the same base-4 Burnside sum with a `2^(#parts)` fixed-vertex tax; the method tournament ranks transporter quotient and bisection factor above raw next-term chasing. | HYP-2209, HYP-2208, HYP-2074, HYP-2064, HYP-2086, T702, S633, `04-computation/sequence_shadow_lab_s633.py`, `05-knowledge/results/sequence_shadow_lab_s633.out`, `07-reflections/sequence-shadow-recursions-s633.md` |
| **Tournament gap transfer principle** | `H=7` and `H=21` as unavailable evaluations that propagate as coimage certificates | S615 reframes the permanent tournament gaps as forbidden values of `I(Omega,2)` whose transfer value is not scalar numerology but a side-channel test: when another problem's visible quotient demands an unavailable state, the quotient forgot forced completion data. The mechanism tournament ranks coimage side-channel retention above H-gap completion forcing, determinant walls, exponent normalization, carry-residue transfer, and raw scalar matching. The shared `1.014` observation is recorded as a normalization target: compare entropy dividends of retained carriers, not raw A000568/H growth. | HYP-2178, HYP-2176, HYP-2175, HYP-2157, THM-002, THM-029, THM-079, S615, `04-computation/tournament_gap_transfer_exponent_s615.py`, `05-knowledge/results/tournament_gap_transfer_exponent_s615.out` |
| **Tournament obstruction/amplification carrier** | Cross-problem reading of forbidden tournament values, LRC shell arithmetic, Collatz two-block gaps, and the unit-distance `1.014` arithmetic amplification | S615 separates the supported bridge from the open exponent claim and now sits downstream of the transfer principle and strong-component mechanism.  Supported: `7=Phi_3(2)` and `21=3*Phi_3(2)` are permanent tournament `H` gaps; HYP-2180 explains them via multiplicativity over strong components and strong-value gaps; HYP-2179 verifies robust finite-spectrum propagation to round LRC worry-set tournaments; LRC reuses the same small arithmetic through `C=21=3*7`, `C=27=3^3`, incoming HYP-2177's `gcd(n-2,2n-1)=gcd(3,2n-1)` lattice law, and S612b's reflection-vs-doubling single-swap classification; Collatz has the parallel two-block residual `2^E-3^k`.  Open: no public or in-repo tournament theorem currently gives the unit-distance `n^1.014` exponent.  The next proof object should be a constrained carrier tournament on split primes, norm-one generators, shell/carry states, or proof obligations. | HYP-2181, HYP-2180, HYP-2179, HYP-2178, HYP-2177, HYP-2176, HYP-2175, HYP-2167, HYP-2159, T682, S615, `04-computation/tournament_obstruction_amplification_s615.py`, `05-knowledge/results/tournament_obstruction_amplification_s615.out`, `07-reflections/tournament-obstruction-amplification-s615.md` |
| **LRC infinite-game n+2 recursion** | Finite successor automaton on the clean even ladder plus `omega^2` row/column proof rank | S616 turns the `n+2` ladder into a residue automaton inside S599t's partition-function/free-energy frame: `n->n+2` sends `C=2n-1` to `C+4`, so `C mod 3` cycles through three states. With HYP-2177, the doubled AP sporadic `V*` is alive exactly on the `3|C` phase: `n=8,14,20,26,...` (`n==2 mod 6`). HYP-1881/HYP-1891 become a well-founded denominator proof rank `rho(N)=omega*((odd(N)-1)/2)+v_2(N)`, where rows discharge inside odd-root columns and columns move by a finite automaton. This makes n=14 the `C=27=3^3` alive phase plus a HYP-2166/HYP-2167 lift/carry side-channel problem, not a one-off sporadic. | HYP-2183, HYP-2182, HYP-2181, HYP-2177, HYP-2091, HYP-1881, HYP-1891, HYP-2166, HYP-2167, T683, S616, `04-computation/infinite_game_nplus2_recursion_s616.py`, `05-knowledge/results/infinite_game_nplus2_recursion_s616.out`, `07-reflections/infinite-game-nplus2-recursion-s616.md` |
| **Infinite-Go boundary partition recursion** | Finite ladder partition function with retained boundary state/pole order | S616b isolates the boundary-state face of the infinite-game analogy: for `r=liberties-2`, finite cutoffs give `Z_{r,K}(q)=(1+...+q^K)^r`, while the open-game limit reads serial fuel depth as `omega*r`. The fixed-path tournament `n->n+2` boundary exposes `binom(n+1,2)-binom(n-1,2)=2n-1` variables, matching the LRC floor shell `C=2n-1`; scalarizing to `Z(1)` forgets exactly the side-channel the proof route needs. | HYP-2185, HYP-2182, HYP-2183, HYP-2161, THM-401, T684, S616b, `04-computation/infinite_go_partition_recursion_s616.py`, `05-knowledge/results/infinite_go_partition_recursion_s616.out`, `07-reflections/infinite-go-nplus2-partition-functions-s616.md` |
| **Equinumerosity/equidecomposability fiber bridge** | Distinguish cardinal equality from predicate-preserving scissors quotients | S617 audits the older Royle/even-graph and Hilbert-third tournament threads together, complementing HYP-2186/S599v's strong-component scissors-volume synthesis. Degree-even graph classes are not equinumerous with tournaments after `n=3` (`n=6: 16 vs 56`), and the naive even-order Burnside complement is `graphs - tournaments`, not Royle-even. On the tournament side, scalar `H` is only volume: `H=9` splits by `beta1` at `n=5`, five `H` values split by `beta1` at `n=6`, and disjoint directed-3-cycle packet data refines `24` `(H,beta1)` classes to `29`. The proof target is a fiber functor preserving `H`, `beta1`, and packet side-channels, not a bare A000568-sized bijection. | HYP-2187, HYP-2186, HYP-1592, T685, S617, S599v, S90, S92, S259, S260, `04-computation/equinum_equidecomp_bridge_s617.py`, `05-knowledge/results/equinum_equidecomp_bridge_s617.out`, `07-reflections/equinumerosity-equidecomposability-fiber-bridge-s617.md` |
| **LRC even clean-lane fixed fibre** | Self-converse fixed-node residual inside the round/converse seam for even LRC n | S577 merges HYP-2092 with Opus HYP-2093, Oracle HYP-2094, Opus HYP-2095, and the HYP-2096 unit-spine audit: n=14 compresses from all runner tournaments `48,542,114,686,912` to round `316`, converse-merged `190`, and then the actual fixed-boundary target `64`. The next proof artifact should be a 64-row fixed-round certificate table with unblocked-small-pair, n-clock, pair-sum, D/U/N private-pivot, THM-397 endpoint-owner, unit-spine, and gcd-3/gcd-9 nonunit descent labels; the other 126 merged nodes are controls. | HYP-2097, HYP-2096, HYP-2095, HYP-2094, HYP-2093, HYP-2092, HYP-2091, HYP-2090, THM-397, S577, `04-computation/lrc_even_clean_lane_fibre_collapse_s577.py` |
| **LRC fixed-boundary owner scaffold** | Two-layer n=14 certificate architecture: fixed round boundary plus speed-owner unit-spine labels | HYP-2098 shows the 64 self-converse classes overcount tightness and the real tight speed family appears to be `{AP,V*}` on a mod-7 tie-wall. Opus S570 covers the transversal flip-lattice (`8191/8191` lonely and unblocked, AP unique tight). S578 covers the composite unit-spine fibre (`820` labelled anti-witness vectors -> `64` groups; `63/64` strong/unique-axis). S609 computes exact HP counts on all 64 fixed classes: all odd and distinct, but regularity and HP-scarcity gauges flip `2005/2016` edges, so parity is scaffold and owner labels carry certificates. S610 reframes this as the reattachment leg of the quotient tower: proof now targets lift/CRT conservativity, not another fixed-class scalar. | HYP-2166, HYP-2165, HYP-2164, HYP-2163, HYP-2162, HYP-2099, HYP-2098, HYP-2097, HYP-2096, HYP-2095, HYP-2094, THM-397, THM-407, S610, S609, S578, Opus S570, `04-computation/lrc_n14_res27_quotient_tower_s610.py`, `04-computation/lrc_n14_res27_fixed_bridge_s609.py`, `04-computation/lrc_n14_fixed_round_certificate_s578.py` |
| **LRC apex-lift certificate sheaf** | Local cheap-pair certificates glued over the mod-7 tight tie-wall and unit-shell lift site | S579 reframes the HYP-2099 bridge lemma as a sheaf over six mod-7 antipodal lanes plus the self-antipodal apex, then tests the unit-shell lift site. Local sections are HYP-2095 cheap-pair certificates, side-denominator charts, positive-measure controls, and ledger-failure patches carrying shield/anchor, D/U/N, endpoint-owner, and unit-spine labels. In the bounded audit, `13169` one-lift full covers split into `7943` denominator-14 apex-chart sections and `5226` side-chart sections with zero restriction residuals; named AP/`V*`/open-gap two-lift sites add union-only `{ledger,ledger}->cheap` rows. HYP-2102 supplies the C-prime multiple-of-n positive-measure branch. | HYP-2102, HYP-2101, HYP-2100, HYP-2099, HYP-2098, HYP-2097, HYP-2096, HYP-2095, HYP-2088, HYP-2083, HYP-2063, HYP-2045, THM-397, S579, Opus S571, `04-computation/lrc_n14_apex_lift_sheaf_s579.py` |
| **LRC n=14 certificate calculus** | Formal gate calculus whose vertices are proof obligations rather than runners or classes | S580 wraps THM-398/HYP-2104/HYP-2103, HYP-2101, HYP-2100, HYP-2099, and HYP-2095 into gates `SpeedRow -> CertificateSection or NamedResidual`, with sections `witness_1_over_n`, `cheap_pair`, `positive_measure`, `ledger_failure`, and `residual`. The proof-obligation tournament is a transitive ledger; its top extension is endpoint-cover-circuit positivity, because the all-short Cprime/Vitali AP-alignment residual and the failed apex-sheaf gluing residual both ask for owner circuits that peel to positive measure. | HYP-2111, HYP-2108, HYP-2105, HYP-2104, HYP-2103, THM-398, HYP-2102, HYP-2101, HYP-2100, HYP-2099, HYP-2095, S580, `04-computation/lrc_n14_certificate_calculus_s580.py` |
| **LRC endpoint-cover circuit positivity** | Exact component-cover criterion plus a scalar positivity functional for the Cprime multiple branch | S575/HYP-2108 proves the endpoint-cover criterion: a component `C_i=(a_i,b_i)` is covered by one `v=nw` danger arc iff both endpoints have the same near-integer `v`-phase, equivalently `||v m_i|| <= 1/n - (v/2)l_i`. The functional `P(S)=max_i(||v m_i||+(v/2)l_i-1/n)` satisfies `P>0` iff loose, and the summed circuit corollary says tightness would force average midpoint `v`-phase `<1/n`, far below the sampled generic `~0.245`. | HYP-2108, HYP-2105, HYP-2104, HYP-2103, THM-398, S575, `04-computation/lrc_endpoint_cover_circuit_positivity_s575.py` |
| **LRC small-owner endpoint descent** | Endpoint-cover positivity sieve for the Cprime multiple branch | S574/HYP-2105 turns all-short cover geometry into endpoint-owner congruences and proves Lemma C: both small owners make a component uncoverable. S581/HYP-2110 pushes this further after HYP-2108: one small owner pins a `v=nw` arc centre exactly; off-lattice pins peel by Lemma E and pinned components longer than a half-radius peel by Lemma F. In sampled multiple-of-n rows, Bprime/C/E/F close all rows for n=10,12,14 and leave only large-owner or mixed tiny-window residuals at n=6,8. | HYP-2110, HYP-2108, HYP-2107, HYP-2105, HYP-2104, HYP-2103, THM-398, S581, S574, `04-computation/lrc_small_owner_descent_s581.py` |
| **LRC paper-sieve certificate compression** | Methodology comparison for total n=11,12,13: exact prime fibers versus certificate clocks | S580 quantifies the gap between the 2026 paper verifier and the repo quotient. Paper raw `sum p^k` scales are `10^27.63`, `10^31.27`, `10^35.31`; full terminal ansatz scales are `10^38.05`, `10^55.01`, `10^48.68`; the certificate model using `n`-clock, `C=2n-1` shells, pair pinches, and D/U/N owners is `10^4.86`, `10^5.26`, `10^5.61`. The arithmetic rotates: paper likes odd-prime `n`, while the repo shell clock makes total `n=12` clean because `C=23` is prime; THM-398/HYP-2103 and S573/HYP-2104 supply the interval-Vitali backbone, reducing the multiple branch to all-short AP arc-alignment components of length at most `2/(n^2 w)` for `v=nw`. | HYP-2106, HYP-2111, HYP-2108, HYP-2107, HYP-2105, HYP-2104, HYP-2103, HYP-2102, HYP-2101, HYP-2100, HYP-2095, HYP-2083, S580, `04-computation/lrc_paper_vs_certificate_efficiency_s580.py` |
| **LRC large-owner residue automaton** | Post-S574 residual algorithm: compile large-owner cover congruences into bounded multiplier residue states | S581 extends the cover-to-congruence translator after B' and Lemma C. Each all-short large-owner component becomes allowed residue classes for `w`; component languages are CRT-intersected and pruned by the dominance cutoff. In sampled `n=11,12,13,14` rows, B' or Lemma C covers `94.4%,95.2%,95.8%,97.5%`; exact-classified large-owner rows have zero live bounded `w` states, with remaining rows being modulus overflows. Next algorithm: prime-power CRT/ZDD plus minimal incompatible component certificates. | HYP-2107, HYP-2108, HYP-2106, HYP-2105, HYP-2104, HYP-2103, THM-398, S581, `04-computation/lrc_large_owner_residue_sieve_s581.py` |
| **Three-state tournament automata** | Pair event words over `L/M/R` whose terminal state projects to a tournament edge | S582 reframes a tournament as the binary shadow of local pair automata. `L` and `R` are side/owner states; `M` is the wall, tie, midpoint corridor, live residue, or gluing seam. The wall rule forbids direct `L<->R` flips, so edge changes must pass through auditable middle states before a tie Hamiltonian path resolves terminal `M`. For LRC, the likely cells are components, endpoint owners, cover arcs, residue states, and proof obligations; HYP-2108/HYP-2110 attack all-middle endpoint circuits, while HYP-2107 attacks middle residue languages by CRT/dominance. | HYP-2109, HYP-2110, HYP-2108, HYP-2107, HYP-2105, S582, `04-computation/tournament_three_state_automaton_s582.py` |
| **LRC tournament speedup stack** | Layered use of tournament computation to avoid exact LRC search except on labelled residuals | S583 ranks reusable engines: Cprime/certificate gate calculus, exact observer-source marked reachability, SCC/good-cut protection cores, threshold-decorated quotient stacks, `L/M/R` middle automata, Burnside/A000568 restricted orbit counts, `A^2` isomorphism fingerprints, cheap TDA prefilters, and near-transitive perturbation updates. Safe use requires keeping marks, thresholds, endpoint-owner labels, or middle-state labels; raw phase classes are only heuristics because they mix lonely and non-lonely states. After S576/HYP-2112, the top Cprime gate should call the exact `Phi` gap functional before any raw residual search. | HYP-2113, HYP-2112, THM-381, THM-398, THM-354, THM-395, HYP-2109, HYP-2108, HYP-2107, S583, `04-computation/lrc_tournament_speedup_roadmap_s583.py` |
| **LRC fold-sieve/randomness bridge** | Split the multiple-branch proof route by 3-term additive circuits | S583 uses HYP-2112's exact `G=Phi(C)` value target while alternating Lemma A and Lemma B. Oracle S578o names the AP/shifted-AP phenomenon: 4-term energy is the translation-invariant shadow of folds. At n=14, AP and shifted AP have the same 4-term structure, but AP is tight with many folds while shifted AP is 3-free and safely above the floor. Fold rows are closer to the floor, and deletion-shadow witnesses are often blocked by reattaching the folded shield `c=a+b`; Opus S577 independently confirms fold-as-shield at the `1/(a+b)` pinch clock. The fold proof should route through another modulus, owner pin, or prime-fibre certificate, then feed HYP-2113's speedup stack. | HYP-2114, HYP-2113, HYP-2112, HYP-2111, HYP-2109, HYP-2108, HYP-2106, HYP-2105, HYP-2075, S583, S578o, S577, `04-computation/lrc_fold_sieve_random_bridge_s583.py`, `04-computation/lrc_three_four_term_energy_encoding_s578.py` |
| **LRC hidden four-term folds** | Pair-sum collisions without visible triples, completed by virtual sum nodes | S584 complements HYP-2114: a 4-term relation `a+b=c+d` is two hidden 3-term folds through virtual `s=a+b=c+d`. Fixed hidden-sum deformation families with no visible triples stay safely above `delta`; random 4-rich circuit-free rows through k=10 also stay safe, but adding the virtual sum often drops the augmented margin and the hidden runner is frequently at distance `0` at the original witness. Raw additive energy is therefore only a shadow; the certificate payload is original speeds, pair edges, hidden sum nodes, virtual fold pressure, and deformation labels. | HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2109, S584, S577, `04-computation/lrc_four_term_hidden_fold_s584.py` |
| **LRC index-stratified extremal fractal** | Translated APs as indexed faces of one recursive summand-completion family | S585 complements Opus HYP-2116's 2-adic self-affine tower and HYP-2117's binary IFS/doubling-map view by making the additive deformation index explicit. For `I(k,q)={q,...,q+k-1}`, the 4-term count is translation-invariant while visible folds are the clipped prefix `#{i<j:i+j<=k-1-q}` of the pair-sum tent. At `k=13`, fold counts are `36,30,25,20,16,12,9,6,4,2,1,0,0`, the formula has zero mismatches, 4-term count stays `125`, and `M/delta` rises `1.000 -> 4.789`; the maximal hidden node keeps multiplicity `6` while its shell walks outward. Summand completion has interval widths `13 -> 23 -> 43 -> 83 -> ...` by `w_{r+1}=2w_r-3`, so the proof address should include parity lane, `C`-gcd shell, visible fold clip, hidden sum shell, dyadic debt depth, and owner labels. | HYP-2118, HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2092, HYP-2091, HYP-2079, S585, `04-computation/lrc_index_stratified_extremal_fractal_s585.py` |
| **LRC extremality order** | Rank which "extreme" features carry proof power | S586 compares AP, unit-shift AP, far-shift AP, `V*`, doubled-apex stress, and hidden-fold controls using exact `M`, safe measure, saturation curves, 3-term folds, ordered energy, hidden virtual pressure, and clock binders. The resulting proof order is `critical_sleeve_saturation > Phi_endpoint_gap_terminal > visible_low_3fold_delta_clock > circuit_free_margin_floor > hidden_virtual_fold_pressure > dyadic_fractal_binder_profile > raw_4term_energy`. At n=14, AP and `V*` are critical floor rows, while far-shift AP has the same energy but `M=5/14`, and unit-shift AP has many folds but positive safe measure. THM-400 sharpens the rule: observer-coupled nonzero augmentation beats balanced support. | HYP-2119, THM-400, HYP-2118, HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2108, HYP-2107, S586, `04-computation/lrc_extremality_order_s586.py` |
| **LRC fold-divisibility sieve** | Lemma B upgraded from visible folds to pair-denominator shields | S587 proves the local algebraic gate: for `D=a+b` and `t=m/D`, every speed `v` with `D|v` is at an integer and shields the whole `D`-pinch family. Thus visible `a+b=c` is the primitive shield `v=D`, while `V*` uses `24` to shield the missing `D=12` fold. At n=14, AP and `V*` kill all appearing low `D<n` but leave `D=n` unshielded and stay at `M=1/14`; unit-shift AP shields `D=n` with `14` and becomes loose. HYP-2120/HYP-2121 add the rooted/source-perspective interpretation: denominator shields live in the observer-coupled source slice and carry incident threshold payload, while balanced energy is unrooted background. The proof route is: killed low denominators + unshielded `D=n` -> delta-clock; killed `D=n` -> Cprime/Phi/endpoint; no denominator scaffold -> Lemma A discrepancy. | HYP-2122, HYP-2121, HYP-2120, THM-400, HYP-2119, HYP-2118, HYP-2115, HYP-2114, HYP-2112, HYP-2110, HYP-2108, HYP-2107, S587, `04-computation/lrc_fold_divisibility_sieve_s587.py` |
| **LRC pincer calculus** | Abstract bidirectional certificate machine above pinch clocks | S588 searched the repo and web for "pincer" and found the durable repo concept is named `pinch`; the broader pincer is two monotone fronts closing around a marked observer/basepoint, with a meeting theorem and labelled escape ledger. LRC bottom jaw = S557 pair-pinch residues; top jaw = THM-396/397/HYP-2095/HYP-2122 shield-anchor-endpoint covers; grip labels = observer/source, side, threshold, residue, endpoint owner, and L/M/R middle state. Web analogues sharpen the design: Pincer-Search gives bottom-up/top-down pruning, bidirectional search gives meet-in-the-middle termination, double envelopment gives escape paths, pincer grasp gives force labels, and squeeze theorem gives the bound-meeting proof shape. | HYP-2123, HYP-2122, HYP-2121, HYP-2120, HYP-2109, HYP-2095, HYP-2082, HYP-2075, HYP-2059, THM-397, THM-396, S588, `07-reflections/lrc-pincer-calculus-pinches-envelopments-s588.md` |
| **LRC rigidity cascade** | Labelled fixed-point propagation across source, fold, endpoint, and pincer closures | S589 treats rigidity as a marked local fixed point plus labels that propagate naturally under isomorphism. Exact enumeration through tournament n=6 shows source roots are rigid: source-rooted n-classes equal U(n-1), deleting the source is collision-free, and sources are fixed by automorphisms. Generic shadows leak badly: 296 rooted 6-tournament classes collapse to 36 side split profiles, 56 unrooted classes, 22 score sequences, or 12 delete-root parents; even deletion decks collide. Complements HYP-2124's AP unit-clock rigidity, HYP-2125's VT trienerment taxonomy, HYP-2126's orbit-rigidity taxonomy, HYP-2127's sheaf monodromy, HYP-2131's prime-2 tangent atlas, HYP-2132/THM-401's pair-sum modulus, HYP-2129's parity sector, HYP-2130's exp/helix cyclotomic frame, and HYP-2134's orbit-functor atlas by asking which fixed-point labels survive quotienting. For LRC, the rigid labels are observer-source state, augmentation grade, pair denominator shields, endpoint owner pins, pincer escape ledgers, and L/M/R middle states. | HYP-2133, HYP-2134, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, HYP-2110, HYP-2109, HYP-2108, THM-401, THM-400, THM-385, THM-381, S589, `04-computation/tournament_rigidity_cascade_s589.py` |
| **LRC orbit-functor rigidity** | Rigidity classified by action, retained labels, and projection/monodromy defects | S590 reframes rigidity as `orbit functor + retained label + closure rule`. The same AP unit witness set is statically one orbit under `(Z/n)*`, dynamically fragmented by doubling at even `n`, reflectively paired by `j <-> -j`, and factor-localized by CRT. At n=14, witnesses `[1,3,5,9,11,13]` form one unit orbit, six singleton doubling fragments, three pincer pairs, and a CRT split where mod-7 doubling remains ordinary while mod-2 collapses. This upgrades the local/global split into static, dynamical, reflective, CRT, quotient-transport, monodromy, isostatic, spectral-character, fiber-stiffness, and source-functor rigidities, and absorbs HYP-2127's orbit-sheaf monodromy audit plus HYP-2129's parity-sector, HYP-2130's cyclotomic signals, HYP-2131's prime-2 localization, HYP-2132/THM-401's additive odd-shell modulus, and HYP-2133's labelled cascade. | HYP-2134, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, HYP-2118, HYP-2117, HYP-2116, HYP-1783, THM-401, THM-400, THM-385, THM-381, S590, `04-computation/orbit_rigidity_functor_atlas_s590.py` |
| **LRC round interval beat** | LRC-accessible regular comparisons are interval circulants, not Paley/QR beats | S591 validates the additive-beat picture: the half-turn comparator has contiguous out-neighbourhoods, so LRC lives in round tournaments. The interval circulant `R_m` is round and has `chi=2`; Paley/QR circulants split from the interval at `m>=7`, are non-round in the checked cases, and have higher chi at Paley_7. This turns THM-401's additive modulus into a tournament-level restriction: the additive interval is the beat, while unit multiplication is witness symmetry. | HYP-2141, HYP-2136, THM-401, HYP-2132, HYP-2124, S591, `04-computation/lrc_interval_vs_paley_s591.py` |
| **LRC sumset-support calculus** | Labelled support ledger after the THM-401 modulus identity | S591 formalizes the residual at `C=2n-1` as five support layers rather than raw `V+V`: speed-shell support, pair-sum shell support, actual pinch denominators, divisibility shields, and unit/nonunit/lift visibility. AP rows have perfect speed shells but pair-shell support misses shell `1`, with only `D=n` unshielded. S573 open-gap rows show why this label stack is necessary: n=7 has perfect speed shells and complete pair shells but exact denominator `33`, while n=8 misses only nonunit shell `6` and witnesses at denominator `23`. At n=14, `V*` is floor with a nonunit hole, doubled-apex gives the unit-hole `2/27` edge, and unit-shift AP shields `D=14` and becomes loose. | HYP-2135, HYP-2138, HYP-2141, HYP-2137, HYP-2132, HYP-2134, HYP-2133, HYP-2122, HYP-2118, HYP-2088, HYP-2084, HYP-2083, THM-401, THM-400, THM-397, S591, `04-computation/lrc_sumset_support_calculus_s591.py` |
| **LRC doubled-prime parity gate** | Lemoine bridge cell `N=p+2q` as one odd prime atom plus one first-even prime-core atom | The doubled prime `2q` spends exactly one dyadic factor while preserving odd core `q`; in the n=4 character frame it can break the all-odd `t=1/4` witness while remaining pairwise `chi4`-quiet, so it is boundary-active but pairwise quiet. Use bridge vertices `(N,p,q)`, not raw runners, when testing this layer. | HYP-2051, HYP-2049, HYP-1984, HYP-1983, HYP-2040, THM-391 |
| **LRC proof-route currencies** | Literature/proof-route atlas for Lonely Runner: finite product-sieve, spectral/kernel, zonotope, endpoint-pressure, and denominator-specific aperture proof currencies | The 2025-2026 fixed-case frontier is product-sieve plus finite checking through thirteen total runners; kernel and mixed-threshold tools feed endpoint pressure; S518 adds an n=145-specific CRT aperture route around unit walls `a/145`, zero-residue embryos, and observer-score descent | HYP-1989, HYP-1969, S518, S505, HYP-1812, HYP-1950, HYP-1961 |
| **LRC endpoint-aware tournament gauges** | Arc rules that insert the LRC threshold into Tournament Analysis before computing `H`, score sequence, cycles, and SCCs | Status dominance is the strongest witness detector; theta-close and antipodal-open are the best shape-only gauges; wall/slack gauges are boundary-pressure diagnostics | HYP-1974, S507, HYP-1951, HYP-1967, HYP-1968 |
| **Additive basis normal form** | Representation hypergraph with atoms, arity, local obstruction product, and carry debt | Goldbach is sparse abundance, Helfgott is ternary smoothing, Fermat polygonal is bounded cover, Zeckendorf is unique carry normal form | HYP-1953, HYP-1962, HYP-1963, HYP-1965, HYP-1966, S494, S495, S501, S502 |
| **Cauldron additive-coloring game** | Online placement of naturals into cauldrons/colors; a cauldron boils when it contains an additive witness such as `A+B=C` | Literal first-boil is weak Schur because summands are distinct (`k=3`: safe `23`, force `24`); repeated-summand comparison is classical Schur (`safe 13`, force `14`); all-boiled removal is a separate active-state sacrifice dynamic (`k=3` last boil `27` literal, `20` repeated, `25` finite-sums); finite-sums tightens first-boil to force `22` in the S618 scout | HYP-2189, T687, S618, `04-computation/cauldron_game_s618.py` |
| **Cauldron complexity/variant atlas** | Classification of which cauldron rules import Schur, A000568, LRC, unit-distance, Collatz, or game-complexity structure | S619 says base fixed-`k` last-boil is Schur-frontier hard but not A000568/LRC/unit-distance hard until extra side channels are added. Weak `k=4` first-boil frontiers exceed `250000` canonical states at `n=13`; attack-only adversarial play splits odd/even attack streams, giving B a two-term residue advantage while three/finite terms let A use `1+3+5=9` in the `2v1` case. | HYP-2190, HYP-2189, T688, T687, S619, `04-computation/cauldron_complexity_variants_s619.py`, `05-knowledge/results/cauldron_complexity_variants_s619.out`, `07-reflections/cauldron-complexity-and-variant-atlas-s619.md` |
| **Cauldron two-block adversarial schedules** | Attack-only cauldron schedules as additive residue channels | S620 compares parity, one-handicap two-block `A,BB,AA,...`, and strict `AA,BB,...`. The weak two-term advantage is controlled by whether a player's stream `S` has pair-sum self-hit `S+S` back into future targets. Parity gives A odds with no self-hit and B evens with full self-hit; one-handicap two-block gives A `{0,1}` with full self-hit and B `{2,3}` with a missing `3` lane, flipping weak `2v1` from B-at-14 to A-at-5 and weak `2v2` from B-at-14 to A-at-13. | HYP-2191, HYP-2190, T689, T688, S620, `04-computation/cauldron_two_block_adversarial_s620.py`, `05-knowledge/results/cauldron_two_block_adversarial_s620.out`, `07-reflections/cauldron-two-block-adversarial-schedules-s620.md` |
| **Cauldron block-turn minimax** | Shared-pool adversarial version of the all-boiled removal game with alternating control blocks | S621 models Delayer vs Spoiler over one active cauldron pool. Under the user's first-turn normalization `D,S,S,D,D,...`, weak distinct `k=3` drops from one-player last boil `27` to minimax `13`; `D,D,S,S,...` gives `14`. Repeated/classical gives `20 -> 11/10`, finite-sums gives `25 -> 13/14`. This is structurally distinct from the separate-pool S619/S620 attack-only games. | HYP-2192, HYP-2191, HYP-2190, HYP-2189, T690, S621, `04-computation/cauldron_block_turn_minimax_s621.py`, `05-knowledge/results/cauldron_block_turn_minimax_s621.out`, `07-reflections/cauldron-block-turn-minimax-s621.md` |
| **Prime-pair plane / pair lens** | Midpoint coordinates `(m-h,m+h)` for additive pairs, prime pairs, and pairwise tournament observables | Twin primes are the fixed ray `h=1`; Goldbach fixes `m=N/2`; Hardy-Littlewood is the local pair product; Zeckendorf is carry debt on the edge | HYP-1966, HYP-1965, HYP-1964, HYP-1962, HYP-1953, S495, S502 |
| **Opposite tournament** T^op | Reverse all arcs: T^op(u,v) = T(v,u) | H(T) = H(T^op) for Hamiltonian paths | THM-002, THM-030 |
| **Self-complementary (SC) tournament** | T ≅ T^op via some permutation | Exists only at n ≡ 0,1 (mod 4); |Aut| always odd (Moon) | THM-024, INV-043, T019 |
| **Self-converse tournament** | T has an anti-automorphism | SC ⊂ Self-converse; all circulant tournaments are SC | THM-052 |
| **Paley tournament** T_p | p ≡ 3 (mod 4) prime; i→j iff j−i is QR mod p | Vertex-transitive, self-complementary, max H(T) | T019, T053, INV-042 |
| **Regular tournament** | Every vertex has out-degree (n−1)/2 | Only at odd n; includes all Paley | THM-027, BIBD |
| **Doubly regular tournament (DRT)** | Regular + every pair of vertices has same # common out-neighbors | Unique at n=3,5,7; multiple at n≥11 | INV-068 |
| **Transitive tournament** | Acyclic; unique up to labeling | H = 1; base case for induction | OCF base case |
| **Cycle-rich tournament** | Every vertex in a directed 3-cycle | ⟹ no source/sink (Lemma Q); key for H=21 proof | THM-079, Part R |
| **R-cone** | One vertex beats/loses to all others | Rajkumar et al.; every tournament is flip-class of R-cone | T046, INV-015 |
| **Locally transitive tournament** | rank 2; sub-neighborhoods are transitive | Still has 5-cycles (67% at n=5) | T046 |
| **Vertex-transitive (VT) tournament** | Aut(T) acts transitively on vertices | Includes Paley; M = (H/n)I iff self-converse | THM-052 |
| **Frobenius tournament** F_21 | Non-circulant VT at n=21 via Frobenius group | First VT that is NOT self-converse | MISTAKE-013 |

### Paths and Cycles
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **H(T)** | # directed Hamiltonian paths | Always odd (Rédei); H = I(Ω(T), 2) (OCF) | Central |
| **Directed k-cycle** | k vertices forming directed cycle | Odd cycles only matter for OCF | THM-002 |
| **c_k(T)** | # directed k-cycles in T | c_3 = sum C(s_i, 2) formula; divisibility by C(p,k) | T036, Savchenko |
| **3-cycle matching number** mm(T) | Max # pairwise vertex-disjoint 3-cycles | mm ≤ 2 ⟹ poisoning graph argument | THM-079 Part R |
| **Base path** P_0 | n → n−1 → ⋯ → 1 | Fixed reference path for tiling model | definitions.md |
| **Oriented Hamiltonian path** | Ham path with prescribed arc directions | El Sahili: resilient under arc deletion for n≥8 | INV-099 (new) |

### Conflict Graph and Independence Polynomial
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Conflict graph** Ω(T) | Vertices = odd cycles of T; edge iff share vertex | Central to OCF | THM-002 |
| **Independence polynomial** I(G,x) | Σ_k α_k x^k; α_k = # independent sets of size k | I(Ω(T), 2) = H(T) | THM-002 (OCF) |
| **α_k** | # independent sets of size k in Ω(T) | H = 1 + 2α₁ + 4α₂ + 8α₃ + ... | OCF decomposition |
| **μ(C)** | I(Ω(T−v)|_{avoid C\{v}}, 2) | Weights in Claim A RHS | definitions.md |
| **Ω_3(T)** | Subgraph of Ω restricted to 3-cycles only | Quasi-regular; degree ≈ Johnson graph J(n,3) | INV-041 |
| **Hard-core lattice gas** | I(G,x) = partition function at fugacity x | λ=2 outside uniqueness for Ω; special structure | T006, T050 |
| **Typed independence polynomial** | I_typed(Ω; y_3, y_5, ...) | Separates cycle lengths; specializes to I(Ω,x) at y_k=x | THM-068/069/070 |

### Tiling Model
| Concept | Definition | Key Properties | Where Used |
|---------|-----------|---------------|------------|
| **Tile** (a,b) | a ≥ b+2; bit 0 = forward arc, bit 1 = backward | m = C(n−1, 2) tiles total | definitions.md |
| **Pin grid** Grid(n) | (r,c): r≥1, c≥1, r+c≤n−1 | Isomorphic to δ_{n−2} (staircase diagram) | T009 |
| **Strip** Str(k) | {(r,c): r+c=k} | k−1 tiles per strip | definitions.md |
| **Hamming weight** | # backward arcs in tiling | Bell-curve correlation with H | S15 |

---

## II. KEY THEOREMS AND RESULTS

### Proved Theorems
| ID | Name | Statement (brief) | Proof Method |
|----|------|-------------------|-------------|
| **THM-002** | OCF | H(T) = I(Ω(T), 2) | Grinberg-Stanley (arXiv:2307.05569) + specialization |
| **THM-003** | Claim B | I(Ω(T),2) − I(Ω(T−v),2) = 2Σ_{C∋v} μ(C) | A-clique argument (internal) |
| **CONJ-001** | Claim A | H(T) − H(T−v) = 2Σ_{C∋v} μ(C) | OCF + Claim B |
| **THM-016/017** | Even-odd split | Σ_{|S| even} Δ(S,R) = Σ_{|S| odd} Δ(S,R) | Internal induction |
| **THM-020** | Real roots n≤8 | I(Ω(T),x) all real negative roots for n≤8 | Claw-free + Chudnovsky-Seymour |
| **THM-024** | SC involution anti-aut | Every SC tournament has involution anti-aut | Moon + Cauchy's theorem |
| **THM-025** | Real roots FAIL n=9 | Counterexample: score [1,1,3,4,4,4,6,6,7] | Explicit construction |
| **THM-027** | BIBD maximizes H | BIBD arrangement maximizes α₁ (total directed cycles) | Counting + Jensen |
| **THM-029** | H=7 impossible | α₁=3 with i₂=0 forces c₅≥1 ⟹ α₁≥4 | Combinatorial forcing |
| **THM-030** | Transfer matrix symmetry | M[a,b] = M[b,a] for all tournaments | Induction on |W| via Walsh analysis |
| **THM-052** | Scalar M for SC VT | M = (H/n)I for self-converse VT tournaments | Aut+Anti transitivity |
| **THM-079** | H=21 impossible | H(T) ≠ 21 for ALL tournaments on ALL n | Dichotomy + base cases |

### Key Lemmas in H=21 Proof
| Part | Name | Statement | Method |
|------|------|-----------|--------|
| **Part C** | 3-matching bound | 3 disjoint 3-cycles ⟹ α₁+2α₂+4α₃ ≥ 13 | Independence counting |
| **Part G** | Base case | H=21 absent for n≤8 | Exhaustive computation |
| **Part J** | Non-cyclic removal | v not in any 3-cycle ⟹ v not in any cycle | Structural (shortest cycle argument) |
| **Part Q** | Cycle-rich no source/sink | Every vertex in 3-cycle ⟹ no source/sink | Trivial contradiction |
| **Part R** | Dichotomy Theorem | Cycle-rich n≥9: mm≥3 OR safe deletion exists | Poisoning graph DAG |

### Disproved Conjectures
| ID | Name | What Failed | Where |
|----|------|------------|-------|
| **CONJ-002** | Paley Ham paths | H(T_p) ≠ expected at p=11 | H=95095 ≠ 4455 |
| **THM-025** | Real roots all n | Fails at n=9 | Score [1,1,3,4,4,4,6,6,7] |
| **MISTAKE-010** | Hereditary maximizer | Only regular maximizers at odd n | n=5: 24/64 hereditary |
| **MISTAKE-013** | VT ⟹ SC | False at n=21 (Frobenius F₂₁) | McKay database |

---

## III. PROOF TECHNIQUES

### Internal Techniques
| Technique | Description | Used In |
|-----------|------------|---------|
| **A-clique argument** | Through-v cycles form clique in Ω; deletion = clique removal | THM-003 (Claim B) |
| **Poisoning graph** | DAG on outer vertices R; w→v iff all w's 3-cycles contain v | THM-079 Part R |
| **Swap involution** | Swap i,j positions in Ham path; matches adj(i,j) with adj'(j,i) | T031, THM-014 |
| **Bracket structure** | 4-way vertex classification (M+, M−, Z0, Z1) under arc flip | T047 |
| **Fourier decomposition** | OCF decomposes into degree-homogeneous identities | INV-050 |
| **Even-odd split** | Alternating-sum vanishing of subset contributions | THM-016/017, T040 |
| **Transfer matrix** | M[a,b] = Σ_S (−1)^|S| E_a(S)·B_b(R\S) | THM-030 |
| **Walsh spectrum analysis** | Expand M[a,b] in Walsh basis; show only even r-powers | THM-030, THM-080 |
| **Pin grid / tiling encoding** | Tournament ↔ binary string via staircase Young diagram | definitions.md |

### External Techniques Referenced
| Technique | Source | How Used |
|-----------|--------|---------|
| **Chudnovsky-Seymour (2007)** | Real roots for claw-free I.P. | THM-020: proves real roots n≤8 |
| **Heilmann-Lieb (1972)** | Matching polynomial real roots | T054: Ω(T) line graph ⟹ real roots |
| **Lindström-Gessel-Viennot (LGV)** | Lattice path determinants | T046: potential bijection route |
| **RSK correspondence** | Tableaux ↔ permutations | T009: pin grid = Young diagram |
| **Stanley-Stembridge conjecture** | e-positivity of chromatic SF | T056: via Mitrovic-Stojadinovic bridge |
| **Lovász Local Lemma** | Independence polynomial zero-free ↔ LLL | Hard-core lattice gas connection |
| **Bethe ansatz** | Statistical mechanics transfer matrices | T006: potential import |
| **Frankl matching bound** | Max k-sets with no (s+1)-matching | INV-100: bounds on Ω with bounded mm |
| **Lichiardopol's conjecture** | min outdeg ≥ (q−1)k−1 ⟹ k disjoint q-cycles | INV-098: proved for q=3 |

---

## IV. NAMED MATHEMATICIANS AND THEIR CONTRIBUTIONS

| Person(s) | Contribution | Connection to Project |
|-----------|-------------|----------------------|
| **Rédei (1934)** | Every tournament has odd # Hamiltonian paths | The foundational theorem we refine |
| **Berge** | Extended Rédei to digraphs (via permanent) | Grinberg-Stanley generalizes |
| **Grinberg & Stanley (2023)** | Rédei-Berge symmetric function; OCF via specialization | Proves THM-002 (OCF) |
| **Forcade (1973)** | F₂-invariance of H for k-block decompositions | We give first combinatorial proof |
| **Moon** | |Aut(T)| is odd for all tournaments | Used in THM-024 |
| **Chudnovsky & Seymour (2007)** | I(G,x) real roots for claw-free G | Proves THM-020 |
| **Jerrum & Patel (JLMS 2026)** | Zero-free regions for H-free graph classes | Extends real roots to subdivided-claw-free |
| **Lichiardopol** | Conjecture: min outdeg condition for disjoint cycles | Proved for q=3; context for H=21 |
| **Bang-Jensen, Bessy, Thomassé** | Proved Lichiardopol's conjecture for q=3 | INV-098 |
| **Chen & Chang (2024)** | Disjoint cycles in tournaments (JGT) | INV-099 |
| **Frankl** | Erdős matching conjecture for k=3 | INV-100 |
| **El Sahili & Ghazo Hanna (2023)** | Oriented Ham path types in tournaments | T and T^op same distribution |
| **El Sahili & El Zein (2025)** | Ham paths stable under arc deletion for n≥8 | NEW: arXiv:2512.09332 |
| **Mitrovic (2025)** | Noncommuting Rédei-Berge; deletion-contraction | arXiv:2504.20968 |
| **Mitrovic & Stojadinovic (2025)** | Chromatic ↔ Rédei-Berge bridge | arXiv:2506.08841 |
| **Savchenko (2016-2024)** | Exact c_k formulas for regular tournaments | Phase transition at n=39 |
| **Herman (2024)** | Terwilliger algebras classify DRTs to n=23 | 237 non-iso DRTs at n=27 |
| **Feng (2025)** | Dual Burnside process; Q=AB factorization | arXiv:2510.25202 |
| **Schweser, Stiebitz, Toft (2025)** | Rédei's theorem revisited | arXiv:2510.10659 |
| **Irving & Omar (2024)** | Walk generating function; det/per formula | arXiv:2412.10572 |
| **Rajkumar et al.** | Flip classes; R-cones; sign rank | arXiv:2110.05188 |
| **Hetyei (2017)** | Alternation acyclic tournaments; Genocchi numbers | arXiv:1704.07245 |
| **Alon (1990)** | max H(T) = Θ(n!/2^{n−1}) via random regular | Upper bound theory |
| **Szele** | max H(T) ≥ n!/2^{n−1} | Lower bound on max H |
| **Guo, Gutin et al. (2026)** | Forward arc maximization in tournament generalizations | NEW: arXiv:2602.10713 |
| **Bencs & Buys (2025)** | Zero-free regions for hypergraph I.P. | Random Struct. Algorithms |
| **Galvin, McKinley, Perkins, Sarantis, Tetali** | Zeroes of hypergraph independence polynomials | CPC 2024 |
| **Stembridge** | Self-evacuating SYT count | T008 connection |
| **Chapman (2001)** | ASM ↔ monotone triangles | INV-020 |
| **Striker (2011)** | Unifying poset perspective (ASM, PP, Catalan) | INV-020 |
| **Leake & Ryder** | Same-phase stability for claw-free | THM-078 connection |

---

## V. CROSS-FIELD CONNECTIONS

### Established Connections
| Field | Connection | Strength | Reference |
|-------|-----------|----------|-----------|
| **Acyclic orientations / Stanley reciprocity** | Tournaments are orientations of K_n; χ_{K_n}(−1) = (−1)^n·n!; U_T connects to chromatic SF | Strong | Stanley 1973, Mitrovic-Stojadinovic |
| **Statistical mechanics** | I(G,x) = hard-core partition function at fugacity x | Strong | T006, T050 |
| **Representation theory** | Pin grid = staircase Young diagram δ_{n−2}; hook lengths all odd | Medium | T009 |
| **Algebraic combinatorics** | U_T = Rédei-Berge symmetric function | Strong | Grinberg-Stanley |
| **Chromatic polynomial theory** | Chromatic SF ≈ Rédei-Berge at poset level | Strong | Mitrovic-Stojadinovic |
| **Hopf algebras / Deletion-contraction** | W_X = W_{X\e} - W_{X/e}↑ for noncommuting Redei-Berge | Strong | Mitrovic arXiv:2504.20968 |
| **Number theory** | Paley tournaments via quadratic residues; 1729 appearance | Medium | T019, T025 |
| **2-adic analysis** | H(T) mod 2^k tower from I(Ω,2^k) | Speculative | T007 |
| **Plane partitions / ASMs** | 2^{m²} count; TSSCPP connection | Weak | T008 |
| **Matching theory** | 3-cycle matchings ↔ hypergraph matching | Strong | Frankl, Lichiardopol |
| **Spectral graph theory** | Gauss sums in Paley tournament eigenvalues | Medium | T024 |
| **Design theory (BIBD)** | Regular tournament 3-cycles as block designs | Strong | THM-027 |
| **Walsh/Fourier analysis** | Spectral OCF: hat{H}[S] = Σ cycle Walsh terms; counting identity | Strong | THM-077, THM-081 |
| **GLMY path homology** | Digraph homology; Tang-Yau Fourier on circulant digraphs | Medium | arXiv:2602.04140 |
| **Root polytope / h*-polynomial** | Deletion-contraction monotonicity for digraph polytopes | Speculative | Kálmán-Tóthmérész 2024 |
| **Oriented matroid circuits** | Tournament cycles = circuits of oriented matroid on K_n | Medium | arXiv:2501.00108 |

### Speculative / Novel Connections (kind-pasteur-S34)

#### 1. GEOMETRIC MODEL: Tournament as Flow on Simplex
**Idea:** A tournament on n vertices defines an orientation of the complete graph K_n. The complete graph K_n is the 1-skeleton of the (n−1)-simplex Δ^{n−1}. An orientation of K_n is equivalent to a discrete vector field on the 1-skeleton. Hamiltonian paths are maximal gradient-like trajectories. The parity constraint (Rédei: H odd) becomes a topological statement about the Euler characteristic of the space of such trajectories.

**Connection to Morse theory:** A discrete Morse function on a simplicial complex assigns values to faces such that the resulting discrete gradient has special properties. A tournament is a discrete Morse function on the 0-skeleton extended to the 1-skeleton. H(T) counts the number of gradient paths of maximum length.

**Potential:** Could import discrete Morse theory (Forman) machinery. The independence polynomial I(Ω,x) might have a topological interpretation as a Poincaré polynomial of some associated complex.

#### 2. TROPICAL GEOMETRY: Permanent as Tropical Determinant
**Idea:** The permanent of a matrix is the tropical determinant (max-plus algebra determinant). For tournaments, H(T) involves counting paths, which is related to permanents of {0,1} matrices. In tropical geometry, the permanent computes the weight of the optimal assignment problem. The tournament arc matrix A defines a tropical hypersurface whose Newton polygon encodes the cycle structure.

**Connection:** Irving-Omar's formula ham(D) = Σ_S det(Ā[S])·per(A[S^c]) mixes determinants and permanents. This is exactly the structure that appears in tropical geometry when studying mixed volumes.

#### 3. HOMOLOGICAL ALGEBRA: Cycle Complex of a Tournament
**Idea:** Define a chain complex C_*(T) where C_k = free abelian group on directed k-cycles, with boundary ∂ = alternating sum of "face" maps (remove one vertex from cycle). The homology H_*(C_*(T)) could encode the independence structure of Ω(T).

**Specific conjecture:** H_0(C_*(T)) = Z (connected), and the Euler characteristic χ(C_*) = Σ (−1)^k rank(C_k) relates to H(T) mod 2.

**Connection to persistent homology:** As we vary a threshold on arc "strength" (e.g., score difference), the cycle complex changes. The persistence diagram could reveal the hierarchical structure of tournament cycles.

#### 4. QUANTUM GROUPS: Tournament as R-matrix
**Idea:** A tournament T defines a solution to the Yang-Baxter equation if we set R(i,j) = q when i→j and R(i,j) = q^{-1} when j→i. The trace of the resulting braid representation counts something related to Hamiltonian paths. The quantum group U_q(sl_2) representation theory at q = root of unity could connect to the 2-adic tower (T007).

**Connection to knot invariants:** Tournaments on n vertices define braids on n strands. The Jones polynomial of the resulting closure might encode H(T). Since the Jones polynomial has connections to statistical mechanics (Potts model), this could unify T006 (hard-core lattice gas) with the tournament structure.

#### 5. ERGODIC THEORY: Tournament Shift Dynamics
**Idea:** The arc-flip operation (reverse one arc) defines a random walk on the space of tournaments. The mixing time of this walk, the spectral gap of the transition matrix, and the stationary distribution all encode information about H(T). The OCF identity H(T) = I(Ω(T), 2) might be a fixed-point equation for some natural dynamics.

**Connection to Markov chain Monte Carlo:** The Jerrum-Patel results on zero-free regions of I(G,x) are motivated by MCMC sampling. Their algorithms use the polynomial method to sample independent sets. For tournament conflict graphs, the specific structure of Ω(T) might allow faster sampling.

#### 6. MODULAR FORMS: Paley Tournament L-functions
**Idea:** The Paley tournament T_p is defined via the Legendre symbol χ_p. The adjacency matrix eigenvalues involve Gauss sums (T024). Define L_T(s) = Σ_{n≥1} a_n n^{−s} where a_n encodes cycle counts. For Paley tournaments, L_{T_p}(s) might factor through Dirichlet L-functions L(s, χ_p).

**Connection:** The factorization H(T_p) = 55 × 1729 at p=11 could reflect the factorization of L(1, χ_11) or related special values. The sequence H(T_p)/|Aut| = 1, 9, 1729 should be checked against modular form coefficient tables.

#### 7. ALGEBRAIC K-THEORY: Tournament Grothendieck Group
**Idea:** Define K_0(Tour_n) as the free abelian group on isomorphism classes of n-tournaments modulo the relation [T] = [T'] + [T''] when T is obtained by "gluing" T' and T'' at a vertex. The OCF H(T) = I(Ω(T), 2) defines a ring homomorphism K_0(Tour_n) → Z.

**Connection:** The vertex deletion formula H(T) − H(T−v) = 2Σμ(C) is a Euler characteristic relation in this Grothendieck group.

#### 8. INDEPENDENCE COMPLEX TOPOLOGY: Ω(T) as Simplicial Complex
**Idea:** The independence complex Ind(Ω(T)) is the simplicial complex whose faces are independent sets of Ω(T). For claw-free graphs (n≤8), Ind(G) is homotopy equivalent to a wedge of spheres (Engström's work). The reduced Euler characteristic of Ind(Ω(T)) equals (−1)^d · I(Ω(T), −1) where d is the dimension. At x=2, I(Ω,2) = H(T) counts 2-colored independent sets, which has a topological interpretation as the number of simplices weighted by a 2-coloring.

**Connection to persistent homology:** Define a filtration on Ind(Ω(T)) by the size of cycles (3-cycles first, then 5-cycles, etc.). The resulting persistence diagram tracks how the topology changes as longer cycles are added. This could reveal WHY certain H values are impossible (gaps correspond to topological obstructions).

**Concrete test:** Compute the homology of Ind(Ω(T)) at n=5,6,7. Check if the Betti numbers relate to α_k coefficients. For the H=21 gap: does the corresponding independence complex have a topological obstruction?

#### 9. SIGNED HP PERMANENT AS PFAFFIAN VARIANT
**Idea (from opus S35c):** S(T) = Σ_P Π B[P_i][P_{i+1}] is a "path permanent" — distinct from the standard permanent (which sums over cycle covers). The standard perm(B) = 0 at odd n for skew-symmetric B, but S(T) ≠ 0. S(T) captures Hamiltonian PATH information invisible to both det(B) and perm(B).

**Key discovery:** S(T) mod 2^{n−1} depends ONLY on n (Universal Congruence Theorem, THM-H). This means the "residual" c_0 = S(T)/2^{n−1} has universal fractional part determined by n alone. At n=5: c_0 ∈ Z (integer). At n=7: c_0 has fractional part 3/4 always. This 2-adic structure connects to the 2-adic tower (T007).

**Connection to knot theory:** The signed adjacency matrix B is skew-symmetric, like the Seifert matrix of a knot. The Pfaffian of B gives the Alexander polynomial value. S(T) is a different "Pfaffian-like" invariant that uses paths instead of matchings. Could S(T) be a tournament analogue of a knot invariant?

#### 10. RANDOM MATRIX THEORY: Tournament Eigenvalue Distribution
**Idea:** The skew-adjacency matrix S_T (with S_{ij} = 1 if i→j, S_{ij} = −1 if j→i, S_{ii} = 0) has purely imaginary eigenvalues. For random tournaments, the empirical spectral distribution converges to the semicircle law scaled by √n. The Pfaffian Pf(S_T) relates to the number of perfect matchings of the underlying graph, which connects to the cycle structure.

**Connection to H(T):** H(T) = permanent of a certain {0,1} matrix. The permanent-determinant gap (permanent is #P-hard, determinant is polynomial) is the central problem in algebraic complexity. For tournament matrices, the permanent might have special structure due to the constraint T(i,j) + T(j,i) = 1.

#### 11. INFORMATION THEORY: Tournament Entropy
**Idea:** Define the entropy of a tournament as H_ent(T) = −Σ_P (1/H(T)) log(1/H(T)) = log(H(T)) (uniform distribution over Hamiltonian paths). The OCF H(T) = I(Ω(T), 2) then gives log H(T) = log I(Ω(T), 2), connecting tournament entropy to the free energy of the hard-core model. The 2-adic valuation v_2(H(T)) measures how "2-divisible" the path count is — a kind of 2-adic entropy.

**Connection to coding theory:** Tournament tilings (binary strings of length m = C(n−1,2)) form a code where each codeword represents a tournament. The Hamming weight structure of this code correlates with H(T) (bell curve, T053). Error-correcting properties of the "tournament code" could import coding theory machinery.

#### 12. CLUSTER ALGEBRAS: Tournament Mutations
**Idea:** In cluster algebra theory, mutations on quivers (directed graphs) produce new quivers. A tournament IS a quiver (complete, no 2-cycles). Cluster mutations preserve certain algebraic structure. The OCF identity H(T) = I(Ω(T), 2) could be a cluster variable identity, with arc-flip = cluster mutation.

**Connection:** Fomin-Zelevinsky cluster algebras on complete quivers have been studied. The exchange graph (quivers related by mutations) for tournaments could be the arc-flip graph. The Laurent phenomenon (cluster variables are Laurent polynomials) might constrain H(T) values, potentially explaining gaps.

**Concrete test:** Check if the cluster variable associated to a tournament quiver equals H(T) or I(Ω(T), x) for some specialization.

**COMPUTATIONAL RESULT (kind-pasteur-S34):** Quiver mutation at vertex k preserves H(T) whenever the result is a tournament! Verified exhaustive n=4 (64/64), n=5 (640/640), sampled n=7 (1126/1126). Zero violations. The mutation only produces a tournament when k is a source (deg 0) or sink (deg n−1), or at specific score-class positions. At n=5, 62.5% of H-preserving mutations produce NON-isomorphic tournaments — this is a genuinely new symmetry, not just relabeling!

**MECHANISM IDENTIFIED:** Quiver mutation at vertex k produces a tournament iff k is a source (deg 0) or sink (deg n−1). The mutation simply reverses ALL arcs incident to k. This preserves H because source/sink vertices are forced to start/end every Ham path, and the reversal creates a bijection between start-paths and end-paths. The cluster algebra connection provides an algebraic framework for this symmetry.

**Scripts:** `tournament_quiver_mutation.py`, `quiver_mutation_h_invariance.py`, `quiver_mutation_mechanism.py`

#### 13. MATROID THEORY: Tournament Matroid
**Idea:** Define a matroid on the arcs of a tournament where independent sets are "acyclic subsets" (subsets of arcs that contain no directed cycle). The rank function r(S) = |S| − (# cycles in S) could encode the cycle structure. The Tutte polynomial T_M(x,y) of this matroid at specific evaluations could give H(T) or cycle counts.

**Connection to Rédei-Berge:** The Rédei-Berge symmetric function U_T is a "symmetric function analogue" of the Tutte polynomial. The chromatic-Rédei-Berge bridge (Mitrovic-Stojadinovic) makes this explicit: the chromatic symmetric function IS the Tutte polynomial in disguise.

#### 14. OPERAD THEORY: Composition of Tournaments
**Idea:** Tournaments form an operad under substitution: given T on n vertices and T_1,...,T_n, the substitution T(T_1,...,T_n) replaces each vertex i by T_i. This operad structure governs how H(T) decomposes. The vertex-deletion formula H(T) − H(T−v) = 2Σμ(C) is a derivation in this operad.

**Connection:** The operad of associahedra (Stasheff) governs compositions of operations. Tournament substitution is the "directed" analogue. The free operad generated by tournaments might have a presentation that encodes the OCF identity.

#### 15. MIRROR SYMMETRY: Tournament Duality
**Idea:** H(T) = H(T^op) is a "mirror symmetry" for tournaments — the path count is invariant under arc reversal. In mirror symmetry, the Hodge numbers h^{p,q} of a Calabi-Yau manifold X equal h^{n−p,q} of the mirror X^∨. For tournaments, the "Hodge numbers" could be the coefficients of the W-polynomial, and the mirror is T^op.

**Connection to THM-030:** The transfer matrix symmetry M[a,b] = M[b,a] is analogous to the symmetry of the Hodge diamond. The W-polynomial W(r) having only even powers of r is analogous to the vanishing of odd Hodge numbers.

#### 16. PERSISTENT HOMOLOGY OF TOURNAMENT FILTRATIONS
**Idea:** Order tournaments by "complexity": T_0 (transitive) → T_1 (one arc flip) → ⋯ → T_m (all arcs flipped). Each step changes Ω(T), and the persistence barcode of this filtration tracks which cycles are born and die. The OCF identity H(T) = I(Ω(T), 2) says the partition function is preserved through the filtration (since H is determined by Ω).

**Novel observation:** The arc-flip graph on tournaments is the Boolean lattice {0,1}^m where m = C(n−1,2). The OCF is a function on this lattice. Its "topological complexity" (number of critical points in the Morse sense) could measure the difficulty of proving OCF.

#### 17. REPRESENTATION STABILITY: Tournament Invariants as n→∞
**Idea:** Church-Ellenberg-Farb representation stability: for families of S_n-representations {V_n}, the decomposition into irreducibles stabilizes. The space of tournament invariants (functions T ↦ f(T) invariant under relabeling) is an S_n-representation. Does it stabilize?

**Connection:** The Rédei-Berge symmetric function U_T lives in Λ_n (symmetric functions of degree n). As n grows, the relevant pieces of Λ_n could stabilize, giving eventually-polynomial formulas for H(T) statistics.

#### 18. PHYSICS: TOURNAMENT AS CAUSAL STRUCTURE
**Idea:** A tournament defines a total preorder on pairs — exactly one of i→j or j→i holds, like a causal relation. In physics, causal sets are partially ordered sets modeling spacetime. A tournament is a "maximally connected" causal set. Hamiltonian paths are "world lines" visiting every event exactly once. The parity of H(T) (always odd, Rédei) is a "topological charge" of the causal structure.

**Connection to loop quantum gravity:** The spin foam models of loop quantum gravity use labeled graphs and their amplitudes. A tournament with weighted arcs defines a spin foam amplitude. The independence polynomial I(Ω,x) at x=2 could be a partition function in this framework.

#### DELETION-CONTRACTION FRAMEWORK (Mitrovic 2025 — KEY LEAD)
**Paper:** arXiv:2504.20968 "The Redei-Berge function in noncommuting variables"

**Core result:** The noncommuting Redei-Berge function W_X satisfies:
- **Deletion-contraction:** W_X = W_{X\e} - W_{X/e}↑ for any edge e
- **Cycle decomposition (Thm 3.16):** For edges e₁,...,e_k forming a cycle: W_X = Σ_{S⊆{e₁,...,e_k}, S≠∅} (-1)^{|S|-1} W_{X\S}
- **Tournament formula (Cor 3.12):** W_X = Σ_{σ: all cycles odd} 2^{ψ(σ)} p_{Type(σ)} — THIS IS OCF!
- **Contraction:** For e=(u,v), X/e merges u,v into vertex e with edges (w,e) iff (w,u)∈E, (e,w) iff (v,w)∈E

**Why this matters:** Deletion-contraction gives an INDUCTIVE proof framework for OCF. The ordinary (commutative) Redei-Berge function does NOT satisfy deletion-contraction, but the noncommuting version does. Specializing to commutative variables at the end recovers all tournament results.

**Connection to transfer matrix:** The contraction X/e produces a digraph on n-1 vertices. If we can track how M[a,b] changes under contraction, this could prove THM-030 (transfer matrix symmetry) inductively.

**COMPUTATIONAL VERIFICATION (kind-pasteur-S34):** H(T) = H(T\e) + H(T/e) holds UNIVERSALLY — 100% at n=4 (384/384 edge tests) and n=5 (10240/10240). The contraction convention: w inherits IN-edges from tail u, OUT-edges from head v. This is the commutative specialization of Mitrovic's W_X = W_{X\e} - W_{X/e}↑ (sign change because W counts signed paths). The reverse convention (IN from head, OUT from tail) does NOT work.

**Implications:** This gives an inductive proof framework: H(T) on n vertices reduces to Ham path counts on (n-1)-vertex digraphs. Combined with base cases, this could yield a new proof of OCF and potentially the H=21 impossibility.

**Status:** VERIFIED COMPUTATIONALLY. Ready for algebraic proof (should follow from Mitrovic's theorem by specialization).

#### 19. ACYCLIC COMPLEX OF TOURNAMENT ARCS
**Idea:** The acyclic complex Acyc(T) = {arc subsets containing no directed cycle} is a simplicial complex. Its f-vector encodes the cycle structure.

**Computational results (n=4):**
- H=1 (transitive): #acyclic = 2^6 = 64, χ(Acyc) = 0
- H=3 (one 3-cycle): #acyclic = 56 = 8×7, χ(Acyc) = 0
- H=5 (all 3-cycles): #acyclic = 49 = 7², χ(Acyc) = -1

**Pattern:** #acyclic = (8-c₃)² at n=4 where c₃ = number of 3-cycles. Verifiable: c₃=0→64, c₃=1→49, c₃=4→16... wait, H=3 has c₃=1 and #acyclic=56≠49. So the pattern is more subtle. The Euler characteristic χ(Acyc) appears to distinguish tournament classes.

**Connection to Rédei:** The acyclic complex dimension = n-1 (acyclic tournaments = spanning trees). The top-dimensional acyclic subsets are exactly the acyclic tournaments on subsets, connecting to linear extensions and Rédei's theorem.

#### 20. EHRHART THEORY: Independence Polytope of Omega(T)
**Idea:** The independence polytope P_{Omega} of the conflict graph has vertices = characteristic vectors of independent sets. The Ehrhart polynomial of P_{Omega} is exactly I(Omega(T), x) evaluated at non-negative integers. So H(T) = I(Omega, 2) = |2P_{Omega} ∩ Z^n| counts lattice points in the second dilate. The Ehrhart reciprocity theorem: I(Omega, -x) = (-1)^dim * I_interior(Omega, x) connects I(Omega, -1) to interior lattice points.

**Connection to parity:** H(T) = I(Omega, 2) is always odd. By Ehrhart theory, the h*-polynomial of P_{Omega} encodes the lattice-point distribution. Real-rootedness of I(Omega, x) (proved n<=8) implies unimodality of the h*-vector.

#### 20. BIRKHOFF POLYTOPE: Tournament-Doubly Stochastic Connection
**Idea:** A doubly stochastic matrix is a convex combination of non-identity permutation matrices iff its inner product with each generalized transitive tournament matrix is >= 1 (arXiv:2406.16284). This connects tournament matrices to the face structure of the Birkhoff polytope B_n. The permanent of a doubly stochastic matrix (van der Waerden conjecture, proved) gives lower bounds; tournament adjacency matrices A satisfy A + A^T = J - I, placing them in a specific affine slice of the Birkhoff polytope.

#### 21. GRUJIC-STOJADINOVIC: Redei-Berge Hopf Algebra (arXiv:2402.07606)
**Established connection (not speculative):** The Redei-Berge symmetric function U_T generates a combinatorial Hopf algebra on digraphs. The induced Redei-Berge polynomial satisfies deletion-contraction (like chromatic polynomial). Berge's theorem on Hamiltonian path parity is a consequence of the reciprocity formula. Tournaments are closed under restriction and products, generating a Hopf subalgebra. The forward-edge polynomial F(T,x) is the X-descent polynomial in this framework. Palindrome F_k = F_{n-1-k} follows from path reversal symmetry.

#### LEE-YANG ZEROS OF F(T,x) (opus S44 — MAJOR)
**Discovery:** The forward-edge polynomial F(T,x) = Σ_P x^{fwd(P)} has zeros that come in RECIPROCAL PAIRS (palindrome: if r is a root, so is 1/r). Key findings:
- Zeros cluster at angles ±2π/3 on the unit circle
- H=9 at n=5: ALL zeros on the unit circle ("Lee-Yang tournament")
- F(T,ω) is REAL at n=7 where ω = e^{2πi/3} (palindrome forces S₁=S₂)
- F(T,i) is PURE IMAGINARY at n=7 (palindrome forces S₀=S₂)
- F(T,ω) ≡ 0 mod 9 at n=7 (universal divisibility)
- F(T,i) ≡ 0 mod 16i at n=7 (universal divisibility)
- **Palindrome phase theorem:** F(ζ_k) lies on a fixed ray in C for all tournaments

**Connection to statistical mechanics:** The Lee-Yang theorem says zeros of ferromagnetic Ising partition functions lie on |z|=1. F(T,x) is a partition function with "activity" x weighting forward edges. The zeros-on-unit-circle property at H=9 is EXACTLY Lee-Yang behavior. This suggests F(T,x) is a partition function of a ferromagnetic-like model.

**Scripts:** f_poly_zeros_leeyang.py, f_poly_roots_of_unity.py

#### 22. GLMY PATH HOMOLOGY OF TOURNAMENTS (Grigoryan-Lin-Muranov-Yau)
**Established theory (2012+):** GLMY path homology is a homology theory specifically for digraphs. Defines a path complex where chains are formal sums of directed paths. Has Künneth formula, homotopy invariance, Mayer-Vietoris. Tang-Yau (2026, arXiv:2602.04140) compute path homology of circulant digraphs using Fourier decomposition — directly applicable to Paley tournaments T_p (which are circulant). Their "symbol-matrix recipe" uses the shift automorphism τ, and results depend on whether n is prime (exactly the Paley case!).

**Connection to OCF:** Path homology detects cycles in directed paths. The boundary map ∂ on 2-chains maps directed triangles → edges. The kernel and image of ∂ encode exactly the cycle structure that OCF uses. H_1(T) (first path homology) should relate to the cycle space of T — and thus to Ω(T) and α_1.

**Testable conjecture:** rank(H_1(T)) = α_1(T) = # independent odd cycles? Or dim(H_1) = c_3(T)?

**Status:** STRONG LEAD — path homology is the natural homological framework for directed cycles.

#### 23. EXTENDED ROOT POLYTOPE AND DELETION-CONTRACTION (Kálmán-Tóthmérész 2024)
**Paper:** arXiv:2409.18902. The extended root polytope of a directed graph has h*-polynomial coefficients that are MONOTONE under deletion and contraction: they never increase. This is a digraph version of the matroid minor monotonicity.

**Connection:** The Mitrovic deletion-contraction for W_X (connection 19/DC framework) operates on the Rédei-Berge function. The Kálmán-Tóthmérész deletion-contraction operates on the h*-polynomial of the root polytope. If these two deletion-contraction relations are COMPATIBLE, the h*-polynomial of the root polytope of T should encode information about H(T).

**Testable:** Compute the extended root polytope of small tournaments and check if its h*-polynomial relates to I(Ω(T), x).

#### 24. CONVERSE INVARIANT DIGRAPH POLYNOMIALS (Ai-Gutin-Lei-Yeo-Zhou 2024)
**Paper:** arXiv:2407.17051. Introduces a new digraph polynomial to give necessary conditions for an oriented graph to be converse invariant (# copies of D in T = # copies of -D for all tournaments T). They characterize which tree orientations are converse invariant.

**Connection:** H(T) = H(T^op) means Hamiltonian paths are converse invariant. The authors' polynomial tests this for all substructures simultaneously. Their polynomial could encode the FULL structure that determines H, not just the path count.

#### 25. SYMMETRIC EDGE POLYTOPES OF TOURNAMENTS
**Background:** The symmetric edge polytope (SEP) of a graph G has vertices {e_i - e_j : ij ∈ E(G)} ∪ {e_j - e_i : ij ∈ E(G)}. For a tournament T, every edge appears in exactly ONE direction, so the SEP simplifies: vertices = {e_i - e_j : i→j in T}. This is a subset of the full K_n SEP.

**Connection to Ehrhart:** Davis-Higashitani-Ohsugi (2024, arXiv:2401.03383) study generalized SEPs via regular matroids. The h*-polynomial is symmetric. The Ohsugi-Tsuchiya conjecture (γ-nonnegativity for graph SEPs) could constrain tournament arc structure.

**Novel idea:** The volume of the tournament SEP might encode H(T). The symmetric edge polytope of K_n has volume n^{n-2} (Cayley formula analog). For tournament subpolytopes, the volume might factor through OCF.

#### 26. ORIENTED MATROID CIRCUIT POLYTOPES
**Background:** arXiv:2501.00108 (2025) studies polytopes from signed circuits of oriented matroids. A tournament defines an oriented matroid on the complete graph edges. The circuit polytope captures the cycle structure algebraically.

**Connection:** The circuits of the tournament oriented matroid are exactly the DIRECTED CYCLES — the vertices of Ω(T). So the circuit polytope is a geometric realization of the cycle data that I(Ω,x) encodes combinatorially.

**Testable:** Compute the circuit polytope for small tournament oriented matroids and check if Ehrhart polynomial = I(Ω, x) or a transform thereof.

---

## VI. POLYNOMIAL OBJECTS

| Polynomial | Definition | Key Properties | Where |
|-----------|-----------|---------------|-------|
| **I(Ω(T), x)** | Independence polynomial of conflict graph | I(Ω,2) = H(T); real roots n≤8 | THM-002, THM-020 |
| **U_T** | Rédei-Berge symmetric function | U_T ∈ Λ_n; p-positive for tournaments | Grinberg-Stanley |
| **G_T(t, x)** | Two-variable generating function | t^m P(u,x), u = t+1/t; P(2,x) = n! | THM-064, S28 |
| **P(u, x)** | Reduced polynomial from G_T | p_m(x) = I(Ω(T), x) as leading coeff | S28 |
| **Q_T(w)** | Size-weighted I.P.: u_T(√w)/√w | All real roots n≤8; fails n≥9 with I(Ω) | THM-078 |
| **W(z)** | Walsh generating function | W(i/2) = (−2)^m p_0(2) | THM-064 |
| **M[a,b]** (transfer matrix entry) | Σ_S (−1)^|S| E_a(S)·B_b(R\S) | Symmetric; Walsh spectrum has 2^s factors; det(M)≠0 at n=5 | THM-030, THM-080 |
| **hat{t_k}[S]** (Walsh cycle spectrum) | (1/2^k) Σ_{k-cycles C⊇S} (−1)^{asc(S,C)} | Each cycle contributes ±1/2^k; sign from traversal geometry | THM-081 |
| **hat{H}[S]** (Walsh of H) | (−1)^{asc(S)} · 2^r · (n−d)!/2^{n−1} | Spectral OCF: equals 2·hat{t3}[S] + 2·hat{t5}[S] + 4·hat{α₂}[S] + ... | THM-077, THM-081 |
| **E_T^perm(t)** | All-permutation forward-edge polynomial | G_T(t; 2,2,...) specialization | S28 |
| **S(T)** (signed HP permanent) | Σ_P Π B[P_i][P_{i+1}], B=2A−1 | S=0 at even n; S mod 2^{n−1} universal | THM-A through THM-H |
| **W(r)** (W-polynomial) | tr(M(r)); only even r-powers | c_{n−1}=n!, c_{n−3} depends on t₃ | signed-hp-permanent-skeleton |
| **D_k** | Σ_P C(forward(P), k) | D_k mod 2^{n−1−k} is universal | THM-H |
| **k-fold partition identity** | Σ_{perms of 2k vertices} Π A[P_{2i-1}][P_{2i}] = (2k)!/2^k | PROVED (THM-I); key to universal congruence | opus-S35c11 |
| **THM-J** (S universality criterion) | S mod 2^{n−1} is T-independent iff n−3 is 0 or 2^k | Via Legendre: v₂((n−3)!) = (n−3)−s₂(n−3); need s₂(n−3) ≤ 1 | opus-S35c11 |
| **F(T,x)** (forward-edge polynomial) | Σ_P x^{fwd(P)} over all n! permutations | PALINDROMIC: F_k = F_{n−1−k} via path reversal | opus-S43b |
| **Palindrome consequences** | H = F_0 = F_{n−1}; S = (−1)^{n−1}F(−1); D_1 = (n−1)/2·n! | Even n: F(−1)=0 ⟹ S=0 (Rédei). G(−1) = F(0) = H | opus-S43b |
| **Eulerian polynomial connection** | F(T,x) = X-descent polynomial in Grujic-Stojadinovic Hopf | Connects to Redei-Berge Hopf algebra (arXiv:2402.07606) | opus-S43b |

---

## VII. GRAPH-THEORETIC HIERARCHY OF Ω(T)

| Property | Holds for n ≤ | Fails at n = | Reference |
|----------|-------------|-------------|-----------|
| Complete | 5 | 6 | T028 |
| Interval | 5 | 6 (13.9%) | T055 |
| Chordal | 5 | 6 | T055 |
| Line graph | 5 | 6 (K₅−e found, 45%) | INV-032 |
| Comparability | 6 | 7 (1%) | T049 |
| Quasi-line | 7 | 8 (49%) | INV-032 |
| Perfect | 7 | 8 (53.8%) | INV-032 |
| Claw-free | 8 | 9 (86%) | INV-032 |
| S_{1,1,1}-free | 11 | 12 | INV-032 |
| S_{2,1,1}-free | 9 | 10 (92%) | INV-032 |
| Real roots of I(Ω,x) | 8 (proved) | 9 (extremely rare failure) | THM-020, THM-025 |

---

## VIII. PERMANENT GAPS IN H-SPECTRUM

| Value | Status | Method |
|-------|--------|--------|
| H = 7 | IMPOSSIBLE for all n | THM-029: α₁=3 with i₂=0 forces extra cycles |
| H = 21 | IMPOSSIBLE for all n | THM-079: Dichotomy + poisoning graph + base cases n≤8 exhaustive |
| H = 23 | ACHIEVABLE at n=6 | alpha=(1,11,0,0) or (1,9,1,0) — NOT a gap |
| H even | IMPOSSIBLE at odd n | OCF: I(Ω,2) always odd when n odd |
| H = 2 (mod 4) | OPEN at even n | Need systematic enumeration |

### Alpha-Decomposition Constraints (kind-pasteur-S34)
H = 1 + 2α₁ + 4α₂ + 8α₃ + ... where α_k = # independent sets of size k in Ω(T).

| H | Required α-decompositions | Why impossible/possible |
|---|--------------------------|------------------------|
| 7 | (1,3,0,0) or (1,1,1,0) | THM-029: α₁=3 with i₂=0 ⟹ c₅≥1 ⟹ α₁≥4; (1,1,1,0) needs 1 cycle + 1 disjoint pair from 1 cycle = impossible |
| 21 | 12 decompositions | All ruled out by THM-079 dichotomy |
| 23 | (1,11,0,0), (1,9,1,0), etc. | ACHIEVABLE at n=6 |
| n=5 | α₂ = 0 always (Ω complete) | I(Ω,2) = 1 + 2α₁ only |
| n=6 | α₂ ∈ {0,1,2,4}; α₃ = 0 | First α₂ contribution appears |
| n=7 | α₂ ∈ {0,1,...,7}; α₃ ∈ {0,1} | Full 3-level structure |

---

## IX. 2-ADIC STRUCTURE AND SIGNED HP PERMANENT

The signed HP permanent S(T) = Σ_P Π B[P_i][P_{i+1}] (B = 2A − 1 skew-symmetric) reveals deep 2-adic structure:

| Property | Statement | Where |
|----------|-----------|-------|
| **S(T) = 0 at even n** | Reversal involution: (-1)^{n-1} = -1 pairs cancel | THM-A |
| **S mod 2^{n-1} universal** | Depends only on n (not T), proved via D_k decomposition | THM-H |
| **S universality criterion** | Universal iff s₂(n−3) ≤ 1 (n−3 is 0 or power of 2) | THM-J |
| **Universal n values** | n ∈ {3, 5, 7, 11, 19, 35, 67, 131, ...} | THM-J |
| **Non-universal n values** | n ∈ {9, 13, 15, 17, 21, 23, ...} — S mod 2^{n-1} depends on t₃ parity | THM-J |
| **S = 0 possible** | Only at n ≡ 1 mod 4 (verified n=5,9) | THM-H pattern |
| **c₀ = S/2^{n-1}** | Fractional part: 0 at n=5, 3/4 at n=7, {0,1/2} at n=9 | opus-S35c11 |
| **n mod 4 chain** | C(n,2) parity → Σ in(v)·out(v) parity → D₂ parity → S parity | Structural |
| **D_k linearity in t₃** | D₂, D₃ linear in t₃; D₄+ depend on higher invariants (t₅, bc₃₃) | opus-S35c11 |
| **k-fold partition identity** | Non-adjacent D_S = n!/2^k (master identity for universality) | THM-I |

### Key Insight: Legendre's Formula Connection
The 2-adic valuation v₂((n-3)!) = (n-3) − s₂(n-3) where s₂ is binary digit sum. For S universality we need v₂(2·(n-3)!) ≥ n-3, i.e., s₂(n-3) ≤ 1. This connects tournament theory to the binary representation of n — an unexpected number-theoretic input.

### Open 2-Adic Questions
- Does S universality extend to higher 2-adic levels (S mod 2^n, 2^{n+1}, ...)?
- At non-universal n, does S depend on invariants beyond t₃?
- Is there a p-adic analogue for odd primes p?

### p-adic Structure for Odd Primes (kind-pasteur-S34)
**H mod 3:** From OCF, H mod 3 = (1 + 2α₁ + α₂ + 2α₃ + ...) mod 3 (since 2^k mod 3 cycles 2,1,2,1,...).
- At n=4: H mod 3 = (1 + 2c₃) mod 3 — UNIQUE per c₃. Formula: c₃=0→1, c₃=1→0, c₃=2→2
- At n=5: H mod 3 = (1 + 2α₁) mod 3 — still unique since α₂=0 (Omega complete)
- At n=6+: H mod 3 depends on both α₁ and α₂

**H mod 7 = 0:** IMPOSSIBLE at n≤6 (all achievable H are {1,3,5,...,45}, none divisible by 7).
Achievable at n=7: H=35, 49, 77, 91, 105, 133, 147, 175, 189 (~10.8% of tournaments).
Note 189 = max H at n=7 = 27×7 (Paley). This parallels H=7 impossibility but is NOT a permanent gap.

---

## X. OEIS CONNECTIONS

| Sequence | What | Values | Status |
|----------|------|--------|--------|
| **A038375** | Max H(T) for n-vertex tournaments | 1,1,3,5,45,189,661,... | Paley achieves max at prime n |
| **A000213** | Tribonacci numbers | Related to transitive tournaments | Connection unclear |
| **H(T_p)/|Aut|** | 1, 9, 1729 | NOT in OEIS | Candidate new sequence |
| **H(T_p)** | 3, 189, 95095 | NOT in OEIS | Candidate new sequence |
| **Tangent numbers** | P_n(0,0) = 2^{(n−1)/2} T_n | Connected via EGF | INV-093 |

---

## X. MISTAKES LOG (DO NOT REPEAT)

| ID | What Went Wrong | Correct Statement |
|----|----------------|-------------------|
| MISTAKE-001 | ind_poly_at_2_restricted() bug | Never use old scripts |
| MISTAKE-003 | Per-path identity extends to all n | FAILS at n≥6 |
| MISTAKE-008 | Even-odd split ≡ OCF | NOT equivalent |
| MISTAKE-010 | Hereditary maximizer for all maximizers | Only REGULAR at odd n |
| MISTAKE-011 | Paley at p ≡ 1 mod 4 | Only p ≡ 3 mod 4 |
| MISTAKE-013 | VT ⟹ SC | FALSE at n=21 (Frobenius) |
| MISTAKE-014 | Scalar M for all VT | Only SC VT at odd n |

---

## XI. WALSH/FOURIER ANALYSIS (opus S35c series — MAJOR)

### Spectral OCF (THM-081)
The Walsh spectrum of directed k-cycle count hat{t_k}[S] = (1/2^k) Σ_{C⊇S} (-1)^{asc(S,C)} where asc(S,C) counts edges in S traversed small→large by C. This transforms OCF into a SPECTRAL IDENTITY:

hat{H}[S] = 2·hat{t3}[S] + 2·hat{t5}[S] + ... + 4·hat{α₂}[S] + 8·hat{α₃}[S] + ...

**Counting identity (new proof path):** Equating THM-077 and THM-081:
Σ_k (1/2^{k-1}) Σ_{C⊇S} (-1)^{asc(S,C)} = (-1)^{asc(S)} · 2^r · (n-d)!/2^{n-1}

This is a PURELY COMBINATORIAL identity. Proving it algebraically gives a new, independent proof of OCF.

**Key uniformity (n=5):** For every P2 monomial, signed 3-cycle count N₃ = -2 and signed 5-cycle count N₅ = -4. This uniformity suggests structural explanation.

**Product independence:** Vertex-disjoint cycle k-tuples FACTOR in Walsh domain: E[I_{C₁}·I_{C₂}·χ_S] = E[I_{C₁}·χ_{S₁}]·E[I_{C₂}·χ_{S₂}]. Key for α₂ terms (bc33).

### Transfer Matrix Walsh Spectrum (THM-080)
hat{M[a,b]}[S] = (-1)^{asc_root} · 2^s · (n-2-d)!/2^{n-2}
- s = number of unrooted even-length components (bug fix: was missing 2^s factor)
- Parity: |S| ≡ n (mod 2) — odd support at odd n, even support at even n
- Row sum cancellation: Σ_b hat{M[a,b]}[S] = 0 for monomials not touching a
- det(M) ≠ 0 CONJECTURED (verified exhaustive n=5)
- M[a,b](T^op) = (-1)^n · M[a,b](T)

### Cycle-Rich Tournament Analysis
| n | Min H (cycle-rich) | Growth | Method |
|---|-------------------|--------|--------|
| 8 | 25 | — | Exhaustive |
| 9 | 45 | +20 | 153k sample |
| 10 | 75 | +30 | 192M sample |

**Definition (novel — not in literature):** T is *cycle-rich* if every vertex is in some directed 3-cycle (equivalently: no vertex with in-degree 0 or n-1 to all others). Key lemma: vertex in no 3-cycle ⟹ in no cycle of any odd length (layered structure).

### H-Gap Conjecture (opus S43)
**CONJECTURE:** H=7 and H=21 are the ONLY permanent gaps in the odd-integer H-spectrum.

Evidence:
- All n=7 gaps (63, 107, 119, 149, 161-169, 173) fill at n=8
- n=8 partial exhaustive (18.6%): only H=7 and H=21 missing in [1, 300]
- n=9 sampling (2M): only H=7 and H=21 missing in [1, 200]
- For w ≥ 13, the number of graph-feasible alpha-decompositions (20+) makes blocking all of them seem impossible
- Grinberg-Stanley Theorem 7.1: H(T) ≡ 1 + 2·(# odd cycles) mod 4

---

## XII. CYCLE STRUCTURE AND H-DETERMINATION

### Key Computational Finding (kind-pasteur-S34)
| n | H determined by (c₃, c₅, ...)? | Explanation |
|---|--------------------------------|-------------|
| 5 | YES — H = 1 + 2(c₃ + c₅) | Omega(T) complete at n=5; all cycles share a vertex; α₁ = c₃ + c₅, α₂ = 0 |
| 6 | NO — same (c₃, c₅) can give different H | α₂ depends on cycle PLACEMENT, not just counts; e.g. (c₃=2, c₅=0) gives H=5 or H=9 |
| 7+ | NO — increasingly non-determined | Higher α_k encode cycle independence structure beyond count vector |

**Fundamental insight:** H(T) = I(Omega(T), 2) encodes the **independence structure** of cycle placement, not just cycle counts. Two tournaments with identical (c₃, c₅, c₇) can have different H because the vertex-disjoint cycle pair count α₂ differs.

### Cycle Zeta Function (Speculative — Novel Model)
Define Z_T(x) = Π_{[C] prime cycle} (1 - x^|C|)^{-1} (Ihara-type). For tournaments, all cycles are odd. The connection:
- log Z_T(x) = Σ_k (c_k/k) x^k (formal power series)
- I(Omega, x) is NOT simply related to Z_T(x) because independence structure is finer than count structure
- The "cycle Euler product" Π_C (1 + λ^{|C|}) would equal I(Omega, λ) IF cycles were pairwise independent — but they're not!
- The gap between Π_C (1 + 2) and I(Omega, 2) measures "cycle interaction"

---

## XIII. DELETION-CONTRACTION THEOREM (kind-pasteur-S34 — NEW)

### Theorem: H(T) = H(T\e) + H(T/e)
For any tournament T and any directed edge e = (u→v):
- **Deletion** T\e: remove edge e (digraph on same n vertices, one arc missing)
- **Contraction** T/e: merge u,v into vertex w; w inherits IN-edges from u (tail), OUT-edges from v (head)
- **H(T) = H(T\e) + H(T/e)** holds universally

**Verified:** 100% at n=4 (384 edge tests) and n=5 (10240 edge tests).

**Proof sketch:** A Hamiltonian path in T either uses edge e or doesn't.
- Paths NOT using e: exactly the Hamiltonian paths of T\e
- Paths using e (say ...→u→v→...): contracting e maps these bijectively to Ham paths of T/e
  - The contraction preserves: in-edges to u become in-edges to w, out-edges from v become out-edges from w
  - Every Ham path through e becomes a Ham path in the contracted digraph, and vice versa

**Connection to Mitrovic:** This is the commutative specialization of W_X = W_{X\e} - W_{X/e}↑ from arXiv:2504.20968. The sign difference is because W counts SIGNED paths (via noncommuting variables).

**Implications:**
1. H(T) satisfies deletion-contraction like the chromatic polynomial — it's a "Tutte-type" invariant
2. Induction on number of edges: every tournament reduces to smaller digraphs
3. Could yield a new proof of Rédei's theorem (H always odd) by tracking parity through D-C
4. Could yield a new proof of OCF by tracking independence polynomial through D-C

### p-adic Structure Beyond 2 (Speculative)
| Prime p | H(T) mod p | Known structure |
|---------|------------|----------------|
| 2 | Always odd (Rédei) | OCF: I(Ω,2) ≡ 1 mod 2 |
| 3 | H mod 3 depends on c₃ | H ≡ 1 + 2c₃ mod 3 at n=5; need investigation |
| 5 | Unknown | Check if v₅(H-1) has universal bounds |

### Sandpile Group / Chip-Firing (Speculative)
The sandpile group K(T) of a tournament T is the cokernel of the Laplacian L = D_out - A. |K(T)| = # spanning arborescences (Matrix-Tree theorem). The structure of K(T) as an abelian group could encode more than just its order. Connection: chip-firing dynamics on tournaments relate to Hamiltonian path enumeration through the recurrent configurations.

### Ihara Zeta Function of Tournaments (Novel Connection)
Define ζ_T(u) = Π_{[C] prime directed cycle} (1 - u^{|C|})^{-1}. By the Ihara-Bass formula:
ζ_T(u)^{-1} = (1 - u²)^{r-1} det(I - Au + (D-I)u²)
where A = adjacency matrix, D = degree matrix, r = |E| - |V| + 1.

For tournaments: |E| = n(n-1)/2, D = diag(out-degrees), A = tournament adjacency.
The Ihara zeta function DIRECTLY encodes the cycle counts c_k that feed into OCF via Ω(T).

**Key question:** Does zeta_T(u) evaluated at specific u relate to I(Omega,x)?

**COMPUTATIONAL RESULTS (kind-pasteur-S34):**
- z_inv(1/2) = det(I - A/2 + (D-I)/4) is STRONGLY correlated with H (r = -0.95 at n=5)
- BUT z_inv is NOT uniquely determined by H (multiple values for same H when multiple score sequences exist)
- At n=4: z_inv IS determined by H for transitive and regular, NOT for H=3
- Adjacency real-part spectrum is mostly unique per H at n=5 (except H=9, H=15 with 2 spectra)
- Conclusion: Ihara zeta CONSTRAINS H but does not DETERMINE it — consistent with cycle counts constraining but not determining independence structure

### Stanley-Stembridge Resolution (2024) — Implications
Hikita proved the Stanley-Stembridge conjecture (2024): chromatic symmetric functions of (3+1)-free posets are e-positive. Via Mitrovic-Stojadinovic, the Rédei-Berge function U_T connects to chromatic symmetric functions. **If** the poset structure of tournament arc orderings is (3+1)-free, then U_T would inherit e-positivity, constraining the symmetric function decomposition and potentially the H(T) distribution.

**Status:** Connection needs investigation. Tournament posets are NOT necessarily (3+1)-free, so this may not directly apply. But the machinery (Hessenberg varieties, LLT polynomials) from Hikita's proof might still inform tournament theory.

### Categorification of OCF (Speculative — Novel Model)
**Idea:** Construct a bigraded vector space V_{i,j}(T) such that:
- dim(V_{*,0}) = 1 (empty independent set)
- dim(V_{*,k}) = α_k (independent sets of size k)
- The "Euler characteristic" Σ_k 2^k α_k = H(T)
- The differential d: V_{i,k} → V_{i+1,k-1} has homology encoding more refined invariants

This would be a categorification of the identity H = Σ 2^k α_k, lifting the numerical equality to a chain complex. The homology groups would be new tournament invariants.

**Connection to Khovanov:** Khovanov homology categorifies the Jones polynomial. OCF gives H(T) as a "Jones-like" evaluation of I(Ω,x) at x=2. A Khovanov-type construction for I(Ω,x) would categorify H(T).

---

## XIV. SOFTWARE AND DATA

### Key Scripts
| Script | Purpose | Status |
|--------|---------|--------|
| `symbolic_proof.py` | Exhaustive OCF verification | n≤7 complete |
| `fourier_homogeneity.py` | Fourier decomposition of OCF | Proved for n=5,7 |
| `symmetry_check.py` | Transfer matrix M[a,b] = M[b,a] | Verified n=4,...,8 |
| `h21_dichotomy_proof.py` | Dichotomy verification at n=9 | 106,424 tests, 0 failures |
| `h21_poisoning_graph.py` | Poisoning graph structure analysis | 51,280 mm=2 cases |
| `interlacing_verify.py` | Clique-deletion interlacing | 0 failures n=5,...,8 |
| `paley_deletion_test.py` | Paley hereditary maximizer | Verified p=3,7,11 |
| `deletion_contraction_test.py` | H(T)=H(T\e)+H(T/e) | 100% at n=4,5 |
| `f_poly_zeros_leeyang.py` | F(T,x) zeros in complex plane | Lee-Yang behavior |
| `f_poly_roots_of_unity.py` | F(T,ω), F(T,i) evaluations | Universal divisibility |

### External Databases
| Database | What | Used For |
|----------|------|---------|
| McKay tournament database | All tournaments to n=10; DRTs to n=23 | Verification |
| OEIS A038375 | Max Hamiltonian paths | Paley maximizer conjecture |
| arXiv | Papers | Literature connections |
| `ihara_zeta_tournament.py` | Ihara-Bass determinant vs H | r=-0.95 correlation |
| `padic_beyond_2_test.py` | H mod 3,5,7 structure | H mod 7=0 impossible at n≤6 |

---

## XV. NOVEL MODEL PROPOSALS (kind-pasteur-S34 brainstorm)

These are speculative high-level connections proposed for future investigation. Sorted by estimated promise.

### Tier 1: Strong theoretical motivation, testable

| # | Model | Idea | Test |
|---|-------|------|------|
| 27 | **Tensor network contraction** | Tournament as tensor network: each vertex = rank-(n-1) tensor, edges = contractions. Contraction gives scalar ∝ H(T). Positive bias (every edge present) makes contraction tractable (Jiang-Chen-Schuch-Hangleiter 2024). | Compute tensor network value for small tournaments, compare to H |
| 28 | **Kazhdan-Lusztig polynomials** | Ham paths are permutations. The Bruhat intervals [e, σ] for σ∈HamPaths(T) have K-L polynomials. Sum of K-L polynomials at q=1 over all Ham paths of T could equal H(T) or relate to I(Ω,x). | Compute K-L polynomials of Bruhat intervals for n=4,5 Ham paths |
| 29 | **F_1 geometry** | Tournament as "variety over F_1" — the scheme of complete oriented graphs. Point counting: #T(F_q) might be a polynomial in q whose evaluation at q=1 gives H(T). Tits building of complete graph = apartment structure. | Check if any q-analog of H(T) specializes correctly |
| 30 | **Galois group of I(Ω,x)** | The splitting field of I(Ω(T), x) over Q. For real-rooted (n≤8), all roots real, Galois group acts trivially on root order. At n=9 (complex roots), Galois group becomes interesting. | Compute Galois groups of I(Ω,x) for specific tournaments |
| 31 | **Dimers on tournament graph** | The dimer partition function of a bipartite graph = permanent. Tournament adjacency matrix A has per(A) = # cycle covers. Irving-Omar: H(D) = Σ det(Ā[S])·per(A[S^c]). This det×per decomposition is a "mixed dimer" model. | Compute per(A) vs H(T) correlation |

### Tier 2: Interesting analogy, needs more development

| # | Model | Idea |
|---|-------|------|
| 32 | **Chern-Simons on tournament** | Tournament as flat connection on K_n; H(T) as holonomy. The Witten-Reshetikhin-Turaev invariant at level k=2 might reproduce OCF. |
| 33 | **Floer homology** | Tournament flow polytope as symplectic manifold; pseudo-holomorphic curves = Ham paths; Floer differential counts paths with signs. |
| 34 | **Langlands correspondence** | Tournament L-function L_T(s) via cycle data. For Paley T_p, L_{T_p}(s) might factor through Dirichlet L(s,χ_p). Functional equation from palindrome F(T,x) = x^{n-1}F(T,1/x). |
| 35 | **Quandle structure** | Tournament defines a quandle-like operation: i◁j = "winner of i,j". Not a proper quandle (fails idempotence) but a "tournament rack". Counting invariants of the rack might encode H. |
| 36 | **Hochschild homology** | Tournament algebra A_T = kQ/I where Q = tournament quiver, I = path relations. HH_*(A_T) captures cycle structure. Bar resolution gives chain complex whose Euler characteristic relates to H. |
| 37 | **Schubert calculus** | Hamiltonian paths as Schubert cells in flag variety Fl(n). Tournament structure selects which cells contribute. H(T) = # cells in a Schubert subvariety defined by T. |
| 38 | **Motivic measure** | Tournament as measure on permutation space. The motivic integral ∫ L^{-ord_T} dμ might compute I(Ω,x) or a motivic lift. |

### Tier 3: Highly speculative, long-shot connections

| # | Model | Idea |
|---|-------|------|
| 39 | **Quantum error correction** | Tournament codes (binary strings of length m encoding tournaments) as error-correcting codes. Weight enumerator relates to Hamming distribution of H. |
| 40 | **Fractional cascades** | Tournament arc-flip graph as cascade lattice. The number of "downward paths" in the cascade = H(T). |
| 41 | **Renormalization group** | Tournament coarsening (contract cycles into super-vertices) as RG flow. H(T) is an RG invariant. Fixed points = regular tournaments. |
| 42 | **Random matrix universality** | Skew-adjacency eigenvalue distribution → GUE for random tournaments. Correlation functions of eigenvalues relate to cycle counts. Phase transition at eigenvalue gap ↔ H-gap. |
