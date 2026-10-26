# Research board, random recovery and connection tests

**Current status:** local Theorems A-G in [README.md](README.md) are proved with
the cited structural prerequisites; receipts give bounded exact checks.
The global prime-tournament and Hamilton-complement classifications remain open.
Base inherited commit: `fc3768ac5f`. The packet makes no claim of priority.

## Portfolio and six live concepts

Anchor: actual child types and quotient stabilizers in prime composition.
Niche: arithmetic period strata, orbit masses and exact quotient operators.
Wildcard: changing size through pointed embeddings and boundary data.

| Live concept | Current mechanism | Decisive observation and effect on the other lanes |
|---|---|---|
| Prime quotient | A canonical indecomposable tournament, with its permutation action | `Q_9` is prime but has group `C_9`: primes must be classified as objects, not as simple groups. This constrains decoration, counts, operators and growth. |
| Actual child type | Colour each quotient slot by a full recursive isomorphism signature | `A+` and `A-` have the same `C3` but are different colours. One change breaks rotations, changes group mass and path-orbit count, but keeps the observer. |
| Observer fibre | Shared labelled cyclic triples and strong modules | The complete 128-element fibre has four classes. Its size alone neither fixes a group nor supplies a move law or growth rule. |
| Stabilizer and divisor | Period `d` gives subgroup of order `m/d`; general groups use a cycle index | The `m=9` middle stratum is real. The Paley-seven control replaces a cyclic group by order 21 and still obeys the same general reduction. |
| Physical move ports | An equivariant move before quotienting, with multiplicity retained | The four-vertex quotient path has row sums 3 only with port weights. Detailed balance and the Walsh spectrum depend on these weights. |
| Pointed growth | A marked gap and an induced embedding | One parent necklace gives two children if its mark is erased. A 44-to-52 embedding loses symmetry; size is not a monotone group invariant. |

After the 81/27 computation, the board changed from an invariant-only atlas
to a decorated-tree classification: every lane needs the child-type stabilizer.
After the full fibre and necklace calculations, the arithmetic and operator
lanes became exact orbit problems rather than analogies. After the growth
hostile, attachment orbits became necessary data for any cross-size claim.
The noncyclic control tests that this reduction survives beyond carousels.

The integrated commit `d761e4cf0f` supplies a stronger theorem on the same
base: [prime-specific onset at 4p](../tournament_composition_classification_20261026/fibre/RESULT.md).
It strengthens the observer/stabilizer lanes without changing our receipts.
Its [Sylow-adapted towers](../tournament_composition_classification_20261026/arithmetic/RESULT.md)
classify a different locus, maximal `p`-valuation, whereas our necklaces
classify a fixed observer fibre. The connecting map is the actual group
order in (2): at `m=p`, changing a constant word to a single marked bit
removes exactly one factor `p`. This preserves the claimed prime loss but
forgets the position and full action; neither locus is identified with the
other. The incoming full path-cover theorem strengthens the count survivor.
The incoming blue-half-edge hostile makes the boundary-data lane sharper.

The view generator used here is explicit:
objects = tournaments / hypergraphs / pointed substitutions;
representations = arcs / strong trees / binary words / Walsh functions;
invariants = `C3`, `H`, `Aut`, period, orbit mass, spectrum;
operations = converse / block cut / substitution / induced insertion;
quotients = labels / rotations / general quotient automorphisms;
scales = 4-vertex blocks, the sharp order-12 boundary, arbitrary odd `m`,
and the `m -> m+2` embedding. Each successful move occupies a concrete cell
of this product; none identifies all these observables or operations.

## Random recovery, with frozen provenance

The user asked for random searches of earlier work. At the inherited commit,
take all **1,915 top-level Markdown files** in `05-knowledge/results`, exclude
names containing `20261025` or `20261026`, sort with Windows `pathlib.Path`
ordering (case-insensitive), and sample 12 without replacement with Python
`random.Random(20261026)`. This is a bounded sample, not a corpus reread.
[random-sample.json](random-sample.json) pins its universe digest, selected
names, normalized source hashes and seed. To reproduce on another platform,
sort these filenames by `str.casefold` before sampling.

| Sampled file (under `05-knowledge/results`) | Disposition |
|---|---|
| `collatz_height_envelope_20261011.md` | Screened opening: source-specific height bounds; no theorem imported. |
| `procgen_tension_20260926_ranks_christoffel_duality.md` | Recovered circulation/potential boundary and primitive-necklace language; used below. |
| `continuing11_20260908_quadratic_family.md` | Screened opening: unit-cyclicity parameter wall; no algebraic-geometric transfer asserted. |
| `planar_jc48_sep07_higher_odd.md` | Screened opening: restricted cusp classification; no transfer asserted. |
| `collatz_mod6_20260917_zsigmondy_triad.md` | Recovered primitive-divisor warning; the elementary length-six hostile below is checked directly. |
| `sc_partition_formula.md` | Historical/superseded note; followed its correction to THM-283 instead of importing its old partition rule. |
| `lrc14_composite_clock_capacity_thm4447.md` | Recovered cyclic-action gcd strata; geometric capacity remains a separate predicate. |
| `gold_prefix_side_poset_20261019.md` | Recovered boundary flags and restriction fibres; used as a design constraint, not a tournament theorem. |
| `metagraph_spectrum_torus_session16_20261020.md` | Recovered exact distinction between `H`, `H/Aut` and a specified return graph. |
| `collatz_directions_20260926.md` | Screened opening: source debt and precision; no convergence claim imported. |
| `amm12592-allR-transient-theorem-boxeph.md` | Screened opening: exact transient skeleton; no transfer asserted. |
| `shared_frontiers_and_reset_two_20261004.md` | Recovered requirement to retain the source and rooted suffix under transport; no Collatz result inferred. |

## Connections with explicit maps and failure boundaries

**1. Burnside parity -> the sharp order-12 threshold.** Source:
[THM-283, corrected SC Burnside formula](../../../01-canon/theorems/THM-283-sc-burnside-formula.md).
Target: the action of a quotient automorphism on child slots. The map sends
a permutation to its cycle decomposition and then to its action on unordered
pairs. Preserved predicate: consistent orientation of each pair orbit.
Lost information: the actual arcs and child-type colouring. Required sidecar:
those arcs/colours when lifting back to an automorphism. The antipodal pair
in an even cycle is the cheapest obstruction; it proves all automorphism
cycles odd, giving `3*4=12` in Theorem B. For anti-automorphisms the corrected
filter is parts `2 mod 4` and at most one fixed point, not arbitrary even
parts. The historical sample is not authority for the latter false rule.

**2. Primitive periods -> divisor inversion, but not primitive primes.**
Source: the sampled [Zsigmondy triad](../collatz_mod6_20260917_zsigmondy_triad.md).
Target: binary word orbits. The actual shared map is restriction from periods
dividing `n` to exact periods, inverted on the divisor poset; in our setting
it acts on `2^n` words and gives (3). Preserved predicate: a specified period
divides another. Lost information: multiplicative orders of primes and
valuations of integer factors. That arithmetic would need its own sidecar;
no word-to-prime map has been constructed. The decisive hostile is length
six: there are nine primitive binary necklaces, but no prime has
`ord_p(2)=6`, since `2^6-1=63=3^2*7` and the two prime orders are 2 and 3.
Even length six is used only for this word/arithmetic control, not as a
carousel tournament parameter. The random hit therefore adds both an
arithmetic connection and a precise limit on transporting it.

**3. Composite clocks -> rotation strata.** Source:
[THM-4447, composite-clock capacity](../../../01-canon/theorems/THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction.md).
Target: the rotation action on subsets of a cyclic index set. The map is
literally `j -> j+t` on `Z/q`: it has `gcd(q,t)` cycles, each of length
`q/gcd(q,t)`. Preserved predicate: invariance under that permutation.
For example, rotation by three on nine positions has three cycles of length
three, exactly the term used in (6). Lost information: the clock's phase,
strict interval endpoints and placement of bad labels. Those are required
to recover a capacity inequality. Orbit arithmetic alone makes no progress
on the remaining LRC(14) entry problem; LRC(14) stays **OPEN**.

**4. Potential/circulation duality -> metagraph mass balance.** Source:
[the sampled tension note, Theorem D](../procgen_tension_20260926_ranks_christoffel_duality.md),
routed there from [THM-4474, strategy-cube cycle densities](../../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md).
Target: `M_m`. The concrete map pushes the uniform directed edge flow on the
cube to necklace orbits. It preserves zero divergence; its edge weight is
`d(u)c(u,v)`, symmetric by (5). It loses the identity of the flipped position;
the multiplicity sidecar restores total flow. A strict descending potential
for *both* directions of every allowed flip is impossible, already on the
two-step walk `000 -> 100 -> 000`. Orienting only increasing-weight edges
would create a different operator. No Collatz weights or convergence
conclusion transfer through this map.

**5. Boundary flags -> pointed extension, a comparison with a tested limit.**
Source: [Gold prefix-side posets](../gold_prefix_side_poset_20261019.md), where
boundary first/last flags control restriction weights. Target: (7).
The common operation is forgetting a marked boundary under isomorphism;
it preserves the unpointed parent and destroys its attachment position.
For the tournament target the exact map is the orbit projection
`{pointed binary words} -> {necklaces}`, and the needed sidecar is a gap
orbit under the word stabilizer. The hostile `011,101 -> 00011,01001` shows
extension does not factor through that projection. This is a recovered
boundary-data principle; no numerical Gold inequality is transported and
no functor from these posets to tournaments is claimed.

**6. Path presentations -> symmetry-sensitive fibres.** Source:
[session16, path spectra and return geometry](../metagraph_spectrum_torus_session16_20261020.md).
Target: the emblem. Map: ordered Hamilton paths to their automorphism orbits.
The action is free, so the quotient has exactly `H/|Aut|` elements. It
preserves the path-orbit object, loses the vertex labels, and needs the
actual group action. The 81/27 pair is the decisive test: equal `C3` and `H`
but presentation masses `5715` and `17145`. This identifies a genuine
obstruction to replacing all metagraph port data by cyclic triples, without
asserting that any particular black/blue adjacency differs.

## Stopping boundary and next precise questions

The positive signals were pursued through a sharp universal threshold,
a full infinite fibre classification, a general prime-quotient reduction,
an exact spectrum, and an explicit growth obstruction. No larger census is
needed to prove these claims. The finite tests are controls, not evidence
that unresolved all-order classification problems have disappeared.

Next: for a specified prime quotient, compute its action on actual child
types and the orbit of an attachment port; for a specified Hamilton-complement
move, transport full endpoint/path-cover responses through that tree. The
question is which response admits a finite compositional signature, not
whether a single scalar can substitute for missing action data. Arbitrary
prime quotients with arbitrary children and the cross-size compatibility of
those responses remain **OPEN** here.
