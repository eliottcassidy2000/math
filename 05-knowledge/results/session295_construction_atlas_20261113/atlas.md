# A construction atlas: what composes, and what does not

Chan–Pak remains open. This inventory classifies proof interfaces, not just
sequences of Hamiltonian path totals. References are paths relative to
`05-knowledge/results`; all quoted bounds retain the cited scope.

| Construction and source | Retained witness | Valid composition | Failure boundary |
|---|---|---|---|
| O/Z/J rooted extension; [S227 boundary escape](../session227_chan_pak_constructions_20261108/boundary_escape/PROOF.md) | Actual marked graph and four-coordinate profile; extension letters constrain terminal faces | Reuse an output root only under its named conditions; a two-vertex gadget can escape all three terminal faces when both roots are nonsinks | A scalar orbit does not certify an inverse or new root; escaping faces is not itself count coverage |
| Rooted gadget grafts; [S227 graft images](../session227_chan_pak_constructions_20261108/graft_images/PROOF.md) | Four state coordinates plus a physically marked module for image recognition | Forward positive transfer and integer inverse on the image plane | Two strong six-vertex graphs share state (37,7,43,9), only one is the graft image; even the exact cone/lattice contain a nonphysical point |
| Multi-root grafts; [S227 multiroot](../session227_chan_pak_constructions_20261108/multiroot/PROOF.md) | Marked root cut and r+1 gadget coordinates | Positive path-cover expansion; adaptive grafts not covered by fixed-menu bounds | Fixed k-type menus give only O((log X)^k) distinct values; complete rooted state can identify binomially many nonisomorphic graphs |
| Scalar gap switches; [S227 propagation](../session227_chan_pak_constructions_20261108/adaptive_arithmetic/PROOF.md) | Pair with equal E,U,d and distinct D | A common continuation multiplies a pure D gap by a positive odd Hamilton count | A graft exposes the gap in U,d; arbitrarily mixing the two pairs does not give independent selectable digits in one final fibre |
| Forced shuffle; [S214](../session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md) | Prescribed pair of Hamilton words on the two parts of a marked cut | A canonical injection preserves both subsequences; nonneutral induced growth is at least 3/2 of old H | It is an inequality, not a total-count recursion; positive induced increment 2k forces old H ≤4k−1 |
| Seven-letter central bridge; [central modular session](../central_modular_chanpak_20261110/README.md) | (H−E,E,C) and exact word | Universal exact content multipliers form 5^b 9^k 21^c; normalize by a unimodular matrix to keep all modular phase | Neither 3 nor 7 is an exact universal content multiplier in this grammar; modular phase retention does not give a physical inverse |
| Five-core compiler; [web](../five_core_web_20261111/README.md) | Nonnegative seven-coordinate transfer; actual core mark and tail | Finite profiles and shell transport; a golden H13 unit-switch family with deletion floor13 | Two H9 cores share unpointed full path-cover profile but differ in marked response; interval hull and modular fullness do not imply occupancy |
| One-port endpoint inflation; [coverage compiler](../endpoint_coverage_compiler_20261111/README.md) | Marked selected set, compatible response and actual transitive inflation | H=delta+2^t U; infinite mixed 2/3 families and 99·2^k+23 with literal witnesses | The gcd-one extension lattice does not supply the requisite endpoint responses on all rays |
| Uniform 141 return and fixed defect strata; [return metagraph](../uniform_return_metagraph_20261111/README.md) | Named persistent vertices, deletion/replacement and pointed transitive module | Literal 425 can be reset to a 1-ready 425 and then reach 1277; a ready chart renews while changing absolute target | Equal-H old sources may not admit induced425; any fixed bounded branching-defect bank supplies only O(log X) distinct values |
| Physical receiver banks; [S281](../session281_chanpak_20261112/receiver/PROOF.md) | Renewable common interface across multiple absolute values | Six-wide bank renewable, twelve-wide intervals occur infinitely often | One wide observed interval is not an induction corridor; residue completeness can remain power-sparse |
| S292 boundary and sparse guards; [S292](../session292_triple_bridges_20261113/README.md) | Labeled cyclic triples + specified quadratic port constraints | Partial guard fibers factor into prime reversals and linear extensions; design port observations with known surviving choices | Interior H survives all such replacements, but attached H need not; full guards leave no hidden compensating freedom |
| S293 cycle opening, feedback and module ideal; [S293](../session293_cycle_returns_20261113/README.md) | First strong SCC, marked entrance and endpoint column, renewed seam | Neutral induced return indefinitely with moving entrance; exact feedback and uniform-context cover response | Fixed entrance cannot repeat; active uniform feedback grows exponentially; uniform ±2 repairs blocked by any nonprimitive context ideal |
| S294 asymmetric content obstruction; [S294](../session294_context_27_20261113/README.md) | Full odd-cycle support polynomial under uniform context | Ordinal composition multiplies packing polynomials; exact valuation laws on mixed exceptional cores and triangles | Mod27 collapses four distinct integral profiles; asymmetry does not imply primitive context ideal |

## Comparing proof types

1. **Forward composition is not physical inversion.** A matrix inverse, or a
   CRT candidate in the intersection of arithmetic cones, does not specify
   the named induced module. An inverse proof can avoid reconstruction only
   by constructing its predecessor independently and applying a certified
   forward move whose actual guards the predecessor has.
2. **Observation equality depends on the next operation.** Unpointed cover
   equality is complete for uniform exteriors, not for a pointed graft or
   selective attachment. A rooted state is likewise not a complete
   fingerprint of inverse graph images. A useful atlas should record a
   continuation class alongside each quotient.
3. **A gap is not an anchored selector.** Reusable switches require absolute
   baseline, sign, compatibility, and a returned interface, not just a gcd
   or difference. The design target is a coverage loop on physical states
   of the form (absolute target, endpoint guards, move capabilities), not a
   large transition graph on scalar values.

## New selective-port filtration

Let C be a named nonempty subtournament. An exterior vertex is *exceptional*
when it has both a win and a loss against C. Assume at most q exterior
vertices are exceptional. For each k, let N_k count full Hamiltonian paths
whose intersection with C consists of k maximal nonempty runs. Then

    max(k−1−q,0)! · max(k−q,0)!  divides N_k.

**Proof.** Fix a path. Its exterior pieces are a prefix, k−1 nonempty
internal excursions between C-runs, and a suffix. Mark an exterior piece
if it contains an exceptional vertex. Let a count marked internal
excursions and b count marked exterior end pieces (prefix/suffix).
Since the exterior pieces are disjoint, a+b≤q.

Join consecutive C-run words linked by a marked excursion, keeping each
linked segment intact as one compound word. There are k−a compound
words, separated by k−1−a unmarked excursions. A marked prefix or suffix
anchors its adjacent compound word; let c≤b count distinct anchored
compounds. Permute the unmarked excursions and permute the unanchored
compound words independently. Both types of word have their internal
arcs unchanged. At any altered junction between two of them, the
excursion endpoint is uniform against C, so its original directed
incidence still works. At any altered junction with an exterior end
piece the same argument applies, since marked end pieces anchor their
adjacent compound. Empty end pieces require no arcs.

This partitions paths into freely permuted sectors: the locations of
anchored compounds, the identities of all words, and the fact that
each compound contains only marked internal excursions are invariant
under both permutations. Both actions commute and are free, since
word vertex sets are pairwise disjoint and nonempty. Thus each orbit
has size (k−1−a)!(k−a−c)!, divisible by the claimed product. Summing
proves the result.

This gives a prime-power cutoff even after nonuniform attachments. For
prime p and precision p^s, all k with
v_p(max(k−1−q,0)!)+v_p(max(k−q,0)!)≥s are automatic. For one exceptional
port at p=3 and s=3 the cutoff is k≥7: 5! and 6! supply 3 powers of
three, while k=6 only guarantees 2. This bounds the
potentially relevant number of runs, not the number of necessary endpoint
profiles, and gives no assertion that the truncated sectors are freely
selectable.

For q=0 this recovers the universal product k!(k−1)! forced by S293.
For q=1 it genuinely differs: take S293's asymmetric core with rows
(56,117,129,22,100,76,141,59), and add a vertex beating only vertex
5 (zero-based) and losing to the other seven. A Held–Karp subset recurrence
gives H=1147≡1 mod3, whereas every uniform context containing the core
has count divisible by27. Thus one selective port can kill the entire
obstruction, even as the higher-run sectors remain subject to a universal
factorial cutoff. The next constructive question is whether one or two
ports can be changed with equal control of the few low-run sectors *and*
with a renewable output interface.

## Cross-field lead: two ways a sector can disappear

The [square-sum transfer](../session53_g27_collatz_cooper_20261024/README.md)
for adding 27 uses both endpoints 9 and 22 on distinct old paths. Only
32 of 1116 spanning two-path covers permit the actual join. This is not
a direct equality of tournament and square-sum cover invariants: the
square-sum condition can make a sector empty because endpoints are
incompatible. The selective-port theorem shows a complementary operation:
even when the sector is nonempty, unconstrained excursions give a free
group action and make the contribution vanish modulo a prime power.
The correct shared state distinguishes **which sectors exist** from
**which sectors are invisible to a given modulus**. Neither the scalar
minimum cover nor a context gcd supplies that distinction.

**Sharp one-port cutoff at modulus 27.** Take a six-vertex transitive C,
eight ordinary exterior vertices in the order x1,y1,...,x4,y4, with
C→xi and yi→C, and a ninth exterior vertex z beating exactly one
vertex of C. Its arcs to the eight exterior vertices are immaterial.
In the six-run sector all C-runs are singletons, and there must be
five internal excursions. Only z can be a singleton excursion, so
the nine exterior vertices force z alone in one excursion and the
other four to be pairs, with empty prefix/suffix. Each pair is xi,yj
with i≤j, so the unique matching is i=j. In the compound C→z→C,
the second C vertex is forced and the first has five choices.
The five resulting compounds can be permuted in 5! ways and the
four ordinary excursions in 4! ways. Thus N_6=5·5!·4! has valuation
exactly two at 3, proving that seven is the best universal cutoff for
one exceptional port at modulus 27.
