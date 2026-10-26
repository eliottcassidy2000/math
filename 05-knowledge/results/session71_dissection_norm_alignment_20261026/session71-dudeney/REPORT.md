# Session71: Dudeney proof dependency audit and exact length-transfer sidecar

## Status and immediately useful findings

- **CITED, not independently re-proved in full:** Demaine–Kamata–Uehara, *Dudeney's Dissection is Optimal*, requested [arXiv:2412.03865v3 HTML](https://arxiv.org/html/2412.03865v3), 31 July 2025. Theorem 1 in HTML / 1.1 in the arXiv PDF excludes at most three connected polygonal pieces under translations and rotations, with no piece flipping. The current published primary reference is [ITCS 2026, LIPIcs 362, Article 47, pp.47:1–47:22](https://doi.org/10.4230/LIPIcs.ITCS.2026.47), published 23 January 2026; its Theorem 1 has the same scope. The existing four-piece construction supplies the matching upper bound. This is not an optimality claim for arbitrary curved pieces or allowed reflections.
- The proof uses coupled edge and vertex matching diagrams, angle balance **inside each connected component**, orientation-sensitive endpoint transport, and metric restrictions. Fermat's integer two-square theorem is **not a cited or actual explicit dependency** of the inspected proof. The text does not invoke general rational independence either; it uses concrete lengths and inequalities.
- **FINITE-EXACT transcription check:** the displayed cut-graph lists have 44 square types and five triangle types. Their Cartesian product has 220 pairs, cc-diff leaves 64, the further angular filter leaves 46, partitioned as Case A:14, B:6, C:26. These are counts of types/listed pairs, not counts of pieces, congruence realizations, or a mathematical isomorphism with another occurrence of 44 or 5.
- **PROVED here as an elementary deduction, no priority claim:** a well-behaved length-splitting chain has an exact partial-affine carrier `(sign, offset, admissible interval)`. The interval cannot be dropped. This is an operation-relative quotient for the stated length-chain experiment, not a complete dissection carrier.
- **PROVED here from the paper's angle balance:** for a matching component between regular n- and m-gons, its corner counts a,b satisfy `2a/n - 2b/m ∈ Z`. When gcd(n,m)=1 this is equivalent to `n/gcd(n,2) | a` and `m/gcd(m,2) | b`. Thus every odd regular n-gon versus a square has all-or-none n-gon corners and an even number of square corners in each component. This does not extend the paper's piece-count theorem to pentagons.

## 1. Provenance and verification boundary

Authors: Erik D. Demaine, Tonan Kamata, Ryuhei Uehara. Source license: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), linked by the arXiv version record and HTML. This report is an independently written explanatory adaptation with additional deductions and audit notes; the authors do not endorse it. Mathematical results are attributed below. No figures or long passages are reproduced in this report.

A material version warning must remain visible. [The v4 record](https://arxiv.org/abs/2412.03865v4), dated 1 August 2025, says that v3 compiled an outdated file and v4 supplies the intended corrected/refined analysis. I did not silently substitute v4 for the user's v3. On this retrieval (26 October 2026), the directly downloaded v3 and v4 PDFs had different SHA256 values, but normalized extracted text differed only in the version/date stamp; extracted text after page one was byte-identical. The HTML text similarly differed only in the stamp and website CSS. Page 18 rendered identically. This observation does **not** explain the authors' correction notice or establish all visual/source-file identity. Do not claim a substantive v3→v4 delta was verified from those served artifacts.

A narrow cross-check against the [published ITCS 2026 text](https://drops.dagstuhl.de/storage/00lipics/lipics-vol362-itcs2026/html/LIPIcs.ITCS.2026.47/LIPIcs.ITCS.2026.47.html) confirms the same theorem scope and the same prose cut-type/pair inventories (44/5,64,46; A14/B6/C26). Its numbering is sequential in both HTML and PDF. The publication is 22 pages rather than 26: Case B-1 refers to Section 4 and Case B-2's analogous proof is explicitly omitted, so it is not an expanded independent case certificate. Rendered published page 47:17, Figure 23, retains the label discrepancies seen in arXiv Figure 24 (for example square B3 versus prose A3; triangle E1 versus prose F1; square F2/F4 versus prose E2/E4). The degree-three tree mistakenly described as ED in the statement is still correctly treated as VD in its proof, and no definition of “monochromatic” was located. Cite the publication as current primary literature, but retain the requested version and these limited audit observations. Do not infer venue from stale SoCG 2025 template metadata in the preprint.

Evidence retained locally: `v3.html`, `v4.html`, `v3.pdf`, `v4.pdf`, extracted text, `download-receipt.json`, `version-observation.json`, and selected rendered pages. The final repository packet need not duplicate the full source PDFs/HTML; URL/version/hash pins plus this explanation suffice.

Read scope: all v3 HTML proof text, source PDF text, and rendered proof figures on pages 5–6,8–25. These include the cut-type lists and final case diagrams. This is a dependency/typing/consistency audit. It is **not** an independent generator of all planar cut topologies or a rerun/reproof of every geometric case. Some displayed labels have apparent presentation discrepancies: Figure 24's labels do not consistently agree with its case-list prose; the Case B statement calls its degree-three tree a tree of ED although ED has maximum degree two, whereas the proof correctly works in VD. The adjective “monochromatic” is used without an explicit definition located in the inspected text. These are reasons to avoid presenting our reading as a full formal certification, not claims that the optimality theorem is false.

Repository inheritance: clean base `c2039c32f9726b09ed0e09d5f5c126fc72975897`; AGENTS, START-HERE, relevant CURRENT-FRONTIER, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, CORE-PAPERS and Session70 synthesis read. The read-only processor command could not run because this checkout has no `.machine-id`; none was fabricated. Parent owns repository coordination/publication.

## 2. What is actually finite, and what matching diagrams preserve

A small piece count does not bound the number of cut segments. Consequently, enumerating a bounded collection of coordinate polygons would not establish the theorem. The paper first suppresses arbitrary strings of paired internal vertices and retains the cut graph's incidence structure and vertex types. Its equivalence is type-preserving subdivision equivalence, **not** an assertion that all metric or assembly information has disappeared harmlessly.

The geometric reduction (HTML Lemmas14–16; PDF5.1–5.3) uses:

1. Normalize triangle side `τ=2`, square side `σ=3^(1/4)`, both areas `sqrt(3)`.
2. The square diameter is `sqrt(2 sqrt(3)) < 2`. No single piece contains two triangle corners. Thus fewer than three pieces are impossible; with three, each piece contains one triangle corner.
3. A further geometric argument excludes a piece carrying diagonal square corners, and a reduction allows simple pieces.
4. The abstract disk cuts give six topological families. Restoring boundary attachments and vertex types yields the displayed five triangle and 44 square classes.
5. Convex/reflex angle differences are invariant under inserting paired vertices, so they filter class pairs before metric matching.

The matching diagrams are not ordinary graphs of neighboring pieces. They record **which intrinsic side/corner occurrence of which same piece occurs in both target assemblies**:

- ED has nodes corresponding to appropriate cut edges/virtual whole-side nodes, and one link per intrinsic piece side. Flat vertices require explicit virtual-node treatment; a subdivided target cut edge is not always an entire piece side. Every ED node has degree one or two, so each component is a path or cycle. Along an ED path, corresponding piece-side lengths agree.
- VD has target cut-vertices as nodes and one link per intrinsic piece corner. Vertex weights are the sums of the actual corner angles incident there. Every component has equal total weight on its two assembly shores.
- ED and VD share the intrinsic piece's side–corner incidence relation. Without it, separate isomorphism types of ED and VD are not enough to perform the proof's endpoint continuation.
- With no flips, crossing a common-side contact exchanges clockwise and counterclockwise endpoint roles. Composing these contacts retains a two-state endpoint transport parity. The paper uses this to force or rule out cycles in VD. Forgetting orientation would erase this test.
- Along a well-behaved VD path, equality runs in ED and divide-in-two events give alternating length identities. A flat-node/three-way boundary exception cannot simply be counted as an ordinary two-way split.

## 3. Dependency map, with HTML/PDF numbering

| Obligation | Source mechanism | Dependency/role |
|---|---|---|
| ≤2 pieces impossible | HTML14 / PDF5.1 | Diameter comparison; not integer arithmetic |
| Topological candidate atlas | HTML14–17 / PDF5.1–5.4, Figs17–20 | Geometric exclusions, simple-piece reduction, typed subdivisions, cc/tri angle differences |
| Equal-length propagation | HTML Definition2 + Observation3 / PDF3.1–3.2 | ED fiber degrees1/2 and congruence of the same piece side |
| Angle complement/equality propagation | HTML Observation5 / PDF3.4 | Paired-vertex total angle2π |
| Component conservation | HTML Observation6 / PDF3.5 | Same piece-corner links counted on both shores |
| Joint side/corner continuation | HTML7–9 / PDF3.6–3.8 | Intrinsic adjacency and well-behaved splits |
| No-flip endpoint tests | HTML10–12 / PDF3.9–3.11 | Clockwise/counterclockwise role transport; incidence and cycle structure |
| Triangle/square corner restriction | HTML18 / PDF5.5 | Component balance moduloπ; triangle angleπ/3 versus square angleπ/2 |
| U-shaped boundary restriction | HTML20 / PDF5.7 | Corner path structure plus perpendicularity/intersection geometry |
| Exhaustion into A/B/C | HTML21 / PDF5.8 | Listed atlas and structural case analysis; not independently generated here |
| Case A exclusion | HTML22 / PDF5.9 | More required U-boundary path terminals than exceptional flat/middle terminals |
| Case B exclusion | HTML23 / PDF5.10 | Degree-three VD tree; alternating identities imply impossible equalities involvingσ andτ; some subcases force cycles in a tree |
| Case C exclusion | HTML24 / PDF5.11 | Well-behaved square-corner path; equality/alternating length contradictions and orientation/cycle restrictions |
| Four-piece optimality | Theorem 1 /1.1 plus existing construction | All above; no reflection or curved-piece extension supplied |

This table is a dependency explanation, not a replacement proof of the figure-dependent last rows.

## 4. The arithmetic input: real metric, not Fermat two squares

The area relation is `σ²=(sqrt(3)/4)τ²`; hence at `τ=2`, `σ^4=3`. The diameter test is `2sqrt(3)<4`. Later contradictions include `σ=τ` or an alleged sum that would force `τ` to be one of incompatible multiples ofσ. Strict-shorter boundary pieces also matter. These are exact real Euclidean statements. The paper neither assumes rational/integer coordinates nor factors Gaussian integers.

Fermat's two-square theorem concerns representations by the **integer** norm `a²+b²` (and its general extension characterizes prime exponents). The ordinary squared-distance formula uses the same polynomial over **real** coordinates, but that is not a theorem-preserving identification of the problem. Arbitrary real rotation/dissection vertices are permitted. Restricting to a square lattice would change the quantified universe.

A useful but weaker arithmetic sidecar is that no nondegenerate equilateral triangle has all vertices in `Q²`: its squared side length would be positive rational and its area both rational (determinant) and that rational times `sqrt(3)/4`. This is an elementary rational-coordinate obstruction, not a lower bound on dissection pieces. The actual four-piece dissection is an immediate hostile to any argument asserting that the arithmetic field mismatch prevents finite dissection altogether.

The 44/5 pair counts are likewise a combinatorial bookkeeping fact. A relation to representations of the integers44 and5 needs its own map and preserved predicate; the paper's atlas labels supply no such map.

## 5. Derived lemma: angle-charge denominator obstruction

**Claim (all finite polygonal dissections, any piece count, flips allowed).** Suppose target shapes are regular n- and m-gons. Let a VD component contain a original n-gon corners and b original m-gon corners. Then

`2a/n - 2b/m ∈ Z`.

**Proof.** The n-gon's corner weight divided byπ is `1-2/n`; the m-gon's is `1-2/m`. Every noncorner VD node has weight an integer multiple ofπ: eitherπ or2π in the paper's convention, which omits the straight angle occupied by a side at a flat node. Divide component weight balance byπ and reduce moduloZ. Integer terms disappear and the displayed relation remains. Reflections do not alter the interior angles or the two-shore count. QED.

If gcd(n,m)=1, multiply bynm and reduce modulo each denominator: `n|2a` and `m|2b`. Conversely those two divisibilities make each of `2a/n` and `2b/m` integral, so the condition is equivalent to

`a ≡ 0 mod n/gcd(n,2)`, `b ≡ 0 mod m/gcd(m,2)`.

For odd n and m=4, `0≤a≤n` forces `a∈{0,n}`, and `4|2b` is exactly `b even`. In particular the **same angle-charge constraint** occurs for an equilateral triangle and for a regular pentagon against a square. It does not make their cut atlases or metric obstructions equivalent. A regular pentagon, a Penrose kite and a Penrose dart are also different target shapes; apply their actual angle weights rather than transferring a regular-pentagon formula to every fivefold object.

More generally the additive target of this argument is the subgroup of `R/Z` generated by normalized target-corner angles. For these two regular polygons it is cyclic of order `lcm(n/gcd(n,2), m/gcd(m,2))`: each corner charge has the indicated reduced denominator, and two cyclic subgroups of `Q/Z` generate the subgroup with their least-common-multiple denominator. Triangle/square therefore has charge group `Z/6`, whereas pentagon/square has `Z/10`. The shared all-or-none/even-corner constraint is not an isomorphism of these charge groups. This is a componentwise charge test. It is only necessary: satisfying the congruence does not provide lengths, embeddings, nonoverlap, or rigid assembly transformations.

## 6. Derived lemma: exact partial-affine carrier of split chains

### Experiment and theorem

Let a chain consist of equality transport and positive split operations. An equality preserves a positive side-length variable. A split at a total side length `L>0` sends `x` to `R_L(x)=L-x` and is legal exactly when `0<x<L`. Start with a positive input x; every intermediate length must remain positive. We observe the final real length or a distinguished failure symbol, and allow future chains composed with the current chain. This is an abstract length-constraint experiment; not every such chain is asserted realizable by a dissection.

Every chain induces either the empty partial function or

`F(x)=εx+c,   x∈I`,

where `ε∈{+1,-1}` and I is a nonempty open interval within `(0,∞)` (possibly unbounded for the empty word). A split updates `(ε,c)` to `(-ε,L-c)` and intersects I with `{x:0<εx+c<L}`. This set is an interval; intersecting all prefix constraints proves the statement by induction.

For `G(y)=δy+d` on J, composition is exactly

`G∘F=(δε, δc+d, I∩F⁻¹(J))`,

with the empty result identified as failure. This gives a closed carrier and an associative composition law because it is actual partial-function composition. It is the exact minimal quotient for the stated named-input observation: different domains are separated by an input in their symmetric difference; with the same nonempty open domain, two distinct affine maps cannot agree everywhere. Thus only equal triples, or two empty partial maps, are observationally identical. The same quotient is exact even for Boolean legality observations when **all future positive real split lengths** are allowed: if two legal outputs satisfy u<v, append a split with total L strictly between u and v. It accepts the smaller output and rejects the larger. A domain difference is already detected by present legality. This stronger conclusion depends on the unrestricted future threshold alphabet; it is not claimed for a fixed finite menu of split lengths.

For split lengths `L1,...,Lr`,

`x_r=(-1)^r x_0 + Σ_{j=1}^r (-1)^(r-j)Lj`,

provided all intermediate inequalities hold. This recovers the alternating-length identity used in the paper. A closure equation `x_r=x_0` requires offset0 when r is even, or `x_0=offset/2` when r is odd, and still requires `x_0∈I`.

### Conditional closed-chain coefficient balance

If the permitted split totals `λ1,...,λd` are linearly independent over Q, an even closed chain can satisfy its formal equation `x_r=x_0` only when each total type occurs equally often in odd and even positions. Indeed the offset is an integer linear combination of the λi, with its coefficient the even-minus-odd occurrence count; independence forces every coefficient to vanish. Conversely those equal counts make the **formal affine map** the identity. They do not establish nonempty admissibility, much less a polygonal realization. An odd chain instead fixes its candidate input at half its offset.

For the paper's normalized lengths, `σ/τ=3^(1/4)/2` is irrational: a rational value would makeσ rational, but in `a^4=3b^4` the exponent of3 is divisible by4 on the left and congruent to1 on the right. Thusσ andτ are Q-linearly independent. This supplies a valid conditional arithmetic sidecar for the abstract chain model, **not an assertion that v3 uses this argument** and not an application of Fermat's two-square theorem.

A particularly decisive positive-guard hostile is the word `(σ,τ)^3(τ,σ)^3`, in the written left-to-right application order. Its first six reflections induce `x↦x+3(τ−σ)` and the last six induce `x↦x−3(τ−σ)`, so the full formal map is identity. Each type appears three times in odd positions and three in even positions. Nevertheless the fifth output would be `3σ−2τ−x=3σ−4−x<0` for every positive input: `3σ<4` follows by fourth powers from `243<256`. The legal domain is therefore empty. Coefficient balance, exact cancellation, and the positivity guard are separate obligations.

### Cheap hostile examples

- **Order loss:** totals `(3,4)` and `(4,3)` have the same unordered multiset, but at input `3/2` return `5/2` and `1/2`, respectively. Every intermediate split is strictly valid. An endpoint target guard distinguishes them.
- **Domain loss:** `(2,2)` and `(3,3)` both induce the formal map x↦x, but their domains are `(0,2)` and `(0,3)`. At input `5/2`, the first fails and the second succeeds. Sign+offset alone is therefore not compositional feasibility data.
- **Parity conflation:** the sign here counts split reflections. The paper's clockwise/counterclockwise endpoint permutation counts common-side contact transitions. Those are two different event counts. Do not identify the two parities without a proved relation in the specific path.
- **Disconnected transfer:** lengths and local endpoint maps do not determine that all pieces form one intended target region. Global component connectivity remains separate, just as Petersen's local run options do not alone imply one spanning cycle.

## 7. Session70 connection contract

Source: the paper's coupled matching diagrams and the abstract split chains above.
Target: exact named-interface continuation semantics, as in `00-navigation/S70.md` and `05-knowledge/results/session70_attachment_guard_alignment_20261026/session70-root/SYNTHESIS.md` at the pinned base.
Map: replace each equality run by its endpoint transfer, each legal split by its partial reflection, and retain named boundary incidences and the correspondence between the side and corner diagrams.
Preserved predicate: the resulting numerical boundary length response and split legality; with endpoint permutations retained, the specified named local endpoint correspondence as well.
Lost information: piece interiors, embedding, cyclic corner order if omitted, nonoverlap, component connectivity, allowed Euclidean motions, and whether the abstract chain is geometrically realizable.
Required sidecar for larger tasks: joint named port pairing/transport, intrinsic side–corner incidence, metric/angle data, admissibility domains, and a global geometric acceptance predicate.
Cheapest tests: reverse totals3/4; compare equal formal identities on different domains; exchange named endpoints with equal lengths; compare local run covers against actual global connectivity.

The connection is an exact common **operation-relative carrier mechanism**, not an isomorphism between the set of all dissections and tournaments. Session70's all-conjunction rank theorem concerns laws on its native alignment states. It cannot be claimed for Dudeney states merely because some endpoint transports are bits: one must exhibit realizable independent dissection states and native distinguishing contexts first. Graceful overwrite's permitted sign erasure similarly does not license erasing no-flip endpoint orientation here; the observers differ.

## 8. Reproduction and scope of finite checks

Run:

```bash
python3 verify_transport.py --out check.normal.json
python3 -O verify_transport.py --out check.optimized.json
```

The script checks all words of length0..6 over total-length alphabet{1,2,3,4}, exact rational input/output comparisons, every split-position composition, the two hostile examples, and the transcribed pair-list partition. It also checks the coprime charge congruence equivalence for3≤n,m≤40 and all possible corner counts. None of these computations is a search for dissections or an independent exhaustion of the paper's topological atlas. The general affine and charge statements are established by the proofs above, not extrapolated from the bounded tests.

### Live board and outcome

1. Dudeney matching diagrams: inherit exact incidence, not scalar counts; whole theorem remains CITED.
2. Session70 named guards: continuation-specific loss ledger transfers exactly; native all-order alignment claim does not automatically transfer.
3. Fermat/Gaussian norms: arithmetic sidecar only; arbitrary real coordinates block the proposed theorem dependency.
4. Regular odd-gon/square: component charge obstruction transfers; topology and metric exclusion do not.
5. Penrose/ports: angle and edge matching invite a typed comparison, but the actual tile language and local rules must be separately specified.

Meta-patterns used: type every analogy/implication; keep the missing coordinate; retain local incidence through counting; test composition after quotienting. No new universal meta-pattern or priority claim is made.
