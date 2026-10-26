# Independent audit of the unprescribed six-colour boundary packet

2026-10-26. **PASS.** Independent proof review and finite structural probe by
the information lane. No mathematical correction or producer-file edit is
required. Parent integration controls status promotion.

Audited [six_color/RESULT.md](../six_color/RESULT.md) and
[six_color/check.py](../six_color/check.py). The checked Python source SHA256 is
`04e865fc9972d292fe109efe5ada8e5fcc14ed90534906661d2ca8bd6dac10da`.
The producer replay exactly matches its saved receipt. The parent's separate
[colour-class recurrence](../audits/check_six_color_independent.py) supplies
an independent check of the numerical extension weights; this audit focuses
on unrestricted proofs, chart restrictions, and the gluing hypotheses.

## 1. Unrestricted symmetry and chart scope

The single-frame statement is valid for every finite graph, including graphs
with no six-colouring. Applying one global palette permutation to an entire
extension is a bijection. Consequently the three-port count depends only on
the equality partition, and all eight disjoint-pair endpoint tuples occupy
one orbit. There is no assumption that an arbitrary graph has this pairing
as a graph-enforced domain. The statement concerns the restriction of its
ordinary boundary table to that specified slice.

For two frames with the same corresponding pairs, a simultaneous endpoint
flip in any pair is a legitimate global palette permutation. It acts on both
frames, proving the difference-kernel law, for counts as well as support.
On the even completion only even simultaneous flips stay in the completion;
these are the regular additive action of the four-element group E. Its
centralizer among deterministic maps consists of translations: putting x=0
in T(x+b)=T(x)+b gives the formula immediately. The 16 Boolean kernels and
four deterministic permutations are therefore exact necessary classifications
of invariant tables, not graph realizability claims.

The port order and corresponding pair assignment matter. With the three port
positions held fixed, the setwise stabilizer in S6 of the eight-tuple endpoint
chart has order 8, and that of its even half has order 4. The independent
probe enumerates all 720 permutations and verifies both. The larger abstract
S4 action does not thereby become a palette action on this fixed chart.
The result explicitly permits changed charts, reordered ports, and larger
retained boundaries as different typed problems. It makes no unrestricted
prohibition on implementing a four-state relation with additional palette
reference vertices.

## 2. The fresh-colour proof and the chromatic thresholds

The bound k+t+max(0,s−1) is valid for any finite boundary and any proper desired
boundary partition, not just for independent three-port sets. Each nonsingleton
block is independent by propriety and receives its own colour absent from the
old colouring. A repeated old colour among singleton blocks is replaced by a
different unused colour; the first singleton of each old colour can remain.
No newly equal pair lies across an edge, and no old interior colour is changed.
The final palette permutation is legal because the boundary colours are
distinct precisely on the prescribed blocks. When the target palette has at
least the displayed bound, the injection of used colours into it extends to
a permutation of that palette.

For three ports the fresh-colour costs are 1 for each nondiscrete partition
and 2 for the all-distinct partition. This proves all stated four- and
five-chromatic conclusions. If an independent boundary's all-equal assignment
does not extend at six, then its exterior cannot be five-colourable: otherwise
one fresh sixth colour colours all three ports. If the whole graph is
six-colourable, this actually makes the exterior six-chromatic, so the claimed
necessary lower bound is valid and sharp.

The independent test does not repeat the producer's recolouring algorithm.
It contracts each target boundary block and makes the distinct block vertices
a clique. The resulting chromatic number is the least palette size permitting
that equality pattern. It checks every simple graph on 3, 4, and 5 labelled
vertices, with the first three as ports: 1,096 graphs and 2,877 proper partition
cases, without a chromatic or geometric filter. All satisfy the bound.

The two abstract sharp examples have minimum extension-palette vectors
respectively (5,6,6,6,7) and (7,6,6,6,6), independently computed by exact
chromatic numbers. Their whole-graph chromatic numbers are 5 and 6. The NAE
example has a six-chromatic exterior. Their K4 subgraphs prevent unit-distance
embeddings in the plane: the three based displacement vectors would have the
positive-definite rank-three Gram matrix specified in the proof. Hence these
examples establish abstract sharpness only, as the packet states.

## 3. Gluing keeps the necessary graph decomposition

The corollary is valid with its stated hypotheses. First colour the entire
retained port graph H. Every fragment sees at most three ports, and this
restriction is proper because every boundary edge has been retained in H.
Apply the four-colour fragment extension independently. The extensions agree
on every overlap, since overlaps are ports whose values were fixed first.
Private interiors are disjoint, and every edge is either in H or in a fragment;
there is therefore no unexamined cross-interior edge whose colouring could
invalidate the combination.

This argument requires that H itself be six-colourable. It also requires the
literal graph decomposition: private interiors, at most three ports per
fragment, and no missing cross-unit edges. The wording explicitly states these
requirements. It does not say that an arbitrary union of low-chromatic graphs
is six-colourable, or that a conditional frame loop supplies such a decomposition.
The conclusion covers support existence, not constant extension counts.

## 4. Geometry, the existential pair sector, and unpinning

The seven coordinates have exactly the induced graph C3+C3+K1. The conditional
assignment of one colour pair to each entire triangle is inconsistent, since
each triangle cannot be two-coloured. Expressing this as three gains 110 is
correct and parity-preserving. Its inconsistency is a property of that fixed
pair sector. The genuine graph is three-colourable, and the packet gives an
explicit colouring, so there is no implication to a six-colour obstruction.

The independent ternary calculation is direct: each triangle has 18 pair maps
using exactly two pair names, giving two components in the same-pair subgraph,
and 6 maps using all three names, giving three components. The isolate has
three choices. Thus the component histogram is

    (5:18²·3, 6:2·18·6·3, 7:6²·3) = (5:972, 6:648, 7:108).

This gives 1,728 successful sectors and 86,400 six-colourings after weighting
by 2 to the number of components. It checks the quantified pair assignment,
not just the failed prescribed one.

The leaf factors in the old gadgets are also legitimate. The probe checks
directly in the prior coordinate receipt that every formerly pinned vertex
has degree one and its unique neighbour is retained, rather than another
formerly pinned leaf. Each therefore has exactly five independent choices
after the retained graph is coloured. This rules out a hidden isolated
leaf-pair component, which would require a different factor. The exact metric
recomputation and retained-core checks in the producer then justify its
conversion from prescribed tables to ordinary extension counts.

The final target scope is accurate: these results establish necessary gates
and exact counterexamples to a proposed compiler. They produce no six-colour
equality gadget, no non-six-colourable unit-distance graph, and no new
Hadwiger–Nelson bound.

## 5. Receipts

Run from the repository root:

```powershell
python 05-knowledge/results/hn_metagraph_information_20261026/information/check_six_color_audit.py
python -O 05-knowledge/results/hn_metagraph_information_20261026/information/check_six_color_audit.py
```

[Independent source](check_six_color_audit.py) uses no producer function and
has 6,364 always-active checks. Normal and optimized receipts agree:
[normal](six_color_audit_receipt.json),
[optimized](six_color_audit_receipt_optimized.json), SHA256
`f389d1b6959b485d0d3ed834f4220ce39586beb3a90b29765070435e8af44f4e`.
The [producer replay](six_color_producer_replay.json) agrees byte-for-byte with
the producer receipt, SHA256
`f478b578c06f4908e280407f5282d18c632d3eb1937ed4c2caa1787395a13b2e`.

All-order acceptance rests on the reviewed proofs. The finite probes cover
the explicitly listed universes, not all finite graphs or all unit-distance
realizations. **Recommendation: promote the unrestricted statements in their
stated scopes; retain FINITE-EXACT labels on the numerical censuses.**
