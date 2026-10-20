# Validation and scope

**Status: exact reproducible computations and scoped proof audits.**
The session used Python 3.13.14 and NetworkX 3.5 on Windows. Only the
Petersen checker needs NetworkX, for graph isomorphism and diagnostics.
Colour extensions and full-graph Tait counts use separate exact routines.

## Mathematical universes and controls

| Experiment | Declared universe | Independent or hostile controls |
|---|---|---|
| [Budget transitions](budget_check.py) | Complete capped states at transitions \(0\to1\to2\); every arc of every tournament class of orders 2–6 | Equality with the earlier endpoint compiler; 991 literal full-vertex flips; weighted-feedback checks; heavy-location and endpoint-offset hostiles |
| [Deficiency strata](strata_check.py) | All canonical states of deficiencies 0, 1, 2 at \(k=1,\ldots,4\), generated from the proved finite menu | Complete endpoint-compiler equality at \(k=1,2\); all five core-order sectors of the local collision formula; all three proved coefficients; normalization-loss and duplicate tracking |
| [Blue rerouting](blue_check.py) | Literal order-five certificates; all 187 singleton/triangle outer words of weights 1–12; blue presentations and inherited line tables through order seven; canonical-SCC nonstrong blue presentations at orders eight and nine | 561 trivalent hub routes; 608 padded mixed-center controls; regular rotation failures; fixed-endpoint hostiles and shortest black BFS; unrealized wiring kept separate |
| [Petersen reductions](petersen_check.py) | Six specified donor cuts; all 75 independent-edge pairs of labelled Petersen; all three such pairs in \(K_4\); 100 successive insertions | Labelled colour backtracking versus matching/even-cycle full-graph counts; explicit rooted paths; two authenticated Blanuša banks; empty-support and wrong-pairing hostiles |

The normal and optimized versions of all four experiments produce identical
transcripts. The sources use explicit exceptions, so `python -O` does not
remove the checks. Saved outputs sit beside the scripts.

The budget and stratum compilers share the inherited clone normalizer and
weighted-core canonicalizer. Their generation routes are independent, not
their entire code stacks. The direct arc controls compare the bounded
refinement with actual full-tournament reversal. Proofs, rather than a
finite census, establish the all-budget statements.

The blue census replay reads the committed THM-4731 TSV tables and does not
claim a newly generated order-eight census. The all-order rerouting proof
uses the independently repaired pure-blue classification and its stated
rigidity dependencies. An additional agent checked the six-vertex shortest
detour with fresh adjacency matrices and NetworkX isomorphism, independently
of the inherited tournament routines.

The Blanuša names refer to the two explicit edge banks inherited from
`04-computation/experiments/snark_boundary_transfer_20261018a.py`.
The checker copies their edge lists, not their algorithms. The 75-case
identification does not claim to classify all snarks of order eighteen.
The hundred-step subdivision audit checks the implementation; the proof
establishes certificate composition at arbitrary depth.

## Proof audit

The root independently checked the cap-lift quantifiers, the two-singleton
split construction, and completeness by reversing an optimal backward arc.
A separate agent audited the deficiency identity, exceptional-motif bound,
and completeness of its four deficiency-two presentation types.
The third-coefficient proof was separately audited by the root and a second
agent. Both checked the five possible core-order sectors, the unique
three-description triangle collision, and the exhaustive clone merger
cases. The all-budget promotion relies on this local proof, not a fit.

The root and the Petersen agent independently audited the blue rotation
degree inequality, the common-hub argument, and the padded-center proof.
The six-vertex counterexample refutes only the proposed fixed-boundary
shortcut. Its blue edge has a verified black detour; it is not a bridge.
The shared-edge cut repair and the larger fixed-endpoint hostile family
retain their explicit boundary conditions. The incoming endpoint-current
gate is tested against the abstract rewiring and rejects it.

The root checked the four-port Kempe-switch bijection, rooted wire extraction,
three-pole claw reduction, and composition of subdivision certificates.
The support theorem itself is inherited from **THM-4116**,
[ordered boundary gluing](../../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md).
The primary *Factorisation of Snarks* paper was opened and its relevant
section and theorem numbers checked during this session; the dot-product
and colour-type language is prior work, not a novelty claim.

## Reproduction

From the repository root:

```text
python 05-knowledge/results/defect_transitions_20261020/verify.py
python 05-knowledge/results/defect_transitions_20261020/verify.py --replay
python agents/check_docs.py
git diff --check
```

The first command checks the packet and computational dependencies against
the normalized SHA256 manifest and checks source syntax. The replay runs
each experiment in normal and optimized Python and compares its complete
stdout with the committed transcript. CRLF and LF line endings are treated
equally. The manifest is provenance, not an independent mathematical proof.

The shared-state compression and global graph claims are scoped separately
in [README.md](README.md). No finite check here proves general black-core
connectivity or a universal reduction of arbitrary snarks.
