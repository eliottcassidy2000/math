# Prime powers and natural bijections: publication ledger

20 October 2026 UTC. Baseline: [geometry/arithmetic session close](https://github.com/eliottcassidy2000/math/commit/aba3e8fe3140aabd1f1eb23037400d17d5a3e6fe). All session checkpoints are ordinary single-parent commits directly on `main`, with nonforced expected-head updates. No branch or credential was created.

## Research and status checkpoints

1. [research: prove uniform equivariance obstruction for all orders at least five](https://github.com/eliottcassidy2000/math/commit/d782e9417d8458e39c98635002892d647a0c554d) — 11 files, 110,797 bytes.

2. [research: prove switching branch integrality and delimit pairing carriers](https://github.com/eliottcassidy2000/math/commit/23088940185dceeb6d4b099d440f2e4a9c0d63e4) — 21 files, 108,810 bytes.

3. [research: publish pronormal quotient symmetry bijections](https://github.com/eliottcassidy2000/math/commit/a800d5da4c8d82b4c731c6d7cc393dc1a9beb4cf) — 20 files, 124,967 bytes.

4. [research: classify generator-free circulant naturality and definability scope](https://github.com/eliottcassidy2000/math/commit/85c9e1318267042194e5c078a9e0d161908aaa8f) — 19 files, 120,360 bytes.

5. [docs: reconcile equivariance and branch statuses with scoped symmetry extensions](https://github.com/eliottcassidy2000/math/commit/cd7648240e718ae224592ec5a109af61a7c65102) — 8 files, 1,138,482 bytes.

6. [research: prove unbounded Hall deficits from joint-prime bipartite families](https://github.com/eliottcassidy2000/math/commit/25752f463a749e64383e485490bd4aec932b912c) — 11 files, 118,085 bytes.

These six checkpoints contain **90 file versions on 90 distinct paths, totaling 1,721,501 bytes**. Every committed version and current file has been checked against its approved frozen bytes. Totals include full replacement theorem/router files, not just changed lines. Mathematical proof reviews and finite-control receipts are retained in each packet; their scope is not replaced by the publication checks.

## Closing checkpoint

The seventh session checkpoint contains the unchanged approved [map](README.md), this ledger, and the two compact startup-pointer replacements. The complete session deliverable inventory is **94 file versions on 94 distinct paths, totaling 1,785,388 bytes**. No deliverable path is replaced twice within this session, so published-version and latest-file totals coincide. The two incoming commits listed below are additional preserved main history, not counted as our seven checkpoints or as our deliverable bytes. The closing commit deliberately does not embed its own hash.

Staged documentation checks pass with startup **179,723 / 180,000 bytes**, frontier **450 / 450 lines**, and guardrails **190 / 190 maintained lines**. The pointer replacements have zero byte and line growth. The final procedure repeats exact Git-blob/SHA256 verification, relative-link checks, ordinary-parent lineage, the documentation gate, and clean local/origin/GitHub equality. Final commit and pinned map URL are reported separately. GitHub statuses are checked directly; an empty status list is not a passing CI run. No Lean rebuild is claimed.

## Mathematical status and scope

- HYP-9172(a) is PROVED for every n>=5 in its precise formulation: no S_n-equivariant map on labelled objects induces a type bijection from tournaments to sign-even graphs, or from switching classes to sign-even Euler graphs. This does not rule out arbitrary or differently natural bijections between unlabelled types, nor settle every interpretation of the broader literature request.
- The quantitative append proves unbounded Hall deficits for its stated joint-prime complete-bipartite families. The original uniform proof is unchanged; the append has its own HALL_GROWTH_SHA256SUMS inventory and independent review.
- HYP-9172(b)'s numerical branch integrality is PROVED for every n>=3, with the n=2 half-integer exception. A structural decorated pairing remains OPEN. Twenty-seven even-Aut types at n=8 cannot be paired by a fixed-point-free involution, and N_lev(8)=15 is not half their count.
- HYP-9340's unrestricted odd-group equality remains REFUTED. All odd p-group cases survive; the new theorem extends to odd pronormal bases and one-prime normalizer-quotient extensions. Quotient q-group symmetry does not require its full inverse image to be a q-group.
- Circulant claims distinguish the existing all-prime algorithm, abstract multiplicative-group naturality, field naturality, and first-order definability. The parity projection's nondefinability does not exclude every prime-field rule. The prime-power Frobenius obstruction is confined to the stated multiplier/necklace model; it is not silently transferred to fully sign-even graph types or rings Z/p^r.
- Current HYP-9172, HYP-9340, THM-4531, THM-479, THM-4721 and direct routers are synchronized with those scopes. Old finite certificates and HYP-9340 sequence data remain. The originally unclosed HYP-9172 metadata now has a closing YAML delimiter before its new body.

## Incoming work preserved

Two concurrent commits arrived during publication:

1. [cc266a5fdde7](https://github.com/eliottcassidy2000/math/commit/cc266a5fdde738269705785552528b9097500f28), session `mac-mini-2026-10-19-hyp9172`: reserved THM-4722, THM-4723 and HYP-9341 as empty stubs and adjusted the hypothesis-index budget. The owner's reservation/index changes were preserved; no new IDs were reserved here.
2. [b9aac7f88921](https://github.com/eliottcassidy2000/math/commit/b9aac7f88921f2bf18f060fccc8fc5ed80db56a5): supplied a separate pre-audit THM-4722 Jordan/Hall proof for tournaments at all n>=5 and an odd-n Euler-parity route for switching classes, with scripts, certificates and a note. Its theorem and proof artifacts remain untouched. The then-current HYP-9172 status is explicitly recorded as concurrent provenance; the independently audited uniform leaf controls the stronger all-n current status.

The exact-baseline check detected the incoming HYP-9172 change before uploading or committing the stale status overlay. Publication paused, compared the incoming work, preserved its contribution, and resumed only from a freshly reviewed manifest against the new baseline. No stale overwrite or forced ref update occurred. The two large router blobs were each uploaded once and hash-verified. No uncertain or pending operation was duplicated.

Future incoming proofs and reserved-node promotions remain their owner's work. Final reconciliation checks for additional incoming commits before reporting the completed head rather than overwriting them.

## Prior session maps

- [Geometry, arithmetic content and actual-source constraints](../geometric_arithmetic_session_20261020/README.md)
- [Earlier maps and publication chain](../geometric_arithmetic_session_20261020/PUBLICATION_LEDGER.md#prior-session-maps)

[Current prime-power/bijection map](README.md)
