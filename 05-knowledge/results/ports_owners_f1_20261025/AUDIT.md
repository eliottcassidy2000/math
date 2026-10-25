# Audit and reproduction ledger

Baseline: origin/main `2b96df316e`, fetched 2026-10-25 in the clean isolated
tournament-metagraph worktree. The later refresh before integration found
no incoming commits. Source status was recovered from proofs and explicit
corrections, not historical summaries. Targeted processor inspection was
unavailable because this worktree has no `.machine-id`; none was invented.

## Proof reviews

- **Anchor:** the investigator derived the fixed-rank augmentation and
  parent forest. The root separately checked cross-rank transitivity,
  finite comparison delay, existence of a rank-lowering deletion, exact
  owner quotient, finite-column Ramsey bound, and the growth obstruction
  to P-recursiveness. The review added the explicit G-to-L recanonicalization
  between successive transfers; paired canonical words, not nondeterministic
  witness paths, are counted. General deterministic machines are not built.
- **Gap:** the investigator derived the band formulas. The root checked
  the deficit identity, exhaustive matching/P3 alternatives, intrinsic
  weighted-vertex order, unit-cell ownership, external-partner guard for
  strong T4, clone rejection terms, and J=zI'-I+I^2. The final review
  includes delta=2k-1. The independent quotient code does not reuse the
  producer's cell owner, but does reuse its structural record generator;
  this limitation is explicit.
- **Endpoint:** the investigator supplied the complementary ledgers and
  prefix extension. Both root and the separate gap investigator reviewed
  actual joins, path-edge XOR, all three anti-reversal contradictions,
  the separate n=11 witness, and the n>=2h+8 guard. No assumption about
  x versus y was added. Actual vertices remain present throughout.
- **F1:** the root supplied the conic/nested-set construction and scoped
  unification map. The anchor investigator independently checked the
  Veronese charts, integral power formulas, composition, primitive scale,
  characteristic-two bad model, exceptional-set orbit counts and growth,
  and Euler-product identity. It reread Connes--Consani section4 directly
  and checked the stated moment normalization. It also independently ran
  the exact checker. No RH positivity is proved.

## Corrections and hostile controls

The F1 review corrected three prepublication type errors: d=1 alone does
not put a cardinality in the old Boolean family; arbitrary integer
coordinates give a rational projective point before content normalization;
and polynomial bounded-proof search needs a unary length cap. The formulas
survived unchanged. Mechanisms and repaired scopes are recorded in
`01-canon/MISTAKES.md`.

The endpoint inheritance pass found that THM-4734's title and body still
said n>=6 after MISTAKE-646 had already repaired the universal threshold
to n>=7. We corrected that stale text and added the literal pure-blue-target
witness, preserving the original pre-audit receipts/status. This is a
scope repair, not a new order-six result or a promotion of the old census.

Positive and hostile controls include the regular-pentagon projection
failure, deletion-created C3 clone pair, wrong-direction endpoint ledgers,
pure-blue order-six target, the 45/42 missing-gap rejection, the two
strong-T4 endpoint roles, the first conic power (63,16,65), and the bad
characteristic-two projective base point.

## Computational scope

| Checker | Universe and independence |
|---|---|
|anchor_check.py|All tournament classes through order7, with imported pinned class generation; rank independently by transitive subsets and triangle hitting;525 positive-rank parents; actual port permutations, owner normalization, clone disagreements; complete alpha3 forest, HP by permutations and subset DP|
|endpoint_check.py|All16,384 relevant n11 second-port sources;1,168 at n8..10; n7 diagnostics;1,120 deterministic larger controls;1,156 prefix controls including34 exact-threshold cases. Arc-set and bit-row evaluators; no producer imports, SCC/mixedness/isomorphism filter|
|gap_check.py|Complete declared structural record layers including k3 delta4 and k4 delta5,6; all literal clone tests, ownership and strongness; k2 independent weight/degree quotient and exact feedback DP;7,739,819 explicit gates|
|gap_owner_audit.py|Every k3 delta4,5 record, plus all k4 delta7 missing-gap/unit-gap records. Imports generation but independently uses full weight/outdegree permutation cells, literal clone rows and reachability;46,725 explicit gates|
|f1_check.py|Projective conic bijections over six prime fields,602 primitive parameter/power cases,18 explicit nested-set orbit universes, coefficient Frobenius controls and local zeta recursion; no producer imports|

All five checkers use explicit exceptions; no proof gate disappears under
Python -O. Normal and optimized stdout agree with the saved receipts.
Larger random controls use fixed seeds and are labelled samples. Neither
those samples nor finite numerical prefixes establish all-order theorems.

The exploratory endpoint scripts are retained for provenance and are not
part of the theorem replay. Their results are not used as a substitute for
the independent final evaluator or the symbolic proof.

## Reproduction and sources

From the repository root:

```
python 05-knowledge/results/ports_owners_f1_20261025/verify.py --replay
python agents/check_docs.py
```

`manifest.json` pins normalized LF bytes of all packet text/code/receipts
and the named inherited inputs. `--record` is a maintainer operation after
review; normal verification refuses changed pins. The required replay list
is explicit, checks each producer in normal and optimized Python, and
compares stdout byte-for-byte after newline normalization.

External primary sources were read on 2026-10-25: Borger, *Lambda-rings and
the field with one element*, sections1,2.4; Lorscheid, *Blueprints—towards
absolute arithmetic?*, its stated RH research motivation; Connes--Consani,
*Geometry of the arithmetic site*, abstract, and *The Scaling Hamiltonian*,
sections3–4; Cook's official *P versus NP* description, definitions of
polynomial reductions and verifiers. Exact URLs and imported scopes are
in [f1_bridge.md](f1_bridge.md). The elementary conic and finite-set proofs
are written out locally; no new priority claim is made for standard
toric/lambda constructions or elementary growth bounds.

Outstanding: unrestricted counting compression, all endpoint-four cases,
all-order black-core connectivity, and the proposed RH-facing positivity
representation. The classification of proof duties is not a cross-problem
equivalence theorem.
