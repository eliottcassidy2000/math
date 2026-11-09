# Session235: reciprocal arcs, missing pairs and the tournament boundary

**Status: PROVED (pair-state dictionary); FINITE-EXACT (complete oriented
census through six, with independently checked controls).** 2026-11-09.
This packet continues [Session234](../session234_digraphs_20261109/README.md).
The present checkpoint publishes the exact finite baseline; all-order
extensions and final independent review will be added after verification.

The key distinction is **tournament = no reciprocal pairs and no missing
pairs**. Forbidding reciprocal arcs alone produces oriented graphs, which
are subgraphs of tournaments. Missing pairs remove the path-count gaps:
seven paths first occurs at order five, and twenty-one at order six.
Exactly one missing pair suffices in each case.

Treating raw 11 as empty is a valid alternative encoding, not an ordinary
adjacency rule. Since the question does not specify 00, we keep three
explicit alternatives: forbid 00; also interpret it as empty; or swap its
meaning with 11 so 00 represents both directions.

- [Pair-state dictionary and exact encoding counts](PAIR_MODELS.md)
- [Complete oriented census and witnesses](census/REPORT.md)
- [All 408 H7/H21 classes through six, explicit paths and completions](census/witnesses.json)
- [Independent small model verification](models/validation.json)

The complete labelled oriented counts are 1,440 H7 graphs at order five,
and 255,600 H7 / 29,520 H21 graphs at order six. Their isomorphism counts are
12, 355 and 41 respectively. These are finite censuses, not an all-order
isomorphism classification.

One particularly revealing raw five-vertex array has H=21 under ordinary
reciprocal semantics but H=7 when 11 means empty. Its four tournament
completions have H=11,13,13,15. Full witnesses and the complete joint
histograms are retained in the census, rather than conflating the two H's.

No conclusion here settles the Chan–Pak tournament spectrum conjecture,
which remains OPEN. Changing the graph class changes the question.
