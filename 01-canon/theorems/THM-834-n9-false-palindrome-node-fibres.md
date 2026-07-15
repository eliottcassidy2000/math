---
id: THM-834
title: The 58 n=9 false-palindrome pairs occupy 53 converse-merged tournament nodes
status: RESERVED + FINITE-EXACT FORWARD MAP COMPLETE; metagraph-flow audit in progress
source: codex-2026-07-15-S13b
depends_on: [THM-810, THM-828, THM-832]
related: [THM-796, THM-809, HYP-3809, HYP-6880]
planned_verification:
  - 04-computation/n9_false_palindrome_node_forward_map_codex_S13b.py
  - 05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.out
  - 05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.json
---

# THM-834 — reservation: false palindromes inside node fibres

This number reserves the exact size-nine tournament forward map for THM-828's
116 endpoints.  They occupy 53 converse-merged nodes.  Forty-eight nodes
carry one collision pair and five carry two; all five repeated nodes lie in
the dominant difference sector `4c41818`.  Forty-nine of the 53 nodes are
self-converse, and 54 of the 58 pairs lie inside one self-converse ordinary
class; four bridge two distinct converse classes.

The completed theorem will include canonical-code verification, node/sector
fibres, score-axis and Hamiltonian-path placement, and the blue/black flow
guardrail.  A collision sector is not asserted to be an isomorphism node.
