---
id: THM-843
title: The n=8 merged-metagraph flow and the radius-one shell of the n=9 defect support
status: RESERVED — exact full-atlas computation in progress; no theorem claim yet
source: codex-2026-07-15-S13e
depends_on: [THM-781, THM-785, THM-790, THM-834]
related: [THM-811, THM-830, THM-842, HYP-6825, HYP-6880]
planned_verification:
  - 04-computation/n8_metagraph_flow_defect_shell_codex_S13e.py
  - 05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.out
  - 05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.json
---

# THM-843 — reserved n=8 flow and defect-shell census

This namespace is reserved for an exact replay over the full `2^21` size-eight
tiling atlas.  The computation will retain literal complement-line
multiplicity, projected merged-node support, blue/black colour, cyclic-triangle
flow, node phase, and the radius-one neighbourhood of THM-834's 155 face
nodes.  It will test monotone and two-phase reachability from the transitive
node and quantify where blue incidence first re-enters the black-only selected
line bank.

Known only from a preliminary read-only probe: the black radius-one shell
meets 3,308 merged nodes, and blue incidence touches exactly two of the 155
support nodes.  These figures remain provisional until the stored verifier and
assertions are committed.
