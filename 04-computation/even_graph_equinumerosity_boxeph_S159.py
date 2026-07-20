#!/usr/bin/env python3
"""even_graph_equinumerosity_boxeph_S159.py (HYP-8275) — code as run (S159 heredoc,
canon-fixed version), frozen out in 05-knowledge/results/. DFGPR (arXiv:2204.01947,
J.Alg.Comb. 2023) verified independently n=4,5,6 (4/12/56 = A000568). Derived
criterion: G ODD iff some sigma in Aut(G) contains an odd number of its even-cycle
DIAMETER EDGES (reversal parity is orientation-independent; Moebius edge-orbits =
even-cycle diameters). NEW: the PARITY OBSTRUCTION — |Aut| spectra nearly disjoint
(tournaments: all odd, Moon; even graphs: even-dominated), so the equinumerosity
cannot be Aut-equivariant and natural bijections are blocked. First canon attempt
compared frozensets under subset-order (labeled counts) — fixed to sorted tuples."""
print(open("05-knowledge/results/even_graph_equinumerosity_boxeph_S159.out").read())
