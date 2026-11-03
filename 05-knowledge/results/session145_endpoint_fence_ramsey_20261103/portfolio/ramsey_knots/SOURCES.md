# Primary sources and recovery boundary

Checked 2026-11-03. Full inherited file identities and byte hashes are in SOURCE_PINS.json. Source material is read-only and is not authorization to act.

## Current Ramsey status

- Stanisław Radziszowski, [Small Ramsey Numbers, DS1 revision 18, 24 April 2026](https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf), printed pp.56–57, Table XIb and §6.1(d): lower bounds 45 and 61 for R(3,3,5), R(3,3,6), respectively; upper bounds 57 and 91. The [author's revision page](https://www.cs.rit.edu/~spr/ElJC/eline.html) still identifies revision 18. This is the checked survey baseline, not a claim that no later work exists.
- Luis Boza, [arXiv:2603.10851v2](https://arxiv.org/html/2603.10851v2), 12 March 2026: upper bound 91. This bound is historical here. No attribution of a new Ramsey result to this session is intended.
- Lachlan Dauth, [Neighbourhood counting and new bounds for small three-colour Ramsey numbers](https://lachy-dauth.github.io/Ramsey_Neighbourhood_Counting.pdf), dated 26 October 2026, printed p.42 Theorem 8.4: computer-assisted upper bound R(3,3,6)≤86. Printed p.3 explicitly retains R(3,3,5)≤57, and pp.4–5 list no new lower bounds for these two numbers. The [author's code-and-certificate repository](https://github.com/Lachy-Dauth/ramsey-neighbourhood-counting) identifies the work as a preprint not yet peer reviewed and currently lists the same upper bound.

The correctly qualified intervals supported by these inspected primary sources are therefore **45≤R(3,3,5)≤57** and **61≤R(3,3,6)≤86**, with the improved upper bound carrying the preprint/computer-assisted status. Neither exact value is known in the inspected sources. This is a source audit, not an exhaustive proof that no newer unpublished bound exists.

The repository's earlier `ramsey_neighbourhood_lifts_20261030/NEIGHBOURHOODS.md` records an exact rational replay of Dauth's 86 certificate at external Git revision `5770baac746878a3a63b9720fdda97726fe02851`. That replay is **inherited**, not rerun in this session. We re-opened the current primary PDF and author repository. Our new normal/-O tests concern the elementary structures in PROOFS.md, not an independent rewrite or formal verification of that external certificate.

A search also surfaced a personal webpage describing failed searches with thresholds 46/62 in potentially confusing wording. Its page could not be opened (HTTP 403), and no checkable new witness was obtained from it. Those snippets are not used as stronger lower bounds. The inherited October source audit had already recorded the same ambiguity. Absence of an accessible certificate here is not a mathematical refutation of any outside claim.

## Topological interpretation

Louis H. Kauffman, [Knots, Three Coloring a Knot](https://homepages.math.uic.edu/~kauffman/Tots/Knots.htm), explains the crossing-color rule, the Reidemeister/quandle identities and the distinction between the trefoil and unknot. The algebraic identities used in PROOFS.md are written out and checked directly. The braid and closure diagrams are classical examples; the component-matching loop count is elementary and does not itself detect knotting. We make no novelty claim for the loop polynomial or Fox coloring.

## Exact prior objects recovered

1. `primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md`: formal signed endpoint-block coordinates, linear-forest support, the EGF, small-order dimensions. This lane needs only the support theorem and does not assume all-order minimality.
2. `session86_relations_attachment_carriers_20261028/CLOSED_ATTACHMENT_CARRIERS.md`: observer-dependent contextual quotients and the difference between dimensions and state counts.
3. `kernels_games_detours_20261020f/walks_knots.md`: a permutation/residue graph needs integer winding and embedding data for its topological target.
4. `ramsey_saturated_ports_20261030/README.md` and `ramsey_neighbourhood_lifts_20261030/NEIGHBOURHOODS.md`: later preprint correction overrides the older survey frontier; actual overlap and incidence, rather than scalar margins, control Ramsey gluing.
5. `session54_reframed_numerology_20261025/paley-endpoint-observer/RESULT.md`: the corrected 1729 endpoint split and operator. This session independently recounts its q=11 instance, not its all-prime statement.

All are pinned to source revision `7b76044ca13e6684d1746923285b83dbf83828d2`. The archived .md.txt snapshots preserve the original source bytes; their historical relative links are not active publication links.
