# Independent improved-distribution audit

PASS, 2026-11-06. Independently reconstructed the native labelled graph, its40-vertex core and two edge types. Both saved sets are39-edge spanning trees consisting exclusively of actual core edges, with connectivity verified by traversals. Exact original-edge tree-path enumeration gives the stated sums X=(186,36,122,100), Y=(140,96,168,68).

Generated and checked all48 relabelling/converse automorphisms independently. Recounted all paths and canonical loads on all96 image trees. For every one of72 individual core edges, the rational26/83,57/83 mixture has expected stretch exactly267/83; the separate16/55,39/55 mixture has expected canonical load exactly1063/330. All bridge leaves remain exact. Both improvements and the at-most96 support bounds are valid. The prior at-most48 support capacity construction retains that separate advantage.

The proof correctly presents feasible upper bounds only, leaves earlier lower bounds unchanged, does not identify the two distributions, and does not claim search optimality. No preprint theorem or floating calculation is required. independent_check.py is portable, uses only the standard library, reads sibling WITNESS.json and prints independent.out exactly. Frozen base release was not modified.
