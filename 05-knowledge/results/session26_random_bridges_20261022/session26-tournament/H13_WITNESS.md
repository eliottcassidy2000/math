# A five-vertex H13 witness outside the one-apex family

The exact arc list and all13 Hamilton paths appear in H13_WITNESS.json. This is FINITE-EXACT and elementary, not a new spectrum-priority claim. On vertices0..4, the increasing-pair bitmask is41, with bit1 meaning the smaller label points to the larger. Independent subset DP gives13; literal path DFS lists exactly13.

The directed cycle0->4->2->1->3->0 proves strongness. No vertex deletion is transitive: after deleting0, use triangle1->3->2->1; after deleting1 or3, use0->4->2->0; after deleting2 or4, use0->1->3->0. Deleting0 and1 leaves the transitive chain4->3->2, so the feedback vertex number is exactly2.

Thus13 is genuinely attained in the full tournament spectrum but absent from the whole feedback-vertex-at-most-one spectrum by CUT_ADDENDUM.md. The obstruction is structural, not a modular obstruction, and does not depend on the global spectrum conjecture.
