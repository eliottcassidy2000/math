# Three consecutive requested ready targets with one helper

**FINITE-EXACT physical construction:**881 ->1691 ->3311 ->6551. Each of the last three targets has an actual strong unit-ready endpoint. Every intermediate adjacency matrix, physical name map, patch, birth, and deletion is independently recounted in [chain_881_1691_3311_6551.json](chain_881_1691_3311_6551.json). This is a finite program; no all-k return or Chan-Pak proof is claimed.

**Incoming correction to the research target (commit c85ab9a188):** an all-k indegree-one unit-ready405 ray is **REFUTED**. The [complete intrinsic inverse](../../target_directed_cores_20261112/TARGET_INVERSE.md) and its [independent exhaustive audit](../../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md) exclude every tournament order at53,084,231=405*2^17+71. The ordinary Hamilton count is attained, and Chan-Pak remains OPEN. The finite constructions and module formulas here remain valid; future continuation must change the required interface or skip an obstructed target.

## Actual operations and count budget

Start from the ten-vertex A2 source in [PROOF.md](PROOF.md). Add one helper beaten only by old vertex3. The helper has physical name10 and raises H from881 to1589. Thereafter all old-arc repairs use the same two physical vertices0,2, and all target births are beaten only by physical vertex0.

| Stage | Selected profile(H,s,q,F) | Strong? |
|---|---|---|
| Helper above881 | (1589,698,881,1772), marked helper10 | yes |
| First repaired parent | (351,1,3,1341), port0 | yes |
| First requested target | (1691,350,351,1692), newborn11 | yes |
| Second repaired parent | (853,1,5,2459), port0 | yes |
| Second requested target | (3311,852,853,3312), newborn12 | yes |
| Third repaired parent | (495,1,33,6057), port0 | no: one dominated singleton |
| Third requested target | (6551,494,495,6552), newborn13 | no: the same singleton |
| Final endpoint | (6551,494,495,6552), same physical root13 | yes, after deleting old vertex8 |

The three patches contain5,6,12 arcs, respectively:

```
first:  {0,9},{1,2},{2,4},{2,5},{2,9};
second: {0,2},{0,9},{0,4},{0,7},{2,3},{2,11};
third:  {0,2},{0,9},{0,4},{0,6},{0,7},{0,8},
        {0,11},{0,12},{2,9},{2,4},{2,11},{2,12}.
```

Thus the old Hamilton count deliberately falls1589->351,1691->853, and3311->495 before the next target birth. These noninduced decreases are the mechanism that escapes the pure-birth obstruction. The separate intermediate arc reversals are not claimed to preserve H or readiness.

Each target birth has a unique-start parent, so the new root has indegree1 and q-s=F-H=1 exactly. The last old vertex8 is a literal global sink after the third patch and remains one after the birth. Its deletion preserves the full rooted state; it can equivalently be deleted immediately before the last birth. Four births and one deletion take the original order10 to final order13.

The final outgoing masks, with retained names0,1,2,3,4,5,6,7,9,10,11,12,13, are

```
(4166,128,2274,1015,231,195,130,1,247,503,1023,2043,4094),
root array index12 = physical name13.
```

## Named branch transport is not endpoint replacement

The first repaired parent has one Hamilton path beginning at0. Its unique order supplies the canonical-to-actual name map

```
(0,2,3,10,9,4,5,6,1,7,8,11).
```

The replay verifies the whole adjacency relation under this map against the inherited1691 representative from parent(order11,rank303). It then transports the discovered second star patch through this map while retaining the actual graph and vertex names. The third patch uses the same map extended by the second newborn. No step discards the current graph and substitutes a separately enumerated target.

The two-star positive searches stopped on their first hits. Their finite search prefixes do not classify all successful programs. In contrast, the earlier failure claims for the smaller rank385 branch are complete for their declared one-star/two-star menus. Scalar1691 and unit readiness do not determine continuation capability.

## Exact four-letter representation with fixed physical ports

Removing ports{0,2} leaves a transitive tournament at every stage. Immediately after the helper its ordered physical names are

```
(3,10,9,4,5,6,1,7,8).
```

For each backbone vertex u encode the letter([u->0],[u->2]). Keep the separate orientation of the port edge. The backbone order, letters, port edge, actual names and active root reconstruct the entire marked graph. The independent `two_port_trace.py` checks this representation exactly.

| Stage | H | Port edge | Four-letter word along the transitive backbone |
|---|---:|---|---|
| source881 | 881 | 0->2 | 10 01 11 11 10 01 10 10 |
| helper1589 | 1589 | 0->2 | 10 11 01 11 11 10 01 10 10 |
| parent351 | 351 | 0->2 | 10 11 10 10 10 10 00 10 10 |
| ready1691 | 1691 | 0->2 | 01 10 11 10 10 10 10 00 10 10 |
| parent853 | 853 | 2->0 | 00 11 11 00 00 10 10 00 00 10 |
| ready3311 | 3311 | 2->0 | 01 00 11 11 00 00 10 10 00 00 10 |
| parent495 | 495 | 0->2 | 10 11 11 11 11 11 10 00 00 10 00 |
| ready6551_weak | 6551 | 0->2 | 01 10 11 11 11 11 11 10 00 00 10 00 |
| ready6551_strong | 6551 | 0->2 | 01 10 11 11 11 11 11 10 00 00 10 |

Each target birth at0 prefixes a new01 letter and makes that newborn the active root. The first helper inserts11 after the initial backbone vertex3. Every repair changes only the two old-port comparisons and their mutual edge; all remaining old pairwise relations are retained. The final neutral deletion removes the last backbone vertex8, then a00 letter.

The backbone/word representation remains an exact two-port compiler for this finite trace. Its future predicate must change: a selector visiting every405-ray count with unit readiness is impossible because k17 has no such source at any order. A regular selector or all-height accepted language that permits a nonunit defect, another marked interface, or skips the obstructed count remains OPEN. The full-cylinder doubling obstruction alone did not exclude a thin language; the new intrinsic target obstruction excludes this particular all-k unit-ready language even without a word-compiler restriction.

## Reproduction

Run each of `chain_881_1691_3311_6551.py` and `two_port_trace.py` normally and with `-O -B`. Their receipt pairs are byte-identical. `search_replay.py` compiles the C++ searches into a temporary directory, repeats all eight declared searches, and independently recounts all12 positive parent/newborn controls in Python. No executable is included in the repository.

The growing-module family in [PROOF.md](PROOF.md) uses the different183->205 parent branch. Its universal coefficient theorem and positive offset tax should not be attributed to this351->853->495 finite chain.
