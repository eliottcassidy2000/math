# THM-006: F3 — Minimum Type-II Count for L-Cycle Embeddings

**Type:** Lemma
**Certainty:** 5 — PROVED
**Status:** PROVED
**Last reviewed:** SYSTEM-2026-03-05-S1
**Disputes:** none
**Tags:** #type-ii #l-cycle #embedding #minimum

---

## Statement

For any odd L-cycle consecutively embedded in path P' (i.e., the L-cycle's non-v vertices appear in a consecutive window of length L−1 in P'), the signature window s[j..j+L−2] has:
- s[j] = 1 (first vertex of window: v→P[j])
- s[j+L−2] = 0 (last vertex of window: P[j+L−2]→v)

Therefore, the window contains **at least 1 Type-II position**.

The minimum (exactly 1 Type-II position) is achieved by the "monotone" signature pattern 1,1,…,1,0.

---

## Proof / Proof Sketch

The L-cycle uses vertices v, P[j], P[j+1], …, P[j+L−2]. The arcs of the cycle are:
- v → P[j] (so s[j] = 1)
- P[j] → P[j+1] → … → P[j+L−2] (the path arcs, not determined by the cycle orientation alone — these are fixed by P' being a path)
- P[j+L−2] → v (so s[j+L−2] = 0)

Since s[j]=1 and s[j+L−2]=0, the signature must transition from 1 to 0 somewhere in the window, giving at least one 1→0 transition, i.e., at least one Type-II position. □

The monotone pattern (1,…,1,0) achieves exactly one such transition (the final position). □

## Verification Scripts

- `04-computation/per_cycle_identity.py` — Type-II position analysis for cycle embeddings
- `04-computation/ballot_sequence_test.py` — verifies monotone signature patterns achieving minimum

## Notes & History

This lemma establishes that every cycle embedding "witnesses" at least one 3-cycle-detectable Type-II position. It's a key component in understanding how the per-path identity can potentially capture the contributions of longer cycles (though as OPEN-Q-001 shows, this capture is imperfect for n≥6).
