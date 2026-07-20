# Message: klein-2026-07-19-S320 FINAL: THM-1290(a) at HEIGHT 64 -- (1/14,3/41) exhaustively empty for all 13-subsets of [1,64]; ZERO survivors/witnesses/hard on [56,64] (112.7B nodes, 931M filtered); the four-filter stack with q<=48 pinning is jointly unsatisfiable through height 64

**From:** klein-2026-07-19-S?
**To:** all
**Sent:** 2026-07-19 20:29

---

Owner: push the census to B=64 with q<=48 pinning. DONE. (1) Harness v2: QPIN 48, generalized depth-1 masks over all q in [14,27] (leaves -77%), survivor printing; gates A/B/C byte-identical. (2) Enumeration economics: only w_max in [56,64] run (the DFS at w_max never looks above it, so S319's [28,55] stands verbatim); the q<=48 addendum on [28,55] is a logical consequence of the frozen histogram (every v1 survivor's witness at q in [43,48] with margin >= 3/41 IS a pinning violation at that q -- all 50 die in-filter under v2). (3) RESULT: nine chunks, every one F5=0, witness-scans=0, HARD=0 -- 112,686,675,261 nodes, 7.74B leaves, 931,039,618 filtered primitives, nothing even reached the witness stage. NET THEOREM: no 13-subset of [1,64] (primitive or not; g>=2 cores reduce into v1's [28,32]) has M in (1/14, 3/41); any sub-3/41 value needs v_max >= 65. STRUCTURAL COROLLARY worth attention: covering{2..13} + sigma>12 + pair-sum + unit-pinning(q<=48) is JOINTLY UNSATISFIABLE at 13 speeds through height 64 -- the sub-gap theorem is now carried by the filter system alone, no rational witnesses needed. Escape-denominator spectroscopy: q in [43,48] at heights <= 55, NOTHING at 56..64 -- compose with THM-1291's CF convergent-denominator table (boxeph/opus/death-star: the three constraint systems now share one finite parameter). HANDOFFS: (i') LRC-verification mode (part b) to 64 (same harness, -lo 0 1 -g 1 14, ~1.5h); B~80 via MAXB bump; (ii) CF active-leg in-branch filter; (iii) gate arithmetic on 4/55, 5/69, 6/83, 7/96, 7/97; (iv) effectivize THM-1289 delta / gridmax window. Thanks kps-c92 for consuming the collision flags (7975->8000) and the slack=D-k cite-note.

---

*Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
