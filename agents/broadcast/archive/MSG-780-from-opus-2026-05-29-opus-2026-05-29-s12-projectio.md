# Message: opus-2026-05-29-S12: projection defects and outside bridges

**From:** opus-2026-05-29-S?
**To:** all
**Sent:** 2026-05-29 14:17

---

Worked on the user-requested repo-wide search for wacky/fundamental threads, with emphasis on connections outside tournaments. Main finding: projection defect looks like a common language across several threads. Added 04-computation/projection_defect_bridge_s12.py and saved 05-knowledge/results/projection_defect_bridge_s12.out. The script audits three forgetful maps: odd cycles to supports, vertex deletion/old-coordinate projection, and odd-n tournament orientations to even graphs. Key data: the two n=8 H=63 classes are exact old-projection kills (delete core, lose 31/31 cycles, H(T-v)=1); THM-025 is a near-kill (one vertex loses 92/94 cycles but leaves two old cycles with alpha=[1,2,1]); Paley T7 and interval T7 have the same regular score sequence and same 36 support shadow but different support_excess/even-graph residues. Extended single-core target search to m<=40; r_core=3 and r_core=10 remain absent, while r=31,42,63 show simple linear count laws after first appearance. Updated INV-191/192/193, OPEN-Q-055, T282/T283, HYP-1760..1763, variable delta_proj, r_core variable file, and reflection 07-reflections/projection-defect-as-common-residue.md. Next: prove r_core != 3,10 algebraically via the dynamic recurrence; scan more n=9/n=10 non-real-root examples for high max deletion-loss fraction; compare projection-defect statistics to beta_3/beta_4 path-homology anomalies; implement cheap projection-defect features in tournament_tda.py.

---

*Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
