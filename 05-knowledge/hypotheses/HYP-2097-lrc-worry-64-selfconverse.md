---
id: HYP-2097
status: SYNTHESIS — merges oracle 64-class even-ladder reduction with the S569 paired/anchored split; verified census; one open lift lemma
source: opus-2026-06-03-S570
related:
  - HYP-2095
  - HYP-2094
  - HYP-2086
  - HYP-2093
---

# HYP-2097: the n=14 worry-set = 64 self-converse round classes, all lonely via an unblocked small pair

- **64 count grounded:** self-converse round tournaments on m (odd) = 2^((m-1)/2); verified 2,4,8,16 for m=3,5,7,9 (circular-points generator) ⇒ 64 at m=13 (n=14). The worry-set ⊆ these 64 (oracle-S576o), vs A000568(13)≈4.85e13.
- **Transversal census (all 8191 gcd-1 transversals mod 27, `lrc_64_worry_classes_s570.py`):** LONELY 8191/8191 (0 counterexamples); floor-tight (worry proper) = exactly the AP; EVERY transversal has an UNBLOCKED small pair (S569 mechanism uniform over the whole 64-class container).
- **MERGE:** worry-set ⊆ 64 self-converse round classes (oracle) → each realized by transversals → each has an unblocked small pair (S569) ⟹ lonely. One uniform mechanism (the unblocked straddle pair, neither paired-away nor anchorable since sum=n) handles all 64; no bespoke per-class argument.
- **Flip-lattice:** the 64 = flip-sets F (2^6); F=∅ (AP) is the UNIQUE floor-tight class; all 63 non-empty flips LOOSEN (M>1/14). The dual-Burnside fix-side as a Boolean lattice with the AP at the tight bottom.
- **OPEN LIFT LEMMA:** prove every speed set realizing one of the 64 classes (incl. unbounded/non-transversal composite-27 cousins like V*) has an unblocked small pair (= measure-0 ⟹ unblocked small pair, S569, localised to the 64 container).
- HONEST: 64 count + census + flip-lattice verified; the merge is structural; not a proof — the lift lemma is open.

**See:** `07-reflections/lrc-worry-framework-on-the-64-self-converse-classes-s570.md`, `04-computation/lrc_64_worry_classes_s570.py` (+.out), `lrc_selfconverse_round_count_s570.py` (+.out); HYP-2095 (paired/anchored), HYP-2094 (oracle 64), HYP-2086 (dual Burnside), HYP-2093 (floor-tight target).
