---
id: HYP-2095
status: REDUCTION (rigorous given THM-396/397) + verified split n=4..14; one sharp open lemma
source: opus-2026-06-02-S569
related:
  - THM-396
  - THM-397
  - HYP-2059
  - HYP-2093
  - HYP-2091
---

# HYP-2095: paired-or-anchored splits LRC(n) — the worry-set is the cheap route

**DICHOTOMY (THM-396/397, verified 0/3428):** a small pair (a,b), D=a+b, reduced sum s=D/gcd<=n, is BLOCKED only via a PAIRED blocker (sum-multiple shield D|c) or an ANCHORED blocker (endpoint ||c/D||<1/n); for D<=n the anchor window is empty (=> must be a shield). So a small pair is UNBLOCKED <=> no shield & no anchor, and then its pinch is a witness.
**REDUCTION (rigorous):** LRC(n) <== every measure-ZERO speed set has an unblocked small-reduced-sum pair. (positive measure => lonely; measure-0 + unblocked pinch => witness => M=1/n; so no M<1/n counterexample.)
**SPLIT (verified exhaustive small boxes, n=4,6,8,10,12,14, no exceptions):** every WORRY (measure-0) config has an UNBLOCKED small pair (the HARD set = the CHEAP route); every block-all config is positive-measure (the residual = the EASY spread set).
**WHY (the hint precise):** the worry-set's tight extremiser has straddle pair (a,n-a) summing to n (S557 N2); D=n => anchor window empty, and the config has no multiple of n (C'/S564) => no shield => the pair is UNBLOCKED. The private obstruction IS the straddle pair, which can be neither paired-away (no mult of n) nor anchored (D=n).
**OPEN LEMMA (sharp, verified n=4..14):** every measure-0 set has an unblocked small pair; equiv. if every small pair is shielded-or-anchored then positive measure. Sub-cases: no-small-pair (spread, prove positive measure via S563) and all-shielded (covering-design contradiction).
**HONEST:** dichotomy proved (THM-396/397); reduction rigorous; split verified not proved; remaining lemma open but arithmetic/sharp. Genuine simplification: concentrates LRC(n) into one shield/anchor-cover lemma, hands the worry-set to the cheap route.

**See:** `07-reflections/lrc-private-obstruction-paired-or-anchored-the-proof-split-s569.md`, `04-computation/lrc_private_obstruction_split_s569.py` (+results lrc_private_obstruction_*_s569.out, lrc_proof_split_even_n_s569.out); THM-396/397, HYP-2059 (pinch), HYP-2093 (worry-set target), HYP-2091 (parity), S564 (C').
