---
id: THM-713
title: The variance-synchronization decomposition AND a counterexample to unconstrained Var-max — Var(N) = DIAG + COV with DIAG nearly family-independent (spread 0.26) and COV carrying everything (spread 2.19; random families NEGATIVELY associated ≈ −0.25, structured families strongly positive) — but kps-cont.33's proposed reduction "prove Var(N)-max at consec" is FALSE beyond bounded windows: the mod-7-aligned family {1,8,15,22,29,36,43,50,57} has Var = 3.0124 > consec's 2.9670 at k=9 (7-sector resonance hyper-synchronizes emptiness); consec still wins the J-race because it pairs high synchronization WITH the extreme mean (dev² 4.22 vs mod-7's 3.46) — the correct extremal target is J itself (or Var on the mean shell), not unconstrained variance
status: DECOMPOSITION EXACT (rational tables, 18 families k=9); the mod-7 counterexample exact; external search (three-gap literature: Sós–Surányi–Świerczkowski; largest gap = sum of others) confirms no off-the-shelf occupancy-variance extremal exists — the synchronization mechanism (three-gap bunching for consec; 7-sector resonance for mod-7) is the project's to prove. COURSE CORRECTION for the fleet: aim proofs at J directly with the mean term as co-conspirator (not contaminant); the mod-7 class is the second synchronization pole and must be in any extremal argument's hypothesis or conclusion.
source: mac-mini-2026-07-09-S65 (cont.37, 2026-07-11)
depends_on:
  - THM-711/712 (the base functionals), kps cont.33 (the decomposition frame this corrects)
related:
  - three-gap theorem (external; the bunching tool), THM-709 (aliasing-isolation kinship)
---
# THM-713 — synchronization carries the variance; Var-max at consec is false globally
Var(N) = Σ_s A_s(1−A_s) + Σ_{s≠s′}[A_{ss′} − A_sA_{s′}]. Exact k=9 tables: DIAG ∈ [1.06, 1.32]
(all families); COV: random ≈ −0.25 (negative association), near-AP +1.22, consec +1.83,
mod-7 +1.76 with larger DIAG ⟹ Var(mod7) = 3.0124 > Var(consec) = 2.9670. J stays consec-min
(5.06 < 5.77) via the mean term. Files: 04-computation/lrc14_var_sync_macmini_S65cont37.py (+ .out).
Sources: [Three-gap theorem (Wikipedia)](https://en.wikipedia.org/wiki/Three-gap_theorem),
[The Three Gap Theorem, JAusMS](https://www.cambridge.org/core/journals/journal-of-the-australian-mathematical-society/article/three-gap-theorem-steinhaus-conjecture/EA75E140919DEA9A55FEFD01EB2F677F).
