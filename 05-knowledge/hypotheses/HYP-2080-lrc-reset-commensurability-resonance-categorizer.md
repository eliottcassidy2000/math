---
id: HYP-2080
status: SYNTHESIS — dynamical reframe; commensurability=hard regime; resonance-folding categorises (verified correlation)
source: opus-2026-06-02-S563
related:
  - HYP-2059
  - HYP-2075
  - THM-369
---

# HYP-2080: the LRC "total reset" is the hard regime; resonance-folding categorises difficulty

The 13 runners = one point on T^13 along gamma(t)=(v_i t mod 1); loneliness = gamma meets B={x:||x_i||>=1/14}.
- **Reset = orbit closes** at t=1/gcd(v_i). NEVER resets (incommensurate ratios) => gamma DENSE (Weyl) => trivially lonely. RESETS (commensurate=integer) => closed 1-D loop, can miss B => ALL of LRC's difficulty. (Tao's reduction to integers = keep the reset case, discard the easy dense case.)
- **Among resetting (integer) sets, the reset LENGTH is constant (t=1) for primitive sets**, so it does not categorise; the **resonance lattice {m: Σ m_i v_i=0} (orbit FOLDING)** does.
- **Verified (`lrc_reset_orbit_resonance_s563.py`):** more short resonances <=> smaller lonely measure <=> harder. AP: 518 resonances (max), safe-measure 0 (tight), 176 distinct critical moments (MIN — resonances make moments coincide). Random: ~70 resonances, positive measure, ~1050 distinct moments. The hardest set condenses its 4D model the most.
- **The moments are pair-sums** (S557 optimal witness t=m/(v_a+v_b)); resonance = many coincide. Unifies S544/S545 (resonance), S557/S562 (pinch), S553 (tight family) under the time/orbit lens.

**Statement:** categorise speed sets by resonance-folding, not reset length; difficulty is monotone in folding (incommensurate->trivial; generic->easy; AP->tight). LRC@14 = even maximal folding can't fold the closed orbit out of B; the only candidates are the maximally-resonant (AP family).

**Honest:** reframe/organiser, not a new bound; cleanly explains why AP is extremal and incommensurate is free.

**See:** `07-reflections/lrc-time-as-orbit-the-reset-commensurability-categorizer-s563.md`, `04-computation/lrc_reset_orbit_resonance_s563.py` (+.out); HYP-2059 (pinch), HYP-2075 (multi-sieve), THM-369; S544/S545, S553, S557.
