---
id: HYP-2042
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S545
related:
  - HYP-2039
  - HYP-2011
  - HYP-2003
  - HYP-2040
  - HYP-2041
  - THM-395
---

# HYP-2042: LRC is a cascade of conditional clearances; the cycle-exclusion is the Helly-3 shadow (necessary, not sufficient); collar ladder = 1/(k+1)

**Cascade (verified, `lrc_cascade_cycle_exclusion_s545.py`):** loneliness telescopes
as |SAFE| = prod_i c_i, c_i = |S_{<=i}|/|S_{<i}| (conditional clearance of runner i
given the previous). Generic: every c_i>0 (cascade clears). Tight AP: a ZERO factor
(the last runner trapped by the earlier clearances, at the wall).

**Cycle-exclusion = Helly-3 layer (verified):** the user's hidden fact (an arc X->Y
forbids the 3-cycle Z->X->Y->Z) is a TRIPLE condition. Test 'every triple clears =>
all clear': the worst k-subset collar = 1/(k+1) (the k-runner LRC, sub-AP tight), so
the worst TRIPLE collar = 1/4 for AP n=5,6,7 (>1/n), yet the FULL collar = 1/n. So
the cycle-exclusion is NECESSARY but NOT SUFFICIENT; LRC is irreducibly the whole
cascade (order n), not any fixed Helly-k.

**Relation to THM-395/HYP-2040:** THM-395 proves the exact tournament identity behind
this cycle-exclusion: transitivity is equivalent to zero backward-wedge debt, and the
total wedge mass is `3*c_3(T)`.  This hypothesis records the LRC limitation: that
exact 3-cycle ledger is necessary, but it must be lifted through the whole
conditional-clearance ladder before it can certify LRC.

**NEGATIVE (recorded):** 'runner sub-tournament at the loneliest time t* is a 3-cycle
<=> inside-debt active' is FALSE (104/262 ~ chance) -- t* is a wall with antipodal
TIES, so the geometric 3-cycle there is ill-defined. The cycle-exclusion lives in the
CLEARANCE LADDER (Helly/1/(k+1)), not the geometric tournament at t*.

**CLAIMS (open):**
- (A) worst-k-subset collar = 1/(k+1) exactly (every k-runner sub-system clears at
  >=1/(k+1), AP tight). Prove it; relate to a lower bound c_i >= f((i+1)-runner LRC)
  on the cascade's conditional clearances. Then LRC(n) by induction down the cascade
  needs only the LAST clearance to survive (the AP-critical rung).
- (B) the Helly number of the clearance arcs (runner i's safe set = v_i arcs) is the
  obstruction order; data says n-1 (full), not 3. Quantify it.

**Files:** `04-computation/lrc_cascade_cycle_exclusion_s545.py` (+.out). Reflection:
`07-reflections/lrc-cascade-of-conditional-clearances-cycle-exclusion-is-helly3-s545.md`.
