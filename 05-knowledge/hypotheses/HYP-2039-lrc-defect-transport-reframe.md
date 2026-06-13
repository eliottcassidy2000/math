---
id: HYP-2039
status: OPEN
source: oracle-2026-06-01-S544
related:
  - HYP-2038
  - HYP-2008
  - THM-382
---

# HYP-2039: LRC is defect TRANSPORT, not existence — global spread guarantees a hole (pigeonhole); the conjecture transports it to the observer

**Reframe.** The n points (observer + runners) have gaps summing to 1, so the
largest gap G(t) >= 1/n at EVERY instant -- local emptiness (a hole of width >=1/n)
is GUARANTEED by global spread (pigeonhole). LRC is therefore not about creating
emptiness but about TRANSPORTING the guaranteed hole to x=0.

**Split (verified, `lrc_wandering_hole_transport_s544.py`):**
- G(t) = largest gap = the HOLE; min_t G >= 1/n for every system (the guarantee).
- O(t) = observer collar = min_i ||v_i t||; max_t O >= 1/n is LRC (the transport;
  a width->=2/n hole sits around the observer).
- Tight AP (n=5,7,14): min_t G = 1/n AND max_t O = 1/n -- the hole is minimum size
  and only just reaches the observer (critical confinement). Generic: both > 1/n.
- The hole = the apex/source-sink (S530) = the fat collar (S541).

**Three equivalent transport reframes:** (1) defect/quasiparticle worldline
theta_hole(t), LRC = its width->=2/n part covers x=0; (2) non-adaptive linear guards
(budget 2-2/n>1 by measure, but constant-speed) cannot cover every instant -- the
obstruction is the LINEARITY of the motion; (3) hole-sweep coverage: by frame-shift
(S541) which point is observer is free, so LRC = the wandering hole sweeps EVERY
point = a coverage/equidistribution of the defect trajectory. The AP gives the least
sweeping (most rigid driving) = the tight case.

**CLAIM / OPEN:** the existence half is a one-line pigeonhole; ALL difficulty is the
transport, concentrated at the AP (G=O=1/n critical confinement). A proof should
target the transport law of a width->=1/n defect under a linear flow:
- (A) state LRC as { x : exists t with a width->=2/n hole containing x } = whole
  circle; is its complement empty unless AP-like?
- (B) a flow/equidistribution lemma: a measure-(2-2/n) union of n-1 linear tubes
  cannot cover a single rational fiber as a SET (only up to the AP boundary)?

**Files:** `04-computation/lrc_wandering_hole_transport_s544.py` (+.out). Reflection:
`07-reflections/lrc-global-spread-guarantees-local-emptiness-the-defect-transport-reframe-s544.md`.
Builds on S530 (apex), S541 (collar/frame-shift), S543 (critical line), S526.
