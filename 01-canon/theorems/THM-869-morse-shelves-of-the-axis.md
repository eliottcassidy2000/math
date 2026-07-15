---
id: THM-869
title: THE MORSE SHELVES OF THE AXIS — a tournament above the floor has no descending flip iff n is odd, d ∈ {0,±2}, and every d=−2 vertex beats every d=+2 vertex (upset saturation); even n has NO local minima above the floor; shelves sit at x = 8a and the a-th shelf opens at r ≥ a + ⌈(a−1)/2⌉ + 1-ish (a=2 first realized at n=9, constructed) — the descent obstruction is odd-n floor-parity geometry, NOT a quintic threshold
status: PROVED (characterization, both directions; even-n emptiness; overload argument for M ≥ 4) + VERIFIED exhaustively n = 5, 6, 7 (0 mismatches; stuck = 40, 0, 16800) + a=2 shelf at n=9 CONSTRUCTED and machine-checked
source: mac-mini-2026-07-15-S110; refines kind-pasteur cont.22's false-peaks census (their 0/40/368 at n=4/5/6 counts UNIT-step stuckness; with all flip moves the n=6 count is 0 — the obstruction is odd-n, not n≥5)
depends_on:
  - THM-855 F3 (Δx = 4(d_l − d_w) + 8; descending flip ⟺ an arc with d_w − d_l ≥ 4)
  - THM-866 (the upward walk; this file is the downward Morse picture)
related: [kps cont.22 HYP-6945 false peaks, THM-466(iv) (the (1,2,2,2,3) fiber), THM-854(I)]
script: 04-computation/ryser_fiber_shelves_macmini_S110.py -> 05-knowledge/results/ryser_fiber_shelves_macmini_S110.out (+ stuck censuses in truncation_ladder_thm868_macmini_S110.out §7)
---

# THM-869 — the Morse shelves of the axis

Call T **stuck** if x(T) > x_floor and no single arc flip decreases x (an x-local-minimum
above the floor). By F3, a descending flip is an arc w → l with d_w − d_l ≥ 4 (d's have
fixed parity, so "> 2" means "≥ 4").

**Theorem.**
1. **(even n)** No tournament is stuck: descent to the near-regular floor is
   single-flip-unobstructed everywhere.
2. **(odd n, characterization)** T is stuck ⟺ every d_v ∈ {0, +2, −2} (not all 0) and
   **every d = −2 vertex beats every d = +2 vertex** (upset saturation). Hence stuck
   tournaments sit exactly on the levels x = 8a, a = #{d = +2} = #{d = −2} ≥ 1.
3. **(shelf tower)** Shelf a requires r = (n−1)/2 ≥ a + ⌈(a−1)/2⌉ + [a ≥ 2]-slack from
   the score arithmetic (some −2 vertex must carry a cross-wins plus its share of the
   internal round-robin inside its own r−1 score); a = 1 exists for all odd n ≥ 5
   (counts 40 at n = 5, 16800 at n = 7 — the upset-oriented sub-fibers of
   (r−1, r^{n−2}, r+1)); **a = 2 is first realizable at n = 9** (explicit construction,
   machine-verified stuck at x = 16).

**Proof.**
*(≤ 2 forces the shape.)* Suppose T stuck, let M = max d ≥ 2 at vertex v. Every
out-neighbor of v has d ≥ M − 2, and |out(v)| = r + M/2 (odd n; the even case is
identical with half-integer shifts). If M ≥ 4: v plus its out-set give ≥ r + 3 vertices
of deviation ≥ M − 2 ≥ 2, so Σ_{d>0} ≥ 2(r+2) + (M−2) ≥ 2r + 6, forcing some
d ≤ −4; the symmetric argument at the minimum gives ≥ r + 3 vertices of deviation
≤ −2; the two sets are disjoint, so n ≥ 2r + 6 = n + 5 — contradiction. Hence M = 2,
and symmetrically min d = −2: d ∈ {0, ±2}. A +2 → −2 arc has margin 4 and would
descend, so all a² cross-pairs must be upsets (−2 beats +2). Conversely, in an
upset-saturated {0,±2} tournament every arc has d_w − d_l ≤ 2: no descending flip. ∎
*(Even n.)* The same overload argument kills M ≥ 3, leaving all d = ±1 — which IS the
floor. ∎
*(Shelf arithmetic.)* The a lows beat all a highs and play C(a,2) internal games, so some
low has score ≥ a + ⌈(a−1)/2⌉, i.e. r − 1 ≥ a + ⌈(a−1)/2⌉. For a = 2 this needs
r ≥ 4 (n ≥ 9), matching the exhaustive absence at n = 7 and the explicit n = 9
construction (scores (3,3,4,4,4,4,4,5,5), lows beat highs, x = 16, stuck — in the .out). ∎

## Corrections and readings

- **The unit-step/full-move distinction matters.** kps cont.22's counts (0/40/368 at
  n = 4/5/6) are tournaments with no Δx = −8 move; with ALL flip moves allowed the n = 6
  count is **0** (a Δx = −16 move always exists there). The Morse asymmetry of the axis
  is: **upward to the transitive always flows (THM-866); downward to the floor is
  obstructed exactly at odd n by upset-saturated shelves**. The threshold is floor
  parity (regular floors exist only at odd n), not the quintic n = 5.
- What remains genuinely n = 5: THM-466(iv)'s digit-1 fiber escape (see the S110
  monodromy probe: 0 parity-flipping Ryser reversals at n = 4, 240 at n = 5, all inside
  the (1,2,2,2,3) fibers — which are exactly the score class hosting the a = 1 shelf).
- The stuck set is a union of sub-fibers: at a = 1 it is "the gap-2 pair plays its
  upset", fraction 40/280 = 1/7 (n=5), 16800/72240 = 10/43 (n=7) of the near-regular
  fiber — counts data for the shelf-census sequence (backlogged).
