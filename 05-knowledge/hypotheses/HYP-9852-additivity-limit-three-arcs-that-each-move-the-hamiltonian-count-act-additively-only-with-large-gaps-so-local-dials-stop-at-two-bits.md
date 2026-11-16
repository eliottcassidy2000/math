---
id: HYP-9852
title: "Additivity limit for arc flips: if three arcs of a tournament each change H when flipped and their eight joint orientations act additively (every interaction coefficient zero), then at most one of the three half-gaps is below 3. Hence no additive 3-bit dial (half-gaps 1, 2, 4) exists, local dials stop at two bits, and exact placement of a prescribed count needs non-local edits. (The stronger form 'no +-2 arc at all' is REFUTED at order 8 by 16 cubes of half-gaps (1, 11, 34).)"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-16-dials, wave 41 of chanpak_waves_20261113.
  FINITE-EXACT, every 3-subset of arcs of every class representative:
  - order 6: 52 additive cubes, all with an inert arc;
  - order 7: 688 additive cubes; those with three moving arcs have half-gaps (3,3,24), (3,8,24), (15,15,15) only;
  - order 8: 11,340 additive cubes; 192 with three moving arcs, of half-gap types (1,11,34): 16, (3,3,24): 32,
    (3,8,24): 32, (5,5,39): 24, (6,23,23): 32, (6,37,37): 32, (15,15,15): 24. The second-smallest half-gap is always
    at least 3. The type (1,11,34) REFUTES the first draft's stronger form, which said no +-2 arc at all.
  Additive 2-cubes are plentiful, including the complete 2-banks with half-gaps (1,2), so the limit is specific to
  three moving arcs. Not proved.
session: mac-mini-2026-11-16-dials
date: 2026-11-16
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 41)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w41_arc_affine.c (usage: arc_affine n k < class masks)
related:
  - HYP-9842 (no complete 3-arc bank, through order 8), THM-5302, HYP-9762 (no affine counter), unit_shields_wheel_followup_20261116
---

# HYP-9852 - the additivity limit and local dials

**Setting.** Fix a tournament T and arcs e₁, e₂, e₃, and write s_i = 1 if e_i is flipped. The count is multilinear in
the s_i:

    H(s) = c₀ + Σ c_i s_i + Σ c_ij s_i s_j + c₁₂₃ s₁ s₂ s₃.

- The cube is *additive* if c₁₂ = c₁₃ = c₂₃ = c₁₂₃ = 0.
- Arc e_i *moves* if c_i ≠ 0, and its half-gap is |c_i|/2. Every c_i is even, since H is odd.

**Conjecture.** In an additive cube with three moving arcs, at most one half-gap is below 3 (the second-smallest is
at least 3).

The first draft said every half-gap is at least 3. Order 8 refutes that: 16 cubes of half-gaps (1, 11, 34) combine a
±2 arc additively with two large moving arcs.

**Why it matters.**
- An additive cube with half-gaps {1, 2, 4} would be a complete 3-arc bank: a 3-bit dial whose eight settings give
  eight consecutive odd counts.
- More generally, a local k-bit dial with independent digits needs k additive moving arcs with half-gaps
  1, 2, …, 2^(k−1).
- The conjecture forbids this for k ≥ 3, since a dial needs two half-gaps (1 and 2) below 3. Two-bit dials do exist (counts 9, 11, 13, 15 at orders 5 and 6).
- Exact placement therefore cannot be finished by local multi-arc dials. It needs a non-local edit, such as the
  two-apex gap patch whose one arc carries any programmed even gap (unit_shields_wheel_followup_20261116,
  apex/PROGRAM.md).
- The patch's absolute count is not affine in its address (wave 41 §4). Affine structure appears only in gaps, which
  matches HYP-9762.

**Data.**

| order | 3-arc cubes | additive | additive with three moving arcs (half-gaps: count) |
|---:|---:|---:|---|
| 6 | 25,480 | 52 | none |
| 7 | 606,480 | 688 | (3,3,24): 16, (3,8,24): 16, (15,15,15): 12 |
| 8 | 22,538,880 | 11,340 | (1,11,34): 16, (3,3,24): 32, (3,8,24): 32, (5,5,39): 24, (6,23,23): 32, (6,37,37): 32, (15,15,15): 24 |

**What the data settle.** No additive cube with three moving arcs has two half-gaps below 3. This is exact through
order 8, so no additive 3-bit dial exists there.

**The symmetric additive cube lives on the Paley tournament.** The (15, 15, 15) cubes of order 7 are based at the
Paley tournament P₇ (H = 189), the unique maximiser at order 7. Three of its arcs each remove exactly 30 Hamiltonian
paths, independently, so the eight counts are 189, 159, 159, 159, 129, 129, 129, 99. Seen from its other corners
the same cube appears with base counts 159, 129 and 99. The seven triples are the reflection matchings of Z/7: for
each vertex v, the arcs {v − d, v + d}, d = 1, 2, 3. The other additive types, (3, 3, 24) and (3, 8, 24), sit on
counts 83–159.

**Falsifier.** One additive cube with three moving arcs and two half-gaps below 3. The first open order is 9.
