# 7 and 21 as conflict graphs: one defect, a shield level, and the join (11, 11)

## Inherited (cited, not re-derived)

- **Half-index chart.** h(c) = 1 + 2c turns a conflict-graph join into addition and the ordinal sum
  into c ∘ c′ = c + c′ + 2cc′. The triangle shield is c ↦ 3c + 1. Then 7 = h(3), 21 = h(10) = h(3 ∘ 1),
  and 63 = h(31) is the next shield image, which is a count.
  - Sources: wave 43 §5 of
    [chanpak_waves_20261113](../chanpak_waves_20261113/README.md#wave-43-the-carry-odometer-tested-an-exact-triangle-shield-two-apex-gluing-rigidity-and-7-21);
    fresh_unit_prime_gluing_20261116/SPINE_AND_HALF_INDEX.md (the (11, 11) row is 5 + 5 = 10);
    Session 210 and THM-4940 (the same coordinate as the α-sum).
- **Exclusion of 7 and 21.** THM-1370 proves it via the crevasses of the strong floor, and refutes the
  "7·3^k tower" at k = 2.
- **Carry obstruction.** THM-4940 (4): the 21-hole is a carry obstruction, since 1 + 10b is a value
  only in bases 4 and 8. THM-4940 (2) is the block lemma; Session 16's set A is the intercyclic spectrum.
- **Two specific forms.** THM-202: P_4 is never Ω(T). THM-343 / THM-029: three pairwise-intersecting
  odd cycles force a fourth.

## New: every conflict graph of value 7 or 21 (PROVED, finite enumeration)

I(G, 2) = 1 + 2|V| + 4α_2 + 8α_3 + ⋯ ≤ 21 forces |V| ≤ 10 and at most four non-edges.
`omega_forms.py` enumerates every complement with at most four edges on at most ten vertices. For each
resulting G it lists all splits into join factors (unions of co-components, i.e. connected components
of the complement). Output: `omega_forms.out`.

| value | graph G (α_1, α_2) | structure | join splits (I(G_1), I(G_2)) |
|---:|---|---|---|
| 7 | K_3 (3, 0) | join | (3, 5) |
| 21 | K_10 = **h(10)** (10, 0) | join | (3,19), (5,17), (7,15), (9,13), **(11,11)** |
| 21 | K_8 − e (8, 1) | join | (3,19), (5,17), (7,15), (9,13), **(11,11)** |
| 21 | K_6 minus a matching of size 2 (6, 2) | join | (3,19), (5,17), (9,13), **(11,11)** |
| 21 | K_6 minus a path of length 2 = K_3 ∨ (K_2 ⊔ K_1) (6, 2) | join | (3,19), (5,17), (7,15) |
| 21 | K_3 ⊔ K_1 (4, 3) | disjoint union: the shield of 7 | none |
| 21 | P_4 (4, 3) | neither (the unique non-cograph) | none |

So 21 is the place where four mechanisms coincide. Only P_4 is a genuinely new obstruction; the other
five forms are excluded through 7 or through a failed join:
- **the triangle shield of the 7-defect**: one form;
- **the join of two 11-forms**: three forms. Since 11 = I(K_5) = I(P_3), the joins K_5 ∨ K_5 = K_10,
  K_5 ∨ P_3 = K_8 − e and P_3 ∨ P_3 are precisely the (10, 0), (8, 1) and matched (6, 2) forms;
- **the +2 dial on a 19 and the double +2 dial on a 17**: (3, 19) and (5, 17) split **every** join form;
- **P_4**, excluded separately by THM-202.

**The two notes on 7.** 7 has a single form, K_3 = K_1 ∨ K_2: the join of the two smallest prime counts.
And K_3 is never a connected component of any Ω(T).
- *Proof.* Three pairwise-meeting odd cycles span a strongly connected set U.
  - |U| = 3 gives one cycle, and a 4-vertex tournament has at most two directed triangles. So |U| ≥ 5.
  - Moon's theorem gives U a further odd cycle, of length 3 or 5.
  - That cycle lies in U, so it meets one of the three.
  - So the three cycles' component of Ω has at least four vertices. ∎ (This is the form of THM-343 used
    here.)

## The two dials and {7, 21} in one algebra

Read the moves of the radix_three_shields packet in conflict-graph terms:

| move | conflict-graph operation | count |
|---|---|---|
| +2 dial | join with K_1 (a cycle meeting every cycle) | I(G ∨ K_1) = I(G) + 2 |
| triangle shield | disjoint K_1 | I(G ⊔ K_1) = 3I(G) |
| +6 dial | a +2 dial applied beneath one shield | I((G ∨ K_1) ⊔ K_1) = 3I(G) + 6 |

So the 24-vertex bank's two dials are **one** conflict-graph move at two shield levels. In its exact form
(G ∨ K_j) ⊔ K_1 ∨ K_i, with i top-level universal cycles and j cycles missing only the shield cycle, the
count is 3I(G) + 6j + 2i.

**The rectangle over the empty base.** Here Ω_ij = (K_j ⊔ K_1) ∨ K_i, with values 3 + 2i + 6j for
0 ≤ i, j ≤ 3: thirteen consecutive odds, 3 to 27. Its forbidden cells are exactly the cells where K_3
appears:
- (i, j) = (2, 0): K_3, value 7;
- (0, 3): K_3 ⊔ K_1, value 21;
- (3, 2): K_3 ∨ (K_2 ⊔ K_1), value 21.

So {7, 21} is the same defect, K_3, read at the top level and under the shield. Correspondingly
{+2, +6} is the same move read at the two levels.

**Why exactly once (ANALOGY; the pieces are PROVED or cited).**
- **The literal shield.** The triangle shield multiplies a window exactly by 3 once per apex or wrapper
  (wave 43 §3; radix_three_shields_20261116/OBSTRUCTIONS.md, the mod-5 congruence). A second exact ×3
  needs a second, disjoint channel. That channel multiplies steps rather than adding them.
- **The defect.** K_3 ⊔ kK_1 is excluded at the conflict-graph level for every k (the lemma above). Its
  value 7·3^k is excluded only for k ≤ 1: 63 is realized by other conflict graphs at order 8 (THM-1370).
- **A shared reading.** Exact ×3 travels only through the disjoint channel I(G ⊔ K_1) = 3I(G), and that
  channel transports windows and defects alike. The strongly connected ("live") channels perturb windows
  and rescue defects, and at the second level they take over.

## The leg 11

The target map (TARGET_MAP.md) shows that 11 never sits as an induced (J) leg below 125, with any
partner. 13 sits only as 13 + 13 − 1 = 25. So the forbidden join (11, 11) → 21 is the tip of a wider
rigidity of the two smallest binary strong singletons (THM-4940 (8)).
- *Status.* FINITE-EXACT below 125.
- *Question (Q-11).* Does 11 sit as an induced (J) leg at any count?
