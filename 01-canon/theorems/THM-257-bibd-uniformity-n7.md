# THM-257: BIBD Uniformity Theorem for n=7 Regular Tournaments

**Status:** PROVED (exhaustive + representative computation)
**Filed by:** kind-pasteur-2026-03-20-S1

## Statement

Among the 3 isomorphism classes of regular tournaments on 7 vertices:

| Class | c3 | c5_dir | c7_dir | alpha_1 | alpha_2 | S=a1+2a2 | H |
|-------|-----|--------|--------|---------|---------|----------|-----|
| BIBD | 14 | 42 | 24 | 80 | 7 | 94 | 189 |
| Type II | 14 | 28 | 17 | 59 | 14 | 87 | 175 |
| Type III | 14 | 36 | 15 | 65 | 10 | 85 | 171 |

The BIBD class (Paley T_7) is characterized by:

1. **ALL 21 five-vertex sub-tournaments have identical structure:**
   - Same score sequence (1,2,2,2,3)
   - Same H_5 = 13
   - Exactly 2 directed 5-cycles each

2. **No other class has this uniformity property:**
   - Type II: 3 types of sub-tournaments, scores {(1,1,2,3,3), (1,2,2,2,3), (2,2,2,2,2)}
   - Type III: 3 types, scores {(1,1,2,3,3), (1,2,2,2,3), (2,2,2,2,2)}

3. **The uniformity maximizes total directed cycles (alpha_1=80) while
   minimizing disjoint pairs (alpha_2=7, the minimum of the 3 classes).**

## Key Mechanism: Uniformity vs Disjoint Pairs Tradeoff

The BIBD property means each pair of vertices appears in exactly 2 three-cycles.
This creates uniform sub-tournament scores, which maximizes:
- c5_dir = 42 (avg 2.00 per 5-vertex set, ALL sets contribute equally)
- c7_dir = 24 (the full 7-vertex set has the most Hamiltonian cycles)

But it MINIMIZES disjoint 3-cycle pairs because the balanced 3-cycle distribution
means any two 3-cycles tend to share vertices. alpha_2 = 7 = C(7,1) — exactly
one good partition of V\{v} into two 3-cycles per deleted vertex v.

For the vertex-transitive Paley T_7 (Z_7 action), this 1-per-vertex count follows
from the symmetry.

## Average Sub-Tournament H

| Class | Avg H_5 | Explanation |
|-------|---------|-------------|
| BIBD (H=189) | 13.00 | All H_5=13, perfectly uniform |
| Type III (H=171) | 12.43 | Mixed: 3x9, 6x11, 6x13, 6x15 |
| Type II (H=175) | 11.67 | Mixed: 7x9, 7x11, 7x15 |

The average H_5 does NOT monotonically order with H_7! Type III has higher avg H_5
(12.43) than Type II (11.67), yet lower H_7 (171 < 175). This is because alpha_2
matters: Type II has alpha_2=14 vs Type III's 10.

## Formula

At n=7 regular, alpha_3=0 (three disjoint 3-cycles need 9 > 7 vertices), so:

H = 1 + 2*alpha_1 + 4*alpha_2

The weighted sum S = alpha_1 + 2*alpha_2 determines H = 1 + 2*S.

BIBD maximizes S = 94. The BIBD's advantage is:
- +21 in alpha_1 over Type II (80 vs 59)
- -7 in alpha_2 (7 vs 14)
- Net: +21 - 14 = +7 in S (94 vs 87)

The alpha_1 contribution (weight 2) always overwhelms alpha_2 (weight 4) for
large alpha_1 differences.

## Proof

1. Classification: exhaustive enumeration of all 2,097,152 tournaments on 7 vertices
   identifies 2640 regular ones in 3 isomorphism classes (240 + 720 + 1680).
2. Sub-tournament analysis: computed H_5 and score for all C(7,5)=21 five-vertex
   subsets of representatives from each class.
3. Cycle counts: computed directed odd cycles of all lengths using Held-Karp DP.

Scripts: `04-computation/n7_alpha_deep.py`
Results: `05-knowledge/results/n7_alpha_deep.out`

## Connection to THM-255

At n=6 regular, the SC Maximizer has TWO routes to max H:
- Route A (alpha_2-driven): More disjoint pairs, fewer total cycles
- Route B (alpha_1-driven): More total cycles, fewer disjoint pairs

At n=7, only Route B works. The BIBD class IS Route B: maximum alpha_1, minimum alpha_2.
The mechanism shift from n=6 to n=7 is explained by the increasing importance of
higher-length cycles (5-cycles and 7-cycles) as n grows. At n=6, alpha_2 can compensate
because the weight ratio 4/2=2 is significant. At n=7, alpha_1 dominates because the
total cycle count grows faster than the disjoint pair count.

## Related

- THM-255: SC Maximizer Dichotomy at n=6
- THM-256: Paley cycle advantage
- OPEN-Q-016: SC Maximizer conjecture
- THM-027 (S18h): BIBD cycle maximization
