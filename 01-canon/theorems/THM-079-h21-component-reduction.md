# THM-079: H=21 Component Reduction

**Status:** PROVED for all n (structural + computational base case)
**Author:** opus-2026-03-07-S39/S41/S42, kind-pasteur-2026-03-07-S33
**Date:** 2026-03-07
**Dependencies:** THM-029 (H=7 impossibility)

## Statement

### Part A (PROVED): Disconnected Omega Case
If Omega(T) is disconnected, then H(T) ≠ 21.

### Part B (PROVED): P_4 Component Case
Omega(T) cannot have a connected component isomorphic to P_4 (path on 4 vertices).

### Corollary
H(T) = 21 requires connected Omega(T) with I(Omega(T), 2) = 21 and |V(Omega)| >= 6.

## Proof of Part A

By OCF multiplicativity for disjoint union:
I(G1 ∪ G2, 2) = I(G1, 2) · I(G2, 2)

For H = 21 with disconnected Omega: need I(C_i, 2) for components with product = 21.
Since H is always odd, each factor is odd.
Factorizations: 21 = 3·7 = 1·21 = 1·3·7 = 1·1·21 etc.

The only non-trivial factorization is 3·7.
I(component, 2) = 7 requires the component G to be K_3 (3 pairwise conflicting cycles).
By THM-029: 3 pairwise-conflicting cycles in Omega(T) ALWAYS generate a 4th cycle
(either via common vertex forcing 5-cycle, or via tournament arc constraints).
So |V(G)| ≥ 4 and I(G,2) ≥ 1+2·4 = 9 > 7. Contradiction.

Therefore I(component, 2) = 7 is impossible for any component of any Omega(T).
So H = 21 via disconnected Omega is impossible. QED.

## Proof of Part B

I(P_4, 2) = 21 (direct computation: 1 + 4·2 + 3·4 = 21).
A P_4 component of Omega means 4 cycles C1-C2-C3-C4 with adjacent pairs sharing vertices
and non-adjacent pairs vertex-disjoint.

Consider the middle pair C2, C3. They share a vertex (adjacent in P_4).
Their vertex union V(C2) ∪ V(C3) has 5 vertices (3+3-1 from the shared vertex).

**Key lemma:** Two directed 3-cycles sharing a vertex on 5 tournament vertices
always force at least one additional 3-cycle (verified exhaustively at n=5:
when C1={0,1,2} and C2={2,3,4} are both cycles, minimum t_3 = 3, never 2).

The forced additional cycle C* lies on a triple from V(C2) ∪ V(C3).
Since C* shares vertices with C2 or C3, it is adjacent to C2 or C3 in Omega.
But C* ≠ C1 (C1 has vertices outside V(C2)∪V(C3)) and C* ≠ C4 (same reason).
So C* is a 5th vertex in Omega adjacent to the middle of the path, and Omega ≠ P_4. QED.

## Graph Classification for I(G,2) = 21

Connected graphs with I(G,2) = 21 for |V| ≤ 6:
- |V| = 4: P_4 (path on 4 vertices). I = 1+8+12 = 21. [RULED OUT by Part B]
- |V| = 5: No connected graphs with I = 21.
- |V| = 6: K_6 minus 2 edges (2 non-edges). I = 1+12+8 = 21. [OPEN]

## Part C (PROVED): alpha_3 Elimination

All decompositions with alpha_3 >= 1 are impossible.
alpha_3 >= 1 means 3 mutually vertex-disjoint cycles, which forces alpha_2 >= 3
(the 3 pairwise-disjoint pairs). Then 4*alpha_2 >= 12 and 8*alpha_3 >= 8,
giving 4*alpha_2 + 8*alpha_3 >= 20, leaving at most 2*alpha_1 = 0.
But alpha_1 >= 3 (the 3 cycles themselves), so 2*alpha_1 >= 6.
Total >= 6+12+8 = 26 > 20. Contradiction.

## Remaining Viable Decompositions

Only 4 decompositions of alpha_1 + 2*alpha_2 = 10 remain:
- (10, 0): 10 pairwise-conflicting cycles. Needs alpha_1=10, all pairs sharing vertex.
- (8, 1): 8 cycles with exactly 1 disjoint pair.
- (6, 2): 6 cycles with exactly 2 disjoint pairs.
- (4, 3): 4 cycles with exactly 3 disjoint pairs.

## Part D (PROVED): (4,3) Decomposition Fully Eliminated

For (4,3): the only connected graphs on 4 vertices with 3 edges are trees: P_4 and K_{1,3}.
- P_4 ruled out by Part B.
- K_{1,3} ruled out: 3 mutually disjoint leaves force alpha_3>=1, contradicting Part C.
- The only other graph with 3 edges on 4 vertices is C_3+isolated, which is disconnected.
Therefore (4,3) is fully eliminated.

## Part E: Directed Cycle Counting Correction

**Important:** Omega(T) vertices are DIRECTED odd cycles, not vertex sets. A vertex set
of size 5 can support 0, 1, 2, or 3 distinct directed 5-cycles, each a separate vertex
in Omega. Multiple directed cycles on the same vertex set are always adjacent (sharing
all vertices), forming a clique. This means alpha_1 (total cycles) >= vertex-set count,
making I(Omega,2) LARGER. The decomposition analysis is unchanged but harder to achieve.

## Part F: Tournament Forcing of Extra Cycles

**Computational finding (opus-S40):** Even when 6 three-cycle triples have the K_6-2e
conflict pattern and are realizable as directed 3-cycles in a tournament, the tournament
structure ALWAYS forces additional directed 5-cycles. These extra cycles expand Omega
beyond the 6-vertex K_6-2e skeleton.

At n=6 (exhaustive, 2160 tournaments with K_6-2e 3-cycle pattern):
- Every realization has 4, 6, or 8 additional directed 5-cycles
- I(Omega(T), 2) is always 29, 33, or 37 (minimum 29, never 21)
- The sharing lemma alone is insufficient (2790/5040 K_6-2e triple arrangements
  pass the sharing lemma test), but tournament forcing of 5-cycles provides
  the actual obstruction

## Part G (PROVED, computational): H=21 impossible at n<=8

**Exhaustive verification (opus-S40, S41):** H=21 never occurs among:
- All 2,097,152 tournaments on 7 vertices (opus-S40)
- All 268,435,456 tournaments on 8 vertices (opus-S41, h21_exhaustive_n8_v3.c)

The n=8 exhaustive check used three pre-filters:
1. Skip source/sink tournaments (31,719,424 skipped) — induction from n=7
2. Skip t3 > 10 (218,589,056 skipped) — alpha_1 > 10 implies H > 21
3. Skip tournaments with vertex not in any 3-cycle (31,360 skipped) — Key Lemma (Part J)
Only 18,095,616 tournaments (6.7%) required Held-Karp computation.

Combined with n<=6 (also exhaustive), H=21 is impossible for all n <= 8.

Additionally, random sampling at n=9 (2,000,000 tournaments, opus-S41) found H=21: 0.

At n=7, the H-spectrum has exactly two "permanent gap" values below 30: H=7 and H=21.
H=35, H=39 are gaps at n=6 but achieved at n=7. H=63 is a gap at n=7 but achieved at n=8.

## Part H: i_2 Jump Pattern (The Core Obstruction)

**Key discovery (opus-S40):** The achievable (alpha_1, i_2) pairs in tournaments
systematically SKIP the exact values needed for H=21.

For H=21: need alpha_1 + 2*i_2 = 10 (with i_3=0). The four decompositions require:

| Decomposition | alpha_1 | needed i_2 | n=6 achievable i_2 | n=7 achievable i_2 | n=8 achievable i_2 |
|---------------|---------|------------|--------------------|--------------------|---------------------|
| (4, 3)        | 4       | 3          | {0}                | {0}                | {0, 4}              |
| (6, 2)        | 6       | 2          | {0, 1}             | {0, 1}             | {0, 1, 5}           |
| (8, 1)        | 8       | 1          | {0}                | {0}                | {0, 3, 7}           |
| (10, 0)       | 10      | 0          | {2}                | {2}                | {2}                 |

In EVERY case, the needed i_2 value is NOT in the achievable set. The i_2 values
"jump" between discrete values, never hitting the H=21 target. This is verified:
- n=5,6: exhaustive
- n=7: exhaustive (2,097,152 tournaments)
- n=8: EXHAUSTIVE (268,435,456 tournaments, opus-S41)
- n=9: 2,000,000 random samples (opus-S41)

### Structural explanations for each blocking:

**(4,3) blocked:** The two graph structures with i_2=3 on 4 vertices are
complement(P_4) and complement(K_{1,3}). P_4 eliminated by Part B. K_{1,3}
has 3 pairwise-sharing leaf cycles; by THM-029 these force a 4th cycle,
making alpha_1 >= 5. Contradiction. (Part D)

**(6,2) blocked:** K_6-2e structure. Tournament forcing always creates
additional 5-cycles. At n=6 exhaustive: I(Omega,2) >= 29 for all K_6-2e
tournaments. At n=7 exhaustive: K_6-2e tournaments have H >= 29
(1,597,968 K_6-2e tournaments found, min H = 29). (Part F)

**(8,1) blocked:** K_8-e structure. 8 cycles with exactly 1 disjoint pair.
Computationally: at n=7, alpha_1=8 always gives i_2=0 (H=17). At n=8
(detailed sampling, opus-S41), alpha_1=8 gives i_2 in {0, 3, 7} only, never 1.
- i_2=0: composition (4,4,0), with source/sink
- i_2=3: composition (5,3,0), with source/sink, all disjoint pairs are (3,3) type
- i_2=7: composition (5,3,0) or (6,2,0), WITHOUT source/sink, star K_{1,7} pattern
  in disjointness graph (one 3-cycle disjoint from ALL other 7 cycles)
Structural proof OPEN but star pattern is key insight.

**(10,0) blocked:** K_10 structure. 10 pairwise-sharing cycles. At n=6,7,8:
alpha_1=10 always gives i_2=2 (H=29). ALL alpha_1=10 tournaments have a
source or sink vertex (verified in 5M random n=8 samples, opus-S41 — every
alpha_1=10 score sequence contains 0 or 7). So they reduce to n=7 by Part K.
Structural proof OPEN for general n.

## Part J (PROVED): Key Lemma — No 3-cycle implies no odd cycle

**Lemma:** If a vertex v in a tournament T is not contained in any directed 3-cycle,
then v is not contained in any directed cycle of any length.

**Proof:** Suppose v is in no 3-cycle. Partition V \ {v} into:
- N+(v) = {w : v -> w} (out-neighborhood)
- N-(v) = {w : w -> v} (in-neighborhood)

Claim: ALL arcs between N-(v) and N+(v) go from N-(v) to N+(v).

Proof of claim: Suppose for contradiction that some u in N+(v) beats some w in N-(v),
i.e., u -> w. Then w -> v -> u -> w is a directed 3-cycle containing v. Contradiction.

So the tournament has a layered structure: N-(v) -> {v} -> N+(v), and also N-(v) -> N+(v).
Every directed path from v goes into N+(v) and can never return to N-(v) (since all
cross-arcs go N-(v) -> N+(v)). So no directed path from v can return to v. Therefore
v is in no directed cycle of any length. QED.

**Consequence for H=21 induction:** If T has n >= 8 vertices and some vertex v is in
no 3-cycle, then Omega(T) = Omega(T - v) (deleting v doesn't remove any odd cycle).
So H(T) = I(Omega(T), 2) = I(Omega(T-v), 2) = H(T-v). By induction on n
(base case n <= 8 from Part G), H(T-v) != 21, so H(T) != 21.

## Part K: Source/Sink Induction

If T has a source (score 0) or sink (score n-1), that vertex is in no directed cycle
(trivially — it has no in-arcs or no out-arcs). So by Part J's argument, removing it
preserves Omega(T) and H(T) = H(T - v) != 21 by induction.

This is a special case of Part J but was discovered first and is worth stating separately.

## Part L (PROVED, computational): Cycle-rich min-H bound at n=8

**Exhaustive finding (opus-S41):** Among ALL 18,095,616 cycle-rich n=8 tournaments
(no source/sink, every vertex in 3-cycle, t3 <= 10), the minimum H is **25**.

H distribution for small values: H=25 (40,320), H=27 (53,760), H=33 (26,880), ...

Since 25 > 21, this provides an independent proof that H=21 is impossible for
the "hard case" (no source/sink, all vertices cyclic) at n=8. Combined with
Parts J and K for induction, this strengthens the base case.

## Part I: H-Spectrum Gap Analysis

H=7 and H=21 are the ONLY permanent gaps in the achievable H-spectrum:
- H=7: permanent gap at n=5,6,7,8 (proved by THM-029)
- H=21: permanent gap at n=6,7,8 (computational; this theorem)
- All other odd values <= 200 are achieved at some n <= 8

At n=8 (500k sample): the only odd values NOT seen in [1..200] are 7 and 21.

## Part M (SKETCH): (10,0) Star Capacity Proof

**Claim:** alpha_1=10 with i_2=0 is impossible at n=7 (and hence all n by induction).

**Proof sketch (opus-S41 agent):**
At n=7, vertex v has score s_v, creating at most s_v*(7-1-s_v) mixed pairs. Max is 9
(when s_v=3). So a star of 3-cycles through v holds at most 9 cycles.

To reach alpha_1=10 with all pairwise-sharing, need a 10th cycle not through v.
A non-star 3-cycle on {w_i,w_j,w_k} is disjoint from star cycles {v,w_a,w_b}
where both w_a,w_b are in the complement {w_1,...,w_6}\{w_i,w_j,w_k}. There are
C(3,2)=3 such pairs. To maintain i_2=0, ALL 3 must be transitive (not 3-cycles),
limiting the star to at most 6 cycles. Total alpha_1 <= 7 < 10. Contradiction.

For 5-cycles and 7-cycles: at n=7, pigeonhole gives 5+3>7, so any 5-cycle shares
with any 3-cycle. The only disjoint pairs are between pairs of 3-cycles, and the
capacity argument above still applies.

## Part N (PROVED at n=8): (8,1) Impossibility via Cascade Forcing

**Claim:** alpha_1=8 with i_2=1 is impossible at n=8.

**Proof (opus-S41 agent, exhaustive case analysis):**

At n=8, a disjoint pair of odd cycles can only be:
- (a) Two 3-cycles (using 6 of 8 vertices), or
- (b) A 3-cycle + 5-cycle (using all 8 vertices).

**Case (b) impossible:** The 5 vertices of the 5-cycle must contain no 3-cycle
(else that 3-cycle is disjoint from C_a, giving i_2 >= 2). But a tournament with
no 3-cycle is transitive and has no cycles at all — contradicts having a 5-cycle.

**Case (a) proved by sub-tournament analysis:**
Let C_a={a,b,c}, C_b={d,e,f}, bridge vertices {g,h}.
- Sub-tournaments on {a,b,c,g,h} and {d,e,f,g,h} each need exactly 1 odd cycle.
- At n=5: t3=1 forces t5=0 (verified exhaustively). Score sequences: {(0,1,3,3,3),
  (0,2,2,2,4), (1,1,1,3,4)}.
- 18 compatible (g,h)-score combinations: 14 produce source/sink (reduces to n=6
  where alpha_1=8 always gives i_2=0, contradicting the disjoint pair). 4 produce
  no source/sink but exhaustive enumeration of all 512 cross-edge assignments shows
  alpha_1 >= 20, far exceeding 8. Contradiction in all cases.

**Cascade forcing mechanism:** If any disjoint pair exists among 8 cycles, bridge
vertices either (a) create a source/sink collapsing cycles to 6 vertices where
i_2=0 is forced, or (b) generate so many extra cycles that alpha_1 >> 8. No
intermediate regime allows exactly 1 disjoint pair with exactly 8 total cycles.

## Part O (VERIFIED, computational): n=9 Cycle-Rich Structure

**Key findings at n=9 (opus-S42, 10-50M samples):**

1. **alpha_1=10 composition is rigid:** ALL alpha_1=10 n=9 tournaments (no src/sink)
   have t3=6, t5=4, t7=0, t9=0. i_2 is always 9 or 10 (never 0). H is always 57 or 81.
   The (10,0) decomposition requires i_2=0, which NEVER occurs.

2. **3-cycle matching:** Among cycle-rich n=9 tournaments:
   - mm=1 (star-like): 5.0%
   - mm=2 (two disjoint): 71.1%
   - mm=3 (three disjoint): 23.9%
   The matching approach (3 disjoint => Part C) only covers 23.9%.

3. **mm<=2 minimum H = 45:** When max matching <= 2, min H = 45 >> 21. Paradoxically,
   fewer disjoint 3-cycles forces MORE 5-cycles, giving LARGER H.

4. **Deletion+Matching dichotomy:** In 153,444 cycle-rich n=9 tournaments tested,
   EVERY one satisfies at least one of:
   (a) 3 pairwise-disjoint 3-cycles exist (Part C gives H >= 27), or
   (b) A vertex deletion yields cycle-rich n=8 (induction gives H >= 25).
   Zero counterexamples found (50M trials).

5. **I(Omega,2) monotonicity:** For any v, Omega(T-v) is an induced subgraph of
   Omega(T) (cycles not through v). So I(Omega(T),2) >= I(Omega(T-v),2), i.e.,
   H(T) >= H(T-v). Combined with cycle-rich n=8 exhaustive (min H=25), any
   inductive reduction to cycle-rich n=8 gives H >= 25 > 21.

## Part P (PROVED): Complete Inductive Proof Structure

**Theorem:** H(T) != 21 for any tournament T on n vertices.

**Proof by strong induction on n.**

**Base case:** n <= 8. Exhaustive verification (Part G). H=21 never occurs.

**Inductive step (n >= 9):** Assume proved for all tournaments on < n vertices.

**Case 1:** T has a vertex not in any 3-cycle.
By Part J, this vertex v is in no directed cycle of any length.
Omega(T) = Omega(T-v), so H(T) = H(T-v). By induction, H(T-v) != 21.
(Note: source/sink vertices are a special case since they have no in-arcs
or no out-arcs, hence are in no 3-cycle.)

**Case 2:** Every vertex is in a 3-cycle (= "cycle-rich").
Note: this automatically implies no source/sink (Part Q).
By the Dichotomy Theorem (Part R), either:
(a) 3 pairwise-disjoint 3-cycles exist. Part C gives H >= 27 > 21. Done.
(b) Some deletion v yields cycle-rich T-v at n-1. By H(T) >= H(T-v) + 2
    (Part O.5, since v is in at least one cycle), and applying the strong
    inductive hypothesis H(T-v) >= 25 (for cycle-rich at n-1 >= 8), we get
    H(T) >= 27 > 21. Done.

**Strong subsidiary claim (proved simultaneously):** For all cycle-rich T on
n >= 8 vertices: H(T) >= 25.
- Base: n=8 cycle-rich, min H = 25 (exhaustive, Part L).
- Step: Case (a) gives H >= 27 >= 25. Case (b) gives H >= 25 + 2 = 27 >= 25.

## Part Q (PROVED): Cycle-Rich Implies No Source/Sink

**Lemma:** If every vertex of a tournament T is in a directed 3-cycle,
then T has no source (score n-1) and no sink (score 0).

**Proof:** A source v beats all n-1 opponents (score n-1). In any 3-cycle
v -> a -> b -> v, we need b -> v. But v beats b. Contradiction.
A sink v is beaten by all opponents (score 0). In any 3-cycle v -> a -> b -> v,
we need v -> a. But v beats nobody. Contradiction. QED.

Verified exhaustively at n=4-7 (0 counterexamples in all tournaments) and
by random sampling at n=8,9 (0 counterexamples in 500k each).

## Part R (PROVED): Dichotomy Theorem

**Theorem (Dichotomy).** Let T be a cycle-rich tournament on n >= 9 vertices.
Then at least one of:
(a) T contains 3 pairwise vertex-disjoint directed 3-cycles, or
(b) There exists vertex v such that T-v is cycle-rich on n-1 vertices.

**Proof.** If the max 3-cycle matching mm >= 3, conclusion (a) holds.
Assume mm <= 2 (hence mm in {1, 2}). Take maximal matching: disjoint
3-cycles A, B (or just A if mm=1). Let R = V \ (A union B), |R| >= 3.
Every 3-cycle C satisfies C cap (A union B) != empty (else 3-matching).

**Define the poisoning graph P = (R, E_P):** w -> v in P iff EVERY 3-cycle
containing w also contains v.

**Lemma R.1 (Bounded out-degree).** Out-degree <= 1 in P.
*Proof:* If w -> v and w -> v' with v != v', all w's 3-cycles contain both
v and v'. Only possible 3-cycle set is {w,v,v'} subset R. But then
{A, B, {w,v,v'}} is a 3-matching, contradicting mm <= 2. QED.

**Lemma R.2 (Acyclicity).** P is a DAG (directed acyclic graph).
*Proof:* Suppose P has cycle w1 -> w2 -> ... -> wk -> w1. Since wk -> w1:
all of wk's 3-cycles contain w1. So wk has a 3-cycle {wk, w1, x} with
x in A union B (since 3-cycles intersect A union B). This cycle contains
w1 but NOT w2 (since x in A union B, wk != w2, w1 != w2 for k >= 2).
But w1 -> w2 in P means all w1's cycles contain w2. Cycle {wk,w1,x}
contains w1 but not w2: contradiction. QED.

**Lemma R.3 (Safe deletion exists).** There exists v in R such that every
vertex in T-v is in a 3-cycle.
*Proof:* P is a DAG with out-degree <= 1 (Lemmas R.1-R.2), hence a forest
of directed paths. Every non-empty DAG has a source (in-degree 0 vertex).
Let v be such a source. For any w in R \ {v}:
- If w has out-degree 0 in P: w has a 3-cycle C_w within {w} union (A union B).
  Since v not in A union B and v != w: C_w persists in T-v.
- If w -> f(w) in P: since v has in-degree 0, f(w) != v. A 3-cycle of w
  is {w, f(w), x} with x in A union B. None of w, f(w), x equals v
  (f(w) != v, x in A union B, w != v). So this cycle persists in T-v.
For u in A: cycle A persists (v not in A). For u in B: cycle B persists.
So every vertex in T-v is in a 3-cycle. By Part Q: T-v is cycle-rich. QED.

**Computational verification:** 0 counterexamples in:
- 106,424 cycle-rich n=9 tournaments (kind-pasteur-S33)
- 153,444 cycle-rich n=9 tournaments (opus-S42)
- All mm=2 cases: S (vertices with out-degree 0 in P) is always non-empty
  (51,280 tests, 100% S non-empty, kind-pasteur-S33)

## Part S (VERIFIED): Bottleneck Analysis and No-Good-Deletion Structure

**Instance:** opus-2026-03-07-S43

The poisoning graph P from Part R is a special case of the bottleneck function.
We independently analyzed the bottleneck structure at n=9:

**Computational Results (100M random n=9, 93M cycle-rich):**
- Found exactly 20 no-good-deletion cycle-rich tournaments
- ALL 20 share IDENTICAL structure:
  - Score sequence: (1,1,1,4,4,4,7,7,7)
  - t3 = 3, max matching mm = 3, H = 27
  - sum|B(v)| = 9 (EVERY vertex uniquely bottlenecked)
  - Number of bottleneck centers = 6
  - The 3 disjoint 3-cycles partition ALL 9 vertices

**Observation:** All no-good-deletion cases have mm = 3, confirming Part R's
dichotomy is tight: the only obstruction to good deletion is having mm >= 3.

**n >= 10 verification:** 0 no-good-deletion in 192M cycle-rich n=10
tournaments (200M random trials), suggesting good deletion always exists
for cycle-rich n >= 10.

## Status: PROVED

**H(T) != 21 for any tournament T on any number of vertices n.**

The proof combines:
- Part G: Exhaustive base case n <= 8
- Part J: Key Lemma (no 3-cycle => no cycle => removable)
- Part Q: Cycle-rich => no source/sink (automatic)
- Part R: Dichotomy for cycle-rich n >= 9 (poisoning graph DAG argument)
- Part C: 3 disjoint 3-cycles => H >= 27
- Part L: Cycle-rich n=8 minimum H = 25
- Part O.5: H(T) >= H(T-v) + 2 when v is in a cycle (IP monotonicity)

## Scripts

- `04-computation/h21_gap_mechanism.py` — (t3, t5, p33) exhaustive analysis at n=6
- `04-computation/t3_t5_parity.py` — parity obstruction: t3=5 forces t5 even
- `04-computation/t3_t5_parity_law.py` — complete case analysis at n=6
- `04-computation/k6_2e_omega_check.py` — K_6-2e exhaustive check (opus-S40)
- `04-computation/k6_2e_sharing_proof.py` — sharing lemma approach (opus-S40)
- `04-computation/h21_n7_k6_2e.py` — n=7 exhaustive with C (opus-S40)
- `04-computation/h_spectrum_n7.py` — full H-spectrum at n=7 (opus-S40)
- `04-computation/h_spectrum_all_n.py` — H-spectrum gaps at n=3..7 (opus-S40)
- `04-computation/h_vs_alpha_n7.py` — alpha_1 vs H at n=7 (opus-S40)
- `04-computation/h_vs_alpha_n8_sample.py` — alpha_1 vs H at n=8 (opus-S40)
- `04-computation/k8e_focused_n8.py` — K_8-e i_2 analysis at n=8 (opus-S40)
- `04-computation/i2_jump_pattern.py` — i_2 jump pattern at n=5,6 (opus-S40)
- `04-computation/h21_alpha_structure_n7.py` — exhaustive alpha_1=8,10 at n=7 (opus-S41)
- `04-computation/h21_targeted_n8_v2.py` — fixed 5-cycle counting, n=8 sampling (opus-S41)
- `04-computation/h21_exhaustive_n8_v3.c` — fast exhaustive n=8 checker (opus-S41)
- `04-computation/h21_induction_n9.py` — random n=9 sampling (opus-S41)
- `04-computation/h21_matching_n9.c` — max matching of 3-cycle sets at n=9 (opus-S42)
- `04-computation/h21_10_0_n9_check.c` — (10,0) decomposition check at n=9 (opus-S42)
- `04-computation/h21_alpha10_h_values_n9.c` — alpha_1=10 H-value analysis at n=9 (opus-S42)
- `04-computation/h21_no_deletion_no_matching.c` — dichotomy check: both-bad = 0 (opus-S42)
- `04-computation/h21_match2_h_values.c` — min H for mm<=2 cycle-rich (opus-S42)
- `04-computation/h21_inductive_min_h.c` — good deletion check at n=9 (opus-S42)
- `04-computation/h21_poisoning_graph.py` — poisoning graph DAG analysis (kind-pasteur-S33)
- `04-computation/h21_lichiardopol_check.py` — Lichiardopol threshold check (kind-pasteur-S33)
- `04-computation/h21_cycle_rich_auto_no_ss.py` — cycle-rich => no source/sink (kind-pasteur-S33)
- `04-computation/h21_source_sink_avoidance.py` — source/sink avoidance analysis (kind-pasteur-S33)
- `04-computation/h21_dichotomy_proof.py` — dichotomy verification at n=9 (kind-pasteur-S33)
- `04-computation/h21_bottleneck_v2.c` — bottleneck structure at n=9: 20 no-good-del all (1,1,1,4,4,4,7,7,7) (opus-S43)
- `04-computation/h21_no_good_del_n10.c` — 0/192M no-good-deletion at n=10 (opus-S43)
- `04-computation/h21_moon_bound_cycle_rich.py` — Moon's formula t3 lower bounds (opus-S43)
- `04-computation/h21_extreme_low_t3.py` — min H by t3 for cycle-rich n=9 (opus-S43)
- `04-computation/h21_achievable_w_n7.c` — exhaustive w-value spectrum at n=7 (opus-S43)
- `04-computation/h21_gap_pattern_analysis.py` — algebraic structure of w=3,10 gaps (opus-S43)
