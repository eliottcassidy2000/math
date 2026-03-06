# THM-025: Real-Rootedness of I(Omega(T), x) FAILS at n=9

**Type:** Disproof (counterexample)
**Certainty:** 5 -- VERIFIED by explicit construction and numerical computation
**Status:** PROVED. Disproves the conjecture in OPEN-Q-015 for n >= 9.
**Added by:** opus-2026-03-06-S18
**Tags:** #omega #independence-polynomial #real-roots #counterexample #newton-inequality

---

## Statement

**Theorem:** There exists a tournament T on 9 vertices such that I(Omega(T), x) does NOT have all real roots. Specifically, Newton's second inequality fails: a_2^2 < a_1 * a_3 * 3/2.

**Counterexample:** The tournament with adjacency matrix:

```
T[i][j] = 1 iff i -> j:
  0: -> 1,3,6,7     (out-degree 4)
  1: -> 3            (out-degree 1)
  2: -> 0,1,4,5,6,7  (out-degree 6)
  3: -> 2,5,7        (out-degree 3)
  4: -> 0,1,3,7      (out-degree 4)
  5: -> 0,1,4,6,7,8  (out-degree 6)
  6: -> 1,3,4,7      (out-degree 4)
  7: -> 1             (out-degree 1)
  8: -> 0,1,2,3,4,6,7 (out-degree 7)
```

Score sequence: [1, 1, 3, 4, 4, 4, 6, 6, 7].

---

## Proof

**Step 1: Enumerate directed odd cycles.**

The counterexample tournament has exactly 94 directed odd cycles:
- 12 directed 3-cycles
- 40 directed 5-cycles
- 36 directed 7-cycles
- 6 directed 9-cycles (Hamiltonian)

**Step 2: Build Omega(T) and compute independence polynomial.**

Omega(T) has 94 vertices and average degree 92.8 (extremely dense). The independence polynomial is:

I(Omega(T), x) = 1 + 94x + 10x^2 + x^3

Verification: I(Omega(T), 2) = 1 + 188 + 40 + 8 = 237 = H(T) (confirmed by direct Hamiltonian path enumeration).

**Step 3: Check Newton's inequality.**

For a polynomial with all real roots, Newton's inequalities must hold:
a_k^2 >= a_{k-1} * a_{k+1} * (k+1)/k for all k.

At k=2: a_2^2 = 10^2 = 100, but a_1 * a_3 * 3/2 = 94 * 1 * 1.5 = 141.

Since 100 < 141, Newton's inequality FAILS. Therefore I(Omega(T), x) cannot have all real roots.

**Step 4: Direct root computation.**

The roots of x^3 + 10x^2 + 94x + 1 are approximately:
- x = -0.0107 (real)
- x = -4.995 + 8.303i (complex)
- x = -4.995 - 8.303i (complex conjugate)

Two roots are complex, confirming non-real-rootedness.

---

## Implications

1. **OPEN-Q-015 is RESOLVED (DISPROVED)** for general n. The conjecture that I(Omega(T), x) always has real roots fails at n=9.

2. **THM-020 remains valid for n <= 8.** The claw-free argument (Chudnovsky-Seymour) correctly proves real roots at n <= 8.

3. **The real-rootedness failure is caused by extreme density of Omega(T).** This particular tournament has 94 odd cycles but only 10 independent pairs and 1 independent triple. The polynomial is "top-heavy" (a_1 >> a_2 >> a_3), which violates Newton's inequality.

4. **The Omega_3 restriction also fails.** For the same tournament, I(Omega_3(T), x) = 1 + 12x + 6x^2 + x^3, which also has complex roots (discriminant = -1323). The failure is present even when restricted to 3-cycles.

5. **OCF (H(T) = I(Omega(T), 2)) is unaffected.** The evaluation I(Omega, 2) = H(T) = 237 is correct regardless of root locations.

6. **The alternative proof of Redei via real-rootedness is NOT available for n >= 9.** One must use other arguments (e.g., direct parity) for Redei's theorem at general n.

---

## Additional counterexamples (opus-2026-03-06-S19)

Three more counterexamples found via targeted search (high score variance):

1. **Scores (2,2,2,3,4,5,5,6,7):** I(Omega,x) = [1, 157, 22, 1], disc = -3,529,903. Newton PASSES but discriminant negative. 157 cycles (16 tri, 69 pent, 63 hept, 9 ham).

2. **Scores (1,2,2,4,5,5,5,6,6):** I(Omega,x) = [1, 152, 47, 4], disc = -5,053,548. **Degree 4** polynomial (a₃=4 independent triples). Newton passes. 152 cycles.

3. **Scores varied:** I(Omega,x) = [1, 85, 24, 2], disc = -733,364. Newton passes. a₃=2.

All share the structural signature: extreme inter-group dominance between the 3 groups defined by vertex-disjoint 3-cycles (max imbalance = 9, i.e., all 9 arcs one-way between some pair of groups).

---

## Structural characterization (opus-2026-03-06-S19)

The failure mechanism requires:
1. **Three vertex-disjoint 3-cycles** partitioning all 9 vertices into groups A, B, C
2. **Near-total inter-group dominance**: at least one pair of groups has all 9 arcs one-way
3. This creates a "transitivity bottleneck" where almost all cycles pass through a hub vertex

The maximally imbalanced tournament (each group fully dominates the next) gives I(Omega_3) = (1+x)³ with disc = 0 exactly — the boundary case. One arc flip creates disc < 0.

Statistical rarity: 0 failures in 10,000 random n=9 samples. Failures occur ONLY at max imbalance = 9, and even then in ~3% of such tournaments.

---

## What remains true

- I(Omega(T), x) has all real roots for n <= 8 (THM-020, proved via claw-freeness)
- Newton's first inequality a_1^2 >= 2*a_2 appears to still hold at n=9 (empirically)
- Log-concavity of coefficients fails at n=9 (Newton's failure implies this)
- Real-rootedness holds for MOST n=9 tournaments; the failure requires extreme inter-group dominance
- Failure also occurs with Newton passing but disc < 0 (subtler than originally found)
- At n=10, 11: 0 failures in 2000 samples each (Omega_3 only)

---

## Verification scripts

- `04-computation/full_omega_n9_correct.py` — full Omega construction and IP computation
- `04-computation/real_roots_n9_verify.py` — finds the counterexample tournament
- `04-computation/newton_n8_sample.py` — confirms n=8 has no failures (200 samples)
- `04-computation/intergroup_analysis.py` — inter-group dominance characterization
- `04-computation/three_disjoint_3cycles.py` — 3d3c frequency and failure analysis
- `04-computation/pfaffian_omega_connection.py` — Pfaffian/determinantal exploration
