# THM-022: Tiling model structure and independence polynomial relationships

**Status:** Collection of PROVED and COMPUTATIONALLY VERIFIED results
**Certainty:** 5 (proved results), 4 (computational for n <= 7)
**Author:** opus-2026-03-06-S16
**Dependencies:** Tiling model (definitions.md), OCF (H(T)=I(Omega(T),2))

---

## Definitions (from tournament-tiling-explorer.html)

**Base path:** n -> n-1 -> ... -> 1 (fixed Hamiltonian path).

**Tiles:** Non-path edges (a,b) with a >= b+2. Count: m = C(n-1,2).

**Tiling:** Binary vector t in {0,1}^m. Bit 0 = forward (a->b), bit 1 = backward (b->a).

**Flip (F):** Complement all bits. F(t)_i = 1 - t_i. Reverses all non-path arcs.

**Grid transpose (G):** Permute bits by the map (x,y) -> (n+1-y, n+1-x). At the tournament level, this corresponds to the vertex relabeling v -> n+1-v combined with non-path arc reversal. A tiling is **grid-symmetric** if G(t) = t.

**Self-flip member:** A tiling T where flip(T) is in the same isomorphism class as T.

**Blueself:** A self-flip member that is also grid-symmetric.

**Blackself:** A self-flip member that is NOT grid-symmetric.

---

## Theorem 1: Transpose preserves I(Omega(T), x)

**Statement:** For any tournament T, I(Omega(T), x) = I(Omega(T^op), x).

**Proof:** If v_1 -> v_2 -> ... -> v_k -> v_1 is a directed cycle in T, then in T^op all arcs reverse, giving v_k -> v_{k-1} -> ... -> v_1 -> v_k, a directed cycle on the SAME vertex set. The map C -> reverse(C) is a bijection between directed cycles of T and T^op that preserves vertex sets. Since the conflict graph Omega depends only on which pairs of cycles share vertices, Omega(T) = Omega(T^op) as graphs. Hence I(Omega(T), x) = I(Omega(T^op), x). QED

**Verification:** Exhaustive at n = 4 (8), n = 5 (64). Every tournament has identical odd-cycle vertex sets in T and T^op.

---

## Theorem 2: Flip and transpose commute

**Statement:** For any tiling t, F(G(t)) = G(F(t)).

**Proof:** Flip is a pointwise operation (invert each bit independently). Transpose is a permutation of bit positions. Pointwise operations always commute with permutations:

F(G(t))_j = 1 - G(t)_j = 1 - t_{sigma^{-1}(j)}
G(F(t))_j = F(t)_{sigma^{-1}(j)} = 1 - t_{sigma^{-1}(j)}

These are identical. QED

**Verification:** 100% at n = 4 (8/8), n = 5 (64/64), n = 6 (1024/1024).

---

## Theorem 3: Flip preserves grid-symmetry

**Statement:** If a tiling t is grid-symmetric, then F(t) is also grid-symmetric.

**Proof:** Grid-symmetry means G(t) = t. By Theorem 2, G(F(t)) = F(G(t)) = F(t). So F(t) is grid-symmetric. QED

**Corollary:** Grid-symmetry is a flip-invariant. Equivalently, if t is NOT grid-symmetric, then F(t) is also NOT grid-symmetric.

**Corollary:** Blueself tilings come in flip-pairs: if T is blueself, then F(T) is also blueself (grid-symmetric + same isomorphism class).

---

## Theorem 4: Score-sequence structure of grid-symmetric tilings

**Statement:** For a grid-symmetric tiling at any n, the non-path out-degrees of the two path endpoints satisfy k_0 + k_{n-1} = n - 2.

**Proof:** For a grid-symmetric tiling, the non-path arc directions satisfy: a -> b iff (n+1-b) -> (n+1-a). Under the vertex relabeling v -> n+1-v, position 0 (vertex n) maps to position n-1 (vertex 1), and vice versa. The non-path neighbors of position 0 are {positions 2,...,n-1}, which map bijectively to {positions 0,...,n-3} = non-path neighbors of position n-1.

A non-path arc from position 0 to position j corresponds (by grid-symmetry) to a non-path arc from position n-1-j to position n-1. So each non-path win of vertex n corresponds to a non-path loss of vertex 1, giving k_0 + k_{n-1} = n - 2. QED

---

## Theorem 5: No blueself at odd n (Score obstruction) — PROVED FOR ALL ODD n

**Statement:** For odd n, no grid-symmetric tiling has the same sorted score sequence as its flip. Consequently, no blueself tilings exist at any odd n.

**Proof:** By Theorem 4, grid-symmetric tilings satisfy k_0 + k_{n-1} = n-2. The flip maps endpoint scores as:
- s_0 = 1 + k_0 -> s'_0 = n - 1 - k_0
- s_{n-1} = k_{n-1} = n-2-k_0 -> s'_{n-1} = k_0

For T and flip(T) to be isomorphic, they must have the same sorted score sequence. In particular, the multiset {s_0, s_{n-1}} must equal {s'_0, s'_{n-1}}. This requires either:

**Case A:** s_0 = s'_0, i.e., 1+k_0 = n-1-k_0, i.e., k_0 = (n-2)/2. At odd n, (n-2)/2 is not an integer. IMPOSSIBLE.

**Case B:** s_0 = s'_{n-1} and s_{n-1} = s'_0, i.e., 1+k_0 = k_0 and n-2-k_0 = n-1-k_0. The first equation gives 1=0. IMPOSSIBLE.

Both cases lead to contradiction. Therefore, no grid-symmetric tiling at odd n can have the same sorted score sequence as its flip. Since isomorphic tournaments must have the same score sequence, no grid-symmetric tiling is self-flip, i.e., no blueself exists at odd n. QED.

**Verification script:** `04-computation/blueself_odd_n_proof.py`
**Verification:** Exhaustive at n=3 (0/2), n=5 (0/16), n=7 (0/512). Algebraic proof covers all odd n.
**Control:** At even n, score matches DO occur: 2/4 at n=4, 24/64 at n=6.

---

## Theorem 6: Blueself/blackself mutual exclusivity at the class level

**Statement:** For n <= 7, no isomorphism class contains both blueself and blackself tilings.

**Verification record:**

| n | Classes | Self-flip classes | Blueself classes | Blackself classes | Mixed | Status |
|---|---------|-------------------|------------------|-------------------|-------|--------|
| 4 | 4 | 1 | 1 | 0 | 0 | PROVED |
| 5 | 12 | 2 | 0 | 2 | 0 | PROVED |
| 6 | 56 | 8 | 2 | 6 | 0 | PROVED |
| 7 | 456 | 30 | 0 | 30 | 0 | PROVED |

**Proof for odd n:** At odd n, Theorem 5 shows blueself tilings don't exist, so mutual exclusivity is vacuous.

**Proof for even n (n = 4, 6):** Verified exhaustively. Every self-flip class has uniform grid-symmetry among its self-flip members. This follows from Theorem 3 (flip preserves grid-symmetry), which ensures self-flip members come in same-grid-symmetry pairs, combined with the empirical fact that no class has self-flip members of both types at n = 4, 6.

**Open:** Whether mutual exclusivity holds at n = 8 and beyond. At even n, blueself classes exist and the question is non-trivial.

---

## Theorem 7: Blueself existence dichotomy

**Statement (odd part PROVED for all n; even part verified for n <= 6):**
- Even n: Blueself classes exist. At n=4: 1 class (H=5). At n=6: 2 classes (H=37, H=45).
- Odd n: No blueself classes exist. **PROVED for ALL odd n** by Theorem 5.

The blueself classes at even n are characterized by:
- Score sequences close to regular: (2,2,1,1) at n=4, (3,3,3,2,2,2) at n=6
- HIGH H values (among the largest in their n)
- Positive alpha_2 (non-trivial disjoint cycle structure)

---

## Additional observations (not proved in general)

### Flip does NOT preserve I(Omega(T), x)

The flip operation changes the independence polynomial in general. This is because flip(T) is NOT isomorphic to T^op — it reverses only non-path arcs while T^op reverses all arcs. Flip preserves I.P. only for self-flip tilings (trivially, since same isomorphism class).

### flip(T) is NOT always isomorphic to T^op

| n | flip(T) ≅ T^op | Total | Percentage |
|---|-----------------|-------|------------|
| 4 | 2 | 8 | 25% |
| 5 | 8 | 64 | 12.5% |
| 6 | 40 | 1024 | 3.9% |

### I(Omega, x) is NOT a complete isomorphism invariant

Multiple non-isomorphic classes can share the same independence polynomial. Already at n=4: classes with scores (2,2,2,0) and (3,1,1,1) both have I.P. = [1,1].

### Self-flip classes have structured score sequences

At every tested n, self-flip classes have score sequences that are "balanced" — close to the regular tournament score pattern. This connects to the known result that self-converse and near-regular tournaments tend to maximize H(T) (see T070, T079).

---

## Significance

1. **Grid-symmetry as flip-invariant (Theorem 3):** Established via the elegant commutativity of pointwise operations with permutations. This structural result constrains how the flip operation interacts with the tiling model.

2. **Odd-n obstruction (Theorem 5):** The score-sequence argument shows that the parity of n controls whether grid-symmetric tournaments can be self-flip. This is an arithmetic obstruction with no known analogue in general graph theory.

3. **Blueself-blackself dichotomy:** The clean separation into blue (grid-symmetric) and black (non-grid-symmetric) self-flip classes, with even/odd n controlling existence, reveals deep structure in the tiling model that may inform the independence polynomial conjecture.

---

## Verification Scripts

- `04-computation/blueself_odd_n_proof.py` — Theorem 5 (odd-n blueself obstruction)
- `04-computation/blueself_sc_maximizer_connection.py` — blueself vs SC maximizer analysis
- `04-computation/blackself8_investigation.py` — blackself investigation at n=8

## See also

- THM-020 (real-rootedness via claw-freeness)
- THM-021 (real-rootedness via discriminant/Turán)
- T070 (self-converse tournaments maximize H)
- T079 (perpendicular bell curve)
- OPEN-Q-015 (real roots for general n)
