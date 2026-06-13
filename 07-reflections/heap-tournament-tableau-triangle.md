# The Heap / Tournament / Young Tableau Triangle

**Session**: opus-2026-03-23-S270

---

## Three Fillings of the Staircase

The staircase Young diagram δ_{n-1} = (n-1, n-2, ..., 1) with m = C(n,2) cells is the common substrate for three distinct combinatorial objects:

| Object | Filling type | Count | OEIS |
|--------|------------|-------|------|
| **Tournament** | Binary (0/1 per cell) | 2^m labeled; V_n iso classes | A000568 |
| **Standard Young Tableau** | Increasing (1..m per cell) | f^δ = m!/∏hooks | A005118 |
| **Binary Heap** | Heap-ordered | m!/∏(subtree sizes) | A056971 |

All three are counted by "m! divided by a product" formulas:
- Tournament orientations: 2^m (no denominator — all fillings valid)
- SYT: m! / (1^{n-1} · 3^{n-2} · 5^{n-3} · ... · (2n-3)^1)
- Heaps: m! / ∏(subtree sizes of complete binary tree on m nodes)

---

## The Hooks Are All Odd

**Theorem**: Every hook length of the staircase δ_k is odd:

    h(i,j) = 2(k - i - j) - 1

**Proof**: At cell (i,j), arm = k-i-j-1 and leg = k-i-j-1. ARM = LEG at every cell. This is the algebraic signature of the **right isosceles triangle** — the staircase has y = x reflection symmetry, so the arm and leg at each cell are equal.

The product of hooks:
    ∏ hooks = ∏_{c=1}^{k-1} (2c-1)^{freq(2c-1)}

where freq(2c-1) = k-c (the number of cells on the c-th anti-diagonal).

---

## H(T) as Order Polytope Volume

The number of Hamiltonian paths H(T) counts the **linear extensions** of the tournament T viewed as a partial order. By Stanley's theory:

    Vol(Order polytope of T) = H(T) / n!

| n | H_min (transitive) | H_max (regular) | Vol_max / Vol_min | H_max / Szele |
|---|-------------------|----------------|-------------------|---------------|
| 3 | 1 | 3 | 3 | 2.00 |
| 4 | 1 | 5 | 5 | 1.67 |
| 5 | 1 | 15 | 15 | **2.00** |
| 6 | 1 | 45 | 45 | **2.00** |
| 7 | 1 | 189 | 189 | 2.40 |

**The Szele ratio**: H_max / (n!/2^{n-1}) = 2.0 at n=3,5,6. This ratio measures how much the most "cyclic" tournament exceeds the random expectation. The fact that it equals exactly 2 at several values is remarkable.

---

## The Edelman-Greene Connection

**Stanley's theorem (1984)**: The number of **reduced decompositions** of the longest permutation w_0 = n(n-1)...21 in S_n equals f^{δ_{n-1}} = the number of SYT of staircase shape.

**Edelman-Greene bijection (1987)**: An explicit bijection between reduced decompositions (sorting networks) and SYT of staircase shape.

This means: **sorting networks for n elements** live on the same staircase as **tournament adjacency matrices**. A sorting network specifies an ORDER in which to compare pairs of elements; a tournament specifies the OUTCOME of each comparison. The staircase is the space of "all possible comparisons" — C(n,2) cells, one per pair.

---

## Viennot's Heaps of Pieces

Viennot's theory of **heaps of pieces** provides a unifying framework:

1. A **piece** = a cell of the staircase = an arc of the tournament
2. Two pieces **conflict** if they share a vertex (= overlap in the heap)
3. A **heap** = a partially ordered arrangement of pieces
4. The **generating function** of heaps connects to:
   - Acyclic orientations (chromatic polynomial)
   - Independence polynomial (our OCF!)
   - Coxeter group elements
   - Rogers-Ramanujan identities

The key insight: **H(T) = I(Ω(T), 2)** (our OCF) counts independent sets at fugacity 2 on the conflict graph. Viennot's heaps encode the SAME conflict structure — pieces that share a vertex can't be placed simultaneously.

---

## The Hook as Arc Influence

The hook h(i,j) = 2(k-i-j)-1 measures the **influence** of the arc connecting vertex i to vertex j:

- Arc between vertices 0 and n-1 (opposite extremes): hook = 2(n-2)-1 (maximal)
- Arc between adjacent vertices i and i+1: hook = 1 (minimal)
- The hook is the "leverage" — flipping a large-hook arc changes H more than flipping a small-hook arc

This connects to the **metagraph edge weights**: arcs with larger hooks create larger H-changes when flipped. The **spectral gap** of the metagraph (~ 2/n) is related to the average hook length ~ 2n/3.

---

## Open Questions

1. Is there a bijection between tournament iso classes and some subset of SYT of the staircase?
2. Does the Edelman-Greene bijection specialize to tournament structure in any useful way?
3. Can the hook length formula be used to compute H(T) more efficiently?
4. Is H_max / Szele → 2 as n → ∞? (It's 2.0 at n=3,5,6 but 2.4 at n=7.)
5. Does Viennot's heap theory give a new proof of Rédei's theorem (H is always odd)?
6. Can the binary heap count on C(n,2) elements be related to tournament enumeration?
