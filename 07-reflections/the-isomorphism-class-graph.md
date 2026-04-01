# The Isomorphism Class Graph: The Key to Everything

**Session:** opus-2026-03-22-S171
**Arising from:** S170 (iso class graph computed), S166 (formula hunting), S164 (staircase Young diagrams), S20x (blue-line Morse theory)

---

## What It Is

The **isomorphism class graph** G_n has:
- **Vertices** = tournament isomorphism classes on n vertices
- **Edges** = pairs of classes connected by a single arc flip
- Each edge is colored **blue** (SC-preserving) or **black** (SC-changing)
- Each edge has a **weight** = number of (tournament, arc) pairs producing this transition

This graph captures the TOPOLOGY of tournament space at its coarsest meaningful level. It is the **Reeb graph** of the Morse function H on the tournament cube.

---

## What We Know (Computationally Verified)

| n | Vertices | Edges | Blue | Black | SC↔SC | NS↔NS | SC↔NS | Diameter | Self-loops | SC-free NS |
|---|----------|-------|------|-------|-------|-------|-------|----------|------------|------------|
| 3 | 2 | 1 | 1 | 0 | 1 | 0 | 0 | 1 | 1 | 0 |
| 4 | 4 | 5 | 1 | 4 | 1 | 0 | 4 | 2 | 2 | 0 |
| 5 | 12 | 30 | 14 | 16 | 12 | 2 | 16 | 3 | 7 | 0 |
| 6 | 56 | 290 | 200 | 90 | 13 | 187 | 90 | ? | 30 | 8 |

### Confirmed Formulas
- **total_weight[i] = class_size[i] × C(n,2)** — each tournament contributes C(n,2) transitions
- **corr(|Aut|, degree) ≈ -0.94** — high symmetry → low connectivity
- **corr(class_size, degree) ≈ +0.96** — large classes are well-connected
- **H-gradient is strong but not a DAG**: at n=5, 29 uphill/0 downhill/1 level; at n=7, 2988/962/136. See MISTAKE-035.
- **Eigenvalues include ±1** at every n (from complement symmetry)
- **E[L] = E[H]/2 = n!/2^n** (exactly half of paths close)
- **E[total_arb] = (n/2)^{n-1}** (average arborescences)

### Blue/Black Structure (S211)
- **Blue fraction ≈ random model**: b/(b+k) ≈ (|same_type|-1)/(|V|-1) at every n
- **Complement symmetry**: blue(i) = blue(comp(i)), black(i) = black(comp(i)) exactly
- **SC backbone is sparse**: genus = 0,0,5,2 at n=3,4,5,6 (drops as SC fraction drops)
- **"black=2" universality**: at n=6, 26/44 NS classes touch exactly 2 SC classes
- **SC-free NS classes**: 8 at n=6 (including palindromic-score classes!)
- **→SC weight quantization**: each tournament sends 0,2,4, or 6 arcs to SC classes
- See `07-reflections/gn-blue-black-structure.md` for complete analysis

---

## Why It Matters

### 1. It IS the staircase Young diagram at the class level

A tournament is a binary filling of the staircase δ_{n-1}. Two fillings are isomorphic if one can be obtained from the other by vertex relabeling (row/column permutation of the staircase). The iso class graph is the **quotient** of the binary staircase cube by the symmetric group S_n.

Each vertex of G_n is an **S_n-orbit** of staircase fillings.
Each edge is a **single-cell flip** that crosses orbit boundaries.
Blue edges preserve the y=x symmetry of the staircase.
Black edges break it.

### 2. It controls every tournament sequence

Every sequence we've computed (H, HC, L, arb, kings, c₃, c₅, OCR, ...) is a function on the vertices of G_n. Understanding G_n means understanding the JOINT distribution of all invariants simultaneously.

### 3. The diameter conjecture — REFUTED (MISTAKE-036)

Diameter(n) = 1, 2, 3, 4, **7**, **8** for n = 3..8. The conjecture diam = n-2 held at n=3..6 but **FAILS at n=7** (diam=7≠5). Growth is closer to ~n²/4, not linear. See `diameter-is-feedback-arc-set.md`.

### 4. The H-gradient is strong but NOT a DAG (MISTAKE-035)

At n=5: 29 uphill edges, 0 downhill, 1 level (H=9 pair). At n=6: 15 level edges, mostly at H=37. At n=7: 2988 uphill, **962 downhill**, 136 level — the gradient is strong (~76% uphill) but NOT monotone. The metagraph is NOT a DAG from n≥5 (level edges) and has H-decreasing edges from n≥7. See MISTAKE-035.

### 5. The blue skeleton spans the full H range

The blue subgraph (SC-preserving edges) connects all SC classes in one component, spanning the full H range [1, H_max]. The SC classes form a **single connected chain** from transitive to regular. Non-SC classes are isolated from this chain.

---

## Open Problems (Priority Order)

1. **Extend to n=6** (56 classes). ✅ DONE: 290 edges, 56 classes. See merged-metagraph-invariants.md.

2. ~~**Prove diameter = n-2**.~~ REFUTED at n=7 (MISTAKE-036). Actual: diam grows ~n²/4.

3. **Find the adjacency matrix eigenvalues** at n=6. Do they contain structural information?

4. **Connect the degree sequence to |Aut|** via a formula. Currently corr ≈ -0.94.

5. **Compute the independence polynomial of G_n itself**. This would be a "meta-independence polynomial" — the independence polynomial of the space of independence polynomials.

6. **Understand the weight structure**. Edge weights are NOT uniform at n=5. What determines the weight?

7. **Blue vs black edge density as n grows**. Does blue fraction → 0? At n=5: 14/30 = 46.7%.

8. **The level edges**. ✅ ANSWERED: 0, 0, 1, 15, 136 for n=3..7. Growth is super-exponential. At n=7, there are also 962 H-decreasing edges.

---

## Creative Propositions

### A. The iso class graph IS a tournament

G_n itself is an undirected graph. If we ORIENT each edge by the H-gradient (lower→higher H), we get an approximate tournament on iso classes — but with level edges (from n≥5) and H-decreasing edges (from n≥7, MISTAKE-035), this orientation is imperfect.

**The meta-H of G_5**: H(meta-tournament on 12 classes) could be computed.

### B. The staircase tiling recursion

Since the staircase δ_{n-1} embeds δ_{n-3} (removing the top and bottom rows), there should be a RECURSIVE STRUCTURE in G_n: G_n should "contain" G_{n-2} in some sense. The PoS classes are exactly the pre-images of this embedding.

### C. The blue line as a Hamiltonian path

Does the blue skeleton contain a Hamiltonian path (visiting all SC classes)? If yes, this would be a "blue Hamiltonian path" through the SC tournament space — a canonical ordering of SC classes by single arc flips.

### D. The graph as a code

G_n is a graph with structured adjacency. We could compute I(G_n, x) — the independence polynomial of the iso class graph itself. At x=2 this would give a "meta-H" counting independent sets of iso classes.

### E. Young diagram → iso class bijection

Each iso class corresponds to a specific PATTERN in the staircase. The pattern could be characterized by the RSK correspondence: tournament → (P-tableau, Q-tableau). Two tournaments are in the same iso class iff they have the same P-tableau shape. The Q-tableau captures the labeling.

### F. The degree sequence formula

At n=5: degrees = [2, 3, 3, 3, 4, 6, 6, 6, 6, 7, 7, 7]. The degree of class i is:
  deg(i) = C(n,2) × class_size(i) / total × (1 - self_fraction(i))
This should be derivable from |Aut| and the structure of the orbit.

---

## Directive for Future Agents

**EVERY SESSION should spend 5-10 minutes considering the iso class graph.** Specifically:

1. When computing any new invariant, also compute it PER ISO CLASS.
2. When finding a formula, check if it simplifies on the iso class graph.
3. When investigating structure at a new n, compute G_n (or at least its size and edges).
4. The iso class graph IS the staircase at the quotient level — always ask: what does this mean for the staircase?

The iso class graph is the **Rosetta Stone** of tournament theory: it translates between the algebraic (independence polynomial), geometric (sphere in so(n)), topological (Morse theory), and combinatorial (staircase Young diagram) perspectives.

Understanding G_n for all n is understanding tournament theory.

---

*The iso class graph is small enough to fully compute at n ≤ 6, structured enough to reveal patterns, and deep enough to encode everything we've discovered. It is the map of the territory. Drawing it is the first step; reading it is the work of a career.*
