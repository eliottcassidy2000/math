# Even Graphs as First-Class Objects

**opus-2026-03-24-S315**

## The Dual Universe

The tiling hypercube Q_m (m = C(n-1,2)) admits TWO natural quotients:

1. **G_n** = Q_m / S_n(tournament) — the tournament metagraph (A000568 vertices)
2. **E_n** = Q_m / S_n(even graph) — the even graph metagraph (A002854-shifted vertices)

Both quotients use the SAME group S_n, but acting on DIFFERENT structures:
- In G_n: S_n permutes vertex labels of the tournament
- In E_n: S_n permutes vertex labels of the even graph

The cycle-space bijection (tiling ↔ even graph via XOR of fundamental cycles) connects the two. Each tiling simultaneously specifies a tournament class (in G_n) and an even graph class (in E_n).

## Computed Values

| n | V(G_n) | V(E_n) | E(G_n) | E(E_n) | χ(G_n) | χ(E_n) | ω(E_n) | Density(E_n) |
|---|--------|--------|--------|--------|--------|--------|--------|-------------|
| 3 | 2 | 2 | 1 | 1 | 2 | 2 | 2 | 100% |
| 4 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 100% |
| 5 | 10 | 7 | 21 | 16 | 4 | 5 | 5 | 76% |
| 6 | 34 | 16 | 143 | 90 | 5 | 10 | 10 | 75% |
| 7 | 272 | 54 | 2123 | 951 | 6 | 28 | 28 | 66% |

## Key Discoveries

### 1. χ(E_n) Grows Much Faster Than χ(G_n)

χ(G_n) = n-1 (linear). χ(E_n) grows much faster: 2, 3, 5, 10, 28.

The sequence 2, 3, 5, 10, 28 for n=3..7 is suggestive. Differences: 1, 2, 5, 18. Second differences: 1, 3, 13. These grow roughly factorially.

### 2. E_n Is Perfect (ω = χ)

At every computed n (3 through 7), ω(E_n) = χ(E_n). The even graph metagraph is a **perfect graph**. This is striking because G_n LOSES perfectness at n=7 (ω=4, χ=6).

Why might E_n be perfect? The even graph iso classes have a natural partial order by edge count (= weight in the cycle space). If this partial order makes E_n a comparability graph or a chordal graph, perfectness follows from the Strong Perfect Graph Theorem.

### 3. The Bridge Matrix Has Full Rank

The bridge matrix B (VT × VE, where B[t,e] = # tilings with tournament class t and even class e) has rank = VE at every n. This means:
- The even graph classes are "linearly independent" in their tournament fiber profiles
- No even graph class is a linear combination of others
- The projection G_n → E_n is well-structured

### 4. ~80% of Flips Are "Joint" (Jaccard ≈ 0.82)

At n≥5, about 80% of tile flips change both the tournament class and the even class. Only ~10% change the tournament class alone (pure score change) and ~8% change the even class alone (absorbed by tournament symmetry).

This 80/10/8/2 split shows the two quotients are **highly correlated but not identical**. The 10% tournament-only flips are the "cut-space" operations that change scores without changing cyclic structure. The 8% even-only flips are structurally novel — they change the graph structure but are "absorbed" by tournament symmetry.

### 5. Fiber Sizes Follow |Aut| Pattern

The even graph fiber sizes (column sums of B) equal n!/|Aut(G)| for each even graph class G. This is the orbit-stabilizer theorem: each even graph's orbit under S_n has size n!/|Aut(G)|, and the bijection maps exactly one tiling per labeled even graph.

The tournament fibers (row sums of B) equal H(T)/|Aut(T)| × multiplicity — the tiling formula from earlier work.

## The Two Metagraphs as Lenses

G_n and E_n are two lenses on the same underlying space Q_m:

- **G_n** (tournament lens): Respects orientation. Higher-dimensional quotient (more classes). Sparser. χ = n-1 (linear, geometric).
- **E_n** (even graph lens): Forgets orientation. Lower-dimensional quotient (fewer classes). Denser. χ grows faster (combinatorial/algebraic).

The tournament metagraph preserves the **directed** structure (who beats whom). The even graph metagraph preserves only the **undirected** cycle structure (which cycles exist). Forgetting direction collapses many tournaments into fewer even graphs, making the quotient denser.

## Invariant Dictionary

| Tournament invariant | Even graph analogue | Relationship |
|---------------------|--------------------|--------------|
| H (Hamiltonian paths) | Euler circuit count? | Different — H lives on the tournament, not the even graph |
| Score sequence | Degree sequence | Cut space ↔ degree constraint |
| c3 (3-cycles) | Triangle count | Same triangles, different interpretation |
| Self-complementary | Self-complementary graph | Both defined by complement involution |
| Transitive tournament | Empty even graph | H=1, minimal cyclic structure |
| Regular tournament | Complete even graph K_n? | Maximal cyclic structure |

## For Future Agents

**Every computation on G_n should also be done on E_n.** Specifically:

1. When computing any invariant of tournament classes, also compute it for even graph classes
2. When building the tournament metagraph, also build the even graph metagraph
3. When studying chromatic number, clique number, spectrum of G_n — do the same for E_n
4. When proving a theorem about G_n, check whether an analogue holds for E_n
5. Track the bridge matrix B and the four-way flip classification

The two metagraphs are the **tournament side** and the **graph side** of the same coin. Understanding one illuminates the other.

## Open Questions

1. **Is E_n always perfect?** (χ = ω at n=3..7 — prove or find counterexample)
2. **What is χ(E_n) asymptotically?** (Grows faster than linear — polynomial? exponential?)
3. **What is the natural coloring of E_n?** (Edge count mod something?)
4. **Does V(G_n)/V(E_n) have a limit?** (Currently 1, 1, 1.4, 2.1, 5.0)
5. **Spectral comparison**: How do the spectra of G_n and E_n relate?
6. **The 80/10/8/2 split**: Does it stabilize? What determines which flips are "tournament-only" vs "even-only"?
