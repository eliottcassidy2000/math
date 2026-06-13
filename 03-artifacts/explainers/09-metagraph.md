# The Map of All Tournaments

## A Universe of Tournaments

For any fixed number of players n, there is a finite collection of structurally distinct tournaments. For n = 7 players, there are exactly 456 distinct tournament types. For n = 8, over 6 million.

This project builds a **map** of all tournament types — a graph (in the mathematical sense) where:
- Each **dot (vertex)** represents one tournament type
- Two dots are **connected by a line (edge)** if the two tournaments differ by exactly one arrow flip

This map is called the **metagraph** G_n. Understanding its structure tells you how the universe of tournaments is organized.

---

## Moving Through the Map

When you flip one arrow in a tournament — reverse one "player A beats player B" to "player B beats player A" — you might stay in the same structural type (if the reversal didn't change anything that distinguishes the tournament from isomorphic copies), or you might move to a neighboring type.

In the tiling picture (see explainer 8), each arrow corresponds to one square. Flipping that square is a minimal move — one step through the hypercube. In the metagraph, one step = one edge.

The metagraph records which tournament types are "adjacent" — reachable from each other by a single arrow reversal.

---

## The Hamiltonian Path Count as a Height

Each tournament type has a specific H(T) value — the number of Hamiltonian paths. We can use this as a **height** in the map: types with few paths are "low," types with many paths are "high."

Does the metagraph have a consistent direction? Is it always possible to increase H(T) by flipping one arrow?

Mostly yes, but not perfectly:
- The vast majority of edges in the metagraph go from lower H to higher H (or vice versa)
- But there are also **level edges** — connecting tournament types with the same H — and even **downward edges** — where flipping an arrow decreases H
- At n = 7: 136 level edges (same H at both ends), and 962 downward edges

So the metagraph has a strong but imperfect "uphill" tendency. It's not a perfectly ordered hierarchy, but it leans strongly in that direction.

---

## Self-Complementary Classes: The Spine

Among all tournament types, the **self-complementary (SC)** ones are special (see explainer 5). In the metagraph, they form a **spine** — a skeleton that runs from the simplest tournament (the transitive one, with exactly 1 path) up to the richest tournaments.

The SC types are connected to each other by SC-to-SC edges — these form the backbone of the metagraph.

### Ribs and Sea

Attached to the spine are **ribs**: edges connecting SC types to non-SC (NS) types. These always connect across the SC/NS divide — never two SC types together (that would be the spine) and never two NS types (that would be the sea).

The **sea** is everything else: the enormous bulk of NS-to-NS connections. At large n, over 96% of all edges are in the sea. The spine is a sparse skeleton; most of the action is in the sea.

### The Geometry

If you draw the metagraph with H(T) as the vertical axis and the SC/NS distinction as the horizontal axis, you see:

- A vertical spine (SC types, varying H)
- Horizontal ribs (SC-to-NS connections)
- A vast bulk of sea (NS-to-NS connections)

This architecture — spine, ribs, sea — is proved to hold for all n, and understanding it is one of the major structural results of the project.

---

## Three Proved Theorems About the Map

1. **No triangles cross the SC/NS boundary twice**: If you take any triangle in the metagraph (three types all connected to each other), you can't have two edges crossing between SC and NS. Ribs can't form triangles.

2. **The rib subgraph is bipartite**: If you look only at the rib edges (SC to NS), the resulting picture has no odd cycles — it's perfectly two-colored. Bipartite means the SC and NS types separate cleanly.

3. **Odd cycles in the metagraph that cross the SC/NS boundary an odd number of times contribute nothing to the trace**: A technical algebraic statement that has deep implications for what spectral properties the metagraph can have.

---

## The Merged Map

If you take any tournament T and complement it (reverse all arrows), you get a related tournament T^c. For non-SC types, T and T^c are always different types. This creates a **pairing** among the non-SC types.

If you collapse this pairing — treating T and T^c as the same type — you get the **merged metagraph** G_n/Z₂. This has:

- Exactly one vertex per SC type (unchanged)
- Exactly one vertex per NS type-pair {T, T^c}

At n = 7: the merged map has 456 types → 228 SC + 456 NS → 228 SC + 228 pairs = 285 merged vertices.

The merged map is the "natural" object — it factors out the complement symmetry that's always present.

---

## What the Map Reveals

### Connectivity

All tournament types of the same size are connected in the metagraph — you can get from any tournament to any other by a sequence of single arrow flips. (This isn't obvious! It was a theorem to establish.)

### H-gradient

The consistent upward bias in the metagraph tells you that optimizing a tournament (finding the one with the most Hamiltonian paths) can be done by a "greedy" hill-climbing algorithm: flip whichever arrow increases H, repeat. This generally works, though the landscape isn't perfectly "smooth" — there are local plateaus and rare descents.

### The Transitive Tournament as Source

The transitive tournament (complete ranking, H = 1) sits at the very bottom of the metagraph. It has exactly 1 Hamiltonian path — the minimum possible. From this single point, all other tournament types are reachable by flipping arrows.

The transitive tournament is the "seed" of the entire universe of tournaments.

---

## Computational Scale

Building the metagraph requires:
- Generating all tournaments on n players
- Grouping them into isomorphism classes
- Computing which classes are adjacent (differ by one arc flip)

This is feasible up to n ≈ 8 (millions of types). At n = 9 and beyond, it requires clever sampling and approximation.

Most of the large-scale computations in this project are building, analyzing, and understanding different aspects of the metagraph and its neighbors.

---

## Key Words

- **Metagraph G_n**: a graph where vertices are tournament types and edges connect types differing by one arc flip
- **Spine**: the SC-to-SC subgraph — the backbone of the metagraph
- **Ribs**: SC-to-NS edges — connecting the spine to the bulk
- **Sea**: NS-to-NS edges — the dominant bulk of connections at large n
- **Merged metagraph G_n/Z₂**: the metagraph with complement pairs identified — the natural quotient
- **H-gradient**: the tendency of the metagraph to increase H as you move "uphill"
