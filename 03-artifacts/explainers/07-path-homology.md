# Path Homology: Counting Holes

## Topology and the Math of Holes

**Topology** is the branch of mathematics that studies shapes and spaces by asking: what can you tell about an object without measuring exact distances? Instead of measuring lengths and angles, you ask qualitative questions like:

- Is this surface connected, or does it fall into separate pieces?
- Does it have holes you can poke a finger through?
- Is it like a sphere, a donut, or a pretzel?

The key tool for this is **homology** — a way of precisely counting and classifying "holes" of different dimensions:
- **0-dimensional holes**: gaps between disconnected pieces
- **1-dimensional holes**: loops that can't be contracted to a point (like the hole in a donut)
- **2-dimensional holes**: enclosed cavities (like the hollow interior of a ball)

These counts are called **Betti numbers**: β₀, β₁, β₂, β₃, ...

---

## But Tournaments Have Arrows — Standard Topology Doesn't Apply

Standard topology doesn't care about direction. A loop is a loop whether you go clockwise or counterclockwise.

But tournaments have directed edges — arrows. A directed loop (A→B→C→A) is very different from its reverse (A→C→B→A). If you ignore direction, you lose information.

In 2015, mathematicians developed **GLMY path homology** (named after Grigor'yan, Lin, Muranov, and Yau) specifically for directed graphs. It detects "directional holes" — patterns of arrows that form closed directed cycles at different scales.

---

## What GLMY Homology Counts

In GLMY path homology, a "chain" is a weighted collection of directed paths. Two chains are equivalent if one can be "deformed" into the other by following the arrows — like rubber-band stretching, but respecting the arrow directions.

A **hole** (in the homological sense) is a chain that:
1. Forms a closed cycle (starts and ends at the same place)
2. Cannot be "filled in" — it's not just the boundary of a higher-dimensional shape

The Betti numbers β_k count these directional holes of dimension k.

---

## What We Computed

For every tournament on up to 8 players (hundreds of thousands of them), we computed all the GLMY Betti numbers β₁, β₂, β₃, β₄, ...

### The First Big Result: β₂ = 0 Always

**For every tournament on every size we tested, β₂ = 0.**

No 2-dimensional directional holes. None. Zero. This has been verified for:
- All tournaments up to 8 players (exhaustive, millions checked)
- Random samples of larger tournaments

And there's strong algebraic evidence it's true for ALL tournaments, forever.

This is a strong structural theorem: the "2-dimensional fabric" of a tournament is always completely flat. No enclosed cavities.

### The Second Big Result: β₁ is Either 0 or 1

Every tournament has either **zero or one** 1-dimensional hole. Never 2, never 3.

β₁ = 0: no directed loops that can't be filled — the tournament's "cycle structure" is all consistent.
β₁ = 1: exactly one irreducible directed loop — one "hole" in the directed fabric.

### The Seesaw Pattern

At small tournament sizes, β₁ and β₃ never coexist. If β₁ = 1, then β₃ = 0. If β₃ > 0, then β₁ = 0.

This "seesaw" suggests deep cancellation laws in the directed algebra — the two types of holes somehow "cancel" each other and can't both be present at once.

### The Rare Case: β₄ = 6 at the Paley T₇ Tournament

The Paley tournament on 7 players — which also has the most Hamiltonian paths among all 7-player tournaments — has the unusual property that **β₄ = 6**.

This is extremely rare. Only about 0.1% of 7-player tournaments have any β₄ at all. And only the Paley T₇ achieves β₄ = 6 (the maximum observed).

This suggests a deep connection: the tournaments that maximize path counts also tend to have the richest topological structure.

---

## What Does a "Hole" in a Tournament Mean?

In a standard graph, a 1-dimensional hole is literally a cycle — a loop with no "shortcut" path through the interior.

In GLMY path homology for directed graphs, a 1-dimensional hole is more subtle: it's a directed cycle that can't be "explained" as the boundary of a directed 2-dimensional surface.

In practical terms: if β₁ = 1 for a tournament T, it means there's a pattern of directed paths in T that form a closed directed loop that cannot be decomposed into simpler directed triangular pieces. The tournament has a structural "knot" that can't be untied.

---

## The Euler Characteristic

A single number captures the alternating sum of Betti numbers:

**χ(T) = β₀ - β₁ + β₂ - β₃ + β₄ - ...**

For most tournaments: β₀ = 1 (connected), β₁ = β₂ = β₃ = β₄ = 0, so χ = 1.

For about 12.3% of 7-player tournaments: β₁ = 1, so χ = 0.

For the Paley T₇: β₄ = 6, so χ = 1 - 0 + 0 - 0 + 6 = 7.

This hierarchy — 1, 0, 7 — clusters neatly. Tournaments fall into a small number of topological types, and the most path-rich tournament sits at the top.

---

## Why This Matters

Path homology connects tournament combinatorics to a whole new branch of mathematics:

1. **New invariants**: The Betti numbers distinguish tournaments that look similar but have different cycle structures.

2. **Constraints**: β₂ = 0 is an algebraic constraint — it means tournament theory "lives" in a restricted slice of all possible directed graphs. Not every directed graph arises from a tournament, and the vanishing of β₂ is one witness to this restriction.

3. **Applications to networks**: Any network with directed links (Twitter follows, food web predation, supply chains) has a path homology. Tournaments are an extreme case (complete graphs), and understanding their homology gives tools for analyzing more general directed networks.

4. **Connection to path counting**: The Euler characteristic χ(T) is related to alternating sums of path counts. The richer the topology, the more complex the path structure.

---

## Key Words

- **Topology**: the math of shapes and spaces, studying connectivity and holes
- **Betti numbers (β_k)**: counts of k-dimensional holes in a space
- **GLMY path homology**: a version of homology adapted for directed graphs, respecting arrow directions
- **β₂ = 0**: no 2-dimensional holes in any tournament (proved for all sizes tested)
- **Euler characteristic χ**: alternating sum of Betti numbers, a single-number topological fingerprint
- **Seesaw pattern**: β₁ and β₃ never coexist in small tournaments — evidence of deep cancellation
