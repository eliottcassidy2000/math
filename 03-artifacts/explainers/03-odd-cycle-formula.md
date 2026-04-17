# The Odd-Cycle Formula

## The Central Discovery

The most important result in this project is a formula that lets you compute H(T) — the number of Hamiltonian paths — by studying something completely different: the **odd cycles** in the tournament and how they conflict with each other.

This connection is surprising. Paths and cycles seem like different things. The formula that links them is called the **Odd-Cycle Collection Formula (OCF)**.

---

## What is a Cycle?

A **cycle** in a tournament is a closed loop: you can start at a player, follow arrows, and return to where you started.

A **3-cycle** (the smallest possible) looks like this:

```
A → B → C → A
```

A **5-cycle** looks like:

```
A → B → C → D → E → A
```

Cycles with an odd number of players (3, 5, 7, ...) are called **odd cycles**.

---

## What is a Conflict?

Two odd cycles **conflict** if they share at least one player.

Imagine the 3-cycle (A, B, C) and the 3-cycle (B, D, E). They both involve player B, so they conflict. You can't "use" both of them together — they're fighting over B.

Two cycles that share NO players don't conflict. They can coexist peacefully.

---

## The Conflict Graph

Now here's a clever move: build a new graph called the **conflict graph** Ω(T).

- Each **odd cycle** in T becomes a **dot** in Ω(T)
- Two dots are **connected by a line** if the cycles conflict (share a player)

This turns the complex tournament into a simpler object: a plain undirected graph (no arrows, just dots and lines).

---

## Independent Sets

In any graph, an **independent set** is a collection of dots with NO lines between them. In the conflict graph, an independent set is a collection of odd cycles that are all pairwise vertex-disjoint — none of them share any players.

So an independent set of size k in Ω(T) = a way to find k completely non-overlapping odd cycles in T.

---

## The Formula

Here's the formula:

**H(T) = I(Ω(T), 2)**

Where I(G, x) is the **independence polynomial** of the conflict graph G, evaluated at x = 2.

What's the independence polynomial? It's defined as:

I(G, x) = 1 + (# independent sets of size 1)·x + (# of size 2)·x² + (# of size 3)·x³ + ...

When you plug in x = 2, you get:

H(T) = 1·1 + (# size-1 indep sets)·2 + (# size-2 indep sets)·4 + (# size-3 indep sets)·8 + ...

Each collection of k vertex-disjoint odd cycles contributes 2^k to the total path count.

---

## A Concrete Example

For a 3-player tournament with one 3-cycle (A→B→C→A):

- There is 1 odd cycle: the 3-cycle itself
- Ω(T) is a single dot (no lines, since there's only one cycle)
- Independent sets: one empty set (size 0) and one set containing the cycle (size 1)
- I(Ω(T), 2) = 1·(2^0) + 1·(2^1) = 1 + 2 = 3

And indeed, this tournament has 3 Hamiltonian paths. ✓

---

## Why This is Remarkable

Before this formula, computing H(T) meant searching through all possible permutations of the players — an approach that becomes impossibly slow for large tournaments (factorial growth). 

The OCF reduces the problem to studying odd cycles and their conflicts. Even better, it tells us WHAT determines H(T): it's entirely controlled by how many ways you can pack vertex-disjoint odd cycles into the tournament, with each packing of k cycles contributing 2^k to the total.

This is a complete description of where Hamiltonian paths come from — not by listing paths one by one, but by understanding the cycle structure.

---

## Who Proved It?

This formula was proved by mathematicians Grinberg and Stanley in 2023 (and made fully rigorous by Irving and Omar in 2024). The research project here independently discovered it computationally — verifying it by brute force for all tournaments up to 8 players — before learning the external proof existed.

The formula appears in an area of mathematics connected to **statistical physics**, specifically the "hard-core lattice gas model." In physics language, the odd cycles are "particles," and two particles that conflict can't both be present at the same location. The number H(T) is the total count of all possible particle configurations, weighted by how many particles you use. This unexpected connection to physics is one of the surprising payoffs of the work.

---

## Key Words

- **Odd cycle**: a closed loop in the tournament involving 3, 5, 7, ... players
- **Conflict**: two cycles conflict if they share a player
- **Conflict graph Ω(T)**: a graph where cycles are dots and shared-player conflicts are lines
- **Independent set**: a collection of dots in a graph with no lines between them — corresponds to vertex-disjoint cycles
- **Independence polynomial**: counts independent sets by size, evaluated at x=2 to give H(T)
- **OCF (Odd-Cycle Collection Formula)**: H(T) = I(Ω(T), 2)
