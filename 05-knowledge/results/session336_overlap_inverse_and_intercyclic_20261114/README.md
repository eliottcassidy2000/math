# Session 336: Invert the target, then choose the overlap

**Chan–Pak remains open.** The session gives a complete exact search in the two-hub class, three direct arithmetic inverses on unbounded families, and an all-order classification of the no-disjoint-triangles strip.

1. [All-order classification](INTERCYCLIC_CLASSIFICATION.md). For each order $n\ge6$, there are $2^{n-3}+n-3$ strong tournaments with no two disjoint directed triangles, up to isomorphism. The $2^{n-3}-1$ with a one-vertex triangle transversal explain the 127 at order ten. The other $n-2$ form explicit families. At $n=10$, six classes have Hamiltonian count **385** and two have **573**. The finite-state language certificate proves this at every order. There is no genuine triple overlap in this subclass.
2. [Exact target-first selector](invert/TARGET_OVERLAP_INVERSE.md). The backward bridge intervals in any strong two-hub tournament force at least $2^{m-3}$ paths, where $m$ is spine length. This supplies a proven finite search cap for a prescribed integer, not a cutoff extrapolated from coverage. At 385 it certifies maximum triple overlap **2** with witness **1030221**; at 573 maximum **3** with witness **20302212**. In both cases the best double overlap at that triple maximum is 1. These are maxima over all two-hub words, not merely over order ten. Their disjoint triangles are a necessary cost of choosing those overlaps.
3. [Two direct target rails](TARGET_RAILS.md). Stretching the exact witnesses produces two mixed-count families with separately controllable binary and ternary overlap. Each bounded list of ternary lengths gives an exact power-of-two test for the required binary length. The 385-based family leaves the offset-one smooth rail everywhere except at its seed.
4. [The three 1/3 comparisons and a third direct inverse](cubic/BALANCE_CUBIC_OVERLAP.md). A cubic bipartite graph of girth at least six can keep the same 36 vertices and a 12-vertex dominating set while projecting to seven different role-list frames. Domination does not determine their overlap. The literal words witnessing those lists lead to an unbounded mixed-count family with at most three candidate patterns per target, selected by $v_2(N-1)$. The poset analogy is exact inside a locked frame but not for the union of all paths on the original tournament vertices.
5. [Whole-path balanced pairs](GLOBAL_BALANCE.md). The hub-free 385 and 573 families have whole-path probabilities $256/385$ and $382/573=2/3$ on a specific named pair. The two designated hubs in the first family can be arbitrarily biased. Neither type of balance supplies triple overlap.

## Replay

- Run **python replay_intercyclic.py**, **python replay_balance.py**, and **python replay_target_rails.py** here.
- In **invert/**, build **g++ -O3 inverse.cpp -o inverse**, run **./inverse 385**, **./inverse 385 intercyclic**, **./inverse 573**, **./inverse 573 intercyclic**, and **python replay_order_masks.py**.
- In **cubic/**, run **python replay_cubic_balance.py**.

The remaining obligation is not to check more targets. It is to find a structural rule guaranteeing that every admissible absolute target enters one of the direct inverse families, or a new family with the right overlap. Neither local 1/3 balance nor a small undirected dominating set makes that guarantee.
