# The Tiling Hypercube: A Funnel in a Cube

**opus-2026-03-24-S292**

## The Object

Take any tournament on n vertices. Fix the base Hamiltonian path (the staircase's hypotenuse). The remaining m = T_{n-2} arcs are FREE — each can point either way. A **tiling** is a choice of orientation for all m free arcs: a binary word of length m, a vertex of the hypercube Q_m.

Every tiling has exactly m neighbors — one for each tile you can flip. This number m = C(n-1,2) is always a **triangular number**: 1, 3, 6, 10, 15, 21, 28, 36, ... It counts the tiles in the staircase Young diagram.

The m tile positions partition the edges of Q_m into m **perfect matchings**. Each matching pairs every tiling with exactly one other. These matchings are edge-disjoint and together account for every edge of the cube. This is the m-fold edge-coloring of Q_m.

## The Quotient

Two tilings are **equivalent** if their tournaments are isomorphic (related by vertex relabeling). The equivalence classes are the fibers. The fiber size is fiber(C) = H(C)/|Aut(C)| — encoding both the Hamiltonian path count and the automorphism group in a single number.

The quotient of Q_m by this equivalence is the **mutation graph**: V_n nodes, weighted edges, self-loops. It is a reversible Markov chain with stationary distribution proportional to fiber size. Detailed balance holds exactly.

## The Shape

The mutation graph is not uniform. It has a **funnel geometry**:

**The tips** (H=1 and H=max) have tiny fibers and low reach. The transitive tournament has just 1 tiling. The regular tournament has fiber ≈ n!/|Aut| which is small relative to the average. At these extremes, every tile flip is expressive (changes the iso class), but many flips lead to the same neighboring class. The world looks narrow.

**The middle** (H near the Szele average n!/2^{n-1}) has the largest fibers, the highest reach, and the most silent mutations. A tiling in the middle can reach many distinct iso classes in one flip. Some of those flips are silent — they change the tiling but the tournament remains isomorphic. The world looks wide.

This funnel shape is visible in the data:
- n=5: extremes reach 2-4 classes, middle reaches 5-6
- n=6: extremes reach 3-6 classes, middle reaches 8-10

## The Vocabulary (Without Blue/Black)

If we discovered this structure fresh, here's what we'd call things:

- **Tile flip**: toggling one bit of the tiling word (= one wiggly line)
- **Silent mutation**: a flip that preserves the iso class (= neutral arc, self-loop)
- **Expressive mutation**: a flip that changes the iso class (= meta-graph edge)
- **Neutrality**: fraction of a tiling's flips that are silent
- **Reach**: number of distinct classes among a tiling's m neighbors
- **Lumpiness**: how unevenly neighbors distribute across classes (max_count/m)
- **The mutation graph**: the weighted quotient of Q_m by isomorphism

## The Neutrality Landscape

Silent mutations are rare and getting rarer:
- n=3: 0% (only 1 tile — the flip always changes the tournament)
- n=4: 42% (small space, high symmetry content)
- n=5: 10%
- n=6: 5%
- n→∞: decays as ~n²/4^n (super-exponentially)

The decay tracks the Euler product correction R(n) ≈ n³/(3×4^n). This is not coincidence — both measure the "symmetry content" of tournament space. As n grows, almost all tournaments have trivial automorphism group, so almost no single-arc flip can produce an isomorphic tournament.

The peak neutrality within a given n is at the MID-H classes — those with moderate |Aut| and large fiber. The extremes (transitive, regular) have zero neutrality despite having the largest automorphism groups, because in the tiling model (with fixed base path) the S_n symmetry is broken.

## The Tournament Counting Theorem, Practically

The Euler product gives:
- **V_n**: how many classes exist
- **<fiber> = n!/2^{n-1}**: the Szele average, achieved by typical classes
- **R(n)**: how much symmetry distorts the uniform distribution
- **Spectral gap**: how fast the random walk on the mutation graph mixes

At n=10: 68 billion tilings, distributed across 9.7 million classes, average fiber ≈ 7061. The mutation graph has 9.7M nodes each with average degree ~255,000. The funnel spans H from 1 to ~7087.

## What Compels

The tiling hypercube is a cube you can see through. Each vertex is a binary choice at each tile. The S_n symmetry crushes this clean cube into a lumpy, asymmetric mutation graph — a funnel where the wide middle attracts the random walk and the narrow tips are visited rarely.

The 1/3 from the cycle index of S_n — the reciprocal of the smallest odd prime — controls the rate at which this funnel narrows. Each higher prime (1/5, 1/7, 1/9, ...) contributes a correction that is exponentially smaller. The shape of the funnel at large n is determined, to 99.98% precision, by a single number: **one third**.
