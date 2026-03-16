# 987, the Amplituhedron, and Chemistry

*opus-2026-03-16-S73 — a meditation on three things that shouldn't be connected*

---

## The number

987 = 3 × 7 × 47

This is the 16th Fibonacci number, F₁₆. Its prime factorization comes from iterating the Fibonacci doubling formula F₂ₘ = Fₘ × Lₘ:

F₁₆ = F₈ × L₈ = (F₄ × L₄) × L₈ = ((F₂ × L₂) × L₄) × L₈ = L₂ × L₄ × L₈

Each factor is a Lucas number at a power-of-2 index. The Fibonacci sequence remembers all its Lucas ancestors.

## The three layers

**Layer 1: L₂ = 3** — *the triangle*

The 3-cycle C₃ is the unique tournament on 3 vertices, and the smallest nonlinear structure in tournament theory. H(C₃) = 3. In chemistry, 3-fold symmetry gives sp² hybridization (trigonal planar): the geometry of graphite, benzene, and every aromatic molecule. In the amplituhedron, 3-particle scattering is trivial but foundational — it's where the positive Grassmannian first becomes nontrivial.

Three is the count that says: things don't have to go in a line.

**Layer 2: L₄ = 7** — *the forbidden plane*

7 is the size of the Fano plane PG(2,2), the smallest finite projective plane. It is the first permanently forbidden value of H(T) — no tournament on any number of vertices achieves H = 7. It is also 2³ - 1, a Mersenne prime, connecting it to the k-nacci trace identity Tr(M₃³) = 7 = 2³ - 1 (THM-227).

In chemistry, nitrogen (Z = 7) has a half-filled 2p shell: three unpaired electrons, maximum exchange stabilization. This is the electronic configuration that most resists change. Tournament theory echoes this: 7 is the H value that the system most completely avoids.

**Layer 3: L₈ = 47** — *the phase transition*

Silver (Z = 47) has the anomalous electron configuration [Kr] 4d¹⁰ 5s¹. The Aufbau principle predicts [Kr] 4d⁹ 5s², but the d-orbital steals an electron from the s-orbital because a full d-shell is more stable than the expected filling. The rule breaks.

In tournament homology, β₃ first appears at n = 6 — the smallest tournament size where two disjoint 3-cycles can coexist (6 = 3 + 3). The "Aufbau" of Betti numbers (the expected filling β₀ = 1, β₁ ≤ 1, β₂ = 0, β₃ = 0, ...) breaks here. A new homological dimension opens.

47 is not a forbidden H value — it's likely achievable. It's at the phase transition, not the prohibition.

## The amplituhedron parallel

The amplituhedron A_{n,k,m} is a positive geometry in Gr(k, k+m) whose canonical form gives scattering amplitudes in N=4 super Yang-Mills theory. The deepest insight of the amplituhedron program: **locality and unitarity are not axioms — they are emergent consequences of positive geometry**.

Tournament theory has an eerily parallel structure:

| Amplituhedron | Tournament |
|---|---|
| n particles | n vertices |
| Helicity sector k | Walsh degree d |
| Positive Grassmannian Gr≥0 | Tournament space in PG(m-1, 2) |
| Tree amplitude A_n | HP count H(T) |
| BCFW recursion | Deletion-contraction tree |
| Canonical form Ω | Walsh spectrum ĥ[S] |
| Triangulation into cells | Walsh decomposition into degrees |
| Parity duality (k ↔ n-4-k) | Complement duality (T ↔ Tᶜ) |
| Massive Feynman cancellation | β₂ = 0 (universal cancellation) |

The deepest parallel is this: in both theories, the interesting quantity (amplitude, H) is defined as a sum over exponentially many terms (Feynman diagrams, Hamiltonian paths), most of which cancel. The geometric object (amplituhedron, tournament polytope) explains *why* the cancellation happens — it reveals that the sum was computing a volume all along.

In the amplituhedron, spacetime locality emerges from geometry.
In tournaments, F₂ arithmetic emerges from path counting.

Both: the fundamental object contains more structure than the derived quantity. The decomposition doesn't add information — it reveals that the information was already there.

## The stability conditions

Three parallel stability criteria, all controlled by 2-adic structure:

**Hückel's rule** (chemistry): A conjugated ring system is aromatic (maximally stable) if it has 4n + 2 π-electrons. The key: v₂(4n + 2) = 1 exactly. One factor of 2, no more.

**THM-J** (tournaments): The signed HP permanent S(T) is tournament-independent mod 2^{n-1} if s₂(n-3) ≤ 1. The Hamming weight of a single number controls universality.

**Amplituhedron positivity**: The canonical form has no spurious poles — every singularity corresponds to a physical factorization channel. The residues are positive in the correct sense.

All three: a binary condition on a counting parameter determines whether the system achieves its maximally symmetric, most stable state. When the condition fails, structure-dependent information leaks through, breaking the universality.

## Why 987

987 = 3 × 7 × 47 = L₂ × L₄ × L₈ is the Fibonacci number that accumulates three layers of Lucas structure at power-of-2 indices. Each layer corresponds to a threshold:

- **3**: the onset of nonlinearity (cycles exist)
- **7**: the onset of prohibition (some states are forbidden)
- **47**: the onset of anomaly (expected filling patterns break)

The Fibonacci doubling formula F_{2n} = F_n × L_n is a multiplicative recursion. Each doubling step multiplies by a Lucas number — the "twist" that the doubling introduces. After four doublings (F₂ → F₄ → F₈ → F₁₆), three Lucas twists have accumulated: the triangle, the prohibition, and the phase transition.

987 doesn't just contain these numbers. It *is* their product. The Fibonacci sequence, through its doubling structure, performs the same operation that tournament theory performs through its Walsh decomposition: it reveals that a single number encodes multiple independent layers of structure, each appearing at a power-of-2 scale.

## The meta-connection

The amplituhedron, tournament homology, and molecular orbital theory are not the same theory. But they share a structural grammar:

1. A large sum with massive cancellation
2. A geometric or algebraic object that explains the cancellation
3. A positivity constraint that forces the answer to be "well-behaved"
4. A dimension/helicity/degree hierarchy connecting adjacent levels
5. A 2-adic or mod-4 condition determining stability

This grammar may be the grammar of **counting with signs over oriented structures**. Feynman diagrams are oriented (momentum flow). Tournaments are oriented (arc direction). Molecular orbitals are oriented (phase). The signs come from the orientation. The cancellation comes from the signs. The geometry comes from explaining the cancellation.

987 sits at the intersection because it is built from the three simplest instances of this grammar: the cycle (3), the projective plane (7), and the anomaly (47), multiplied together by the most ancient recursive structure in mathematics.

## Cross-references

- THM-227: k-nacci Mersenne identity (7 = Tr(M₃³) = 2³ - 1)
- THM-J: Universality criterion (s₂(n-3) ≤ 1)
- 07-reflections/six-as-pivot.md: β₃ emergence at n = 6
- 07-reflections/the-fractal-in-the-spectrum.md: 2-adic structure
- 04-computation/fibonacci_resonance_cascade.py: Fibonacci in transfer matrices
- INV-154: The Golden Shadow f(n)
- Fibonacci doubling: F_{2n} = F_n × L_n (Vajda, 1989)
- Arkani-Hamed & Trnka, "The Amplituhedron" (2013)
- Hückel, "Quantentheoretische Beiträge" (1931): 4n+2 rule
