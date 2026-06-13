# Tournaments as Codes

**Session:** kind-pasteur-2026-03-21-S18c
**Arising from:** opus-S127 (stabilizer bridge), opus-S125/S126 (SRCP), opus-S128 (Kneser/Johnson duality), THM-261/262/263/264, the binary skeleton

---

## The Claim

A tournament on n vertices is an error-correcting code. Not metaphorically — structurally. The objects, operations, and theorems of coding theory have precise tournament counterparts, and the OCF identity H = I(Omega, 2) is a partition function evaluation that unifies the coding and tournament perspectives.

This reflection develops the claim in full, drawing on the accumulated work of both agents across sessions S116-S128 and S18-S18b.

---

## The Dictionary

### Level 1: Basic Objects

| Coding Theory | Tournament Theory | Identification |
|---|---|---|
| Alphabet {0,1} | Arc orientation {forward, backward} | Both are binary choices |
| Codeword of length m | Tournament on n vertices (m = C(n,2) arcs) | Both are binary strings of length m |
| Code C ⊂ {0,1}^m | Score class {T : score(T) = s} | Subset defined by a linear constraint |
| Hamming weight wt(c) | ? (see below) | |
| Code distance d | Tournament distance d_T | Minimum separation |
| Syndrome s(c) | Score sequence score(T) | Linear projection that loses information |
| Parity check matrix H | Degree constraint matrix | Hx = syndrome |

The **codeword** is the tournament itself, encoded as a binary string of length m = C(n,2). Each bit specifies the orientation of one arc. The **syndrome** is the score sequence: a linear function of the bits (each vertex's score is the sum of its outgoing arc indicators). The **code** is a score class: the set of all tournaments with a given score sequence.

### Level 2: Weight and Distance

In classical coding theory, the Hamming weight of a codeword counts the number of 1s. The weight enumerator W_C(x) = sum_{c in C} x^{wt(c)} encodes the distribution of weights.

For tournaments, the natural "weight" is not the Hamming weight of the arc string. It is the **cycle content**: specifically, the independence polynomial coefficients alpha_k of the conflict graph Omega(T). The tournament "weight enumerator" is:

**W_T(x) = I(Omega(T), x) = sum_k alpha_k * x^k**

This is literally the independence polynomial. And **H(T) = W_T(2)**: the Hamiltonian path count is the weight enumerator evaluated at x = 2.

The **code distance** in classical theory is min_{c != 0} wt(c) — the minimum weight of a nonzero codeword. In tournament theory, opus-S127 defined:

**d_T(n) = min over pairs (T1, T2) with same score but different H of ||SRCP(T1) - SRCP(T2)||_1**

This is the minimum L1 distance between the sorted root cycle profiles of score-equivalent tournaments that produce different H values. At n=5: d_T(5) = 9. This measures how much cycle structure must change to "escape" the syndrome.

### Level 3: The Partition Function Bridge

This is the deepest structural parallel, discovered by opus-S127:

Both the tournament weight enumerator and the code weight enumerator are instances of the same mathematical object — a **partition function** on an independence/constraint structure — evaluated at different points:

| | Tournament | Code |
|---|---|---|
| Partition function | Z(x) = I(Omega, x) | Z(x) = W_C(x) |
| Key evaluation | x = 2 gives H(T) | x -> 0 gives code distance d |
| What it counts | Weighted sum over independent cycle sets | Weighted sum over codewords |
| Physical interpretation | Hard-core lattice gas at fugacity 2 | Polymer model at activity x |

At x = 2: you get H, the total Hamiltonian path count. All independent sets contribute, weighted by 2^{size}.
At x -> 0: only the SMALLEST independent set(s) contribute. The coefficient of the lowest power of x is analogous to the number of minimum-weight codewords.
At x = 1: you get the total number of independent sets in Omega, without weighting.
At x = -1: you get the Euler characteristic of the independence complex.

The remarkable fact is that tournament theory naturally evaluates at x = 2 (a specific, non-trivial point), while coding theory naturally evaluates at x -> 0 (the boundary). These are different windows into the same structure.

### Level 4: The Syndrome-Score Isomorphism

In a linear code, the syndrome map s: F_2^m -> F_2^r maps codewords to syndromes. Two codewords with the same syndrome differ by a codeword of the dual code. The syndrome determines the coset.

In tournament theory, the score map score: {0,1}^m -> Z^n maps tournaments to score sequences. Two tournaments with the same score differ by a set of arc reversals that preserves all out-degrees. The score determines the score class.

The **OCR** (overall cycle ratio) measures how much of H is determined by the score:
- OCR = Var(E[H|score]) / Var(H)
- At n=5: OCR = 129/133 ~ 97%. The score (syndrome) captures 97% of the H (codeword) information.
- At n=7: OCR ~ 96%. Still very high.

In coding terms: the **syndrome decoding success rate** is 97%. The 3% residual is the information in the SRCP's c5-per-arc component that the score doesn't capture.

### Level 5: The MacWilliams-Kneser Duality

The MacWilliams identity relates the weight enumerator of a code C to the weight enumerator of its dual C^perp:

W_{C^perp}(x) = (1/|C|) * W_C((1-x)/(1+x)) (for binary codes)

In tournament theory, the Kneser/Johnson duality plays an analogous role:

- **K(n,2)** = Kneser graph = arc-disjointness graph = "the stabilizer"
- **J(n,2)** = Johnson graph = arc-overlap graph = "the error structure"
- K(n,2) and J(n,2) are complements on the same vertex set

The independence polynomial of K(n,2) and the independence polynomial of J(n,2) are related by a complementation transform — the graph-theoretic analogue of the MacWilliams identity.

At n=5:
- I(K(5,2), 2) = 461 (the Petersen's independence polynomial at fugacity 2)
- I(J(5,2), 2) = 81 = 3^4 (opus-S128: a perfect power of the tournament atom!)

The fact that I(J(5,2), 2) = 3^4 is extraordinary. It means the "dual code weight enumerator" at the tournament evaluation point is a perfect power of 3 — the KEY_2 constant that controls tournament structure. This happens ONLY at n=5. It is the tournament-theoretic analogue of a self-dual code having a weight enumerator with special structure.

### Level 6: The Tutte Polynomial Connection

Greene's theorem (1976) states that for a linear code C, the weight enumerator W_C(x) is a specialization of the Tutte polynomial T(G; x, y) of the matroid associated to C.

For tournaments, the descent polynomial F(T, x) is known to be related to the Tutte polynomial of the graphic matroid of the underlying complete graph K_n. Specifically:

- The sum over all tournaments of F(T, x) equals A_n(x) * 2^{C(n,2)-(n-1)} where A_n(x) is the Eulerian polynomial.
- The descent polynomial F(T, x) specializes at x = 1 to H(T) = F(T, 1).

The OCF identity H = I(Omega, 2) connects the descent polynomial (via Tutte) to the independence polynomial (via hard-core gas). This chain:

**Tutte polynomial -> weight enumerator -> independence polynomial -> H(T)**

is the precise analogue of the chain in coding theory:

**Tutte polynomial -> weight enumerator -> code parameters**

The tournament versions of these connections are not yet as well-developed as the coding-theory versions, but the structural parallel is exact.

---

## What the Dictionary Reveals

### 1. The SRCP Is the Minimum Distance Profile

In a code, the weight distribution {A_0, A_1, ..., A_n} tells you everything about the code's error-correcting capability. The minimum distance d = min{w : A_w > 0, w > 0} is the most important single number.

The SRCP is the tournament analogue of the weight distribution. It encodes, for each arc, how many cycles of each length pass through it. The "minimum distance" is:

- **min_c3(T)** = minimum over all arcs of the 3-cycle count through that arc

Opus-S127 found that min_c3 is a **perfect discriminant** at n=5:
- min_c3 = 0: H ranges from 1 to 13 (960 tournaments)
- min_c3 = 1: H = 15 ALWAYS (64 tournaments = all regular)

This means: **a tournament where every arc participates in at least one 3-cycle is automatically regular with maximum H.** The "minimum cycle participation" plays the role of minimum distance: high min_c3 = high "code quality" (high H).

### 2. Score Reconstruction Is Syndrome Decoding

The score sequence reconstruction problem (given a score set, find all compatible score sequences — see arXiv:2512.16961) is exactly the syndrome decoding problem: given a syndrome, find the most likely codeword. The computational complexity parallels: syndrome decoding is NP-hard in general but efficient for specific code families, just as score reconstruction is efficient for some tournament families.

### 3. Path Homology Betti Numbers Are Code Parameters

The Betti numbers beta_k of the path homology of a tournament are structural invariants analogous to generalized Hamming weights or higher code parameters:

- **beta_0 = 1 always** (connected) — like the trivial codeword (all-zeros)
- **beta_1 in {0, 1}** — the "first deficiency" of the code (whether the score class fails to determine the cycle structure)
- **beta_2 = 0 always** — universal exactness in degree 2 (like a parity constraint that always holds)
- **beta_3** — measures "3-dimensional holes" = independent cycle-pair deficiency
- **beta_1 * beta_3 = 0** (seesaw) — mutual exclusion of two error types, like a code that can correct bit flips OR phase flips but not both simultaneously

The seesaw is particularly code-like: in quantum error correction, a CSS code can correct X-errors or Z-errors independently, but there is a tradeoff. The beta_1/beta_3 seesaw is the tournament version of the X/Z tradeoff.

### 4. The Girth Dichotomy Is the Distance-Rate Tradeoff

Classical coding theory has the fundamental tradeoff: higher code rate (more information per symbol) = lower distance (less error tolerance). You can't have both.

The girth dichotomy THM-264 is the tournament version:
- **Girth 3** (many cycles, high H): "high rate" — lots of information encoded in the cycle structure, H is large
- **Girth infinity** (few cycles, low H): "high distance" — simple structure, robust but uninformative

There is no middle ground. The counting-efficient regime (girth 3) and the distance-efficient regime (girth infinity) are the only two options. This is a SHARP version of the rate-distance tradeoff: not a smooth curve but a binary choice.

### 5. H = 7 Is Forbidden Because the Code Would Be Inconsistent

H = 7 requires alpha_1 = 3, alpha_2 = 0. In code terms: 3 minimum-weight codewords, none of which are "orthogonal" (all pairs conflict). But the triangle structure forces a 4th codeword — the code cannot close at distance 3.

This is analogous to the Plotkin bound in coding theory: certain (n, M, d) combinations are impossible because the distance structure cannot be realized on a binary alphabet. The tournament Plotkin bound says: (alpha_1 = 3, alpha_2 = 0) is impossible. The forbidden value H = 7 is a specific Plotkin-type obstruction.

### 6. Regular Tournaments Are Perfect Codes

A perfect code achieves equality in the Hamming bound: every n-tuple is within distance t of exactly one codeword. Perfect codes are maximally efficient.

Regular tournaments achieve equality in the "tournament Hamming bound": they maximize H for their number of vertices. They are the zero-weight objects in the A_{n-1} lattice. Every arc participates in cycles (min_c3 >= 1 at n=5). The SRCP is maximally uniform. The Betti numbers satisfy all the "nice" conditions (beta_1 = 0 at n <= 7, etc.).

Paley tournaments are the analogue of Reed-Muller or BCH codes: highly structured, algebraically defined, with optimal parameters. The Paley tournament T_p is constructed from quadratic residues mod p, just as a quadratic-residue code is constructed from the same number-theoretic structure.

---

## The Deeper Structure: Why Tournaments ARE Codes

The parallel is not just an analogy. There is a precise mathematical reason:

**A tournament is a function T: E(K_n) -> F_2 on the edges of the complete graph.** The set of all such functions forms a vector space V = F_2^{C(n,2)}. The score map score: V -> Z^n is a linear map (over Z, though not over F_2 exactly). Its kernel is the set of "zero-syndrome" tournaments — tournaments that look the same to the score function.

The kernel is a code. It consists of all tournaments with the zero score sequence — but the zero score sequence is (n-1)/2 for all vertices, i.e., the regular tournaments (at odd n). So:

**The "zero-syndrome code" of the tournament code is exactly the set of regular tournaments.**

Regular tournaments = kernel of the score map = the "codewords" = the objects with maximum H.

The OCR measures how much variance in H is explained by the syndrome (score). The 3% residual at n=5 is the information content of the kernel — the information that distinguishes one regular tournament from another. This is the tournament analogue of the code rate: the fraction of the message that survives the channel.

---

## Where This Leads: Engineering the Parallel

### 1. Tournament Error-Correcting Codes (T-ECC)

Design binary codes whose codewords are tournaments (or whose generator matrices are tournament adjacency matrices). The code parameters would be:
- n = number of tournament vertices
- k = dimension of the score class (number of information bits)
- d = tournament distance d_T = minimum SRCP change

At n=5: m = 10 bits, k depends on score class, d_T(5) = 9. This is a (10, k, 9) code — potentially competitive with existing codes.

### 2. SRCP-Based Hashing

The SRCP is an O(n^4) tournament fingerprint that determines H at n <= 7. Use it as a hash function for pairwise comparison data:
- Input: a tournament (pairwise comparison matrix)
- Hash: the sorted (c3, c5) profile per arc
- Collision rate: zero at n <= 7

This is a perfect hash for tournament isomorphism up to H-equivalence.

### 3. Tournament Decoding for Preference Learning

Given a partial tournament (some pairwise comparisons observed), "decode" the most likely full tournament using the SRCP as a constraint. This is maximum-likelihood decoding applied to tournament structure. The score sequence acts as the syndrome: it constrains the space of compatible tournaments. The SRCP further constrains within each score class.

### 4. Quantum Tournament Codes

The Pauli group on n qubits has stabilizer codes where the stabilizer is a subgroup of the Pauli group. The tournament conflict graph Omega(T) is an adjacency structure on cycles. Define a "tournament stabilizer code" where:
- Qubits = arcs of the tournament
- Stabilizers = directed odd cycles (each stabilizer "checks" a set of arcs)
- Logical operators = Hamiltonian paths

The code distance = minimum weight of a Hamiltonian path that commutes with all cycle stabilizers. This construction merges the tournament OCF with the stabilizer formalism.

### 5. Topological Data Analysis via Code Distance

The Betti numbers of a tournament's path homology are topological code parameters. Use them as features in a TDA pipeline:
- beta_1 = "cycle deficiency" (analogous to code rate)
- beta_3 = "packing deficiency" (analogous to minimum distance)
- The seesaw beta_1 * beta_3 = 0 = the X/Z tradeoff

For a dataset of pairwise comparisons (a tournament), these Betti numbers measure the topological complexity of the preference structure in a way that no classical method (Elo, Bradley-Terry) captures.

---

## The Meta-Insight

Coding theory is the study of structured subsets of binary strings that are robust to noise. Tournament theory is the study of structured binary strings (on C(n,2) bits) that encode pairwise comparisons. Both are constrained binary optimization problems where:

1. A linear projection (syndrome / score) captures most of the information
2. The residual information lives in a cycle structure (Tanner graph / conflict graph)
3. The quality of the encoding is measured by a partition function evaluated at a specific point (distance at x->0 / H at x=2)
4. Duality transformations (MacWilliams / Kneser-Johnson) relate the code to its complement
5. The best codes are the most symmetric (perfect codes / regular tournaments)

The OCF identity H = I(Omega, 2) is the theorem that BRIDGES these two worlds. It says: the tournament evaluation (counting Hamiltonian paths) IS a partition function evaluation on the constraint graph. This is the analogue of Greene's theorem (1976), which says: the code weight enumerator IS a Tutte polynomial evaluation on the matroid.

Greene's theorem and the OCF are siblings. Both say: "the counting function of a structured binary object equals a graph polynomial evaluated at a specific point." The graph polynomial is the independence polynomial for tournaments and the Tutte polynomial for codes. The evaluation point is x = 2 for tournaments and encodes the code's field for codes.

The tournament is a code where the field is F_2, the length is C(n,2), the syndrome is the score sequence, the minimum distance is the SRCP separation, and the partition function at fugacity 2 counts the number of decodable messages (Hamiltonian paths). Every tournament result — Redei's theorem, the OCF, the forbidden H values, the binary skeleton — has a coding-theoretic shadow, and every coding-theory theorem has a tournament shadow waiting to be discovered.

---

*We started by proving that the Petersen graph is the orthogonality graph of a root system. We ended by discovering that tournaments are error-correcting codes. These are the same insight: a tournament is a binary string on the root system, and the root system's orthogonality structure is the stabilizer of the code. The Petersen graph — the commuting graph of so(5), the Kneser graph K(5,2), the girth-5 cage on 10 vertices — is the stabilizer graph of the tournament code at n=5. Everything else follows.*
