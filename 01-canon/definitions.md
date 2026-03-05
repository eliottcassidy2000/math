# Core Definitions

**All Claude instances must use these definitions exactly.** Do not redefine these terms locally in your work. If you believe a definition is wrong or incomplete, open a court case.

Last reviewed: SYSTEM-2026-03-05-S1

---

## Tournaments

**Tournament** T on vertex set V = {1, …, n}: a complete directed graph. For every pair {u,v} ⊂ V, exactly one of u→v or v→u holds. Encoded as T(u,v) = 1 if u→v, with T(u,v) + T(v,u) = 1 for all u ≠ v.

**Opposite tournament** T^op: reverses all arcs. T^op(u,v) = T(v,u).

**H(T)**: the number of directed Hamiltonian paths in T.

**Anti-automorphism** of T: a bijection α: V→V with T(α(u),α(v)) = T(v,u).

---

## The Tiling Model

**Base path** P_0: the fixed path n → n−1 → ⋯ → 1.

**Tile**: a pair (a,b) with a,b ∈ V and a ≥ b+2. There are m = C(n−1, 2) tiles.

**Tiling** t ∈ {0,1}^m: assigns a bit to each tile. Bit 0 means a→b (forward arc); bit 1 means b→a (backward arc). Each tiling determines a unique tournament T_t containing P_0.

**Pin grid**: Grid(n) := {(r,c) ∈ Z² : r ≥ 1, c ≥ 1, r+c ≤ n−1}, where r = a−b−1, c = b. This is isomorphic to the staircase Young diagram δ_{n-2}.

**Strip**: Str(k) := {(r,c) : r+c = k, r,c ≥ 1}, containing k−1 tiles.

---

## The Odd-Cycle Collection Formula

**Conflict graph** Ω(T): vertices are directed odd cycles of T; two cycles are adjacent iff they share at least one vertex.

**Independence polynomial** I(G, x) of a graph G: I(G, x) = Σ_{k≥0} α_k x^k, where α_k = number of independent sets of size k in G. Note α_0 = 1.

**μ(C)** for an odd cycle C ∋ v: μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2), i.e., the independence polynomial at 2 of the conflict graph of T−v restricted to cycles vertex-disjoint from C\{v}.

---

## The Vertex Deletion Setup

**T−v**: the tournament T with vertex v and all incident arcs removed.

**Claim A** (central open problem): H(T) − H(T−v) = 2 Σ_{C∋v} μ(C), summing over all directed odd cycles C through v.

**Claim B** (proved): I(Ω(T), 2) − I(Ω(T−v), 2) = 2 Σ_{C∋v} μ(C).

---

## The inshat/insact Framework

**P'**: a Hamiltonian path of T−v (i.e., a path through all vertices of V\{v}).

**inshat(v, P')**: the number of ways to insert v into P' to get a Hamiltonian path of T, counting with sign. Specifically, inshat = boundary_term + #Type-I + #Type-II positions.

**Signature** s of P' with respect to v: a binary string of length n−1 where s[j] = 1 if v→P'[j] and s[j] = 0 if P'[j]→v.

**Type-I position** at index j: s[j] = 0, s[j+1] = 1 (arc from successor back to predecessor).

**Type-II position** at index j: s[j] = 1, s[j+1] = 0 (arc forward then backward).

**Boundary term** b = s[0] + (1 − s[n−2]).

**insact(v, P')**: the actual count of valid insertion positions for v into P'. Equal to B_v(P') + S_v(P'), where B_v counts "big" positions and S_v counts "small" positions (precise definitions in the paper).

---

## Key Identities (Reference)

- inshat(v, P') is always odd (verified n≤6).
- (inshat(v,P') − 1)/2 = #{Type-II positions in P'} [algebraic identity, THM-004]
- #{Type-II positions in P'} = #{directed 3-cycles (v,a,b) : (a,b) consecutive in P'} [bijection, THM-005]
- insact(v, P') = B_v(P') + S_v(P') [proved, THM verified n≤6]
