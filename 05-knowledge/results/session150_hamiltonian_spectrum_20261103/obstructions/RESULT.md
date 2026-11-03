# No fixed uniform substitution gives universal additive induction

**Proved; no claim of priority.** Let Q be any fixed finite tournament with a specified vertex v. For any nonempty tournament T, write Q[v←T] for the substitution of T at v, with every vertex of T inheriting v's external orientations. No constants a,b with b≠0 satisfy

    H(Q[v←T]) = a H(T) + b

for every finite tournament T. In particular, no such wrapper implements H↦H+2 on every input. The theorem holds with arbitrary real constants a,b, not merely integer constants.

Proof. Strong components of a tournament are totally ordered by domination. Each directed Hamiltonian path traverses its strong components in that order, so its count is the product of their counts.

If {v} is a strong component of Q, replacing it by T inserts the strong components of T into this total order. Thus H(Q[v←T])=H(Q)H(T). An assumed affine identity applied to a one-vertex tournament and to the directed triangle, whose counts are 1 and 3, gives H(Q)=a+b and 3H(Q)=3a+b, hence b=0.

Otherwise v lies in a strong component S with at least three vertices. By the standard directed Hamiltonian-cycle theorem for strong tournaments (Camion), S has a Hamiltonian cycle v,a₁,…,aᵣ,v. Let T be the transitive n-vertex tournament. For each nonempty proper subset A of V(T), list A in its unique directed order, then a₁,…,aᵣ, then V(T)\A in its unique directed order. The two module-to-outside edges are the cycle edges v→a₁ and aᵣ→v; hence these are distinct directed Hamiltonian paths in S[v←T]. Prepend and append fixed Hamiltonian paths in the remaining strong components of Q, respecting their domination order. Therefore H(Q[v←T])≥2ⁿ−2, but H(T)=1 for every n. The proposed affine identity would give the constant a+b, a contradiction. ∎

The strong-cycle input may alternatively be replaced by a direct proof or standard citation; it is not derived from finite computation. The argument applies also when Q has a single vertex. It requires a single fixed quotient with one uniformly connected variable module and a formula valid for *all* tournament inputs, including transitive ones. It does **not** rule out adaptive quotients, rooted or nonuniform interfaces, use of multiple interacting copies of T, restrictions to a specially chosen family, nonlinear profile-dependent updates, or a finite-state induction confined to suitable invariant seed states. In particular, it does not settle the odd spectrum conjecture or disprove an H+2 induction that avoids counts 7 and 21.

## Complementary exact transfer obstruction

Let cₖ(T) count ordered k-tuples of nonempty vertex-disjoint directed paths covering T. If F(T) is the directed-triangle wrapper T→x→y→T, direct word counting gives

    H(F(T)) = 3c₁+c₂,
    c₂(F(T)) = 6c₁+8c₂+2c₃,
    H(F²(T)) = 15c₁+11c₂+2c₃.

The coefficients count possible words of module-runs, x and y and choices of cuts. These identities explain why scalar H alone is not generally a closed transfer state; the first-iterate counterexample and the general formula are in the accompanying transfer packet. The higher coordinate appearing symbolically does not by itself prove that any specific lower-dimensional invariant state is impossible, because profile coordinates can satisfy identities.

## Reproduction and provenance

Run `python3 verify.py`; its independent subset-and-last-vertex dynamic program checks all Q of order ≤4, every chosen v, and every T of order ≤3 in the singleton-component case (1,485 assertions). For each of the 150 non-singleton components it also explicitly finds a Hamiltonian cycle, and checks the bound for transitive T of orders 2–5 (600 assertions). Orient pair (i,j), i<j, from i to j iff the corresponding lexicographic bit is 0. `verify.out` records output. These bounded checks are controls, not a substitute for the all-order proof.

Prior repository work already develops path-cover substitution kernels, including `math/01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`, `math/01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`, and `math/05-knowledge/results/session26_random_bridges_20261022/session26-attachments/proof.md`. The bounded search made here does not establish novelty or literature priority for this obstruction. The separately named THM-2195 concerns arc-reversal distance, not Hamiltonian-path counts, and is not a source for this theorem.
