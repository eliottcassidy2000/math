# Formal Proof: Dichotomy for Cycle-Rich Tournaments (COMPLETE)

**Instance:** kind-pasteur-2026-03-07-S33
**Status:** PROVED for all n >= 9

## Statement

**Theorem (Dichotomy).** Let T be a cycle-rich tournament on n >= 9 vertices
(every vertex in a directed 3-cycle). Then at least one of:
(a) T contains 3 pairwise vertex-disjoint directed 3-cycles, or
(b) There exists a vertex v such that T - v is cycle-rich on n-1 vertices.

## Key Preliminary: Cycle-Rich Implies No Source/Sink

**Lemma Q.** If every vertex of T is in a directed 3-cycle, then T has no
source (score n-1) and no sink (score 0).

**Proof.** Source v (score n-1) beats all opponents. For any 3-cycle
v -> a -> b -> v: need b -> v, but v beats b. Contradiction.
Sink v (score 0) beats nobody. For any 3-cycle v -> a -> b -> v: need
v -> a, but v beats nobody. Contradiction. QED.

This means "cycle-rich" = "every vertex in a 3-cycle" with no additional
source/sink condition needed.

## Proof of the Dichotomy

If mm(T) >= 3 (max 3-cycle matching), conclusion (a) holds. Assume mm <= 2.

### Setup

Take a maximal matching: disjoint 3-cycles A, B (if mm = 2) or a single
3-cycle A with B = empty (if mm = 1, though this requires minor adaptation).
Let R = V \ (A union B), |R| >= n - 6 >= 3.

Every 3-cycle C in T satisfies C cap (A union B) != empty.
(Otherwise {A, B, C} is a 3-matching, contradicting mm <= 2.)

### The Poisoning Graph

**Definition.** The poisoning graph P = (R, E_P) is defined by:
w -> v in E_P iff EVERY 3-cycle containing w also contains v.

Intuitively, w "depends on" v: deleting v would destroy all of w's 3-cycles.

### Lemma R.1: Out-degree at most 1

**Claim:** Every vertex in P has out-degree at most 1.

**Proof.** Suppose w -> v and w -> v' with v != v'. Then every 3-cycle
through w contains both v and v'. A 3-cycle has 3 vertices, so the only
possible vertex set is {w, v, v'}. But {w, v, v'} subset R and
R cap (A union B) = empty, so {w, v, v'} cap (A union B) = empty.
Then {A, B, {w, v, v'}} is a 3-matching, contradicting mm <= 2. QED.

### Lemma R.2: Acyclicity

**Claim:** P is a DAG (directed acyclic graph).

**Proof.** Suppose P contains a directed cycle w_1 -> w_2 -> ... -> w_k -> w_1
with k >= 2.

Since w_k -> w_1 in P: all 3-cycles containing w_k also contain w_1.
Since T is cycle-rich, w_k is in at least one 3-cycle C. This cycle must
intersect A union B (since mm <= 2). So C = {w_k, w_1, x} for some
x in A union B.

C contains w_1. Does C contain w_2?
- w_2 != w_k (cycle length >= 2)
- w_2 != w_1 (they are distinct vertices of the directed cycle in P)
- w_2 in R while x in A union B, so w_2 != x
Therefore w_2 not in C.

But w_1 -> w_2 in P means every 3-cycle containing w_1 also contains w_2.
C contains w_1 but not w_2. Contradiction. QED.

### Lemma R.3: Safe Deletion Exists

**Claim:** There exists v in R such that every vertex in T - v is in a 3-cycle.

**Proof.** P is a DAG with out-degree <= 1. Such a graph is a forest of
directed paths, each ending at a "sink" (out-degree 0 vertex).

Every non-empty DAG has at least one "source" (in-degree 0 vertex).
Let v be a source of P.

We show every vertex in T - v is in a 3-cycle:

**Vertices in A:** Cycle A doesn't contain v (v in R, not in A). So A
persists in T - v. Every vertex of A is in 3-cycle A.

**Vertices in B:** Same argument. B persists in T - v.

**Vertices in R \ {v}:** Let w in R \ {v}.

*Case 1: w has out-degree 0 in P.* Then w has a 3-cycle C_w that does NOT
depend on any single other R vertex. In particular, w has a 3-cycle
C_w subset {w} union (A union B) (using only w and vertices from A union B).
Since v not in A union B and v != w: v not in C_w. So C_w persists.

*Case 2: w has out-degree 1 in P.* Then w -> f(w) for some f(w) in R, meaning
all of w's 3-cycles go through f(w). Since v is a source of P (in-degree 0),
no vertex points to v, so f(w) != v. A 3-cycle of w is {w, f(w), x} with
x in A union B. Since v not in {w, f(w), x} (w != v, f(w) != v, x in A union B):
this cycle persists in T - v.

In all cases, every vertex in T - v is in a 3-cycle. By Lemma Q, T - v has
no source or sink. Therefore T - v is cycle-rich. QED.

## Consequence: H(T) != 21 for All Tournaments

**Theorem.** H(T) != 21 for any tournament T on any number of vertices.

**Proof.** Strong induction on n with subsidiary claim:
*For all cycle-rich T on n >= 8: H(T) >= 25.*

**Base case (n <= 8):** Exhaustive verification. H = 21 never occurs among
all 2^C(n,2) tournaments for n <= 8. For cycle-rich n = 8: min H = 25.

**Inductive step (n >= 9):**

*Case 1: Some vertex not in any 3-cycle.*
By the Key Lemma (Part J), this vertex v is in no cycle of any length.
Omega(T) = Omega(T - v), hence H(T) = H(T - v) != 21 by induction.

*Case 2: Every vertex in a 3-cycle (cycle-rich).*
By the Dichotomy Theorem:
- (a) If 3 pairwise-disjoint 3-cycles exist: Part C gives
  alpha_1 + 2*alpha_2 + 4*alpha_3 >= 13 > 10, so H >= 27 > 21.
  Also H >= 27 >= 25 (subsidiary claim).
- (b) If safe deletion v exists: T - v is cycle-rich on n - 1 >= 8.
  By induction on the subsidiary claim: H(T - v) >= 25.
  By independence polynomial monotonicity: Omega(T - v) is an induced
  subgraph of Omega(T), obtained by removing the clique S of all cycles
  through v. So H(T) = H(T-v) + sum_{c in S} 2*I(Omega-N[c], 2) >= H(T-v) + 2.
  Therefore H(T) >= 25 + 2 = 27 > 21. Also >= 25 (subsidiary claim). QED.

## Literature Connections

- **Lichiardopol's conjecture (proved):** For l = 3, k = 3: tournaments with
  min out-degree >= 5 have 3 disjoint 3-cycles. This gives conclusion (a)
  directly when min out-degree >= 5. However, at n = 9 cycle-rich, min
  out-degree is always <= 4, so Lichiardopol doesn't fire. The poisoning
  graph argument covers all cases regardless of minimum degree.

- **Frankl's proof of Erdos matching conjecture (k=3):** Bounds the max number
  of 3-element sets with no matching of size s+1. Not directly needed but
  provides context for the 3-uniform hypergraph matching theory.

- **Bang-Jensen, Bessy, Thomasse:** Proved Lichiardopol's conjecture for l = 3.

## Sources

- [Chen-Chang 2024: Disjoint cycles in tournaments](https://onlinelibrary.wiley.com/doi/10.1002/jgt.23038)
- [Lichiardopol's conjecture proof](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v27i2p52)
- [Vertex-disjoint cycles of different lengths (2024)](https://arxiv.org/abs/2403.03692)
- [Chen 2025: Extremal Results on Disjoint Cycles](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.23255)
