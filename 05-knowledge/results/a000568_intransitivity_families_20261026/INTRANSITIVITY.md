# Intransitivity is the relation modulo its modular tree

> **CORRECTION (2026-10-26).** Section 3's original all-order `|Aut|` and
> strong-component-isotype claims are retracted. The sharp first counterexample
> has 12 vertices and orders 81/27. The repair retains the child-colour stabilizer;
> the cited fibre theorem, H-invariance and finite counts survive. See
> [the sharp 4p theorem and independent controls](../tournament_composition_classification_20261026/fibre/RESULT.md).

klein, 2026-10-26. The object is a tournament `T` on a finite set `V`, i.e. a complete
antisymmetric binary relation, and its **intransitivity structure**: the 3-uniform hypergraph
`C3(T)` whose edges are the cyclic triples `{a,b,c}` with `a -> b -> c -> a` (the literature's
"C3-structure", Boussaïri–Ille–Lopez–Thomassé 2004). Transitivity fails on exactly these triples,
and every failure of transitivity of a tournament is witnessed by one of them (a tournament is
transitive iff `C3(T)` is empty). The question of the directive — "binary intransitivity in the
relation itself" — is: how much of the relation does its intransitivity structure determine?
The answer is exact, and it is governed by the modular decomposition tree.

## 0. Vocabulary

A **module** of `T` is a set `M` such that every `z` outside `M` is uniform on `M` (`z` beats all
of `M` or loses to all of `M`). Gallai's theorem for tournaments: for `|V| >= 2`, either `T` is
not strong and its strong components `C_1 -> ... -> C_k` (`k >= 2`) partition `V` with
transitive quotient (**linear node** of arity `k`), or `T` is strong and its maximal proper
modules `M_1, ..., M_k` (`k >= 3`) partition `V` with a **prime** quotient `Q` (no non-trivial
module; **prime node**). Recursing gives the modular decomposition tree `𝒯(T)`; leaves are
vertices. `p(T)` = number of prime nodes, and `k_N` = arity of a linear node `N`. A tournament
is **prime** iff `𝒯(T)` is a single prime node with leaf children (`n >= 3`). Note `C_3` is
prime and no 4-tournament is prime (THM-1960's seed census, extended below).

`T*` is the converse (all arcs reversed). Dualizing a prime node means replacing its quotient `Q`
by `Q*` while keeping the children; permuting a linear node means reordering its children.

## 1. Theorem I (the fibre of the intransitivity structure)

**Statement.** For every tournament `T` on `V`, the tournaments `T'` on `V` with
`C3(T') = C3(T)` are exactly those obtained from `T` by, independently, (a) permuting the
children of each linear node of `𝒯(T)` arbitrarily, and (b) dualizing or not each prime node.
These moves commute and act freely, so the number of such `T'` is

    N(T) = 2^{p(T)} · ∏_{linear nodes N} k_N! .

Consequently `2 <= N(T) <= n!` for `n >= 3`, with `N(T) = 2` iff `T` is prime and `N(T) = n!`
iff `T` is transitive.

**Status.** The hard direction is the published theorem: for a prime `T`, `C3(T') = C3(T)`
forces `T' ∈ {T, T*}` (Boussaïri, Ille, Lopez, Thomassé, *The C3-structure of the
tournaments*, Discrete Math. 277 (2004) 29–43; restated as Theorem 24 of Boussaïri, Chergui,
Ille, Zaidi, *3-uniform hypergraphs: modular decomposition and realization by tournaments*,
arXiv:1805.04607), together with the fact that a realizable 3-uniform hypergraph and its
realizations share the same strong modules (same paper, Section 4), whose Corollary 57 is the
counting formula `|R(H)| = 2^{|D_△|} ∏_{X ∈ D_○} |Π(H[X])|!` — identical to `N(T)` above
(`△` = prime nodes, `○` = linear nodes). **CITED**. Here: the easy direction and Lemma C are
proved in full; the formula is **VERIFIED-EXACT** on every isomorphism class for `n <= 7`
(530 classes) by exhaustive depth-first enumeration of all tournaments on the same labelled
vertex set with the same cyclic triples (`scripts/c3class.c`, receipt `receipts/c3class-check.txt`,
0 mismatches), and indirectly for `n <= 10` through the exact counts of §4.

**Proof of the easy direction (moves preserve `C3`).** (a) Let `N` be a linear node with
children `C_1, ..., C_k` (strong components of the corresponding strong-module quotient). A
triple inside one child is untouched. A triple with two vertices in one child and one outside is
transitive before and after, because the outside vertex is uniform on the child. A triple meeting
three children is transitive before and after, because the quotient is transitive in both. (b) Let
`N` be a prime node with blocks `M_1, ..., M_k` and quotient `Q`; dualize `Q`. Triples inside a
block are untouched; a triple with two vertices in one block is transitive before and after (the
third vertex is uniform on the block); a transversal triple is cyclic in `Q[...]` iff cyclic in
`Q`, iff cyclic in `Q*`. (c) The moves at distinct nodes reverse pairwise disjoint sets of arcs:
an arc `{x,y}` belongs to the unique node whose children separate `x` from `y` (their least
common ancestor), so the group `∏_N S_{k_N} × (Z/2)^{p(T)}` acts freely on the labelled
tournaments, which gives the count once the hard direction is granted. ∎

**Proof of the extremes.** Writing `n - 1 = Σ_N (k_N - 1)` over internal nodes and using
`a! · b! <= (a + b - 1)!` for `a, b >= 2`, `∏ k_N! <= n!` with equality iff there is a single
internal node and it is linear (transitive `T`); and for `n >= 3`, `N(T) = 2` forces a single
internal node that is prime (a linear node contributes `k! >= 2`, two nodes contribute `>= 4`). ∎

## 2. Lemma C (the strong components are the hypergraph components)

**Statement.** The connected components of `C3(T)` (a vertex on no cyclic triple is a
singleton component) are exactly the strong components of `T`.

**Proof.** Every cyclic triple lies inside one strong component, so the hypergraph components
refine the strong components. Conversely let `C` be a strong component, `|C| >= 3`, and suppose
`C = A ⊔ B` with `A, B` non-empty and no cyclic triple meeting both. For `a ∈ A` put
`B⁺(a) = {b ∈ B : a -> b}`, `B⁻(a) = B ∖ B⁺(a)`.
(i) *The sets `B⁺(a)` form a chain.* If `b ∈ B⁺(a) ∖ B⁺(a')` and `b' ∈ B⁺(a') ∖ B⁺(a)`, the
transitive triple `{a, b, a'}` with `a -> b -> a'` forces `a -> a'`, and `{a', b', a}` with
`a' -> b' -> a` forces `a' -> a`: contradiction.
(ii) *`B⁻(a) ⇒ B⁺(a)` for every `a`* (every vertex of the first set beats every vertex of the
second): for `b ∈ B⁺(a)`, `b' ∈ B⁻(a)`, the transitive triple `{a, b, b'}` with `a -> b`,
`b' -> a` forces `b' -> b`.
(iii) Let `a₂` have the largest `B⁺(a₂)` and `X = B⁻(a₂) = ∩_a B⁻(a)`. Then `X ⇒ A` by
definition and `X ⇒ B⁺(a)` for every `a` by (ii), i.e. `X ⇒ C ∖ X`. Since `T[C]` is strong and
`C ∖ X ⊇ A ≠ ∅`, `X = ∅`: `a₂ ⇒ B`.
(iv) The same argument with `A` and `B` exchanged yields `b₂ ∈ B` with `b₂ ⇒ A`. Then
`a₂ -> b₂` and `b₂ -> a₂`: contradiction. ∎

(Lemma C is the root case of the "shared strong modules" statement of arXiv:1805.04607; the
proof above is elementary and self-contained.) Verified on every class for `n <= 8`
(`scripts/conn_check.py`, receipt `receipts/conn-normal.txt`).

## 3. Corollary (what the intransitivity structure determines, and what it does not)

The following are functions of `C3(T)`: `c3 = |E(C3(T))|`; the strong-component
vertex sets and their C3-isomorphism types; the common strong-module tree with
prime quotients up to duality; and the entire ordered path-cover profile,
including `H(T)`. The score sequence and automorphism order are **not**.

**H-survivor, with the needed induction strength.** At a linear node, reorder
within each cover path its runs through the children into the new child order.
At a prime node, reverse the sequence of child runs while dualizing the
quotient, keeping the internal direction of each run. Both are bijections on
ordered path covers. A move in an internal module lifts to the whole tournament
by replacing its runs in their existing outer positions; uniform adjacency
outside the module preserves entry and exit arcs. Thus all path-cover counts
are C3-invariant. Details and independent exact profile controls are in the
[repair, F4](../tournament_composition_classification_20261026/fibre/RESULT.md).

**Failed implication and repair.** At a canonical prime node,
`|Aut(Q[T_i])| = |{sigma in Aut(Q): T_sigma(i) ≅ T_i for all i}| product_i |Aut(T_i)|`.
While `Aut(Q*)=Aut(Q)`, a move *inside* one child can change its oriented
isomorphism type and hence this subgroup. For `A=C3→1`, compare `C3[A,A,A]`
and `C3[A,A,A*]`: orders 81 and 27, same 67 labelled cyclic triples, and
`H=462915`. This is the first failure; the p-part is invariant below 4p and
can first change at 4p for every odd prime p (repair, F2–F3). Strong-component
isotypes are not invariant either: these two examples are themselves strong.

The original finite confirmation through order 9 was correct, but its range
cannot see the first failure. Its recorded H, Aut and component-size tests and
counts are preserved. At the class level the valid orbit mass is
`sum_[T] 1/|Aut(T)| = |R(C3)|/|Aut(C3)|`; it does not require equal stabilizers.

For a global converse pair `(T,T*)`, automorphism order and all other
converse-invariant quantities agree. This explains those ties on THM-931's
wild pairs; it does not make every such quantity invariant under all local
C3-fibre moves, nor prove that every wild pair has equal C3 structure.

*Scores plus intransitivity.* Adding the score **function** `V -> N` to `C3(T)` kills every
reordering move (a reordering shifts the scores of a whole block by a positive amount) but not a
dualization at a prime node whose quotient is balanced with respect to the block sizes
(`Σ_{j: i -> j} |M_j| = Σ_{j: j -> i} |M_j|` for every `i`, e.g. a regular quotient on
equal-size blocks). Exhaustively, the classes on which `(C3, score function)` does **not** determine `T` are
**exactly** the classes with a size-balanced prime node — `1, 2, 5, 14, 53` classes for `n = 3..7`
(`scripts/c3score.c`, `scripts/balanced_check.py`, receipts `c3score-n3-n7.txt`,
`balanced-normal.txt`); `C_3` itself is the first. So the only cancellation-free statement that
is finite-exact here is: scores kill every reordering and every unbalanced dualization, and
moves at different nodes never cancel in the score function (`n <= 7`).

## 4. Counting intransitivity structures

Let `c(n)` be the number of isomorphism classes of `C3`-structures of `n`-tournaments (equivalently
of realizable 3-uniform hypergraphs on `n` vertices), `c_s(n)` those of strong tournaments, and
`ℓ(n)` the number of distinct `C3`-structures on the labelled set `[n]`.

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| A000568 | 1 | 1 | 2 | 4 | 12 | 56 | 456 | 6880 | 191536 | 9733056 |
| c(n) | 1 | 1 | 2 | 3 | 9 | 32 | 247 | 3328 | 93612 | 4774417 |
| c_s(n) | 1 | 0 | 1 | 1 | 6 | 22 | 214 | 3074 | 90255 | 4680547 |
| ℓ(n) | 1 | 1 | 2 | 11 | 208 | 8563 | 680426 | 101747483 | 28904318404 | 15824793243139 |
| ℓ(n)/2^{C(n,2)} | 1 | .5 | .25 | .172 | .203 | .261 | .325 | .379 | .421 | .450 |

`ℓ(n) = Σ_{[T]} (n!/|Aut T|)/N(T)` (sum over isomorphism classes; Theorem I). The unlabelled
counts obey a species equation parallel to session56's `T = SEQ(S)`, `S = X + P∘T₊`:

    c = MSET_{>=1}(c_s)             (Lemma C: a C3-class of a non-strong tournament is the multiset
                                     of the C3-classes of its strong components, in no order)
    c_s = X + Σ_{[Q]_± prime up to duality} Z_{Aut±(Q)}[c]

where `Aut±(Q)` is the group of automorphisms and anti-automorphisms of the prime quotient
(`|Aut±| = 2|Aut|` iff `Q` is self-converse) and `Z_G[c]` is Pólya substitution of the series
`c(x), c(x²), ...` into the cycle index of `G`. Both equations were evaluated from the prime
cycle-index census (`data/prime_cycle_types_n9.tsv`, plus the 7,906,161 primes of order 10, of
which 5,525 are self-converse) and agree **exactly** with the direct count of nauty canonical
forms of the hypergraphs for every `n <= 10` (`receipts/species-normal.txt`,
`receipts/species10-normal.txt`, `receipts/census-n10.txt`). Neither `c(n)` nor `c_s(n)` nor
`ℓ(n)` is in the OEIS (searched 2026-10-26).

The distribution of the number of tournament classes sharing one `C3`-structure (`n = 10`):
`1: 8587, 2: 4673145, 3: 176, 4: 89133, 5: 10, 6: 3039, 7: 1, 8: 269, 9: 1, 10: 15, 12: 36,
15: 1, 16: 2, 18: 1, 20: 1` — almost every structure is realized by exactly a chirality pair.

## 5. The information reading (Hartley fibre of the intransitivity observable)

Fix the ambient universe `U_n` = tournaments on `[n]`, the equivalence = equality, the observable
`L = C3`. In the only sense the repo allows after MISTAKE-231 (finite-fibre Hartley information
after fixing universe, relation and observable), the information the intransitivity structure
fails to carry about the relation is

    H₀(T | C3) = log₂ N(T) = p(T) + Σ_N log₂(k_N!)  bits,  1 <= H₀ <= log₂ n!,

minimal exactly on the prime tournaments (one labelled converse choice, distinct
from THM-1420's `A_n` orbit split) and maximal exactly on the transitive relation (the only one whose
intransitivity structure says nothing beyond "a linear order"). Under the uniform law on `U_n`:

| n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| E[log₂ N] | 2.189 | 2.969 | 3.040 | 2.602 | 2.057 | 1.640 | 1.378 | 1.223 |
| P(T prime) | .250 | 0 | .258 | .315 | .464 | .601 | .721 | .815 |

**Proposition.** `E[log₂ N(T)] -> 1`. *Proof.* `N(T) = 2` unless `T` has a non-trivial module,
and `log₂ N(T) <= log₂ n! <= n log₂ n`; the expected number of non-trivial modules is at most
`ε_n = Σ_{m=2}^{n-1} C(n,m) 2^{-(m-1)(n-m)}`, whose two extreme terms are `C(n,2) 2^{2-n}` and
`n 2^{2-n}` and whose interior terms are smaller, so `ε_n n log₂ n -> 0`, and
`E[log₂ N] - 1 <= ε_n · n log₂ n -> 0`. ∎ (At `n = 10`: `ε_10 = 0.225`, observed `P(not prime) = 0.185`.)

So a typical binary relation of this kind is determined by its intransitivity structure up to one
bit; at the class level, `c(n)/A000568(n) = 0.4906` at `n = 10` and tends to `1/2` for the same
reason (almost all classes are prime and not self-converse, and the crude bound
`#non-prime classes <= #non-prime labelled / (n!/3^{(n-1)/2})` suffices since `(√3/2)^n n² -> 0`).
Level-resolved (by `c3`): the mean of `log₂ N` falls monotonically from `log₂ n!` at level 0 to
about 1.2 at the top levels, and the fraction of prime classes rises from 0 (levels `< n - 2`) to
above 0.8 (`receipts/analyze-normal.txt`; FAMILIES.md §E). Minimum number of cyclic triples of
a prime tournament: `1, 3, 4, 5, 6, 7, 8` for `n = 3, 5, 6, 7, 8, 9, 10`, i.e. Moon's bound
`n - 2` for strong tournaments is attained by primes from `n = 5` on; the extremal strong
tournaments are the `2^{n-4}` chain-reversed transitive tournaments and the extremal primes the
`F_{n-4}` of them with all chain gaps `<= 2` (FAMILIES.md §E, HYP-9413).
