# Labelled free-arc flips force the entire tiling state

Status: PROVED all n>=3; FINITE-EXACT independent partition/refinement certificate n=3,4,5. Separate addendum: the audited main NOTE.md, check.py, and results.json are unchanged.

## Scope and inheritance

Use the source X={0,1}^m, m=binom(n-1,2), fixed directed Hamiltonian path, and converse-merged node observation o from NOTE.md. All free-edge labels remain pinned to path positions. The allowed operations now include each individual labelled free-edge flip g_k(x)=x XOR 2^k, not merely their product tau. The output remains only the current merged tournament node; color is unnecessary.

Inherited scaffolding: THM-840-hamming-five-continuation-congruence-boundary.md's operation-kernel criterion and THM-796-three-sorted-recursive-tiling-line-node-incidence.md's observation/operation-relative continuation semantics. Targeted searches for singleton transitive fibres, coordinate-flip continuation, translating observations, and full Boolean continuation found no exact prior theorem in the relevant theorem/results snapshot. The general mechanism below is elementary group-action separation, not a claimed novel theorem of automata theory.

## 1. Singleton-anchor separation lemma

Let a group G act transitively on a finite set X, and suppose an observation o:X->Y has a singleton fibre {a}. Then equality of all continuation outputs o(gx), g in G, forces equality of source states.

Proof: for x!=y, choose g with gx=a. Since g is a bijection, gy!=a, so o(gx)!=o(gy). Thus no distinct pair is continuation-equivalent. This applies to an operation monoid whenever its allowed words include the required group action.

For fixed-path tilings, the all-forward mask a=2^m-1 is the unique tiling with transitive tournament node. A transitive tournament has a unique linear order compatible with its arcs. The prescribed directed Hamiltonian path 0->1->...->n-1 forces that order, hence every free edge points forward. Converse merging changes nothing: the converse of a transitive tournament is isomorphic to it, and no nontransitive tournament is isomorphic to either.

An alternative elementary witness: every non-all-forward tiling has a shortest backward free edge j->i. With k=i+1, minimality makes i->k and k->j forward, giving a directed triangle i->k->j->i. Thus it cannot belong to the transitive node. The checker records one such triangle for every nonanchor at n<=5.

The labelled flips generate the regular translation action of F_2^m. For x!=y the word toggling the coordinates of x XOR a sends x to a and y to a XOR x XOR y !=a. This word has length at most m. Hence the coarsest deterministic continuation-complete quotient is equality on X, with exactly 2^m states.

## 2. The linear closure also becomes full

In H=R^X, the current node-indicator space contains delta_a because the transitive fibre is singleton. Its translates under all flip words are the delta_x for every x in X. Consequently the observable linear closure is all H, dimension 2^m. The unobservable subspace is zero, so any exact linear continuation state for arbitrary signed source vectors must have dimension at least 2^m. The literal tiling space attains it.

This is stronger than the tau-only result precisely because the operation alphabet is larger. At n=5 the tau-only cyclic observer has dimension19, but all labelled flips have dimension64. No appeal to colors, line multiplicities, hidden phase, or recursive faces is needed.

## 3. Positive stochastic realization has the same full lower bound

Suppose every pure tiling x is represented by a probability distribution p_x over h hidden states. Let O be a fixed column-stochastic node-output matrix and T_k fixed column-stochastic transitions for the labelled flips. Require exact deterministic outputs for every word w:

    O T_w p_x = e_(o(g_w x)).

If x!=y, the singleton-anchor lemma gives a common word w for which the two outputs differ. Extremality of the output simplex forces every hidden state in supp(p_x) to have output e_(o(g_w x)) under O T_w, and every state in supp(p_y) to have the different output e_(o(g_w y)). Thus supp(p_x) and supp(p_y) are disjoint. All 2^m nonempty source supports are pairwise disjoint, so h>=2^m. The literal tiling permutation model attains equality.

The result allows nonlinear x->p_x and does not require full-state intertwining. But it requires all pure tiling initial conditions and all operation words; it does not apply automatically to a restricted family of averaged initial conditions. Fixed transition labels are essential: unlabelled "flip some arc" is a different operation semantics.

## 4. Tiny certificate and comparison

Run after the main checker:

    python3 metagraph-memory/check_full_flips.py
    python3 -O metagraph-memory/check_full_flips.py

The addendum reads the independently computed actual node map from results.json and records its SHA256. It checks the singleton transitive fibre, explicit triangle witnesses for all other masks, iterative continuation partition refinement for tau versus all coordinate generators, every distinct pair's explicit anchor-separating word, and full rank of the translated anchor delta matrix. Checks are explicit require/raise and survive optimized Python. The receipt is full_flips_receipt.json.

| n | tilings | tau deterministic states | tau linear dimension | full-flip deterministic/positive/linear minimum |
|---|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 2 |
| 4 | 8 | 5 | 5 | 8 |
| 5 | 64 | 38 | 19 | 64 |

## 5. What this obstruction does not say

This is an exact state-size obstruction for a continuation-complete carrier under every labelled free-coordinate flip. It is NOT a runtime or memory lower bound for static enumeration of tournament isomorphism classes (A000568), computing a node census, using an implicit representation, or a particular metagraph algorithm. Static counting need not retain an independently addressable state for every possible future labelled flip. The theorem also does not say that the full tiling state remains necessary under a restricted operation alphabet; the tau-only main result is a strict counterexample.

Source -> target: raw tilings -> any operation-complete quotient. Preserved predicate: all current-node outputs after labelled flip words. Lost information permitted: none. Decisive obstruction: a singleton observable anchor plus transitive reversible operations. This identifies the exact boundary between a successful operation-specific compression and an impossible all-flips quotient.
