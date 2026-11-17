# The extension graph at order 10, and a re-arranged toothpick tree: the heaviest parent

Owner: *"keep pushing the extension graph and toothpick tree themes but also re-arrangements and extensions of those
ideas."* The previous packet (`toothpick_grammar_creators_20261116`, THM-5266 (D)) stopped at order 9; here order 10
(`extgraph10.c`, nauty `labelg` in chunks, 86 s) and a second tree.

## 1. Order 10 ([receipts/extgraph10.out](receipts/extgraph10.out))

| n | classes | contacts | mean distinct parents (of n) | children per class of order n−1 (min, mean, max of 2^{n−1}) | `M_{n−1}` children | least-string tree: max branch, childless |
|---|---|---|---|---|---|---|
| 8 | 6,880 | 50,357 | 7.32 | 11, 110.4, 121 of 128 | 106 | 106, 215 (3.1 %) |
| 9 | 191,536 | 1,656,619 | 8.65 | 22, 240.8, 248 of 256 | 227 | 227, 2,977 (1.6 %) |
| 10 | 9,733,056 | 95,441,218 | 9.81 | 23, 498.3, 503 of 512 | 475 | 475, 78,498 (41.0 %) |

- 8,150,261 of the order-10 classes (83.7 %) have ten pairwise non-isomorphic deletions; only 3 have a single deletion
  class (the transitive tournament among them). The most fertile order-9 class realises 503 of its 512 extension
  patterns as distinct classes; the least fertile 23.
- The transitive tournament keeps its two roles: not the most fertile (475 < 503) but the largest branch of the
  least-string tree, because it is the least canonical string. The childless share of that tree at order 9 is 41.0 %
  (43 % at order 8); S-6's 45–65 % was refuted.

## 2. The heaviest-parent tree ([receipts/heavy_tree.out](receipts/heavy_tree.out), [receipts/heavy_ties_ic.out](receipts/heavy_ties_ic.out), [receipts/heavy_trunk.out](receipts/heavy_trunk.out))

Attach every class to its deletion class of largest count (ties by least string). **Deletion never increases the count**:
every Hamiltonian path of `T − v` extends by inserting `v` at the front, the end, or before the first vertex it beats, so
`H(T − v) ≤ H(T)` (0 violations in 198,944 classes; a classical insertion argument). Hence every class is the heaviest
parent of its own sink extension: the tree has **no childless node** at any order (against 41–43 % in the least-string
tree). Branching: max 3, 8, 13, 44, 115, 247 at orders 4–9; the largest branch at 8 → 9 is no longer the transitive
tournament but the class of maximum count: the order-8 class with `H = 661` and score sequence `3, 3, 3, 3, 4, 4, 4, 4` (247 children), and the next four branches also carry `H = 661` or `659` — the heaviest lineage is the maximal-count lineage, against the transitive trunk (`H = 1`, 227 children) of the least-string tree.

**Every heaviest deletion of a first creator is a first creator** — for all first creators of orders 4–9 (1, 6, 33, 328,
5,457, 163,676 of them), checking every tied heaviest deletion (5,064 ties at order 9), with the single exception of
the converse pair creating 63 (their heaviest deletions have count 37, born at order 6). Compare the least-string tree,
where a first creator's parent is a first creator for 0 %, 17 %, 27 %, 44 %, 61 %, 73 % of them at orders 4–9. So the
first creators form a tree under the heaviest-deletion map, rooted at the 3-cycle and the single vertex, with 63 the
only break: the owner's "relationship between these constructions and themselves". Conjectured for all orders —
HYP-9863 (a); THM-4922's genealogy ("some deletion is a first creator") is its weak form.

## 3. Toothpick reading

The least-string tree draws the growth from the transitive trunk outward (previous packet's figures); the heaviest-parent
tree draws it from the heaviest lineage — at each order the toothpick with the most children is a class of large count,
and every toothpick sprouts at least one child (its sink padding). The two trees are the two extreme re-arrangements of
the same contacts: lightest parent (the least string is usually the transitive-most deletion) versus heaviest parent.
