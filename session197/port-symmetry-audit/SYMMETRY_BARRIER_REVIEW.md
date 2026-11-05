# Independent review: persistent two-port symmetry barrier

Verdict: **PASS**, 2026-11-05. No correction required.

Frozen files reviewed in `/workspace/scratch/07d5532f2fce/session197-hybrid-constructions/symmetry-barrier/`:

- PROOF.md SHA256 `4091d1e1ca9bcfd692a6fd0070fda4faa148b9edaa0d61cde17b56e912da9aff`
- verify.py SHA256 `7c5b04ef8f79675dcb94d3a91450eb7b7e72de756fe7747a3069e3cb6223af3d`

The exact persistent-port language is correctly restricted: no initially unexposed or forgotten old vertex may be re-exposed. For either birth orientation the new adjacency function is constant on all nonports, and the current ports are individually fixed by the surviving subgroup. Extending that subgroup to fix the new vertex therefore preserves every old and new edge. Any forgetting choice retains only fixed vertices. The induction does not require the action-selection strategy to be equivariant, and allows arbitrary graph- and target-dependent choices.

The action on complete Hamilton vertex words is free, because a fixed word fixes every vertex individually. Thus subgroup order divides every descendant count. Empty words and zero births are covered by the initial subgroup hypothesis.

The displayed nine-vertex graph is exactly the strong cyclic substitution C3[C3,C3,C3]. The three independent internal rotations are actual disjoint commuting order-three permutations. Every pair of ports leaves at least one module untouched; a same-module pair leaves two untouched rotations and a subgroup of order9. The conclusions H divisible by3 or9 are therefore valid at arbitrary depth, and exclude both literal nonzero ternary children from that input and every descendant. No global absence of9479 or9481, or universal failure from every small seed, is claimed.

Independent replay of the supplied verifier normally and with `python3 -O` produced byte-identical JSON receipts. The independent permutation and subset-DP path counts both give3159. All72 initial ordered pairs,1728 first births,41472 second births,72 selected Hamilton divisibility controls, and the full7! automorphism test of the rigid H49 contrast pass. Frozen hashes remained unchanged.

Output files beside this report: `symmetrybarrier-review-normal.json` and `symmetrybarrier-review-optimized.json`. No repository or worktree was changed.

Separate possible extension, not needed for the frozen claim: for k continuously persistent ports and m disjoint directed-triangle modules in any strong quotient, at least m-k modules remain untouched, giving a divisor3^max(m-k,0) forever. This is a quantitative initial-port lower bound, with the same no-re-exposure restriction.
