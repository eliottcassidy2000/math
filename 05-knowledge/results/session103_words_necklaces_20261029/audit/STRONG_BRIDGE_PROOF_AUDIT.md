# Audit of the all-order forbidden-factor / strong-bridge equivalence

**Verdict:** the combinatorial proof is valid after one load-bearing definition
repair: explicitly require both blocks to have at least three vertices, or say
“nontrivial strongly connected induced blocks.” With the unqualified standard
meaning of strongly connected, the original statement is false.

This is an audit of a combinatorial equivalence for the exact recovered grammar
at pin `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`. It does not establish either
implication involving observer singularity.

## Required definition repair and minimal hostile

Singleton digraphs are strongly connected under the standard convention. In the
four-vertex seed, take R={2}, Q={0,1,3}. Q is the directed triangle
0→1→3→0. The only cross arc R→Q is 2→3. Thus the seed satisfies the unqualified
“both induced blocks strong, exactly one reverse cross arc” definition, while
its empty word contains neither forbidden factor.

The sentence “Both blocks necessarily have at least three vertices” must be a
hypothesis in the definition, not an asserted consequence. The original
computational checker already imposes this hypothesis (`len(V)<3` is rejected
and cuts have sizes 3 through n−3); its finite data need no recomputation for
this correction. Every occurrence of the intrinsic observer conjecture should
refer to this corrected definition.

## All-order proof audit after repair

1. **Newest maximum.** The indegree update proves the newest vertex is the sole
   indegree-one vertex or the second vertex of the pair, exactly as claimed.
   This uses generated strong tournaments, so there are no indegree-zero
   vertices. Strongness is preserved by insertion: old vertices reach the
   chosen parent then the new vertex, and the new vertex reaches the old
   tournament. The seed is strong.
2. **Cut intervals.** A Hamilton cycle must traverse its unique R→Q edge once.
   It therefore has exactly one run of each block. This observation applies
   equally when the high block is nonstrong, provided the Hamilton cycle
   exists; no strongness is needed for the interval conclusion itself.
3. **Persistence.** The size-three lower bounds ensure every R vertex has
   indegree at least two. An indegree-one maximum and its predecessor are in
   Q: the only possible entry vertex has both the bridge incoming and an
   internal incoming arc. All legal parents consequently lie in Q; adding the
   new vertex there preserves strongness and the same unique bridge.
4. **First appearance.** The newest vertex cannot lie in R and its sole
   predecessor must lie in Q. For nonstrong Q0=Q−{x}, the parent must lie in its
   final SCC. If |Q0|≥3, that SCC contains neither an indegree-one vertex nor
   the cycle predecessor of one. An indegree-one vertex can only be in C1 or
   in a singleton C2 after a singleton C1. In the latter case its predecessor
   is C1, while in the former its predecessor is within C1 or R. None is the
   final SCC when |Q0|≥3. This exhausts the legal-parent list without assuming
   any unspecified classification theorem.
5. **The two-vertex remainder.** Hence Q0={u,p} with u→p. The unique bridge must
   enter u, since the marked cycle visits both consecutively. Both are the
   indegree-one maximum pair; the final insertion is e. The previous letter
   is b,d,or e. A previous d would leave the former second maximum as a source
   in R after deleting u and p; this contradicts nontrivial strongness of R.
   Thus the previous letter is b or e. No missing seed exception exists:
   |R|≥3 and |Q0|=2 force the parent to have at least five vertices.
6. **Creation.** A final be or ee creates the triangle u→y→x→u. Deleting its
   older selected maximum u recovers a preceding generated strong tournament,
   because u is the newest vertex in all nonseed cases. The sole seed case is
   e,e, and deleting seed maximum 1 leaves 0→2→3→0. Therefore R is nontrivial
   strong in every case. Persistence then handles arbitrary legal suffixes.

The proof works with any one chosen cut; it neither assumes uniqueness of the
cut nor assumes different cuts are disjoint. The unique marked Hamilton cycle
is provided by the grammar, but the interval and SCC arguments only need that
marked cycle to exist. No quotient of histories is taken.

## Independent finite control

`check_strong_bridge_proof.py` imports no production grammar or bridge code.
It reconstructs legal operations using indegrees and enumerates **all vertex
subsets**, not just cycle intervals. All 232 legal histories through order nine
pass the corrected equivalence. It also checks every cut is cyclically
consecutive and every first-appearance cut has a three-vertex high block and a
word ending in be or ee. Nineteen histories are first-appearance histories in
this bounded search. The singleton seed hostile is checked using the standard
strong-connectivity convention.

```sh
python3 -B -O session103_words_necklaces/audit/check_strong_bridge_proof.py
```

Output: `strong_bridge_audit_receipt.json`. Explicit exception checks remain
active under optimization. No spectral implication, all-order nullity theorem,
or pointwise kernel persistence is promoted by this audit.
