# Independent audit: canonical carrier and strong Boolean gate

Status: VERIFIED. All principal all-order arguments in the submitted NOTE.md are valid. Verification is an independent mathematical audit plus finite controls; it is not a literature-priority claim or a spectrum solution.

## Proof audit

- Canonical carrier: the endpoint identities are correct with Q_uv defined by first path ending at u and second beginning at v. Substituting them into the one-apex response gives exactly K_uv/2 on the original arc u→v, with zero on its reverse. Parity makes every such weight integral.
- Deletion floor: the slot argument counts legal endpoint placements by |R_i minus sigma(R_(i+1))|. Its lower bound telescopes to two. Deleting the designated endpoints recovers the original path and split, so the injection is valid, including the empty deletion minor. Thus c_uv >= H(A−{u,v}) >= 1.
- Equality: positivity makes H preservation equivalent to an empty outgoing cut in A. Those sets are exactly terminal unions of the linearly ordered strongly connected components, yielding k+1.
- Minimum cuts: orienting capacity on the original arcs agrees with the z-neighbor convention. Capacity M > total original capacity enforces each required source/sink vertex. Strong extension is equivalent to beating some vertex of the first SCC and losing to some vertex of the last SCC. Distinct endpoints are necessary when the SCCs coincide.
- Polynomial support: each permutation contributes a product of indicators of its consecutive edges. Each expanded monomial is therefore supported on a linear forest. Incidences in the exceptional set prove the 2s−r refinement. Integer Möbius inversion justifies coefficientwise division by two for always-even observer differences.
- Group minimum support: the difference-of-slices proof works for every abelian group, without division or field structure. A tiny induction-boundary clarification is suggested below.
- Literal gate and wrapper: the four rows, H and F(2), marked observer formula, and all-strong assertion are correct. The wrapper splits exactly into paths with a immediately before b and paths with b first/a last. Ordered slot counts multiply under ordinal sum.
- Sharpness and valuations: factoring out 3^(s−k) works for k<s; the parity padding handles k=s. Multiplicities are binom(s,k)3^(s−k), giving (1+3w)^s and precisely one nonzero residue mod 3. The limitation concerning unequal inactive marked observers is accurate and important.

## Source comparison

I retrieved and read both cited S152 files through the GitHub connector at commit 7c371e94aae8bd4994c2569ea92b60b0c56cd822:

- https://github.com/eliottcassidy2000/math/blob/7c371e94aae8bd4994c2569ea92b60b0c56cd822/05-knowledge/results/session152_sequences_20261103/cut_induction/README.md
- https://github.com/eliottcassidy2000/math/blob/7c371e94aae8bd4994c2569ea92b60b0c56cd822/05-knowledge/results/session152_sequences_20261103/recursions/RESULT.md

The credit is accurate. S152 establishes the quadratic response, endpoint covers, submodularity and integral max-flow carrier. Its carrier statement explicitly does not constrain orientation to follow the tournament. The reciprocity document uses the opposite apex-neighbor convention; the submitted NOTE correctly translates it. No comprehensive novelty search was performed.

## Reproduction

`python3 independent.py` and `python3 -O independent.py` use no imports from the submitted `core.py` or `check.py`.

- Exhaust all 1,099 labeled tournaments of orders one through five. Compute Q by direct permutation plus split enumeration. Check the canonical carrier, every deletion floor, all child counts, SCC equality, and minimum strong extension.
- Cross-check all child counts for old order at most four by direct permutations as well as an independently implemented memoized suffix recursion.
- Test every disjoint forced-source/forced-sink assignment on every old tournament through order four and four selected order-six tournaments, comparing an independent augmenting-path max-flow implementation against complete interface enumeration.
- Check forbidden Möbius supports on the full order-four orientation cube and on a seven-variable order-five cube whose variable edges all touch a designated two-vertex set, including its internal edge refinement.
- Test scalar modular controls for moduli 3,4,6,9,15, and genuinely vector-valued controls for moduli 3,4,6.
- Count gate H, F(2), and marked observer by direct permutations. Count every s=1 and s=2 literal wrapped family using an independent path recursion, and verify strong connectivity.
- Check valuation formula arithmetically for all 1<=s<=40 and 0<=k<=s.

The original `check.py` also passed both normal and optimized runs.

## Small fixes identified in the submitted version

These do not invalidate the main theorems and were reported to the publishing worker.

1. `F2(empty)` in `check.py` returns two, while the ordered two-empty-slot cover count is one. This branch is unused by the original tests; return one explicitly for the empty tournament.
2. The original reproduction paragraph claims a two-exception support check that was absent. The independent suite supplies it, or the paragraph can be narrowed.
3. The original test comment claims vector-group controls, but the original loop is scalar only. The independent suite supplies genuine vector-group checks.
4. State that the minimum strong-extension problem is infeasible for a one-vertex parent, where there are no distinct candidate endpoints.
5. In the minimum-support induction, handle d=N immediately before reducing an independent slice to N−1 variables, or explicitly use the slice's actual degree. The result itself is correct.
