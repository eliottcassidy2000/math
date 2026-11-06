# Independent audit: constructive maximal-prefix signed repair

**Verdict: PASS in the stated conditional scope.** The frozen proof and README correctly establish exact target-dependent C1 section repair with maximal retention of a supplied chronological prefix. The two explicit families have the stated optimal deletion counts and, at those optimal retained prefixes, optimal replacement lengths. This is **not** universal C1/C2/C3 occupancy. No mathematical defect or implementation correctness defect was found.

## Frozen inputs and preservation

Audited proof SHA-256:

`3ac0eefd2815134456101cd6ba7841b756bc13667ab6e149417b529b7c2281d9`

Audited implementation SHA-256:

`2868922f90f5ba2cc8d6122408d87df692bfd64d1afa5747b814ca4f070af8dc`

Both author receipts have SHA-256:

`bb10aecd4f968a2868737cd88525d0c5b374ff1119e7e0e21e3fe4197dec7da1`

The frozen dependencies are accepted within their stated independently audited scope, with their inherited physical source grammar and C1 interpretation:

- APPEND-CLASSIFICATION.md: `e8f6ecb95718870e334854cdd0284e4f5e0b089bbab3753610047f3b084ece3a`.
- SIGNED-SECTION.md: `92b300d272d529bb21a504a914226b1245a631249d30d22d0e9362625a698182`.

All packet files were read without alteration. The complete original input packet, including the original README, is archived byte-for-byte in `frozen-inputs/` so later status promotion does not change the audited object. The audit script checks the frozen proof, implementation, and author receipt hashes before and after its computation. The audit is stored separately and makes no Git or GitHub changes.

## Mathematical review

### 1. Finite target table and its completeness

For a genuine source, nonnegative physical coordinates give `U >= 2d >= 2`, `D >= 2E`, and `D-E >= d`. Consequently U strictly increases under every O/Z/J move and D never decreases. A bounded-U profile inventory is finite, and every predecessor of a retained state satisfies both retained bounds.

If a descendant first reaches a larger target D=n, stripping its terminal O moves leaves a nonempty extension ending in Z or J. A Z output satisfies `2D-U >= 3d >= 3`. A J output satisfies `2D-U = D_old-2E_old+U_old >= 3`, apart from U_old=2, which is precisely the seed and gives equality 3. Hence the resulting core has U<=2n-3. Every prefix on the path to it remains inside the target table. This proves the reverse direction of table membership; the forward direction is an actual suffix retained by the search.

The equality branch D_s=n must precede the U cutoff, because arbitrary terminal O moves can push U above that cutoff without changing D. The author implementation does this correctly. The n=1 branch also behaves correctly: D-tail membership is exactly equality, while deletion/short-completion routes are still applied separately.

Merging identical numerical profiles in the table is safe here. Every O/Z/J suffix remains physically executable at every source having that profile, its future profile is determined by that profile, and both the section test and target are profile predicates. The algorithm does not claim to preserve a particular graph isomorphism class after replacement.

### 2. Maximal retained prefix

Nested descendant sets give nested section-target sets. Thus feasible prefix indices form an initial interval, including the possibilities that none or all are feasible. The implementation first checks the root, then binary-searches with a feasible lower endpoint and an infeasible sentinel one past the final prefix. It returns exactly the largest feasible index.

The optimization is the number of terminal letters removed from the **supplied presentation**. It is not an optimization over alternative presentations, arbitrary vertex deletions, total edit count, or isomorphism classes. The theorem and README preserve this distinction.

Root feasibility is precisely `n in {1} union C1`. The author implementation reconstructs an unsigned root witness directly from the D table. For n>1 its last non-O birth has a physical parent whose J closure has the same count. Replacing that final birth by J is the correct C1 reconstruction. A failed C1 query says nothing about C2 or C3.

### 3. Two infinite families and their optimality

The identity b.OJ=3b and linearity prove `t.Z.(OJ)^r = 3^r b`. This is a physically presented scaled family, not an inference that arbitrary coordinate division preserves physicality.

- For target 3*3^r, the old OZ suffix has f>0 and target below D, so zero deletions are impossible. The OJJ output is in the section with the desired target. Retaining b.O leaves no zero- or one-letter section completion, while JJ succeeds. Deletion count one and replacement length two are therefore optimal.
- For target 33*3^r, failure of target 33 from b.O excludes every prefix longer than the retained b. The complete U<=63 inventory contains 55 profiles. Exactly 34 additionally satisfy D<=33, and their D-values are `{3,9,15,19,27,31}`. Thus 33 is absent. The core bound makes this an all-height exclusion, not a depth guess.
- Exhaustion of all 1,093 words of length <=6 from b finds no section output at target 33. At length seven the two successful words are `ZZZZOJJ` and `ZZZZJZJ`; in particular the claimed word is valid. Thus the optimal retained prefix requires a replacement of length seven.

All these suffix tests transfer to every scaled family member because applying the **same literal word** to 3^r b multiplies its profile by 3^r, and section membership and target scale homogeneously.

### 4. Important scope of suffix optimization

The general implementation minimizes retained-prefix loss; it does not minimize its appended suffix over all physical suffixes. For example, `--repair-word ZOOZ` correctly retains `Z`, deletes three letters, and returns the nine-letter suffix `ZOZJOJJJJ`, producing profile `(573,202,674,303)` at target 33. This is a valid result. The explicit optimal seven-letter family repair is a separately established result.

This distinction is already stated in Section 2 of the proof. A future interface that advertises shortest suffixes would need the optional bounded breadth-first search described there, or an explicit family specialization. There is no need to change the theorem or current implementation to satisfy their actual claims.

## Independent computational checks

`check_independent.py` uses positive coordinates `(a,f,e,p)` throughout its reference arithmetic rather than the author's raw `(D,E,U,d)` recurrence. It constructs its own forward descendant inventories and uses a separate adjacency-matrix graph builder with recursive subset path counts. Author arithmetic and graph-counting functions are not reused as reference calculations. The author module is imported only as the system under test.

Exact results:

- 17,536 D-tail queries: all positive odd targets <=63 on all 548 source profiles with U<=127. Each result agrees with independently enumerated descendants under the complete target bound, and every returned suffix is replayed.
- Additional equality-branch checks after 80 O moves, well beyond the target table height.
- 1,024 direct section comparisons: all odd targets <=31 on all 64 source profiles with U<=31. Forward inventories use the inherited complete section height bound and total 81,460 profile visits. This is not merely comparing the author's classification with itself.
- All original words through length seven whose positive signed target is <=301: 1,237 cases, comprising 1,126 successes and 111 exact no-C1 cases. The independent reference checks 9,025 prefix flags, interval structure, maximum cut, returned section witness, unsigned count witness, and the C1 recipe.
- 56 explicit-target cases cover empty words, target 1, target 7, already-section sources, unchanged-prefix repairs, and large O-only presentations.
- Infinite-family identities checked at r=0,...,100. The all-r proof is linearity, not this finite range.
- 1,093 short-word exclusion checks and all 2,187 words of length seven for the target-33 minimum.
- 13 literal rooted tournaments through order 12. All four response coordinates and the unique root-start Hamilton path are independently counted. Relevant C1 witnesses also have exactly one Hamilton cycle.

Selected literal profiles:

| Word from t | Order | (D,E,U,d) | Hamilton cycles |
|---|---:|---|---:|
| ZOZ | 5 | (9,1,11,1) | 1 |
| ZOOZ | 6 | (15,1,23,1) | 1 |
| ZOJJ | 6 | (15,6,18,9) | 1 |
| ZZZZZOJJ | 10 | (165,66,198,99) | 1 |
| Z | 3 | (3,1,3,1) | 1 |
| ZZZZZ | 7 | (33,1,33,1) | 1 |
| ZOJZZZZOJJ | 12 | (495,198,594,297) | 1 |
| JZOZJ | 7 | (33,14,34,15) | 1 |

## Reproduction

From the common workspace parent:

```bash
python3 session201-prefix-repair-audit/check_independent.py > session201-prefix-repair-audit/independent-normal.json
python3 -O session201-prefix-repair-audit/check_independent.py > session201-prefix-repair-audit/independent-optimized.json
cmp session201-prefix-repair-audit/independent-normal.json session201-prefix-repair-audit/independent-optimized.json
python3 session201-prefix-repair-audit/frozen-inputs/verify.py > session201-prefix-repair-audit/author-replay-normal.json
python3 -O session201-prefix-repair-audit/frozen-inputs/verify.py > session201-prefix-repair-audit/author-replay-optimized.json
cmp session201-prefix-repair-audit/frozen-inputs/VALIDATION-normal.json session201-prefix-repair-audit/author-replay-normal.json
cmp session201-prefix-repair-audit/author-replay-normal.json session201-prefix-repair-audit/author-replay-optimized.json
```

All reproduction commands and comparisons passed. Checks use explicit exceptions and survive optimized Python. Both independent receipts are byte-identical; both author replays are byte-identical to the supplied frozen receipts.

## Editorial note, not a correction requirement

For clarity, the sentence “The search has55 profiles” could say “The U<=63 inventory has 55 profiles, of which 34 have D<=33.” The existing sentence is mathematically consistent, but explicitly naming the unfiltered inventory prevents confusion with the doubly filtered target table in Section 2. This audit did not edit the frozen proof.

The conditional existence qualification and the remaining all-odd C1/C2/C3 occupancy problem are correctly retained throughout the packet.
