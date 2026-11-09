# Independent audit: THM-4942 first-creator generation

Source revision: `f870c2abf10f08ac8648923e313ca417e5e797f3`.
Audit date: 2026-11-09.

## Current-status update after the baseline review

The review below audits the exact baseline f870c2abf. Concurrent commit
`eea6f39cc` subsequently supplied the original [full independent audit](../../first_creator_generation_20261108/audit/AUDIT.md), applied the orientation and global-scope repairs, and established class-level creator genealogy at order eleven by an exhaustive deletion-count test. The missing-audit and untested-class-genealogy findings below are therefore historical, not current blockers. Session230 did not rerun that large computation. Its own97-witness audit remains independently reproducible here.

The same incoming audit recovered prior THM-4811(B): the semidegree/connectivity threshold below is a restatement and sharp-boundary elaboration of its existing cut argument, not a new theorem. We retain our independent proof and checks while correcting the inheritance credit. The original Busch-bound list also required restriction to S(11), or all odd h<=80,451, rather than every h<=95,095.

## Verdict

The stated finite five-case decomposition has a coherent, noncircular mechanism. Its order-11 exceptional witnesses independently check. This review does **not** independently rerun the roughly ten-billion-extension exhaustive census, so the completeness premise remains inherited from the recorded computation. The promoted canon links `audit/AUDIT.md`, but no such file exists in the source Git tree. This missing advertised audit must not be represented as a completed audit.

Several scope/wording repairs are needed; none presently refutes the five-case finite decomposition itself. Calling the first-creator inventory a **minimal generating set** is already false at order four under its own operations.

## Independent reproducible checks

Run from this directory:

```
python3 check.py > check.json
python3 check_witnesses.py > check_witnesses.json
```

`check.py` has no original packet code/data dependencies. It enumerates every labelled tournament through order six (2, 8, 64, 1,024, 32,768 at orders 2..6), uses direct subset Hamiltonian-path DP, all directed cuts, and exhaustive module subsets. It checks:

* The five-case property on this entire small universe.
* The sharp general degree/connectivity implication below.
* All 33,864 attachment evaluations for parents of orders 2..5 against separately recomputed extended-tournament DP.
* Sharp strong counterexamples at orders 10, 14 and 18 to extending the degree bound beyond its threshold.

Spectra recovered are S(3)={1,3}, S(4)={1,3,5}, S(5)={1,3,5,9,11,13,15}, and S(6)={1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}. Connectivity floors recovered: f1(3..6)=3,5,9,15; f2(5)=15 and f2(6)=41.

`witnesses_source.txt` is an exact copy of the tracked `lambda3_witnesses_order11.txt`, extracted with `git show HEAD:<path>`. `check_witnesses.py` independently checks each of its 97 entries using the same direct DP and all cuts, and groups by a separate degree-filtered partial-bijection backtracking isomorphism algorithm. It does not use `canon11.c` or ordered colour refinement.

Results: 97 valid entries, exactly 11 isomorphism classes. The seven entries-of-classes labelled `inS10` all have nontrivial proper modules. Their H/class multiplicities are 11475 (one), 15213 (two), 15273 (two), 15561 (two). The remaining four classes have H=14787 (one), 15637 (one), 15703 (two). All have lambda>=3. This verifies witness soundness and the class grouping, **not** completeness of the witness list or the supplied membership labels in S(10).

## Exact scope and code audit

Inspected `PROOF.md`, `README.md`, `gen.c`, `lam11.c`, `canon11.c`, the canon theorem and HYP-9532. Sparse checkout initially hid the packet; Git-object reads resolved that, so packet absence from the filesystem is not an evidence defect. The absent advertised audit is a genuine Git-tree absence.

`gen.c` classifies FC first, then EXT, SUM, SUB, REST. FC is membership complement in S(n-1); EXT tests deletion counts against S(n-2). Modules are detected by closing a pair under outside vertices that distinguish it, a valid test. It tests minimum directed-cut size directly on REST. Its binary path trusts supplied Hamiltonian/deletion counts and strongness from `dump10.c`; the text path recalculates counts. Inputs, not merely the classifier, must therefore remain in the reproducibility chain.

`lam11.c`'s degree>=2 parent filter is safe: deleting a vertex of any child with all in/outdegrees>=3 loses at most one from each degree. All attachment patterns on an isomorphism representative reach every extension isomorphism class. Its child degree checks use the correct bounds [3,7], and its full cut test safely rejects cuts of size <=2. Restricting to H<=15745 is sufficient to find all non-first-creators because 15745 is max S(10). New values below this maximum must still be separated by membership in S(10), as the packet does.

The explicit attachment polynomial is mathematically correct: classify Hamilton paths by the new vertex's position. Coefficients include split path covers without requiring an old joining arc. Hence it avoids the invalid historical assumption that deleting a vertex from a Hamilton path always leaves one.

`canon11.c` uses invariant ordered colour refinement and exhaustive individualization; its bytewise integer comparison ordering is platform-dependent but still a fixed total order on the small nonnegative signatures, not an isomorphism defect on a fixed platform. Our independent grouping avoids reliance on it.

The reported first-creator value genealogy is **existential per H**: each new h has at least one extension representation with a first-creator parent. It does not say every class realizing h has such a deletion. It does not establish that a chosen single representative per value generates all other first-creator classes.

## Repairs and smallest witnesses

1. **Not a minimal generating set.** Let C3 have one vertex replaced by TT2 (a two-vertex transitive tournament). The result is the strong order-four tournament, with H=5; 5 is new because S(3)={1,3}. Yet this first creator is generated by substitution from C3, TT2 and singletons, and TT2 is an ordinal sum of singletons. Thus the complete first-creator inventory contains redundant generators already at order four. Describe it as a distinguished seed inventory, not a proved minimal basis. Specify a singleton seed (or explicit order-one convention for H=1); seeds only of order>=3 cannot construct orders one and two under order-nondecreasing operations.

2. **Neutral-pattern orientation.** PROOF section 1 defines x_a=[a->v], but then calls its neutral support a condensation down-set, importing the convention in THM-1900 for P={a:v->a}. These sets are complements. With component order C1->...->Ct, the x=1 supports are prefixes, while P supports are suffixes. Smallest witness: old edge a->b, x={a}; the insertion gives a->v->b and H remains 1, although {a} is not a suffix. Replace the sentence with “the complement of the x-support is a condensation suffix.”

3. **Moon-chain overstatement.** PROOF section 2 explicitly disclaims class-level first-creator genealogy at order eleven, while the section-5 analogy table suggests chains inside the seed inventory through order eleven. Restrict that conclusion to order ten, with the stated 63 exception, or use value-level existential language at eleven. Ordinary Moon strong-subtournament chains do not imply first-creator chains.

4. **Global strong-first-birth inference.** PROOF section 4.1 concludes globally that first creators are the members of Strong(h) at least order, although section 5 identifies first-creator strongness only as finite-exact through eleven and points to HYP-9481(b). Keep the strong-fiber antichain theorem (valid), but qualify the first-birth identification by the strong-first-birth hypothesis, or restrict it to the verified range. The h=1 singleton/transitive convention also needs an exception to the displayed [3,B(h)] window.

5. **Unnecessary canonization requirement.** Class-level genealogy does not require canonical order-eleven forms. For each child reached by the exhaustive order-ten class/attachment scan, compute all eleven deletion H-values and test whether any lies outside S(9). A zero-counterexample result proves the universal class property despite duplicate representations. Canonization is needed to count counterexample classes, not prove their absence.

## Why the generation claim is not circular

Use a precise grammar: singleton and all designated FC classes as leaves; ordinal sums of smaller objects; proper substitution (quotient size and every nontrivial module size strictly below the result); arbitrary one-vertex attachment **only when the immediate parent is a designated FC**; and a **coupled** operation that takes a proper ordinal sum and reverses at most two cross arcs. For each non-FC terminal object, select its certified decomposition and recurse on smaller blocks. This decreases vertex count. In particular, the arc-repair operation does not recursively request an arbitrary same-size input.

Unrestricted substitution with a quotient equal to the target and singleton modules would be tautological. Unrestricted vertex insertion would generate everything from the singleton and make first creators irrelevant. Neither relaxation is justified by the intended census. Explicitly prohibit both.

Assuming the recorded finite census, induction with this precise grammar proves generation through eleven. The source's “first-creator” branch treats **every** first-creator class as available, not one representative per new value. The finite theorem is consequently weaker than the user's possible minimal-one-object-per-value interpretation.

## All-order structural extension: sharp semidegree threshold

**PROVED.** Let delta0(T) be minimum indegree or outdegree, and let lambda(T) be the minimum number of arcs from one nonempty side of a cut to the other. For any integer k>=1 and 2<=n<=4k+1,

    lambda(T)>=k  iff  delta0(T)>=k.

Singleton cuts prove necessity. For sufficiency suppose a cut X,Y has r<=k-1 arcs Y->X, with a=|X|. Summing indegrees over X gives

    ak <= a(a-1)/2 + r.

If 1<=a<=2k, then ak-a(a-1)/2=a(2k+1-a)/2>=k, contradiction. Thus a>=2k+1. Summing outdegrees over Y gives |Y|>=2k+1, forcing n>=4k+2, contradiction.

For every k>=2 the order threshold is sharp even among strong tournaments. Take two cyclic regular tournaments A,B of order q=2k+1, orient every cross arc A->B, and reverse one cross arc. Every vertex has both degrees at least k, with minimum exactly k; the reverse cut has one arc, so lambda<=1. The result is strong: each block is strong, the unreversed cross edges give A-to-B access and the single reversed edge gives B-to-A access. Hence lambda=1<k at n=4k+2. For k=1 the nonstrong ordinal sum of two C3 blocks gives the sharp obstruction at n=6; no strong obstruction exists since strong is equivalent to lambda>=1.

This recovers n<=13 for k=3 and proves sharpness at n=14. It is a structural all-order generalization, not a claimed new literature result.

A separate all-order observation clarifies the floor argument: for any n,k, if f_k(n)>max S(n-1), every non-first-creator has lambda<k and is at most k-1 cross-arc reversals from a proper ordinal sum. This follows directly from f_k's definition and needs **no** n<=4k+1 restriction. The degree threshold only lets one replace the cut test by score bounds.

## Open obligations after this review

* Reproduce/manifest the complete large parent-class universe and all census slice ranges, or retain its inherited finite-exact status explicitly.
* Supply the missing original audit or replace its dead link with a real audit; do not claim this small independent check reran the full census.
* Test order-eleven class genealogy with the direct universal deletion test.
* Define the operation-relative meaning of “minimal”; the present complete FC inventory is demonstrably redundant.
* No counterexample to the stated five-case property through order eleven was found. Nothing here establishes that property beyond eleven.
