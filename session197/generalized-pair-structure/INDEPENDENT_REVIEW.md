# Independent review: intrinsic generalized-pair structure

**PASS. No correction.** This report binds to `PROOF.md` SHA256 f234993840c00fae0130a2c1b9efc1b19a10601d7f2a6b8fe8fb952273a1fd87, `structure.py` 71e361848b5b2b6b1e6942fbc6ed8ed0619a64eb6be6acce3076a238d448a899, `verify.py` 749d17d10dde6449ed38e3e165ac8cb6479e1749b8114765faa6380a01b2dc52, and both validation receipts 8a99679075710c6e0b86379856108ffe9965122ff350ac468a35b7d0895fb46d. The complete twelve-file source manifest is preserved in `source/RELEASE_FILES.sha256`. Original source files were not changed.

## All-order proof review

The start-two SCC decomposition is used with the correct intrinsic initial component, not a vertex interval. Its tail is a chain of singleton SCCs. The two-entry case lies in C1 and has exactly the two possible pair types from the C1 block theorem. The triangle exception is correctly assigned to the adjacent case.

For adjacent entries, the source obtained by deleting the sole-out-neighbor vertex retains exactly one starting path at the remaining entry. The endpoint argument is reversible. Empty tail would give C1, so the N restriction v>=2 is necessary and sufficient. For a genuine nonadjacent pair, the two forced final vertices of every relevant path are deleted; the source has at least two vertices. The stated map of M's three added vertices to b,r,u preserves every edge. For one entry, the inherited marked C1/start-two clone lemma gives a unique smaller pair, and deleting both clone vertices leaves the required nontrivial start-one source. The Q map and singleton-gadget boundary are correct.

The P=Q(J) identity is an actual rooted graph identity, including the singleton V boundary, where J(V) is the two-vertex transitive source. It does not rely on equal Hamilton-path counts. These maps prove necessity, sufficiency and uniqueness of the N/M/Q presentation for a specified qualifying root. No unrooted uniqueness is asserted.

The canonical rooted generating function uses one rooted source and one binary tail, with gate increments 2,3,3. Its simplification to z^5(1-z)(3-2z)/[(1-3z+z^2)(1-2z)] is correct. P is redundant. Counting qualifying roots is not the same as counting unrooted Gamma graphs. In the finite C2 bank the only possible nontrivial automorphism exception is regular five, which has no start-two root; thus its presence does not invalidate the root counts.

The H77 and H289 exclusions follow from literal root/deletion-SCC information plus the proved criterion. H289 is outside Gamma in both orientations. Their same-order P repairs are scalar examples and do not undermine graph exclusion. The packet correctly makes no scalar omission or scalar-completeness claim.

## Independent computational evidence

`check.py` imports no author code. It uses recursive endpoint path counts, a separate anchored cycle DP, Kosaraju SCCs with explicit dominance-order checks, and literal appended-vertex orientations to rebuild all N/M/Q inverse maps. Its input `c2-inputs.csv` is the order<=10 C2 slice of the separately audited complete low-cycle census, not the author's start-two grammar regeneration.

The checker independently confirms all 865 C2 graph classes, 2,262 start-two roots, every qualifying source and full labelled graph map, and the following inventories at orders5 through10:

- Complete C2 classes: 2,5,18,60,191,589.
- Gamma graph classes: 1,5,18,57,169,482.
- Canonical qualifying rooted classes: 3,10,31,91,258,715.

It separately recomputes both orientations of the H77/H289 hostiles and the first61 coefficients of the exact rooted generating function. These finite checks corroborate the all-order argument; they do not replace it.

The copied source normal and optimized runs reproduce the frozen receipt exactly, including all352 P=Q(J) tests. The independent normal and optimized runs also match exactly. Reproduction:

    python check.py > independent-normal.json
    python -O check.py > independent-optimized.json
    cd source
    python verify.py > ../source-normal.json
    python -O verify.py > ../source-optimized.json

The source verifier writes its receipt beside itself, hence the copy. The independent checker only prints its receipt.

Exact independent hashes:

- check.py: f07efcb3aa2ef7f9ffbb3cd49c36fd3a89ac1e34ec1fc4a4728816084b48490d.
- c2-inputs.csv: 452fc0ef1773a48c8cbb82b51e3003c114536986d2bed841ff938d0c79d45065.
- independent-normal.json and independent-optimized.json: 904c882ecb68fca76a3d3ba77f5dc46c91f932fc0b3a41429f09f69ce4411d6d.
- source-normal.json and source-optimized.json: 8a99679075710c6e0b86379856108ffe9965122ff350ac468a35b7d0895fb46d.

No source correction or mathematical blocker remains.
