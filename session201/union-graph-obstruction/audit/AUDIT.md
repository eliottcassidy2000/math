# Independent audit: union witness graph obstruction

Verdict: **PASS**, within the precise index-preserving witness-relation scope. No mathematical correction is required. The input remains untouched and retains its own RESERVED / UNREVIEWED label; this external audit is not a source edit or publication.

## Frozen input and independence

The complete seven-file input packet was copied to `input/`. `INPUT-SHA256.txt` freezes every original file, including its manifest. The original manifest also verified successfully. Proof SHA-256: `548be4c5ad0918ccb7a67de1e288504ff50b4af0450b04bef647c00aaca18953`. Menu SHA-256: `01dce613993c7c0dd37ec27d47d6644cb30b43b92abe14e784a99e3a36b2a81d`.

`independent_verify.py` imports no author code. It reconstructs every menu coefficient by applying chronological physical scalar operations to basis vectors and evaluating the scalar bilinear pairing. It independently checks positive-coordinate transformations, reconstructs residual coefficients, constructs actual unit sources, and regenerates finite controls. Basis vectors used for polynomial coefficient extraction are not asserted to be physical sources. `math/AGENTS.md` was read. The dirty/stale checkout was not modified; no GitHub action was taken.

## All-order mathematical audit

1. **Coordinate isomorphism and boundary.** Substitution verifies both coordinate maps are inverse. All displayed O/Z/J actions in (a,f,e,p) are exact and coefficientwise nonnegative. Root size d=e+p is unchanged by O and Z; under J it becomes a+2e+p, an increase of a+e. The seed has p=1 and all other positive coordinates zero. Hence every reachable source has nonnegative coordinates and d>=1, including boundary faces with a=0, f=0, e=0, or p=0. No strict positivity of those four individual coordinates is used.

2. **Unit-source membership, made explicit.** The proof calls tau(n) an actual source. This follows inductively from tau(2)=seed, tau(2h)=O tau(h) for h>=2, and tau(2h-1)=Z tau(h) for h>=2. The given recurrence gives each first coordinate exactly; Z sets E=d=1. Parents are strictly smaller for every n>=3. Thus every tau(n), n>=2, really lies in the same physical cone. The n=1 recurrence base is not a source-domain claim. This fills the implicit boundary justification without adding a hypothesis.

3. **Residual matrices and lower bound.** Independent coefficient extraction reproduces all sixteen entries of each displayed R_O, R_Z, R_J. Each R_c minus the matrix with four 2s in the lower-right block is entrywise nonnegative. Therefore y^T R_c z >= 2(e_y+p_y)(e_z+p_z), including all boundary faces. Since d is nondecreasing along the target word, telescoping proves the stated 2|W|d_s d_t bound. Empty words give equality zero; every nonempty word gives a strict increment. The constant 2 is sharp for a one-letter increment at the seed. The pairing Q is symmetric; transforming either argument from seed consequently also increases K. Since K(seed,seed)=3, all actual-source outputs are positive natural numbers, so no signed-output encoding issue arises.

4. **All twelve forms.** Every frozen menu matrix equals Q times the chronological physical word action. The menu contains twelve distinct words and empty exactly once. No assumption about uniqueness among other outputs is needed. The audit verifies the matrices and their words; dagger-partner metadata is not a dependency of this theorem and is not independently certified here.

5. **Minimum selection and synchronized padding.** For each fixed (u,v), empty is the unique minimizing template and its output is strictly below every nonempty template's output. Treat an automatic relation as the regular language of equal-length binary tracks, allowing common leading zero columns. A canonical synchronized language can be saturated with arbitrary common zero prefixes by regular concatenation; canonical representatives can be recovered by intersecting with the regular canonical-length condition. Lift R to four tracks by inverse letter projection, intersect with the regular comparison m<n, and project away m. Projected witnesses may have extra leading zero columns, so normalize by left quotient with the regular language of all-zero columns (and optionally re-saturate). Regularity is preserved. Complement relative to the regular representation domain, then intersect with R, yielding precisely the fiber minimum. In this minimum test m<n additionally means an m-witness never needs a larger numerical width than the retained n. Fixing u=2 is a regular section. Thus the closure claim holds for both padded and canonical conventions, not only a fixed machine width.

6. **Nonregularity without external dependency.** The fixed-point identity g(2^r+1)=2^r+1 is valid for r>=1. For b=1 in (6), the even parent is 2^a+2, whose half is 2^(a-1)+1; the identity applies because a>=2. For b>=2 the parent is the preceding b-case and is odd, so the defect n-g(n) is unchanged. This proves (6). The word A Z^k B C Z^(m-1) A decodes exactly as stated, with a=k+1 and b=m+1. Its membership is equivalent to 2^k=2^m, hence k=m. The prefix family is pairwise separated by the stated suffixes, giving infinitely many Myhill–Nerode classes. Both tracks start with 1, so the obstruction already lies among canonical encodings; admitting leading zeros does not remove it.

7. **Reduction from h to g.** The relation n=y+v-(v mod 2) is finite-carry regular (equivalently n+(v mod 2)=y+v), with parity recognized from the last bit. If h had an automatic graph, intersecting its lift with this relation and existentially eliminating n recovers exactly g on v>=2. The same zero-prefix normalization handles any eliminated track that is longer than the retained tracks. This contradicts item 6. No cited literature theorem specific to g is needed.

## Computation receipts (FINITE-EXACT)

All four runs exited successfully:

- Author normal and optimized Python: 364 physical states, 30,576 strict-increment checks, 1,024 union fibers. Both receipts exactly match the supplied validation JSON files.
- Independent normal and optimized Python: 1,093 chronological source words of lengths 0 through 6; 213,135 sampled source-pair generator checks; 4,095 constructible unit sources n=2 through 4096; 9,216 union fibers u,v=2 through 97; 600 sparse-formula cases a=2 through 25, b=1 through 25; 1,600 separating-word pairs k,m=1 through 40.
- Normal and optimized outputs are byte-identical for each implementation. Checks use explicit exceptions rather than removable assertions.

The 1,093 source count counts words, not distinct states. Source pairs use every first source and the second-source list sliced with stride 17. All-order claims are justified above by algebra and language closure, not extrapolated from these finite universes.

## Exact scope and non-implications

**PROVED by the audited argument:** the ternary relation of exact (u,v,n) witnesses for this twelve-word union is not binary-automatic. The same proof works for any collection of finite O/Z/J words containing empty, finite or infinite, because every nonempty word strictly increases the same pairing.

**Not concluded:** nonregularity of the projected set of attainable n; failure of all-target occupancy; failure of cofiniteness; an obstruction to unbounded-state algorithms. A nonregular relation can have a regular projection. The twelve-template occupancy assertion remains OPEN. This audit does not re-audit the separate digital-graph companion or its 2-regularity claim.

## Reproduction

Run from this audit directory:

    python3 independent_verify.py > INDEPENDENT-normal.json
    python3 -O independent_verify.py > INDEPENDENT-optimized.json
    python3 input/verify.py > AUTHOR-normal.json
    python3 -O input/verify.py > AUTHOR-optimized.json
    cmp INDEPENDENT-normal.json INDEPENDENT-optimized.json
    cmp AUTHOR-normal.json AUTHOR-optimized.json
    sha256sum -c MANIFEST.sha256

`RECEIPTS.json` records statuses, input integrity, and comparison results. The audit manifest covers all deliverables and the input snapshot, excluding itself and its own verification-output file.
