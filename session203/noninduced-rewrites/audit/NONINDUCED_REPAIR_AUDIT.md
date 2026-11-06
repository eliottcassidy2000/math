# Independent audit of the 775 repair and source-ray diagnostic

6 November 2026. **PASS in the exact stated scopes.** No eventual coverage, threshold, or all-input extractor is established by these results.

## Positive physical witness

The proposed construction takes E(129), reverses c -> p0, then adds two vertices with outgoing masks 11 and 83 in sequence. Independent pair-by-pair graph construction gives precisely the published outgoing rows

    (2,0,4091,3059,2019,4035,1923,3843,3587,3075,2059,83).

Literal recursive enumeration, without memoization or subset dynamic programming, gives exactly 775 Hamiltonian paths. The reversal makes the old tournament transitive, so no undefined sector subtraction is involved. This certifies one physical repair, not a globally minimal edit count for E(129).

## Exact identity and hostile boundary

For the same recipe on E(2^k+1), k >= 4, independently rebuilt literal graphs have H=29,91,267,775,2259,6631,19587 for k=4 through 10. The k=6 and k=8 targets are 391 and 1543, so the repaired masks must not be extrapolated as an extractor.

The exact all-k identity is also sound. Retain the nine-vertex k=4 core. Every extra private vertex is a transitive clone beating exactly core vertices r,c,x,y, whose core indices are 0,1,7,8. The independent checker verifies this embedding explicitly for each control length. The named construction makes the same statement immediate at every larger length.

For each ordering of the core, classify the clone gaps into allowed and mandatory gaps; any wrong core adjacency must be repaired by a nonempty allowed gap. The independent enumeration of all 362,880 core orders gives these only contributing classes:

    (allowed, mandatory) -> multiplicity
    (1,0) -> 1
    (2,0) -> 10
    (2,1) -> 8
    (3,0) -> 18
    (3,1) -> 8

Each of m named, internally transitive clone vertices chooses a gap; the internal order is forced. Inclusion-exclusion for the mandatory gaps therefore gives

    H = 1 + 10*2^m + 8*(2^m-1) + 18*3^m + 8*(3^m-2^m)
      = 26*3^m + 10*2^m - 7,
    m = k-4.

This proves every clone length, rather than interpolating the sampled values. Equality to 6*(2^k+1)+1 reduces to 13*3^m-43*2^m=7. Its values at m=0,1,2,3 are -30,-47,-55,7, and its next difference 26*3^m-43*2^m is strictly positive for m>=2. Hence k=7 is the unique hit.

## Restricted source grammar diagnostic

The separate SOURCE_RAY_BOUND.md argument is valid for its specified one-state affine source family. The O/Z/J matrices agree with the inherited positive-coordinate source maps. With F=(2,1,3,3), the coefficient vectors of 2FM-3F are (0,1,3,1), (2,3,1,1), and (0,1,3,1), respectively. Thus a nonzero nonnegative eigenvector with eigenvalue q requires (3/2)^length <= q.

For bases 3,5,7,9, the complete word-length bounds are 2,3,4,5. A separate checker using integer matrix products and handwritten Fraction-based reduced row echelon form, rather than the author's SymPy nullspaces, examines 12,39,120,363 words. Every nonzero kernel is one-dimensional. The complete positive-ray/matrix inventories match the author's receipt exactly; the maximum distinct matrices per common ray is 1,1,0,2.

For an affine family s(u)=uv+w in the nonnegative cone for every sufficiently large integer u, its nonzero slope v is nonnegative. A branch equation M s(u)=s(qu+d) forces Mv=qv and Mw=w+dv. One actual matrix supports at most one digit because v is nonzero. Fewer than q matrices on a common ray cannot supply a digit from every residue class. The claimed exclusions at these four bases follow. This does not exclude multistate grammars, nonlinear profiles, different parameter dimensions, or physical rewrites.

## Files and reproduction

Independent scripts: audit_repair775.py and audit_source_rays.py. Each has NORMAL and OPTIMIZED JSON receipts; the corresponding pairs are byte-identical. The checks use explicit exceptions, not removable Python assertions. REPAIR775_NORMAL.json retains the full core rows, physical witness, five-class gap histogram, exact coefficients, and all literal DFS controls. SOURCE_RAYS_NORMAL.json retains all exact rays and word classes. NONINDUCED_REPAIR_MANIFEST.sha256 freezes the audit and reviewed inputs.

    python audit_repair775.py
    python -O audit_repair775.py
    python audit_source_rays.py
    python -O audit_source_rays.py

No claim is made that the transitive-clone expansion or source-coordinate maps are new. The useful constructive result is the exact 775 witness; the family calculation identifies its precise failure boundary.
