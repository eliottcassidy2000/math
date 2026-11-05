# Independent factorwise inverse-bridge audit

2026-11-05. **PASS for the finite factorwise theorem and the two-sign
H15 unit-cloning obstruction.**

Frozen hashes:
- FACTOR_BRIDGE.md: `7a07a5f6ab451287ec11640c3caf01ea70ae18d61d13e8b29ab4302c3ee9f000`
- verify_factor_bridge.py: `619d1cb9e67a1f24fd7f507b9d07642f271096aecf5863cbc30bb384cde86fdc`
- receipt: `cfa78323ca8a1f358ad1276a47113abdfebe65c930f33c9a4b03e1d908ee8812`

## Proof audit

The inherited complete star construction preserves the two separate
physical factors `alpha=U*L` and `beta=d*g-f*t`. The new argument does not
replace actual rooted pieces by a numerical determinant solution.

For U>0, endpoint reinsertion gives `d+2f<=2U-1`. A negative unit determinant
forces f>d, hence f>=d+1 by parity and `3d<=2U-3`. For a positive unit
determinant with f<d, writing `1=(d-f)t+d(g-t)` forces t=g=1 and d-f=1.
The only remaining positive-unit case again has f>d. Consequently both
the stated uniform d bound and `f<=U-1` are valid. The transitive-partner
conclusion is correctly scoped to the positive source-offset subcase.

For X, normalization gives an actual core plus j K-vertices, preserving
the whole profile. Its positive amplitude satisfies `U=2^j U_core`, so
`j<=v2(U)`. Monotonicity of k and ell then gives exactly B_X. Cases with
`D_beta(U)<1` are rejected before looking for a positive odd deletion count.

For the partner, applying normalization to the converse uses the endpoint
count `starts_c(Y)<=t`, giving nonnegative slack `F<=L-g<=L-1` and
`U_W=L-t`. If U_W=0, the physical inequalities give g=t=L and the unit
determinant forces L=1. The source-rooted transitive partner may then be
replaced by a singleton without changing (g,t,L). This argument is not
asserted for nonunit determinants. Otherwise t<=L-1 and
`j<=v2(L-t)<=ell(L-1)`. Thus core plus tail has order at most
`L+k(L-1)+2ell(L-1)`, exactly B_Y. Taking the converse back preserves every
quantity needed by the star identity.

Since X and Y share one root and F_1 adds one vertex, the certificate order
is |X|+|Y|, as claimed. The factorwise exhaustive search is finite and
complete for the globally classified positive binary architecture; the
balanced-factor observation is explicitly conditional on an actual
physical balanced factorization.

Zero-amplitude sources are excluded by hypothesis. For example singleton
pieces can yield the constant count-one family, where the positive-source
inequality would not apply. For a count-level test, a nonpositive or
nonintegral `(N-beta)/2` is simply not a candidate positive integer
amplitude. No zero case undermines the theorem.

## H15 and the new slack-six table

At H15, beta=+1 means alpha=7, already excluded by the independently
audited prime-amplitude theorem. For beta=-1, alpha=8 and all factor
possibilities have U=1,2,4,8. The budget removes the first two and leaves
only (U,d)=(4,1),(8,1),(8,3). Deletion-one sources with power-of-two U
have f=0, contradicting the negative determinant. At d=3,U=8, the only
remaining even slack values are 4 and 6. Neither exact table has a base
whose dyadic multiples can equal 8. Both signs are therefore impossible.

The d=3,f=6 enumeration is complete for the same reasons as the audited
primary finite-core theorem: the deletion graph has one triangle SCC,
at most four outneighbor singleton components, and at most two interior
singleton components. Every position and root incidence is covered.
The claim does not exclude H15 itself, isolated clone coincidences, or
non-binary higher-run families.

## Independent computation

`audit.py` uses the prior independent forward-path/boundary-cover engine,
with a separately constructed star adjacency matrix and independently
rebuilt normalized pieces. All 42 distinct unit profile pairs from the
complete rooted labelled banks through order four satisfy both bounds
and preserve all required coordinates. Three actual clone multiplicities
per pair give 126 direct star count checks.

The independent boundary-cover classification gives exactly
`{9,15,27,36,60}` for d=3,f=6, and all five witnesses pass literal
permutation/cut profile checks. It inspects 753 positive-outneighbor
normalized candidates; the source's count 754 additionally includes the
irrelevant zero-outneighbor core with U=f=0. This explains the one-count
bookkeeping difference without any difference in the classification.

As a separate physical control, all four strong H15 classes in the
independently complete order-eight census were inspected at every root.
Every positive binary root, detected by its O-to-I star graph, has nonunit
offset when counted after actual cloning. These controls include the
order-six examples with H15; no false lower bound of 17 is used.

The source verifier passes normally and under `-O` with identical output.
No mathematical or implementation blocker was found. Results and witness
data are in `AUDIT.json`; reproduction uses the sibling independent
prime-atoms audit engine and H67 census. Author files were left untouched.
