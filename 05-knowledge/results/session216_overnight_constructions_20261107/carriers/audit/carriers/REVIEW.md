# Independent audit: common right factors and module probes

Status: APPROVED for PROOF.md Sections1–5 and their specified core certificates. No global completeness theorem, no new ambient coverage range, and no obstruction to asymmetric unbounded constructions is approved. Auditor: Session216 independent audit worker.

## Main all-order theorem

Let a physical paired source have equal (X,e,V,d), starts sA=1+3jA and sB=2+3jB, and adjusted bases X=BA+2jA=BB+2jB. For the same rooted right factor W, q=H(W-root) is a positive integer (H(empty)=1). For nontrivial W the inherited row law gives delta B'=q delta B, equal e',V',d', and multiplied starts. For q not divisible by3, assigning actual residue representatives rA',rB' in {1,2} gives

    delta X'=2(-q-rA'+rB')/3.

For q=1 modulo3 this is2(1-q)/3; for q=2 modulo3 it is-2(q+1)/3. Therefore equality holds exactly at q=1. At q=0 modulo3 neither start has the required nonzero residue. If q=1, old residue labels and j values survive, including the original budget. Larger q need not meet that budget; dropping the budget cannot repair the displayed inequality. Common O/Z continuation preserves the discrepancy because the other three row coordinates agree.

The readout separation scales from2 to2q, with the original A/B order retained. If q=2 modulo3 the residue labels swap; this must not be confused with unchanged output ordering.

For singleton W the right product is the identity and q=1. Its e-coordinate uses the explicit singleton exception, not e'=d eW. For nontrivial W, q=1 means W-root is transitive. Its vertices can be adjoined in transitive order, choosing O or Z according to their root incidences, proving exact physical equivalence to a common O/Z word. Conversely such words give a transitive right deletion. The cited SCC/Camion justification is valid. An optional elementary proof avoids that citation: any nontransitive tournament contains a directed triangle, whose three Hamilton paths inject into paths of the whole tournament by vertex insertion.

This is a theorem about common right factors of an already synchronized physical pair. It says nothing comparable for different factors, replacements, root changes, or initially unequal sidecars.

## Marked-module interface and finite source scope

The run decomposition is valid for nonroot replacement by a nonempty module. Number the path bins, then number visits in traversal order. Ordered module path covers are assigned to those visits without a factorial quotient. Adjacent placeholders in one path are forbidden, but placeholders in different bins are distinct visits. Empty ordinary bins are allowed; the designated rooted bin is never empty.

For transitive m-modules, ck=k!S(m,k), c1=1, and these functions are linearly independent by triangular evaluations at positive integers. Therefore coefficient equality proves equality for all transitive module sizes. Constant start requires the higher coefficients to vanish. This is exact symbolic algebra, not extrapolation from a finite m scan.

The complete tested universe is92 nonroot positions across twelve physical sources, grouped into six fixed pairs. Matching distinct positions between the two members gives exactly (759,1,1), (821,1,1), (4225,9,9), (4357,9,9). The first two marked vertices beat only the root and the latter two are global sinks. Their transitive expansions are exactly Z and O padding respectively. These finite position conclusions do not cover new source pairs, simultaneous modules, changing the root, or unrestricted nontransitive replacements.

## Independent implementation and literal controls

`independent.py` uses a different symbolic recurrence from the producer's placeholder DP: contract each potential module visit between external vertices into a degree-one polynomial edge. Add optional initial/final visits, then perform subset DP over external vertices alone. Actual external arcs contribute degree0 and bridges through the module contribute degree1. The resulting path polynomials are combined for ordered two-bin covers.

All eight profile polynomials at each of the92 positions agree coefficient-by-coefficient with the producer. Exactly the four stated position pairs survive. This tests the full finite source universe, not merely the successful pairs.

A separate forward-extension profile counter checks all32 selected literal substitutions, including directed-triangle modules. In addition160 actual rooted right products on randomly seeded small graphs validate every row coordinate and start formula, including25 singleton identity cases. Normal and optimized independent receipts are byte-identical. The independent checks use explicit exceptions rather than assertions.

Producer verify.py and its optimized invocation reproduce the32 module and2754 arithmetic controls byte-for-byte. The entire symbolic search regenerates matches.json byte-for-byte.

## Asymmetric finite searches

All three C++ searches were rebuilt and replayed. Source/right/candidate/key/match receipts agree exactly:

-12,893 sources,420 right inputs,266,768 products,203,170 keys,0 matches;
- same sources/rights/products with sink reset,203,168 keys,0 matches;
-167,609 prefix presentations,420 right inputs,3,467,984 products,1,588,044 reset keys,0 matches.

The code filters q>1, positive starts at most26, and nonzero residues. It inserts unmultiplied sources before products; later overwrites cannot erase the existence of a nontrivial product because all such products are processed afterward. Its start-preserving O/Z/J prefixes use the inherited nine-coordinate closure. These are searches of the supplied bank and stated operation menu, not exhaustive searches over arbitrary rooted tournaments. Signed64-bit arithmetic has ample headroom for these bounded inputs; no claim of arbitrary-height safe implementation is approved.

The provenance hashes for BANK.txt, inherited.py, ROW_NORMALIZATION.md and SESSION215_CARRIERS.md agree with the canonical files in session216-stage at its source revision. Exact checks and canonical relative paths are recorded in REPLAY.json.

## Freeze, scope exclusions, reproduction

`frozen/` preserves reviewed inputs; FROZEN_SHA256SUMS lists their exact hashes. The incidental ALL_PAIRS.txt and macro_search files were present in the workspace and frozen for traceability, but are not dependencies of the approved proof or included in the audited finite claims. Exclude them from publication unless separately justified.

From this audit directory:

    python independent.py > replay-independent.json
    cmp replay-independent.json INDEPENDENT.json
    python -O independent.py > replay-independent-opt.json
    cmp replay-independent-opt.json INDEPENDENT_OPTIMIZED.json
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

Producer commands should run in a separate copy because probe.py overwrites matches.json. The replay/ directory records exact reproduced receipts. No inherited files were edited and no publishing was performed. Any later revision must be checked against this freeze before being identified with this approval.
