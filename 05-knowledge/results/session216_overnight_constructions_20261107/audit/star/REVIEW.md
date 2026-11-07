# Independent audit: matching support, star contexts, and four choices

Status: APPROVED mathematical proof and stated finite certificates. All replay receipts below passed. No interval-width amplification, eventual spectrum, or all-target coverage theorem is approved.

## Matching support and normal form

The support theorem is valid for every nonempty marked tournament, including a singleton. A j-visit skeleton supplies j-1 vertex-disjoint outside paths from A=N+(v) to B=N-(v); each contains an A-to-B arc, and selecting one per path yields a matching. Conversely a matching of size j-1 gives a directed path alternating j root clones with its matched two-vertex separators. Inserting remaining outside vertices never changes the relative order of that path and hence never destroys any existing separator. The final path has exactly j separated clone visits. This proves every degree1 through nu+1 occurs. Degree0 is H(G-v)>0, including H(empty)=1. The criterion describes support, not coefficient multiplicities.

A bipartite graph of matching number one is a star; the elementary shared-endpoint proof is correct. If its center is in A, every other A-vertex loses to every B-vertex, yielding L=>v=>R and L=>R after removing the center. A center in B becomes this case under converse, which preserves all path-cover coefficients. Empty sides are allowed in the formulas; a positive second coefficient requires an actual out-to-in edge. The assertion that a nontrivial strong core has that edge follows because otherwise its marked vertex is a singleton SCC.

## Factorization and boundaries

Removing the single apex x splits a path into at most two ordered pieces, each respecting L,U,R order. The four before/after-x cases are disjoint, and their root-containing pieces give the stated alpha,beta,tau factors. Empty ordinary pieces are properly included: beta includes the empty first R-piece, alpha requires nonempty second L-piece, and tau requires both L-pieces. Thus c2=alpha beta and c1=2 alpha beta+hL hRx-hR tau.

For nonempty U this gives alpha beta Q_U(2)+(hL hRx-hR tau)H(U). No independent control of these factors or the signed offset follows. Empty L and R, empty cuts, singleton cores in the matching theorem, and nontransitive modules were explicitly tested. The substitution identity retains the nonempty-U restriction.

The transitive formulas for alpha and tau are correct. In tau, endpoint pairs in one order contribute2^(gap-1); the reversed order is possible only for adjacent vertices and contributes1. Their combined expression is even. These remain coupled quantities, not free binary coefficient and free correction.

## Returning append update and physical menu

Appending a sink z to cut-labeled L gives exactly

    (h,alpha,tau,e) -> (h,2alpha+epsilon*h,tau+epsilon*(alpha+e),epsilon*h).

The proof separates z's occupied path, including singleton Q={z}; no case is lost. It also works for empty L. Common h,alpha,e therefore stay common while tau differences are fixed. The four physical seed states are (9,245,150+2j,9), j=0..3. Their arbitrary common binary-word descendants retain four choices with spacing2H(U), specialized to four consecutive odd counts for transitive U. This is a returning source-state construction, not closure on unmarked output graphs or a widening operation.

The supplied initial core cycles are valid. Each newly appended L vertex receives an edge from every old L vertex and points to the root, so it attaches strongly to the existing strong core regardless of its apex incidence. All four seeds have nonempty L. Subsequent nonempty module substitution preserves strongness. The bound k+m+9 follows from maximum initial L order8, one apex, and a k-vertex module.

## Independent checks

The independent coefficient routine uses contracted polynomial transitions and subset DP, not the producer's permutation/gap enumeration. Matching sizes are computed by exhaustive disjoint edge subsets. All5405 rooted labeled tournaments through order5 pass; histogram nu=0,1,2 is1215,3350,840.

The factor audit computes endpoint/cover data by literal permutation enumeration. All5625 cut-labeled pairs L,R through order3 produce16,875 composites using TT1,TT2, and a directed-triangle root module. A separately written C++ pull subset-DP counts all of them. The triangle cases go beyond the producer's11,250 controls.

All11 interval decomposition records were independently reconstructed, including converse, star-center membership, L-to-R edges, and every scalar factor. All255 transitive cut cases through order7 match the two closed formulas.

The menu audit independently recounts all four seed states and cycles, all60 seed/word states through word length3, and150 append updates on arbitrary small cut-labeled L including the empty graph. It builds180 closed graphs with TT1,TT2,C3 modules. All counts and all180 strongness tests pass. Combined literal total is17,055; maximum menu order15, so15! bounds all path counts safely within unsigned64-bit arithmetic.

Normal and optimized independent receipts and optimized/debug literal receipts are compared. The producer's degree, factor, and menu controls are replayed normally and with -O against their frozen exact receipts. See the corresponding replay files.

## Provenance and reproduction

Occupied-gap/path-cover machinery is inherited as stated. No literature-priority conclusion is audited. The finite interval factors use the separately approved literal witnesses copied to dependencies/interval_witnesses.json; bank completeness is irrelevant.

    python independent.py > independent-replay.json
    cmp independent-replay.json INDEPENDENT.json
    python menu_independent.py > menu-replay.json
    cmp menu-replay.json MENU_INDEPENDENT.json
    g++ -O3 -std=c++17 literal.cpp -o /tmp/session216-star-audit
    cat GRAPHS.txt MENU_GRAPHS.txt | /tmp/session216-star-audit
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

frozen/ preserves producer files, and FROZEN_SHA256SUMS pins them. No producer files or inherited repository files were edited. No publishing was performed; the parent remains sole publisher.
