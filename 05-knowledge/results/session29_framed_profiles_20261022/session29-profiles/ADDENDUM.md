# Session29 Phase2: genuine wiring edge counts, finite-state seams, and mirror-adjoint transfer

The Phase1 REPORT.md is frozen. This addendum strengthens its carefully limited incidence claims using a recovered repository theorem, then derives exact Hamiltonian statistics of the strong neighbors.

## 1. Genuine outward edges, not merely incidence budgets

Nearest proved mechanism: THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md, law L3. Its nonstrong-world independence is proved directly and does NOT depend on L4's flagged pure-blue classification gap.

For r>=2, every binary join T_w of Phase1 is nonstrong. Each Hamiltonian path P concatenates its SCCs, so its first vertex beats its last. Complementing all arcs outside P preserves the path and reverses that first-last chord, producing a Hamiltonian cycle in the target. The target is strong and lies outside the family. Thus no complement line has both ends in the family. Complementation is a fixed-point-free involution on tiles, so each family tile is incident with a distinct outward line, even when several lines land at the same strong node. Converse merging does not collapse tiles or parallel lines: it merges their node assignments only.

Consequently for each r>=2 the exact numbers of outward lines from this family to strong nodes are:

- all lines:18^r;
- blue lines:b_(2m)=18^m, b_(2m+1)=4*18^m;
- black lines:18^r-b_r.

Their block-length GFs, excluding r=0,1, are

E(z)=1/(1-18z)-1-18z;
E_blue(z)=(1+4z)/(1-18z^2)-1-4z;
E_black(z)=E(z)-E_blue(z).

These are genuine edges counted with multiplicity, not counts of distinct target isomorphism classes. The theorem supplies no target canonical-label or complete-wiring classification.

## 2. The seamed target and a two-block port formula

Take a source tile word t1...tr, each letter one of the18 spine tiles of A or B. Let X_i be the tournament obtained by complementing the nonspine arcs inside tile t_i. Across distinct blocks the target has all arcs pointing from later to earlier blocks, EXCEPT the one forward seam from the final vertex of block i to the first vertex of block i+1. This is an exact source-target construction, not an approximation.

For two blocks write the target as Y join X, except reverse b->a into a->b, where a is the exit of X and b the entry of Y. Let H be the Hamiltonian count, S_X(a) the number starting at a, E_X(a) the number ending at a, and

Q_X(a)=sum_{nonempty proper U containing a} E_{X[U]}(a) H(X-U).

This counts ordered two-path covers whose first path ends at a and whose second path is nonempty. Define Q_Y^start(b) analogously. Then

H_target=H_X H_Y-E_Y(b) S_X(a)
 +(E_X(a)+Q_X(a))(S_Y(b)+Q_Y^start(b)).

Proof: paths not using the new seam are the old join paths minus those using the deleted boundary b->a. A path using a->b has a unique block pattern among XY,YXY,XYX,YXYX; their four contributions factor as the displayed product. Internal subsets and path orders are disjoint choices, so there is no credit reuse or multiplicity ambiguity.

FINITE-EXACT: all324 two-letter targets satisfy this formula against independent subset DP on10 vertices. They are all closed, and exactly18 are mirror-fixed. Target H takes80 distinct values from221 to1049; the total H over all outward edges is195472, with blue contribution11522. These are edge-weighted statistics, not an unweighted target-class census.

## 3. Constructive finite-state proof, not recurrence fitting

At a prefix cut of a seamed target, the only forward arc to the unprocessed suffix is the next seam, from the current portal a. Restrict a final Hamiltonian path to the processed prefix. Its components are directed paths. At most one component can have an end other than a, because any such end cannot later acquire an outgoing edge and therefore must be the unique final endpoint. Hence at most two components occur.

A sufficient five-state interface is:

S: one path, a is its start;
E: one path, a is its end;
I: one path, a is internal;
DE: two paths, the active path ends at a and is nonsingleton;
DB: two paths, the active path is the singleton a.

For two components, order the active one first; this eliminates artificial double counting. To append a block, contract each old path to one atom. There are at most two old atoms and five new vertices, so every transition is obtained by enumerating directed path covers of at most seven atoms. Old-to-new arcs exist only from the active end to the new entry. New-to-old arcs are all available except entry->a, if a was an old component start. Old-to-old arcs are unavailable because such arcs would already have belonged to the previous prefix cover.

Every new cover contracts to exactly one old cover and one transition. Conversely every allowed atom cover expands uniquely. Pruning covers with two nonportal ends is safe because neither could acquire an outgoing edge later. Therefore the transfer is exact for every word length.

The DE and DB rows agree. With two old components, the active end must use the seam; otherwise both old ends would remain nonportal terminal ends. If that active path is a singleton, the exceptional prohibited arc from the new entry back to a would close a cycle, already excluded from a path cover. Thus the start-status distinction has no effect on any valid transition.

Writing state as (H,S,E,Q), with H=S+E+I and Q=DE+DB, gives four coordinates. In all18 letter matrices, E and Q further enter only through R=E+Q; exact verification is in check_compression.py. Therefore each tile t has an explicit3-by-3 matrix K_t and

v(t1...tr)=K_tr ... K_t1 e0,
H_target(t1...tr)=e0^T v(t1...tr), e0=(1,0,0)^T.

The first-column identity K_t e0=(H_X,S_X(a),E_X(a)+Q_X(a)) is independently checked for all18 tiles. All matrices are included in compressed-matrices.json. This finite certificate plus the path-cover proof establishes the all-word theorem; it is not an extrapolation from a few values of r.

## 4. Exact total edge-endpoint Hamiltonian GF

Summing the18 matrices gives

K = [[172,-24,408], [24,-4,68], [408,-68,706]].

The sum of H_target over all words of length r is g_r=e0^T K^r e0, with g0=1. Thus

G(z)=sum_{r>=1}g_r z^r
 =(172z+45144z^2+50400z^3)/(1-874z-43344z^2-50400z^3).

Initial values are172,195472,178348096,164357443072 for r=1,2,3,4. The genuine outward-edge statistic for r>=2 has GF G(z)-172z. Counts weight each line by H of its strong target, including parallel edges separately.

The recurrence g_r=874g_(r-1)+43344g_(r-2)+50400g_(r-3) holds for r>=4 (also for the h0=1 extension at r=3). The denominator has three distinct real reciprocal poles; its dominant characteristic root is approximately921.115395614955. This gives an algebraic exponential growth rate for this edge-weighted subfamily, not the global tournament enumeration growth rate.

## 5. Mirror becomes a genuine Lorentzian adjoint

Let tau be the mirror involution on the18 tile letters. Its four fixed letters are983,689,757,955 in the documented lexicographic bit encoding. Let J=diag(1,-1,1). The finite explicit matrices satisfy, for every letter,

K_tau(t)=J K_t^T J.

All18 identities are checked exactly, and K itself is J-self-adjoint. Thus the mirror operation on this actual combinatorial observer is an indefinite adjoint, not merely an analogy based on a repeated dimension.

For an even mirror-fixed word with first-half port vector v=(H,S,R), the target Hamiltonian count is

v^T J v=H^2-S^2+R^2.

For an odd mirror-fixed word with middle fixed tile u it is v^T J K_u v. In the even case H>=S>=0 and R>=0 because these are actual path counts, so this indefinite expression is nonnegative on the physically realizable vectors. This does NOT make J positive definite or prove a general self-adjoint spectral realization. Nor does it identify this three-dimensional observer with II_(25,1), a Leech cusp, or a Dehn invariant.

## 6. Exact blue and black H-weighted edge GFs

Let moment vector be (H^2,S^2,R^2,HS,HR,SR). Its sum over all prefixes evolves by an explicit6-by-6 integer matrix L obtained by expanding sum_t (K_t v)(K_t v)^T. The matrix is recorded in lorentzian-results.json and computed independently by symbolic polynomial coefficients. The initial vector is (1,0,0,0,0,0).

The even blue cost is (1,-1,1,0,0,0), and the odd blue cost, summing the four possible middle letters, is (42,2,160,-14,178,-32). Consequently the two GFs are exact entries of (I-wL)^(-1), with common denominator

D(w)=1-46154w-111330238w^2+17053296294w^3
 +1232670401484w^4-21361134804664w^5+55258216851600w^6.

Let

A(w)=1-34632w-94666590w^2-6314078182w^3
 +34746696144w^4+825385790712w^5;

B(w)=42+352060w+872465124w^2+58254110460w^3
 -1135327289608w^4+2091226279632w^5.

Then the blue target-H series over all nonempty mirror-fixed source words is

G_blue(z)=A(z^2)/D(z^2)-1+z B(z^2)/D(z^2).

Its coefficients for r1..5 are42,11522,2290528,548450036,111265364432. For actual family-to-strong outward blue edges, subtract42z. The black target-H outward series is G(z)-G_blue(z)-130z. At r2 the black contribution is183950, and at r3 it is176057568.

These degree6 formulas are certified by the symbolic resolvent of the displayed finite-state construction. A separate higher-dimensional exploratory recurrence recovery exists in blue_transfer.py, but is not the final proof. The clean proof/certificate is lorentzian_transfer.py. All72 blue length3 targets were also independently counted by full subset DP on15 vertices and sum to2290528.

## 7. Validation, lineage, and honest limits

New scripts: seamed.py, transfer.py, check_compression.py, blue_transfer.py, lorentzian_transfer.py. The last script's moment resolvent reproduces the exploratory blue rational functions exactly. The transfer coefficients have two independent implementations: permutation/cut enumeration and subset DP on all90 state-letter rows. All324 length2 targets and six selected general length3 targets match transfer counts. The complete72 blue length3 set matches independently. Normal and python -O runs are available for the base seam transfer and clean Lorentzian certificates; explicit exceptions rather than removable assert enforce validation.

Parent independently proposed the four-state reduction, verified the aggregate3-state matrix, and identified its indefinite-adjoint pattern. This lane checked and strengthened that proposal to every individual letter and computed the colour-resolved GFs.

Preserved structure: path-cover concatenation, source mirror, and complement-target H. Forgotten structure: target isomorphism class, target automorphism group, and how different source lines share a target node. Thus these results close a genuine edge-enumeration and edge-observable obligation for one infinite family, while leaving complete metagraph wiring open.
