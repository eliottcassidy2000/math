# Independent audit: sparse feedback arcs and ribbon resets

Status: PASS WITH MINOR PRESENTATION CORRECTIONS. All independent exact replays completed successfully; see RESULTS.jsonl.
Scope: REPORT.md, RIBBON_RESET.md, SPECTRAL_RESETS.md from session231-packing-carries. This audit writes only its own directory and does not commit or edit the producer.

## 1. Arbitrary-cutwidth kernel: PASS

For an order with r backward arcs, fix y->x of span s. Each other reversed pair has an endpoint outside {x,y}; delete one such endpoint. This uses at most r-1 deletions, preserves the target edge, and eliminates all other reversed edges. The restricted interval has at least (s-1)-(r-1)=s-r interior vertices. When this is positive, its Hamilton count is 1+2^t, t the number of surviving interior vertices. This is because each directed cycle uses the reversed edge, and odd cycles correspond precisely to nonempty odd subsets of interior vertices, numbering 2^(t-1); disjoint nonempty cycle packings cannot coexist. Monotonicity then gives 2^(s-r)<=X-1 and s<=r+floor(log2(X-1)). When s<=r the desired bound is automatic.

Vertices outside the closed reversed-edge intervals are beaten by every earlier vertex and beat every later vertex, even if the rest of the tournament has noninterval strong components. They are singleton ordinal factors. Removing all of them preserves H. The remaining order has at most r(L+r+1) vertices. Padding with sinks to N=k(L+k+1) preserves H and the reversed-edge budget. Counting ALL subsets of at most k pairs is an upper bound on distinct counts, with no uniqueness or reachability assumption. The r=0 case is separately represented by the transitive padded order. Thus the binomial bound, fixed-k exponent 2k, and moving-k density-zero consequence hold.

For moving k, log(eM/k)<=log k+2loglog X+O(1)<=3loglog X+O(1), uniformly for the specified k. Floors and finitely small X do not affect the asymptotic. The base-b bound follows from H_b=1+b*2^(t-1), giving s<=r+1+floor(log2((X-1)/b)); the restriction X>=b+1 makes the small-span case automatic.

Credit: 01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md and its source 05-knowledge/results/spine_grammars_cutwidth_apexes_hankel_seams_20261104.md, section 6, already establish the qualitative bounded-budget obstruction. The audited novelty is the explicit arbitrary-cutwidth vertex kernel and description bound, not the qualitative obstruction. The insertion/packing identities and one-back-arc formulas are inherited. check_kernel.py independently enumerates all labelled orders through order 5 to test the deletion, parity, and span bound; this is a hostile control, not the all-order proof.

## 2. Spectator obstruction: PASS, with boundary notation clarification

H(A=>B)=H(A)H(B) since a Hamilton path cannot return from B to A. A nonzero residual for F commuting with C3 ordinal prepending therefore grows by 3^t. For a finite menu and an even integer radix, each initial residual is an odd integer, hence nonzero, and the minimum residual magnitude over the entire menu grows at least as 3^t. This defeats every bounded digit set even if the chosen block depends on the full prefix.

A C3 is a span-one ribbon: reverse both adjacent edges of an ordered triple. For wider words, ordinal concatenation must clear the out-of-range bits of the first w letters of the right operand, rather than blindly concatenate arbitrary raw word encodings. Those bits had no effect in the original standalone prefix. This normalization gives the same tournament and prevents cross-boundary back edges. With |u|>=w, appending a suffix can reach only the final w vertices of u, never the spectator triples. Thus the displayed ordinal equality and all-width conclusion are valid. Add this one sentence to avoid an encoding ambiguity. A finite menu automatically has bounded maximum length; the obstruction is all lengths, not a finite search result. The domain is implicitly nonempty.

No claim follows for strong-only invariant domains, since ordinal prepending leaves that domain, or for operations reaching arbitrarily far into the prefix.

## 3. Spectral certificate design and exact scope

check.py does not import any producer scan or elimination function. It reads each matrix's explicit mask, validates exactly 2^w distinct masks, and builds complete undeduplicated product levels. Expected totals are 5,460 (span2, all letters, lengths1..6), 584 (span3, all letters, lengths1..3), and 29,523 (span2, masks1,2,3, lengths1..9). Every product is tested at eigenvalues 2,4,8 modulo the verified prime 1,000,003, with a separately written batched pivot-elimination routine. A known tripling block is a positive singular control.

Overflow is certified BEFORE every level multiplication by n*max_abs_letter*max_abs_previous_product<2^63-8, bounding even intermediate dot-product partial sums. Each result records the bound and actual maximum entry. Elimination uses residues <p; products are <p^2 and differences are within signed 64-bit range. This is exact finite arithmetic, not floating-point inference.

The quotient is verified in its actual stored basis, without relying on reproducibility of LLL choices. Start with pairs (h,phi), close under corresponding right multiplication by stored and full-profile letters, select an exact independent set of observable rows C with associated rows D, and set B=C^(-1)D. Then independently verify B*init=y0, h*B=phi, and A_a*B=B*M_a for EVERY letter. These exact rational identities imply correctness for arbitrary prefixes and suffixes by induction. They do not merely compare bounded sample trajectories. The full profile transfer is inherited from THM-4880-ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-has-a-rank-14-integer-representation-irreducible-on-nonempty-words-47-needs-span-3.md. No irreducibility claim is needed: the full 14-coordinate span2 representation includes the empty-prefix coordinate; its nonempty reachable space is 13-dimensional, as corrected in THM-4880 and THM-4881.

Cayley-Hamilton applies to every scalar orbit h M^t y. If it equals c a^t-d/(a-1), c!=0, then p(E) gives c p(a)a^t-(d/(a-1))p(1)=0. Comparing consecutive t gives p(a)=0, contradiction. Therefore each finitely tested block has an ALL-TIME obstruction to nonconstant exact affine growth at the tested radix. The block universe remains finite. This does not exclude switching blocks, longer blocks, scalar fixed points, variable digits, or general strong-state constructions.

## 4. Required presentation repairs

1. Span3 file order is masks [0,1,2,4,3,5,6,7], not [0,1,2,3,4,5,6,7]. The producer's exhaustive scan is unaffected, but index-based words must be decoded with the file mask. Span2 indexing is as claimed.
2. search_resets.py's `found` filters modular singularity through an exact determinant check. Its empty list alone certifies no exact eigenvalue, not no modular candidates. The independent checker explicitly counts modular singulars, and finds zero in every listed universe, validating the stronger certificate claim.
3. Clarify ordinal word boundary normalization as above.
4. Include search_long_resets.py in the reproduction command list, which currently omits it.

## Reproduction

Set MATH_REPO to the repository root if needed. Run:

    PYTHONDONTWRITEBYTECODE=1 python3 check.py > RESULTS.jsonl
    python3 check_kernel.py > KERNEL.json

Input digests are recorded in RESULTS.jsonl. Exact hostile characteristic polynomial and stable ranks are checked there; parity_probe.py's modular discovery rank alone is not a rational-rank proof. check.py computes rational rank of M^14 and binary rank of M^14 directly. Stable ranks are attained by exponent n=14 for any n-dimensional matrix, so those checks are decisive.

## Final replay outcome

PASS: both quotient intertwiner systems, all 35,567 enumerated block occurrences (some span2 scans overlap), every tested determinant at 2,4,8 nonzero modulo 1,000,003, and the exact hostile polynomial x^11(x^3-971x^2-1302x+96). Its stable ranks are 3 over Q and 1 over F2. The eigenvalue-3 positive control is detected. Largest actual entry in the longest span2 universe is 8,015; its largest certified raw dot-product bound is 220,080, far below int64. No all-block extension follows from the finite scan.

The original attempted LLL reconstruction selected a different valid basis on this platform; it was not used as evidence against the stored representation. The final checker instead proves exact stored-basis intertwiners, which all pass. No producer or repository files were edited.

## Optimized-mode hardening and replay

All correctness gates in check.py and check_kernel.py are explicit require calls raising RuntimeError; none depends on Python assert. Both full scripts passed normally and with python3 -O, and their normal/optimized receipt files compare byte-for-byte equal. Input paths in receipts are repository-relative. check_active_gates.py additionally checks the syntax trees contain no assert nodes and deliberately passes False to each require function, verifying rejection in normal and optimized mode; its receipts also compare byte-for-byte equal. The unused intrep import was removed; no inherited source was changed.

Optimized reproduction:

    PYTHONDONTWRITEBYTECODE=1 python3 -O check.py > RESULTS.optimized.jsonl
    python3 -O check_kernel.py > KERNEL.optimized.json
    python3 -O check_active_gates.py > GATES.optimized.json
    cmp RESULTS.jsonl RESULTS.optimized.jsonl
    cmp KERNEL.json KERNEL.optimized.json
    cmp GATES.json GATES.optimized.json
