# Independent audit of Session203 rooted-cover closure

Verdict: PASS within the stated source grammar and inherited graph-classification scope. No theorem/status promotion performed. No repository edits, commits, fetches, or publication performed; the shared dirty checkout was preserved as requested. This is a scratch audit of the frozen candidate, not a fresh audit of the repository frontier.

## Inputs and reproducibility

The complete candidate packet was copied into `frozen/` before testing. Every entry of its SHA256SUMS verified. Replay outputs reproduce NORMAL.json, OPTIMIZED.json, COMPILER_NORMAL.json, COMPILER_OPTIMIZED.json, TARGET137.json, TARGET139.json, and INVENTORY997.json byte-for-byte. The commands are the candidate's documented Python invocations, with output redirected to the corresponding REPLAY files here. `python3 independent.py` and `python3 -O independent.py` likewise agree byte-for-byte.

The additional graph-classification document examined is frozen as ROOTED_DEPENDENCY_REVIEWED.md, copied from session197-lowcycle/ROOTED.md. AUDIT_SHA256SUMS records all audit inputs, independent code, and receipts, excluding interpreter caches.

## Transfer and labels

I independently checked all nine rows of each O/Z/J map. A designated first path is nonempty; ordinary paths are labeled and may be empty. These distinctions are essential:

- The ordinary C_Z convolution has two choices for the bin containing z. Its correction is 2T-2L. Splitting a root-start two-cover on the selected subset supplies a designated path plus two distinct ordinary labeled bins, exactly T rather than an unordered three-cover count.
- For T_Z, the new z cannot use z->v because v already starts its designated path. Appending z to one of three old paths therefore enumerates every case. Precisely the b configurations with singleton designated path v fail, giving 3T-b.
- For T_J, singleton designated path z leaves C covers; a nonsingleton designated path must begin z,v, leaving T covers. Their sum is T+C.
- L_J=C-L includes the initially empty first bin of C, which becomes singleton z. R_J=h+R likewise includes singleton z separately. Deletion gives b_J=C, whereas O and Z give 2b.
- Nontrivial Z preserves s because any root-start path ends at a nonroot vertex, which beats z. At the singleton it fails: s_Z=0. No other singleton formula needs changing.

The four nonnegative cover coordinates and every displayed transfer were independently algebraically substituted. C>=2L uses disjoint endpoint-root classes in the two labeled bins; T>=b uses singleton designated root; b>=L uses deletion of the terminal root, retaining an empty first bin if needed. No realization claim for arbitrary cone points is justified or made.

## Positive growth and finite inverse

The monotonicity proof is valid for nontrivial actual seeds with positive s. R>d follows from a root out-neighbor and at least one path on the remaining vertices. Thus J has strict B increase, while Z increases by at least two. The edge-extension identity B=H(E(U)) is a bijection, so tournament path parity upgrades J's strict increase to at least two. O preserves B and s.

Deleting a terminal O suffix preserves the target. If the remaining word is nonempty, its final Z or J is applied to an n-vertex source. Every labeled bipartition contributes at least one cover, including the empty bins, giving b>=2^(n-1), C>=2^n. The Z bound follows from B_Z=h-s+L+R+b-e>=b. The J bound follows from L_J=C-L>=C/2 and h_J>=s_J. Therefore the resulting order n+1 is at most floor(log2 M)+2. Empty words require the explicit seed-order alternative. The weaker bounded-O-block argument also checks; it uses e=0 only after at least one O.

The BFS compiler uses the rooted-edge seed of order two, so bit_length(M)+1 is exactly the proved order cap. Its all-nine-coordinate merge preserves every continuation, and BFS retains minimum recipe length, hence never loses a continuation fitting the cap. Monotone pruning is sound. State exhaustion explicitly has complete=false and forbids absence claims. The target-mode CLI emits INCOMPLETE when no discovered witness is available. Inventory mode uses complete=false/reason instead of a status string; this is unambiguous but is a minor literal difference from prose saying every resource-limited return is labeled INCOMPLETE. Nonpositive max-state limits also return an incomplete result after retaining the unavoidable seed; they do not produce a false absence. Neither issue affects the mathematics or controls.

## Inherited source scope

Response completeness alone would not justify completeness for the nine-coordinate refinement. The frozen dependency's proof, however, explicitly constructs each entire physical gadget M_L as J followed by its O/Z tail, and says that every nontrivial source comes from the rooted edge. I also read the underlying ROOTED.md Theorem 2: the SCC decomposition of T-root identifies an initial strong/singleton source with unique root-start path and a transitive singleton tail, giving the actual one-seam graph, not merely matching responses. For weak sources, the root is in the initial SCC; positive start one forces all later SCCs to be singletons. These are terminal O extensions. This supplies the required graph-level completeness for the new state.

Recommendation for release packaging: include ROOTED.md or its relevant frozen dependency chain beside PHYSICAL_SOURCE_GRAMMAR.md so this stronger graph-level premise is immediately auditable. This is a provenance improvement, not a mathematical blocker. The audit accepts the inherited classification as a proved dependency after inspecting its relevant proof; it does not re-audit unrelated claims elsewhere in that document.

## Independent controls

`independent.py` counts paths by enumerating permutations and nondecreasing cuts, not by subset dynamic programming. Each permutation and cut vector is one ordered cover, with repeated cuts representing empty bins. It checks every input tournament rooted at zero through order four, plus 24 deterministic order-five and 6 deterministic order-six inputs, and all three literal children: 315 child checks in total. Physical graph addition/relabeling is separately implemented using edge pairs.

An unmerged, unpruned enumeration of every word of lengths 0 through 9 (29,524 words total) yields exactly the compiler's 330 bases through 997. The recurrence is shared for this enumeration, so it is an independent check of search completeness/control flow, not a second proof of transfer laws. The independent permutation counter supplies the latter finite evidence. The supplied canonical strong-grammar comparison through order nine also replays successfully.

B=137 is absent from the complete unit-start support, only in that scope. The B=139 literal source is independently counted and has s=1. Resource checks at limits -1, 0, 1, 2, and 10 are all incomplete, with no absence conclusion. The supplied all-labeled order-five transfer checks, all recipe words through length five, and optimized-interpreter controls all pass.

No result establishes arbitrary-start absence, F(8,137), infinite-tail occupancy, or polynomial-time membership. The all-order exclusion for endpoint one at B=137 follows from the proven cap plus the inherited graph classification, not from a guessed census cutoff.

## Packaging addendum

Before finalizing this audit, the parent added dependencies/ROOTED.md and dependencies/PROVENANCE.json and refreshed the candidate SHA256SUMS. I verified every updated checksum, compared ROOTED.md byte-for-byte with ROOTED_DEPENDENCY_REVIEWED.md, and incorporated both additions into frozen/. Their provenance hashes agree with the actual files. The existing proof, executable code and computation receipts remained hash-identical. The release-packaging recommendation above is therefore satisfied; no outstanding mathematical or packaging blocker remains.
