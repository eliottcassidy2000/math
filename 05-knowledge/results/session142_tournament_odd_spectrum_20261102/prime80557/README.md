# Existing-constructor recovery: two actual prime-count tournaments

**VERIFIED literal witnesses / inherited finite coverage / universal spectrum OPEN.**

This packet corrects an obsolete choice of numerical frontier. The reported
Royle prefix through80555 is not the current repository construction frontier.
The pinned `square_sum_spectrum_overnight_20261031d/README.md`, section6,
reports finite unrestricted coverage of all positive odds through99,999,999
except7,21. Its `VERIFICATION.json` records two independently implemented
scans with matching packed support SHA256
`34530aa67d56eb84e89943e9b525f37e2baaecb2b1ccce39ad094f8f4ef6fa7d`.
Their seven-parent bank uses D and D ordinal-sum T_k ordinal-sum D,
k=0..5, with all16144 cuts, source address x<=56, and arbitrary singleton
suffixes. The raw bank misses7,21,43,93,173,247; four inherited graph witnesses
repair the latter four. **This session did not rerun those billion-pair scans
or independently re-prove their whole finite prefix.** It recovered their
small executable constructor, source tables and source pins.

No novelty or minimum-order claim is made. In particular80557 was already
covered by Session124's million-prefix constructor.100000007 is an isolated
positive control just beyond the inherited finite prefix, not a new interval
or a proof of universal completeness. No random search was started.

## Exact recovery results

| Raw directed Hamiltonian-path count | Order | Independent graph checks |
|---:|---:|---|
|80557|19|Full endpoint subset DP, compiled O0/O2; two-lane block count|
|100000007|29|Two-lane block count; direct transitive-parent/apex formula|

Both graphs are valid tournaments, strongly connected and have prime H.
`witness_*.json` freezes outgoing bitset rows, complete block words and labelled
partitions; `adjacency_*.txt` freezes full binary matrices. Rows and columns
are0-indexed, bit j in row i means i->j. A word consists of right ordinal
blocks S (one vertex) or D (rows8,1,3,6), followed by one apex. Cut-mask bit1
means the apex beats that local vertex. Every edge is therefore explicit.

The inherited point query gives:

*80557: D cut8, x=0, y=10935, actual source(F0,q,g,h)=(11,7,1,5),
9 distinct fixed-source queries and422 inverse nodes.
*100000007: singleton source, x=0, y=165010214, one source query and1745 nodes.

The executable inverse has no arbitrary suffix-length cutoff. Its scope here
is x<=56 in exactly the seven-parent bank. A failed query is bank-relative,
not evidence that an integer is globally unattainable. The frozen literal
witnesses no longer depend on repeating the inverse.

## Why the independent graph count is exact

Delete the apex. The parent is an ordinal sum of the specified blocks. Each
Hamiltonian path splits uniquely at the apex into a possibly empty before
path and a possibly empty after path. Both visit blocks in increasing order.
Within each block, enumerate its partition into before/after vertices and
all directed local paths on each side (block size at most4).

For the before side retain whether it is empty, or whether its last vertex
beats the apex. A nonempty new block portion replaces this endpoint. For the
after side retain whether it is empty, or whether its first vertex is beaten
by the apex. Only its first nonempty portion fixes this endpoint. Cross-block
steps are automatic by ordinal orientation. Multiplying and summing these
local choices counts each complete vertex word once. Accept precisely empty
or correctly connected sides. `two_lane_count` first checks the literal
adjacency is a tournament and really has the claimed ordinal partition.
It imports no inverse or source-profile formulas.

For the order29 graph the parent is transitive in its displayed order0..m-1.
Let b_i=1 when the apex beats i. Directly from the graph,

    H = b_0 + (1-b_(m-1)) + sum_i (1-b_i)b_(i+1)
        + sum_(j<i) (1-b_i)b_j 2^(i-j-1).

The first two terms are paths with an empty side. Otherwise let i be the last
before vertex (the largest there) and j the first after vertex (the smallest
there). If i<j, coverage forces j=i+1 and gives one prefix/suffix split. If
j<i, every vertex strictly between j and i may be assigned freely to either
side; all outside vertices are forced. The apex arcs require b_i=0,b_j=1.
This proves the formula for every nonempty transitive parent, not just this
example. The checker verifies all parent edges before applying it.

The two checks share the elementary apex path partition but use distinct
summation schemes. The order19 full subset DP is a further independent
recurrence on arbitrary adjacency. **No full2^29 subset enumeration is
claimed.** The checker additionally compares the block recurrence with
literal permutations on all87 tested short words and with the direct apex
formula on all511 singleton words of lengths0..8, for598 controls in total.
The empty-parent case is separately1.

Primality is checked by trial division through floor(sqrt(H)); strongness by
reachability from every vertex. These two examples are also asymmetric
without a classification: Aut(T) acts freely on Hamiltonian vertex words,
so its order divides prime H. It also divides n!, whereas H>n in both cases.
Thus the automorphism group is trivial. This does not claim smallest order.

## The exact remaining feasibility obligation

For an actual source(F0,q,g,h), h>=1, q>=0,0<=g<=h, suffix output is

    N=F0+q*y+g*(y mod2)+h*(f(y)-1),
    f(0)=1, f(2a)=f(a)+a+(a mod2), f(2a+1)=f(a).

The inherited inverse sets t=N-F0+h and recursively uses

    even y: (t,q,g,h) -> (t,2q+h,h,h),
    odd  y: (t,q,g,h) -> (t-q-g,2q,0,h).

Terminal y=0 means t=h. For q>0 every candidate obeys
max(0,ceil((t-h-g)/(q+h)))<=y<=floor((t-h)/q). The parity branches require
t-h congruent to0 or q+g modulo gcd(2q,2h). The actual q=0 source has
(F0,g)=(h,0); redundant odd digits can be removed and an even representative
selected. This is exact membership for a fixed actual source, not proof
that every odd target belongs to some source image. The inherited proof
supplies a square-root arithmetic-operation bound, not polynomial bit
complexity. The80498 typed source-query controls were rerun normally and
with Python-O and match their inherited receipt exactly.

Universal acceptance of the seven-parent cap56 bank remains OPEN in its
source. Successful decoding, congruence admissibility, and an executable
finite inverse do not discharge source-image coverage. Nor is the singleton
subfamily universal: its inherited missing value13 is already documented in
Session136 local_closure/PROOF.md. No new exclusion theorem is asserted here.

## Portable reproduction and selection

Requirements: Python3 standard library and a C++17 compiler. No network,
repository checkout, installed packages, or large bitmaps are needed.

    python3 /absolute/path/to/prime80557/run_checks.py --work-dir /tmp/prime-replay
    python3 /absolute/path/to/prime80557/construct.py 100000007 \
      --adjacency /tmp/witness.txt --recipe /tmp/recipe.json

`run_checks.py` replays normal/-O receipts, validates pinned source bytes,
compiles the order19 checker at O0/O2, exercises malformed/hostile cases and
reconstructs both graphs. It writes only to its requested work directory.
`WHITELIST.txt` is the portable source/receipt selector; executables and caches
are excluded. `SHA256SUMS` covers that selector and all selected files except
itself. `sources.json` pins exact repository paths, commit9d5bae74f05a685d523df0c5559d71c3cd53cd6a,
Git blob IDs, SHA256, byte lengths and local snapshot selectors.
`recover_sources.py` is an optional provenance materializer tied to the
research checkout; it is not required for portable verification.

Normal/-O hostile checks reject self-loops, one-sided mutations, incomplete
partitions, nonordinal cross edges, and a wrong claimed count. The independent
Session142 audit separately reconstructs both graphs, verifies native
singleton transitions, reverses an apex edge as a hostile, and replays the
order19 full DP. Its report is stored in the sibling independent packet.
