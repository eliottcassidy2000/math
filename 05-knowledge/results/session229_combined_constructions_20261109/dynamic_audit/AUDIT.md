# Independent audit: context-stable core switches

Verdict: PASS for the stated construction identities, finite exact interfaces,
and exhaustive small-order maximality. No coverage theorem or unbounded odd
interval is certified. Audit performed independently without producer imports
or Git writes. The frozen producer proof is INPUT_PROOF.md.

## Proof review

The Session216 three-core convention is (h,q,A,G,T,S,E), with optional first
path in A, optional second path in G, both nonempty in T, and H(empty)=1.
The checker reconstructs these statistics from literal directed path
permutations, separately from the producer code. Outgoing masks and root bits
have their usual vertex-index meaning; binary padding words read high bit first.

For arbitrary transitive padding, the inherited closed formula gives

    N=(qP+A)(2^(b+1)-rev_b(R))+h(2R+1-f(R))
      -PG-hf(P)-T-(P mod 2)S-E first_b(R).

Equal (h,q,A,G,S,E) therefore leaves exactly minus the T difference, for ALL
padding lengths. Finite tests are not the basis of that universal quantifier.
The stated fixed-h,q necessity argument is also valid on the full algebraic
context language: vary b at P=R=0, then P=2 and P=1 at b=0, then suffix first
bit 1. This is not a necessity theorem for restricted legal languages.

The five-vertex masks 3 and 5 give exactly (15,120,42,42,10/8,6,6).
The three six-vertex choices give exactly (45,432,129,129,20/22/24,15,15).
Every listed core is strong. A selected first prefix vertex, or an empty
prefix and nonempty cut in the strong core, guarantees that x reaches the
initial source component. Every body vertex reaches x through v; hence the
closure is strong. An unselected initial prefix vertex is a genuine failure,
not something oddness of the count repairs.

The ordinal composition formula retains every interface coordinate while a
pure T defect is multiplied by outside Hamilton-path counts. Thus t repeated
H15 switches give step 2*15^(t-1), not step 2. The corresponding H45 step is
2*45^(t-1). Fixed singleton insertions do not change these multipliers.
Literal two-copy H15 counts are 9465,9495,9495,9525. These are a concrete
hostile control against claiming an unbounded consecutive odd interval.

H19's six equal-A pairs, H9's compensated carry, and the H29 four-mask
fiber all have the displayed difference formulas. H19's 19-to-42 pair has
difference -12P-2e, hence gap 2 exactly at P=0,e=1; its shortest graphs
have counts 345 and 343. The other five differences cannot have magnitude 2.
The H9 difference is always at least 558, so equal A alone is insufficient.
The H29 triple at empty padding is 175,177,179, with its advertised restricted
suffix stability. Cyclic cores on 5,7,9 vertices independently reproduce the
stated two-mask T sets, including the gaps at orders 7 and 9.

## Independent computations and their universes

check.py uses no producer imports, no assertions, and no external data.
It computes the five main interfaces by literal directed path permutations;
checks every prefix/suffix word of total length at most 3 for all five states
(245 literal graphs), and compares whole-graph Hamilton subset DP against the
closed formula. It tests strongness independently by reachability, including
illegal leading-zero prefixes as hostile controls. It also checks H19/H9/H29,
a repeated-switch literal graph, and cyclic 5/7/9 two-mask fibers. Normal and
optimized Python outputs are byte-identical JSON receipts.

census.cpp is separately implemented and exhausts every labeled tournament
and every cut at each order 1 through 6. No isomorphism bank, strongness filter,
or pruning is used. It retains an accumulated map from (h,q,A,G,S,E) to all
observed T values across orders. Each output line records order, the number
of rooted labeled graphs at that order, cumulative distinct interfaces, and
maximum consecutive step-two T-run through that order. Expected output:

    1 2 2 1
    2 8 6 1
    3 64 18 1
    4 1024 64 1
    5 32768 342 2
    6 2097152 3202 3

This independently verifies minimality of a fixed-interface gap-two switch
at order 5 and absence of a four-value step-two fiber through order 6. It does
not prove anything about order 7 or a different interface/context language.

## Reproduction

    python3 check.py > /tmp/dynamic-audit.json
    python3 -O check.py > /tmp/dynamic-audit-O.json
    cmp /tmp/dynamic-audit.json CHECK.json
    cmp /tmp/dynamic-audit-O.json CHECK.optimized.json
    g++ -std=c++17 -O2 census.cpp -o /tmp/dynamic-audit-census
    /tmp/dynamic-audit-census > /tmp/dynamic-audit-census.txt
    cmp /tmp/dynamic-audit-census.txt CENSUS.txt
    sha256sum -c SHA256SUMS

## Scope and review notes

The inheritance is the audited Session216 path-cover/star compiler, not a new
proof of arbitrary target occupancy. The fixed-interface switch supplies
context-stable windows of two or three counts. It does not by itself show that
these windows cover an interval or all sufficiently large odd integers.

Release review caught a terminology error missed in the initial
audit: calling the H15/Q120 switch "unsaturated". Since gcd(15,120)=15,
the module is arithmetically saturated. Its output gap is 2 (unit difference
in the half-count coordinate); this does not imply an unsaturated module.
The corrected description is "context-stable pure-T difference". Numeric
identities and receipts are unaffected.

The final frozen proof also incorporates these editorial corrections: state the H19 base masks explicitly;
replace 'A Hamilton path reaches x from every L vertex' by the elementary
reachability statement. Neither changes the mathematics. The proof status line references this independent audit directly. The finite-search
interface counts are cumulative, as explicitly labeled above.
