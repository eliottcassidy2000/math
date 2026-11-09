# Final independent review: Session234 digraph packet

**Verdict: PASS. No required mathematical edits identified.**

Reviewed README.md, PROOFS.md, GENERAL_GRAMMAR.md, census/REPORT.md, census/generator.py, audit/AUDIT.md, and supporting theorem citation. The prior independent full scanner audit supplies census completeness; this review additionally checked mathematical deductions, generator correctness and source references.

## Proof checks

- Diagonal independence and the minimal orders four/seven and five/twenty-one are correct. The missing-arc bound on four vertices excludes six permutations, leaving at most eighteen.
- SCC ordering, multiplicativity and uniqueness are valid for semicomplete digraphs. Automorphisms fix SCC positions, including when repeated component types occur.
- The longest-cycle argument for strong semicomplete Hamiltonicity is sound. The order-two exception is handled; spanning strong tournaments exist at orders at least three.
- Tournament insertion gives an injection in the stated direction. Moon's induced strong subtournament theorem and the independently enumerated minima give the all-order cutoffs four and six. This is not an extrapolation from census frequencies.
- Both strong H3 cores are explicitly included. The doubled triangle's extra arc is unused, so strict increase under arc addition is correctly rejected.
- The 3 H7 and 32 H21 core alphabets, mixed ordered products, singleton gaps, automorphism weights, ordinary/ exponential generating functions and loop-orbit coefficients are mutually consistent.
- The weighted tournament-resolution identity, one-digon H21 example, sharp necessary alternatives, and exact arc-increment identity have the correct directions and conventions.
- DG8 is a bijective generator for labelled off-diagonal digraphs once P is a set and the optional arc ground set is the complement of U. Hypergraph avoidance excludes exactly the additional paths. This does not provide a semicomplete-style strong-core classification for unrestricted digraphs, and the packet does not claim one.
- The general finite grammar is valid. The total core-order bound follows from sum d_i <= product d_i, and eventual degree K(h) is asserted only for nonempty fibres. Labelled counts are correctly normalized by n!, loop-decorated counts by 2^n, and h=0,1 are separately treated.

## Executed checks

`check.py`, `check-output.txt`, and `validation.json` accompany this report.

1. Imported the standalone generator without invoking its file-writing CLI; reran its full loopless validation through order eight.
2. Independently enumerated vertex permutations of all 37 embedded strong cores; each has the required H, is strong, and is distinct up to isomorphism within its target alphabet.
3. Exhausted every loop-decorated generator output through order seven, for both targets. For every output, checked H by subset DP, canonicalized actual SCCs, checked unique ordered canonical signatures, computed decorated automorphism groups, and summed n!/|Aut|. All counts agree with the formulas.

| H | n | Loop-decorated isomorphism classes | Labelled matrices |
|---|---|---|---|
|7|4|48|1,152|
|7|5|192|23,040|
|7|6|576|414,720|
|7|7|1,536|7,741,440|
|21|5|480|56,320|
|21|6|2,944|2,088,960|
|21|7|11,008|54,190,080|

These are exhaustive outputs of the classified generator, not fresh exhaustive scans of every input digraph at these orders. Completeness remains the proved structural reduction plus the independently audited small-core census.

## Links and source checks

All relative file links in the three main documents resolve, including the specifically named THM-1370 tournament file. The two README section anchors match their target headings.

Moon's original paper was opened directly via Cambridge:
https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf

The first page confirms title, author, Canadian Mathematical Bulletin volume 9 number 3 (1966), page 297, and Theorem 2 s(n,k)=n-k+1. Its definition of s gives the needed induced strong (n−1)-vertex subtournament statement. The DOI open endpoint itself was unavailable in this tool; the authoritative PDF was accessible. No citation correction is indicated.

https://arxiv.org/abs/2308.10214v2 confirms Chan–Pak's title, authors and journal metadata TCS 1015 (2024), 114776. https://arxiv.org/html/2308.10214v2 confirms Section 6.5, Conjecture 6.2 (all odd totals except 7 and 21), and its attribution to bof/Gordon Royle, March 2016.

No universal spectrum or novelty claim is established by this work, and the packet correctly says so.

## Exact reciprocal-count corollary

Also checked after the initial review: ordered SCC sums add no reciprocal
pairs between components, and singleton components contribute none.
Thus a semicomplete H7 graph has exactly two reciprocal pairs. A
semicomplete H21 graph has one, two or three: its single core has 1–3,
or its mixed H7/H3 word has 2+0 or 2+1. Every stated value occurs.
Consequently arbitrary binary matrices satisfy H7 => an absent pair or
b=2, and H21 => an absent pair or 1<=b<=3. These are necessary, not
sufficient, and arbitrary loops do not change them. No new census is
required for this direct consequence of the verified catalogue.
