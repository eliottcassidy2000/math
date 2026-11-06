# Equal-profile physical sources

**PROVED and independently audited, with limited scope.** Two strong six-vertex rooted sources share the full retained edge-extension profile while their start counts are 6 and 7. Every common O/Z padding word and every common ordered rooted right amalgamation preserve equality of the profile. This is a physical common-base construction, not independent prescription of the base.

Under a right input W the endpoint choices become 6 H(W-root) and 7 H(W-root). The first is always 0 modulo 3. Consequently this pair and these common right lifts alone do not supply both endpoint residues 1 and 2 needed by F(8,137). The hypothesis remains open, and this packet makes no eventual-spectrum or threshold claim.

The packet also supplies two explicit, separately verified profiles for normalized base 999, with four-core child counts 2999 and 3001. These finite witnesses do not alter the previously published conditional compiler or its frozen missing-certificate controls.

- ROW_NORMALIZATION.md contains the theorem, boundary cases, explicit graphs and inherited formulas.
- verify_rows.cpp and NORMAL.json/DEBUG.json contain the author checks: 31 padding words and all 75 labeled right tournaments through order 4 with root fixed at vertex 0.
- BASE999_WITNESSES.json gives every adjacency row for the two extra finite witnesses.
- audit/ contains the independent proof review and standalone controls, including all 285 marked right inputs through order 4, 568 generic padding checks, and both new witnesses.
- dependencies/ freezes the inherited rooted-amalgamation and four-core readout results and reviews.
- reviewed-originals/ preserves the exact proof that was independently reviewed before status-only promotion.

Reproduce the author controls from this directory:

    g++ -std=c++17 -O3 verify_rows.cpp -o /tmp/verify_rows
    /tmp/verify_rows

Independent reproduction commands are in audit/REVIEW.md. The arbitrary-module source-pair extension is separate and is not included in this release.
