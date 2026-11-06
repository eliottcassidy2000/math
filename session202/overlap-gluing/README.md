# Hamilton paths across an overlap

**PROVED, independently audited in the stated scope.**

For two induced subtournaments covering a tournament, Hamilton paths with the
same order on the shared vertices admit an injective merge. Their compatible
order-sector products give a lower bound for the output count. With zero or
one shared vertex, this is the product of the two full Hamilton-path counts,
regardless of cross-edge orientations.

For an arc u->v, exactly one forward-order Hamilton path forces transitivity;
exactly two force a single directed triangle containing the arc, with only
neutral transitive prefix/suffix padding. A literal six-vertex graph shows
that majority Hamilton-path order need not follow the arc.

These results rule out generic additive constructions from independently
large zero- or one-overlap inputs and identify a necessary order-profile
sidecar for larger overlaps. They are not an exact additive compiler and do
not prove any eventual odd-spectrum threshold.

Reproduce using the Python standard library, from this directory:

    python3 verify.py > /tmp/overlap-author-normal.txt
    python3 -O verify.py > /tmp/overlap-author-optimized.txt
    cmp /tmp/overlap-author-normal.txt verification.txt
    cmp /tmp/overlap-author-normal.txt /tmp/overlap-author-optimized.txt
    python3 audit/check.py > /tmp/overlap-audit-normal.json
    python3 -O audit/check.py > /tmp/overlap-audit-optimized.json
    cmp /tmp/overlap-audit-normal.json audit/NORMAL.json
    cmp /tmp/overlap-audit-normal.json /tmp/overlap-audit-optimized.json
    sha256sum -c RELEASE_FILES.sha256

The independent census includes all 1,100 labeled tournaments through order
five, including the empty tournament, all 254,254 ordered covering pairs, and
877,822 compatible-path merges. It includes empty private parts and containment
boundaries omitted from the author's narrower finite census. It also checks
both structural rigidity equivalences and the exact rigid family through ten
private vertices. These controls support the separate all-order proofs.

The frozen audited proof is retained under originals/. PROMOTION.json records
the status-only promotion. The independent audit is audit/REPORT.md. Frozen
ordered-amalgamation texts under dependencies/ establish inheritance and are
not needed at runtime; the present argument is self-contained. The speculative
clone-sector probe is deliberately excluded from this reviewed release.
