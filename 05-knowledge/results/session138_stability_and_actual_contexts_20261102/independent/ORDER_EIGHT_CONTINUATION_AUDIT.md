# Order-eight continuation audit (excluded from checkpoint one)

**Verdict: ACCEPT.** This is a separate second-checkpoint continuation, not
part of the frozen first Session138 release through order seven.

The complete cycle-cover ledger in `../stability/SHARP_EIGHT_GRAPH_CERTIFICATE.md` was
independently checked against all cycle-length partitions of at most eight.
In particular the full covers4+4,2+6,3+5,2+2+4,2+3+3 and2+2+2+2 have the
correct signs and multiplicities in u,v and both signed margins.

Put P=V-M22+M222-2M2222. The potentially negative4-cycle terms group as

    C4-M24+M34-2M44
      =sum_C w(C)[1-V_(C^c)+tau_(C^c)-C4_(C^c)]
      =sum_C w(C)[u_(C^c)-M22_(C^c)] >= C4/8.

The factor2 on M44 is essential and correct: each unordered pair of disjoint
4-cycles is encountered from either chosen cycle. The lower bound uses
u4>=5/16 and M22 on four vertices<=3/16. For the other signed margin the
negative2M44 term is absent, so inserting it only weakens a lower bound.
Likewise C6-2M26>=C6/2, and M223<=15tau/16 because a complementary five-set
has15 two-edge matchings. All omitted cycle terms are nonnegative. Therefore

    Delta>=P+tau/16+M23+C4/8+C5+C6/2+C7+3M33.

The local triangle variance lemma gives M23>=P21(q). The final graph
polynomial is `F=P+P21-87V/128`, WITHOUT any wedge term. At endpoints,

    512F=41m-32m2+8m3-4m4+8p21.

The separate independent child audit in `n8_audit/FINAL_GRAPH_AUDIT.md`
verifies this stronger polynomial, not merely the earlier+2w variant.
Its seven-vertex orbit BFS uses adjacent transpositions and covers all2^21
labelled graphs in1044 disjoint orbits. Every eighth-vertex graph is represented
by one of128 neighborhoods on one of those cards, proving augmentation
coverage without a full2^28 scan. A vertex-subset matching DP and wedge-first
p21 count check all133632 augmentations. Normal/UBSan receipts match:
minimum0, maximum10808, exactly2 zero augmentations, explicitly empty and
perfect matching. The earlier+2w files are preserved as provenance only.

Multiaffine interpolation proves F>=0 on the whole real q-box. Consequently

    Delta>=87V/128+tau/16+C4/8+C5+C6/2+C7+3M33.

Four fair ordinal pairs have H1,u81/256,|v|1/256,V1, proving coefficient87/128
sharp. Positive-V unreserved equality forces tau0 via the retained positive
reserve; H1 classification then forces precisely four fair ordinal pairs.
No order-nine or all-order optimality claim follows.

Independent graph JSON SHA256:
`043a19e9449fdfefc0acfa0c45f1628598ed37905cbe6beacaccc5ecb2a6eaf5`.
The pinned manuscript/code snapshots are in `ORDER_EIGHT_REVIEWED_SOURCES.json`.
