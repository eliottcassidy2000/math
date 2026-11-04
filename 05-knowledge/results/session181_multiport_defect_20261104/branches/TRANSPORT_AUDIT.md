# Independent audit of positive deleted-root defect transport

2026-11-04. **VERIFIED proof review, replay and independently implemented candidate count.** Audited sibling transport/PROOF.md, verify.py, common.py and candidate_transport.py. No mathematical errors found.

The deleted-root formulas correctly include empty left/right slots and the singleton-root case. The singleton-inneighbor simplification follows because both factors contain the same vertex and their squarefree product vanishes. The arbitrary-birth endpoint updates and defect extraction account correctly for whether the newborn beats the physical root; all old incidences and cross-port terms are retained.

The dominant-prefix formulas (3)/(4) are valid for t>=1, including empty deleted root. Their defect becomes b-S at the first birth, then stays constant, rather than necessarily staying at its original value. The proof and wording distinguish these stages. The author explicitly handles U empty separately; the identity A=S_I+S_O is used only for nonempty U, as required. Strongness of a nontrivial strong source is preserved. A source having b=0 cannot develop negative defect, since (1) already gives b>=S.

The multiport formulas (5) follow from the actual root deletion TT_s->C: all prefix vertices belong to the first slot in the bad-join product, while its second slot starts in W. Consequently b=J is correctly retained. The two-clock argument gives exactly F'=3F and (h'-J')=3(h-J)+d. No condition J'=0 or J'=3J follows. The stated algebraic targets are necessary/sufficient for all-clock scalar matching, not sufficient for the existence or native legality of a physical rewrite.

Replayed verify.py normally and under -O; stdout is identical in transport_normal.out and transport_optimized.out. Counts are285 deleted-root identities,4306 arbitrary-birth identities and855 dominant-prefix checks. These are replays of the author's path-enumeration-versus-DP implementation, not newly independent software.

For the actual17-vertex candidate, audit_transport.cpp independently uses forward extension of path endpoints, computes starting endpoints by reversing all arcs, and extracts both physical cross-joints. Output is transport_candidate_audit.out, generated from transport_candidate.in. It confirms:

- Full H=512425 and E_39(2)=516115.
- Suffix h=18583, E_4(2)=68248, E_21(2)=56135, sum F=124383.
- J(4,21)=1206 and J(21,4)=2484; both diagonal terms vanish.
- Therefore J=3690 and alpha=h-J=14893.

The input's physical names were separately checked against the actual rows: the root's in-neighbor set is precisely{4,21},50->51, and both prefix vertices beat every suffix vertex. The offset arithmetic8308-2076*2^s and next targets F'=373149, alpha'=44681/44683 are correct. The +4 offset occurs uniquely at s=2; +2 is impossible at positive integer clocks. These checks do not duplicate or replace the separate native-guard certificate.

Reproduction:

    g++ -O2 -std=c++17 audit_transport.cpp -o audit_transport
    ./audit_transport < transport_candidate.in
    python ../transport/verify.py
    python -O ../transport/verify.py

At order17 the generous bound18*17! is below2^64, so the independent unsigned64 computation is safe for every two-slot quantity used here. No larger-order arithmetic claim follows.
