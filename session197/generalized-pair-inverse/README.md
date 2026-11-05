# Generalized pair families: reviewed inverse packet

PROOF.md gives four actual C2 graph constructions from arbitrary rooted start-one sources, their exact formulas, and a finite all-order membership algorithm for their union Gamma. It proves source cutoff D<=floor(h/7), preserves the weak-source deletion sidecar and gives finite suffix bounds. It does not claim Gamma is the whole C2 class or that its scalar image is cofinite.

The exact bounded consequence is: Gamma together with the complete C1 scalar image and singleton1 covers every odd value through1861 except7,21,67. SUPPORT.json stores one physical recipe for each of Gamma's883 output values in that range. Every one is independently reconstructed and counted. VALIDATION.json includes explicit repairs533,821,1861.

Source replay, in a copy because receipts are written beside the verifier:

    python verify.py
    python -O verify.py

Both runs must reproduce the included receipt. The independently written C++ checker needs no source imports:

    g++ -std=c++17 -O3 audit/check.cpp -o /tmp/check-generalized-pair
    /tmp/check-generalized-pair /tmp/independent-banks.txt > /tmp/independent.json
    python audit/compare.py . /tmp/independent-banks.txt /tmp/independent.json

A debug check uses -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined instead. The complete independent bank text is regenerated, not required as an input. The audit report explains its separate literal graph counting and source-set comparison.

RELEASE_FILES.sha256 is the unchanged reviewed source manifest. PUBLICATION_FILES.sha256 additionally binds this wrapper and the independent review/checkers/receipts. Mathematical dependencies are copied with their reviewed hashes; no navigation guide or earlier result is modified by this additive packet.

The separately reviewed graph-count-addendum/ proves a necessary scope limitation: at order14 there are at most41,736 generalized-family presentations versus44,394 C2 graph classes, and the graph-image proportion tends to zero. It does not exclude scalar coverage, since different graphs can have the same H. Its independent coefficient checker is graph-count-addendum/audit/check.py. This addendum leaves the audited inverse snapshot and its original manifest unchanged.
