# Reviewed all-order feedback separation at H=533

[The theorem](source/session197-lowcycle/FEEDBACK533.md) classifies exactly six strong tournaments with H=533 and at most two Hamilton cycles. All six have order11, C=2 and feedback vertex number3. Therefore permanently retaining two original ports cannot realize this count within the strong C<=2 class, whereas an explicit changing-two-port word does.

[The independent review](source/session197-bridges/feedback533-audit/REPORT.md) passes without correction. This is an all-order result in its stated class: Busch's strong bound 5^(n-1)<=H^3 gives n<=12, then the complete low-cycle bank supplies the finite classification. The bound is not applied to weak singleton padding, and the theorem excludes no graph with C>=3.

The unchanged reviewed source bytes and all executable dependencies are included. The source-relative directory shape is preserved. The source C3 completeness proof and its named structural proofs are bundled as inherited context. Only its independently established strong-extension census completeness is needed for this separation; no all-order C3 enumeration assumption is used.

Run the self-contained reproduction from this directory:

    python3 reproduce.py

It requires Python3 and g++ with C++17, uses standard libraries only, and stages all output in a temporary directory. It replays the frozen source normally and under -O, independently regenerates the full 27,945-graph bank, recounts all 12,077 C1/C2 graphs with uncapped reverse-endpoint DP, then checks the six cores with walk inclusion-exclusion, exact isomorphism search, every induced subset for feedback number, and every changing-word arc. All regenerated receipts must match the frozen bytes. The large bank CSV is regenerated rather than shipped. FILES.sha256 binds this checkpoint. No navigation or index file is changed.
