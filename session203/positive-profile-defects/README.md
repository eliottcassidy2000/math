# Positive profile defects and exact halving

Status: PROVED / independently audited. This packet concerns physical rooted-cover inequalities and exact search within the complete start-two grammar. F(8,137) and the unconditional eventual odd-spectrum conclusion remain OPEN.

## Results

- Universal forest cancellation proves 2T+2L-C-2R-b+2d >= 0 for every rooted tournament, including the singleton.
- Every nonneutral O/Z/J/K1/K2 grammar step more than doubles its scalar B; O itself is neutral.
- Exact Z-ray formula and strict inverse descent apply to existing representations, not arbitrary scalar targets.
- Every represented start-two pair (B,2) has a representative of order at most max(2, floor(log2 B)+2), one vertex sharper than the inherited cap.

The main proof is PROOF.md. The candidate proof, independent theoretical review, independent permutation-count verifier, and original review hashes are retained under audit/. PROMOTION.json records the exact status-only proof changes and the import-path-only rehoming of the independent verifier. The original audit hash paths name the scratch source locations; audit/SOURCE_PROOF.md preserves the hashed proof, and all other hashed source files are unchanged in the packet.

## Reproduction

Run from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/independent.py
    python3 compiler.py 137
    python3 compiler.py 139
    python3 compiler.py 139 --max-states 1
    sha256sum -c SHA256SUMS

The first two outputs match RESULT.json and OPTIMIZED.json. The independent verifier uses active assertions and is intended for an ordinary, not optimized, Python run. It matches audit/INDEPENDENT.json. Main controls cover 1,099 inputs, 5,490 literal children, and 18,803 forests, including 192 negative forests with 192 distinct cancellation images. Old/new complete grammar support agrees at bounds 139, 997, and 9999. This finite agreement checks the proof; it is not used to infer an infinite tail.

The source dependency proof START_TWO_GATES.md is a verbatim frozen copy, whose dependency references retain their original packet-relative wording. The needed mathematical dependency files are supplied alongside it under dependencies/. PROVENANCE.json records verified source commits and hashes. No navigation or global index is edited.
