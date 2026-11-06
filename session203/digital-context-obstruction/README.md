# Fixed one-port digital-context obstruction

Status: PROVED / independently audited, in the stated restricted scope.

O/Z already encodes every integer at least two in its R coordinate: R maps to 2R or 2R-1. Nevertheless no single fixed one-port closing context can convert the entire orbit to B=alpha R+beta b+gamma e+constant. Its B readout has unavoidable positive h coefficient H(G-anchor), and eight exact orbit states span the hyperplane s=d.

This is not a counterexample to F(8,137). Multiple ports, context choices on different word classes, target-dependent rebuilding, and readouts retaining h or L remain open. Only one fixed context acting on the whole O/Z orbit is excluded.

PROOF.md contains the all-order argument. The original proof and independent review are retained under audit/. PROMOTION.json records the status-only proof promotion and import-and-data-path-only rehoming of the independent script. SOURCE_SHA256SUMS is the original source freeze, whose PROOF.md hash matches audit/SOURCE_PROOF.md after promotion; SHA256SUMS is the release manifest. Verbatim dependency copies preserve their original reference wording, with pinned source commits in PROVENANCE.json.

Run from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/independent.py
    sha256sum -c SHA256SUMS

The first two outputs match NORMAL.json and OPTIMIZED.json. The independent verifier uses assertions and is intended for an ordinary Python run, matching audit/INDEPENDENT.json. Controls include the exact determinant -1152, binary register intervals through word length8, and literal context/input counts. Exploratory finite context searches are not needed for the theorem and are not included in the release.

No navigation or global index is edited.
