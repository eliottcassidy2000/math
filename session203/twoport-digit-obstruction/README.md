# Exact two-port digit obstruction

Status: PROVED / independently audited in the stated restricted architecture.

No actual seed tournament of any order can obey both binary recurrences D0=2D-kappa and D1=2D+6-kappa, D=3B+2s, with one fixed two-retained-port update per digit and one of six fixed elementary closing contexts. The affine shift kappa is arbitrary, and no start-count bound is needed. Peak width is three, including the newborn before a port is forgotten.

The theorem does not exclude other control states, longer digit blocks, more ports, larger or state-dependent closing contexts, or F(8,137). No positive F compiler is claimed.

The proof combines universal forest semantics, ten exact spectral exclusions, and 1,176 rational Farkas certificates. The released verifier requires only Python3 and checks exact arithmetic; it does not call or trust an LP solver. Literal finite tests are controls rather than the source of universality.

Run from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 audit/independent.py
    python3 audit/literal_permutations.py
    python3 audit/closure_rows.py
    sha256sum -c SHA256SUMS

Normal and optimized main outputs match NORMAL.json and OPTIMIZED.json. Independent outputs match the corresponding audit receipts. Some independent scripts use assertions and should run normally, without -O.

The exact reviewed proof, original README, audit source hashes, and original independent scripts are retained under audit/. PROMOTION.json records the status-only proof promotion and import/data-root-only rehoming of independent scripts. The original audit manifest's PROOF.md and README.md entries match audit/SOURCE_PROOF.md and audit/SOURCE_README.md; other source entries are unchanged. SHA256SUMS is the release manifest. Dependency copies retain their source-relative references, with commits pinned in PROVENANCE.json.

No navigation or global index is edited.
