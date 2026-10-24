# Validation and provenance boundaries

The twelve mathematical/audit packets retain their own manifests, exact scopes,
checkers and receipts. Manifest filenames intentionally preserve their existing
case. Standalone replays were run from actual copies containing only the declared
payload, normally and with Python optimization; the compiled order-nine controls
also ran with optimization and NDEBUG. Runtime binaries and Python caches are
not deliverables.

The root integration check verified every selected file's byte hash, rejected
symlinks and generated Python caches, required UTF-8 text, and checked the relative
Markdown links. Final remote publication verification is a separate publisher
step; local mathematical tests alone do not prove a commit was published.

Two provenance distinctions must not be silently changed:

- audit-boundaries contains historical producer-manifest snapshots as reference
  documents. Their complete executable producer packets are the sibling
  identifiability and inertia-counts directories, not an incomplete nested copy.
- audit-order-nine retains an earlier provenance-only SOURCE_PINS.json snapshot.
  The final producer aligned one archived holonomy theorem to its final wording,
  changing that metadata entry, not the mathematical code or results. The auditor
  verified all eleven archived mathematical/input/code/result files, including
  the final order-nine theorem, against final producer bytes. It separately
  verified final producer manifest
  aa3b9f61f0ab862a6072ad7b5b5710e83697be5c3ae0b118017dc988b373fc19
  and its successful42-check final manifest-only replay. The older provenance
  snapshot is not asserted to equal the updated metadata and is not an input to
  the audit's computations.

The latter replay attestation is supplied as
[validation/order-nine-final-replay.json](validation/order-nine-final-replay.json),
outside the producer's hashed payload to avoid a self-referential receipt.

The root synthesis was reviewed separately for growth normalization, the
regular-sector extension and exact order-nine classification. All numerical
claims retain the scope of their actual tests. No fixed-size census is substituted
for an all-size proof, and no known-answer comparison is counted as the predictor's
input-generation method.
