# Lossless receipt storage

Large receipt files are stored as deterministic gzip data, base64-encoded and
split into numbered text chunks. Each `ARCHIVES_CP*.json` records the original
relative filename, exact byte count, SHA-256, Git blob hash, and compressed
byte count/hash. All entries, including named zeros and full response tables,
are preserved; compression does not summarize the mathematical data.

Before replaying a packet that uses an archived file, run from this directory:

    python3 restore_receipts.py ARCHIVES_CP02.json

For any later archive manifest, pass its exact filename instead. Verification
without writing is available with `--verify-only`. The restoration script uses
only the Python standard library and will not overwrite an existing file with
different bytes. Manifest paths are relative to this directory.

Every chunk and ordinary file was retrieved at its immutable publication
commit and compared with the staged exact bytes and Git blob hash. Each
archive was independently reassembled and compared with the complete original
receipt before publication. An absent GitHub CI report is not a CI pass.
