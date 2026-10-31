# Portable replay and payload boundary

The curated payload was copied into a fresh temporary directory and its
absolute `verify_all.py` path was invoked from an unrelated `/tmp` working
directory. All five jobs passed normally and with `-O`, matching the bundled
receipts byte for byte. The full replay summary was identical before and
after relocation; it is stored as `receipts/full_replay.json`.

Every job uses only Python's standard library. No workspace-specific path,
network request, floating solver, symbolic package, exploratory source, or
certificate discovery step is required. Default replay checks the SHA-256
manifest, performs only read operations, and suppresses bytecode writes.

`WHITELIST.txt` lists every curated repository payload file, including itself
and `SHA256SUMS`. The manifest hashes every whitelisted file except itself
(to avoid an impossible self-hash). The replay requires exact agreement
between manifest coverage and the whitelist. SHA-256 provides integrity
consistency, not an independent signature or trust root.

Transport JSON files under `transport/` are delivery material, not part of
the portable repository payload; they are intentionally excluded from the
whitelist. Each contains complete GitHub create-tree arguments with a forty
zero placeholder for `base_tree_sha`. The publisher must replace that
placeholder with the current appropriate base tree SHA and must not submit
it unchanged. Substituting another forty-character SHA leaves the serialized
byte count unchanged. No commits, trees, or refs were changed by packaging.
