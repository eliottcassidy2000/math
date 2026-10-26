# Frozen evidence and reproduction

The packet preserves the explicitly curated producer and auditor payloads
byte-for-byte, including their native manifests, source ledgers, exact reviewed
inputs, and deterministic normal/optimized receipts. `NATIVE_COPY_MANIFEST.json`
records the native-to-packet mapping without rewriting native manifests.
Source snapshots are inert evidence, not additional programs to execute.
Downloaded third-party PDFs, page images, OCR, exploratory scripts and
timing-bearing logs are excluded.

From this packet directory:

```
python3 -B verify_packet.py --replay
python3 -B -O verify_packet.py
```

The first command validates all frozen payload bytes and modes, then replays
every declared current checker in **both normal and optimized Python**, in a
temporary copy outside the repository. It compares actual receipts to immutable
expectations and never updates them. The second command independently checks
the full envelope with Python optimization enabled. The manifest specifies the
exact replay commands and their finite universes.

The independent guard checker emits timing-bearing progress on stdout. Its
`--receipt` output is instead directed to a fresh temporary file outside the
packet, and only that deterministic result is compared with its frozen normal
and optimized receipts. The native guard envelope enforces its own exact
seven-file universe; its generated envelope receipts therefore live outside
that native directory.

`verify_provenance.py` checks native hash lists/manifests, exact reviewed file
hashes and replay equivalences, source SHA-256/Git blobs, and exact commit pins.
During publication, every repository source snapshot was additionally checked
against actual pinned `git show COMMIT:path` bytes and Git blob IDs. Classical
source citation/inspection evidence is preserved in the guard audit; external
source bytes are not republished.

Native mathematical controls include proper-locality parity collisions,
intentionally false bare-K5 lifting, and corrupted tournament reconstruction.
The last two are required to exit nonzero with their specific rejection
messages in both normal and optimized modes. Other native hostiles are caught
explicitly and recorded in deterministic receipts. None uses removable
`assert` statements as its sole rejection gate.

`FILE_MANIFEST.json` is a frozen envelope, not a cryptographic signature. The
publication journal separately tests corruption, deletion, added files, wrong
modes, symlinks, malformed paths, duplicate entries, source pins, and altered
reviewed evidence under normal and optimized Python. Exact remote byte/mode/
parent verification is separate from CI. Empty GitHub check, workflow, or status
lists mean no CI evidence, not a CI pass.
