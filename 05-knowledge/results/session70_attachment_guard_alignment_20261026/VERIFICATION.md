# Frozen evidence and reproduction

This checkpoint preserves 33 curated producer/auditor files byte-for-byte,
including both native hash lists, 18 inert source snapshots, and the precise
Vitali files accepted by the independent audit. `NATIVE_COPY_MANIFEST.json`
records the native-to-packet mapping without rewriting native manifests.
Source snapshots are historical evidence, not additional programs to execute.

From this packet directory:

```
python3 -B verify_packet.py --replay
python3 -B -O verify_packet.py --replay
```

The verifier checks every payload byte and mode, checks the exact reviewed
inputs and replay equivalences, then runs copied current scripts in a temporary
directory outside the repository in both normal and optimized Python. It never
updates an expected receipt. `verify_provenance.py` independently checks the
native hash lists, reviewed hashes, source SHA-256/Git blobs, and exact commit
pin syntax. During publication, every source snapshot was additionally compared
to the actual pinned `git show COMMIT:path` bytes and Git blob ID.

The mathematical controls in the native verifiers include a deliberately false
bare-K5 lifting assertion that must raise an exception, and parity-law pairs
that collide under proper local observables but separate under the full native
guard. These are explicit non-assert checks and survive `-O`.

`FILE_MANIFEST.json` is the frozen envelope, not a cryptographic signature.
The publication journal separately tests payload corruption, deletion, extra
files, wrong modes, symlinks, malformed paths, duplicate entries, source pins,
and altered reviewed evidence under normal and optimized Python. Empty GitHub
check/workflow/status lists mean no CI evidence, not a CI pass.
