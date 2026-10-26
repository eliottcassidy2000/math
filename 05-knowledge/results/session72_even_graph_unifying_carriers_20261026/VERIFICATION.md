# Portable verification and scope

From this directory:

```sh
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
python3 -I -B -O verify_provenance.py
python3 -I -B test_packet_corruption.py .
python3 -I -B test_provenance_corruption.py .
```

`verify_packet.py` binds the complete explicit inventory, bytes, SHA256 hashes,
Git-compatible file modes, native-copy map, reviewed evidence, and declared
normal/optimized receipts. Replay occurs in a temporary isolated copy; no
source-snapshot program runs. Inputs and expected receipts must remain unchanged.

`verify_provenance.py` separately checks native hash lists/manifests and exact
repository source blobs. The original paths were verified against their pinned
Git commits during curation; the portable verifier recomputes the blob identities
without requiring a repository checkout.

The replay certifies only each script's stated finite universe and explicit
controls. It does not convert finite evidence into an all-order proof, promote
an inherited theorem into a new result, or establish a proposed grand theorem.
Mathematical claims and scopes are separated in the README and native reports.

The publisher also tests deliberate corruption of hashes, source pins, inventory,
paths, symlinks, modes, schema and evidence bindings in both Python modes, and
runs the repository documentation and whitespace gates before every checkpoint.
Remote verification compares exact trees, paths, bytes, modes, parent and required
ancestry; CI observations refer to exact remote commit SHAs. Zero CI checks is
an absence of configured evidence, not a successful test result.
