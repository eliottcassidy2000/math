# Independent Session212 audit

Read [AUDIT.md](AUDIT.md) for the PASS verdict, mathematical checks, exact
finite universes, boundary cases and inherited-evidence limits.

This directory is installed beside `cycle-bridge`, `cycle-spectrum`, and
`cycle-levels`. The independent scripts also recognize the original scratch
folder names. They do not import research builders or counters.

From this directory:

```sh
python3 audit_binary_cycles.py > /tmp/binary-normal.json
python3 -O audit_binary_cycles.py > /tmp/binary-optimized.json
cmp /tmp/binary-normal.json /tmp/binary-optimized.json
cmp /tmp/binary-normal.json binary_audit.json

python3 audit_packet.py > /tmp/packet-normal.json
python3 -O audit_packet.py > /tmp/packet-optimized.json
cmp /tmp/packet-normal.json /tmp/packet-optimized.json
cmp /tmp/packet-normal.json packet_audit.json

python3 audit_refinements.py > /tmp/refinements-normal.json
python3 -O audit_refinements.py > /tmp/refinements-optimized.json
cmp /tmp/refinements-normal.json /tmp/refinements-optimized.json
cmp /tmp/refinements-normal.json refinements-normal.json

g++ -std=c++17 -O2 -Wall -Wextra -Werror brute_census.cpp -o /tmp/cycle-audit-brute
/tmp/cycle-audit-brute ../cycle-levels > brute_census.out
python3 verify_brute_histograms.py > /tmp/brute-histograms.json
cmp /tmp/brute-histograms.json brute_histograms.json
```

For full research-census replay use the commands in each sibling packet.
`APPROVED_PACKETS.json` pins the exact reviewed packet manifests and index.
`REPLAY_RECEIPT.json` records independently replayed outputs and source hashes.
No compiled binaries or duplicate working copies are part of this release.
