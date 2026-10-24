# Independent order-nine prediction audit

Read `AUDIT.md` for the independently checked proof, complete numerical result,
finite controls, scope and source boundary. The computed count is 191,536,
from upper trace 1,724,582 and structural error correction 758.

Reproduce: `bash reproduce.sh`. Only a C++17 compiler and Python 3 are needed.
All source inputs and comparison snapshots are local and hash-pinned. No
network, Git, producer code import, known-target input or order-nine
canonicalization is used. `verify_manifest.py` checks the closed packet.
