# Reproduction and scope

Run from the repository root with standard Python3:

    python3 05-knowledge/results/session281_chanpak_20261112/receiver/verify.py
    python3 05-knowledge/results/session281_chanpak_20261112/selection/verify.py
    python3 05-knowledge/results/session281_chanpak_20261112/wildcard/verify.py
    python3 05-knowledge/results/session281_chanpak_20261112/new_operations/verify.py
    python3 05-knowledge/results/session281_chanpak_20261112/audit/independent.py

Repeat with `python3 -O`; checks use explicit exceptions rather than disableable assertions. The selection verifier also constructs residue witnesses:

    python3 05-knowledge/results/session281_chanpak_20261112/selection/verify.py --witness 360 21

This produces an actual count1821 congruent to21 modulo360, not count21.

Optional bounded discovery:

    g++ -O3 -std=c++17 05-knowledge/results/session281_chanpak_20261112/receiver/search.cpp -o /tmp/s281-search
    /tmp/s281-search 14

The discovery census is explicitly limited to its deletion templates. It is not an all-tournament nonexistence theorem. The twelve-count order23 family is certified by exact physical fusion and its established count law, not by claiming a literal2^23 dynamic-programming replay for every output.

Wildcard uses the credited Session275 carrier helper. Independent audit imports no producer modules. Source hashes in audit/SOURCE_HASHES.json identify the restored release files; they are not hashes of the lost runtime filesystem.
