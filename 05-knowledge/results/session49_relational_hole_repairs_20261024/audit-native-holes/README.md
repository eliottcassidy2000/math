# Independent native relational-hole audit

Start with [AUDIT.md](AUDIT.md). The standalone checker reconstructs the original
native records and exhausts seven third-source choices times all16 relation
tables, exact branch boundary costs, and finite route/cost controls.

Run from any directory:

    python3 audit.py
    python3 -O audit.py
    python3 replay_inherited_filtration.py
    python3 -O replay_inherited_filtration.py
    python3 compare_producer.py
    python3 reproduce.py
    python3 verify_manifest.py
    python3 verify_manifest.py --relocate

`reproduce.py` runs both normal and optimized audit and filtration checks and
compares their outputs. `--relocate` copies only manifest-listed payload plus
MANIFEST.json into a fresh temporary directory, checks that closed tree, and runs
the entire reproduction suite there. No original repository, network, Git,
producer installation, absolute input path, or original working directory is
needed. All checks use explicit exceptions and survive Python optimization.

`receipt.json` is the frozen independent result. `audit.py --write` regenerates
it for development; ordinary verification never rewrites it. Producer and earlier
historical proof/code snapshots are preserved as `.txt` archives, not maintained
current proof files. Producer comparison is separate from the fresh audit: no
producer output is an input to `audit.py`.
