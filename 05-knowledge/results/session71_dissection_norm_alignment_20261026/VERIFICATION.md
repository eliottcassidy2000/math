# Verification contract

Run from this packet directory:

```sh
python3 -I -B verify_packet.py --replay
python3 -I -B -O verify_packet.py --replay
python3 -I -B verify_provenance.py
```

Python 3 and SymPy 1.14.0 are the replay dependencies. The norm producer uses
only the standard library; the independent norm/rank audit uses SymPy.

`FILE_MANIFEST.json` freezes every payload file's size, SHA-256 and mode.
`NATIVE_COPY_MANIFEST.json` records byte-preserving curation; producer and audit
manifests remain unchanged. `PROVENANCE_RULES.json` actively checks selected
native manifests, source Git blob identities, permitted exact source pins and
reviewed output bindings. Repository origin paths were additionally compared
against their exact Git commits before publication.

The portable replay runs only declared current scripts in a temporary copy,
both normally and under `python -O`; it compares deterministic output bytes to
native expectations and rejects changed, added or deleted frozen inputs. It
does not execute repository or paper source snapshots. Reports distinguish
proofs, CITED results and finite universes. Mathematical hostile controls are
inside the native checkers; package corruption tests are a separate integrity
check, not mathematical evidence.

The two first-checkpoint mathematical replays contain 824,939 producer checks
and 38,815 independent audit checks, plus five independent golden-template
area controls. They do not independently certify the
paper's entire geometric case analysis. Full optimality remains CITED in the
no-flip polygonal category.

External publication receipts live in the session's durable publication
journal: exact Git tree, modes, bytes, parent and ancestry, expected-head lease,
read-back of main, and exact-commit CI queries. Empty CI is not a pass.
