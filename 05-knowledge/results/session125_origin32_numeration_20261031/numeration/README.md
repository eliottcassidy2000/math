# Session 125: canonical arithmetic from 1, 3, 8

This portable bundle records the classical Chung–Graham numeral system and
the verified derivations used in the origin-at-32 investigation. It provides
integer arithmetic and a bridge to Zeckendorf, Beatty/Wythoff, and quadratic
norm identities. It does **not** prove square-sum Hamiltonicity or identify
actual graph growth with numeral successor.

The exact ceiling-shift formula applies only to canonical words. Shifting
an arbitrary coefficient string and then normalizing need not give the
same result as normalizing first and then shifting.

## Read

- [PROOFS.md](PROOFS.md): the `1,3,8` Diophantine triple, canonical digits,
  uniqueness, terminating carries, minimum summand count, successor, parity,
  and the distinction between graph-order coordinates and vertex labels.
- [BEATTY_BRIDGE.md](BEATTY_BRIDGE.md): exact ceiling shift and division,
  fractional-window dynamics, primitive shift orbits, and the characterization
  of the norm-one sequence `1,3,8,21,...`.
- [Arithmetic audit](audit_port_numeration/PROOF_REVIEW.md) and
  [Beatty/orbit audit](audit_beatty_bridge/PROOF_REVIEW.md): independent proof
  review, oracle checks, and precise scope boundaries.

The numeral system and its literature are explicitly credited; no novelty
claim is made for these established constructions.

## Reproduce

Python 3, standard library only. From this directory:

```sh
python -B verify_bundle.py
python -B test_port_numeration.py
python -B -O test_port_numeration.py
python -B test_beatty_bridge.py
python -B -O test_beatty_bridge.py
python -B audit_port_numeration/independent_audit.py
python -B -O audit_port_numeration/independent_audit.py
python -B audit_beatty_bridge/independent_beatty_audit.py
python -B -O audit_beatty_bridge/independent_beatty_audit.py
python -B verify_bundle.py
```

The scripts resolve their inputs relative to their own files, not the
current working directory. They also accept an optional output-report path;
omit it to leave the checked-in reports unchanged. The arithmetic audit
regenerates its deterministic DFA certificate with the same bytes. `-B`
avoids adding bytecode caches to the strict publication whitelist.

Every substantive test uses explicit failures, so `-O` cannot remove its
checks. Each mode passes:

| Suite | Explicit checks |
|---|---:|
| Original arithmetic producer | 1,359,569 |
| Beatty/orbit producer | 1,983,501 |
| Independent arithmetic audit | 3,639,641 |
| Independent Beatty/orbit audit | 4,676,438 |

## Integrity and portability

- `WHITELIST.txt` is the exact published-file list; no caches or duplicate
  `.out` captures are part of the payload.
- `SHA256SUMS` covers every whitelisted file except itself.
- `PUBLICATION_MANIFEST.json` records payload sizes and hashes; the original
  producer manifests are preserved byte-for-byte as `manifest.json` and
  `manifest.beatty.json`.
- `PACKAGING_VERIFICATION.json` records clean replays, relocation checks and
  intentional corruption/extra-file rejection controls.

The only audit-code adaptations are source-relative import paths required
by nesting the audit folders under this bundle. Their math logic is
unchanged and their normal/optimized and relocated runs are rechecked.

Publication preparation is not proof of a remote commit; any repository
publication receipt is recorded separately by the integrating checkpoint.
