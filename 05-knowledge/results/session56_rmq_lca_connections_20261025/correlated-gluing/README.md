# Correlated overlap repair

Status: elementary proofs plus exact replay; no general constant-time overlap result.

- `THEOREM.md`: original-choice natural join, diagonal-copy law, once-only weights, formal projectors, lossless-join boundary, classical join-tree criterion and zero-safe DP, complexity/loss ledger.
- `verify.py`: self-contained replay using the pinned Session55 arithmetic helper.
- `certificate.json`: exact 112-versus-784 control (128 fabricated distinct arrays), exhaustive tree tests, and adversarial boundaries.
- `SOURCES.md`: inheritance and classical-source audit.

Replay: `python verify.py --verify` and `python -O verify.py --verify`.
