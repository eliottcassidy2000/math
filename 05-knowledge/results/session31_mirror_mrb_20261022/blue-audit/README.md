# Independent colored-forest audit

- `AUDIT.md`: independent derivations, prior-canon reconciliation, proof
  verdict, implementation review, and exact scope.
- `independent_checks.py`: standard-library regression checker. It imports
  no implementation from the audited submission.
- `independent_receipt.json`: exact finite results and hashes of the
  audited source files.
- `independent_checks.out`: captured successful run.

Run `python3 session31-blue-audit/independent_checks.py` from the working
root, or use `python3 independent_checks.py --source-dir PATH` from this
directory. Published packets with sibling `blue-forests` are recognized
automatically. Use ordinary Python, not `python -O`.

No Git mutation was made by this audit.
