# The two-unit corridor and the first actual owner collision

**PROVED scoped, independent audit pending / FINITE-EXACT controls.**

- [Two-unit cells](two_unit_cells.md) extend the exact total-mass formula
  through 0<=delta<=2k-1 for k>=3, preserving every external matching port.
- [The first owner correction](first_owner.md) closes delta=2k for k>=4.
  The only nontrivial owner fibre is one two-empty-gap matching record
  plus two P3/one-unit-gap records representing the same pointed C3.
  The exact correction to the corridor expression is K_k-2H_k, or
  -3(k-1)(2k-3)!! without the strong restriction.
- The incoming `zero_heavy_next_20261020h` already proves deficit k and
  supplies the full-matching two-unit-gap cell argument. This addendum
  inherits it explicitly. Its mixed missing/unit port distinction and
  first-owner correction are the further steps, not a new claim for
  the incoming boundary.
- [Audit and exact finite scopes](AUDIT.md).

The first-owner full small censuses give 114/11 classes at k=2 and
10,692/1,455 at k=3 (all/strong). These exceptional all-unit cases are
FINITE-EXACT extensions of the stated all-k proof range.

Reproduce from this folder, or prefix these paths from the repository root:

    python extension_check.py
    python local_ports_check.py
    python first_owner_check.py
    python owner_four_check.py

Each also passes with python -O and gives the corresponding stored JSON
receipt byte-for-byte. verify_addendum.py authenticates frozen files and
replays all four in both modes. The included zero_heavy_check.py is an
identical, hash-checked helper from the base packet, not an independent
implementation. local_ports_check.py uses no helper imports and independently
canonicalizes the crucial finite marked cells.

This bounded continuation stops at delta=2k. No result or computation for
delta>2k is included. It does not give an unrestricted tournament recurrence.
