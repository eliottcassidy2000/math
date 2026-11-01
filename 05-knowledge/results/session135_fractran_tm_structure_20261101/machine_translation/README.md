# Actual Turing machines for repository FRACTRAN streams

**Result:** explicit literal transition tables, certified checkpoint decoding and exact runtime—not just universality.

- S132's 22-fraction pure-blue stream: **one tape, 76 states, 257 symbols, 10,926 transitions**.
- S133's 43-fraction fixed-H multiplicity tail: **15 tapes, 168 states, 3 symbols, 293 transitions**.
- S134's radix-80 migration: invariant-face reduction from **25 to 13 fractions**, nine to five data registers/tapes, retaining the signed invariant and actual graph seed.
- General finite FRACTRAN → register compiler preserves first-eligible priority by rollback; reverse register → FRACTRAN compiler is cancellation-safe and costs exactly two fractions per register instruction.

For S133's general positive matrix compiler, x→y=Ax at clock n has exact TM cost

    ||x||² + 3||y||² + 5|x|₁ + 9|y|₁ + 6d + 2n + 4.

See [PROOF.md](PROOF.md) for all-input statements, tape convention, domain restrictions, costs, inherited interpretation scope, and loss ledger. No minimum-size, binary-efficiency, novelty, or undecidability claim is made.

`*.tm.json` tables (or their authenticated `*.tm.rows.*.json` plaintext shards, each <=90KB) have rows `[state, read, next_state, write, move]` with no wildcards. The 15-tape fields each contain exactly 15 explicit symbols/moves; inactive heads read/write `#`, not arbitrary symbols. Use `machines.load_tm(path)` to authenticate and assemble sharded tables. `*.ir.json` exposes the finite register intermediate representation. `*.roundtrip.json` contains the reverse-compiled actual fractions. `machines.py` contains the independent literal-table interpreter and compilers; `inherited_linear.py` is the pinned S133 compiler. Build with `python3 build.py`; verify with `python3 verify.py` and `python3 -O verify.py`.
