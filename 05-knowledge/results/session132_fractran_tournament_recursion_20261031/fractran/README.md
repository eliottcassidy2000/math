# Finite FRACTRAN programs, not just computability

- **A000568:** a fixed 5,543-fraction program computes the number of tournament
  isomorphism classes. Input `157*5^n`; output is the terminal exponent of 3.
- **Literal metagraph lines:** a fixed 16-fraction program computes
  `ell_n=2^(binom(n-1,2)-1)`, n>=3. Input `306*13^(n-3)`; terminal state is
  `59*2^ell_n*3^(n-1)`.

[THEORY.md](THEORY.md) gives the complete short fraction list, exact compiler
semantics, invariants, termination, recurrence proofs, and information-loss scope.
[VERIFIED.json](VERIFIED.json) records finite replay evidence and hashes.

Quick actual-integer examples (no compiler invocation):

```bash
python3 run_program.py a000568 3     # output 2; 42606 steps
python3 run_program.py literal16 6  # output 512; 3142 steps
```

Full normal/optimized verification:

```bash
python3 verify.py
cp VERIFIED.json VERIFIED.normal.json
python3 -O verify.py
cmp VERIFIED.normal.json VERIFIED.json
```

Both generated A000568 artifacts are below 100KB. The source and compiler are
standard-library-only; no table of A000568 values is used by the program. The
large program is a transparent correctness construction, not an efficiency or
minimality claim. The line recurrence is inherited carrier geometry compiled
into a new program; the existing mathematical recurrence is not claimed as new.
