# Session134: a global obstruction to bounded additive induced descent

**PROVED, with cited ingredients and exact independent controls. Global raw-H
completeness except7,21 is still OPEN.**

For any induced tournament inclusion A subset T with positive Hamiltonian-count
gain delta, its canonical nonsingleton-SCC core has at most delta+2 vertices.
Intersecting A with that core preserves both H values and gives an actual
bounded pair, yielding an exact finite gain-based decision/reconstruction
procedure. In particular,

    2 H(A) <= (delta+1)!.

Therefore fixed bounded additive gains cannot continue at arbitrarily large
counts, even with arbitrary old witnesses, arbitrarily many inserted vertices,
and neutral padding. Every nested tournament ray has zero-density raw-H image,
with O(X loglog X/log X) distinct values at most X. This does not obstruct a
branching constructor or one that reencodes old arcs/vertices.

The exact endpoint pairs at gains2 and4 are respectively

    +2: (1,3), (3,5);
    +4: (1,5), (5,9).

The +2 result follows directly from the proof. The +4 result uses a complete
strong-order<=6 enumeration justified by an all-order affected-SCC cutoff.
Actual families realize both2^k+1 and2^k+3 for every k>=6, yet no witnesses of
the two counts can nest in that direction. This is an infinite nonvacuous
barrier, not a claim to new finite global coverage.

Read [THEORY.md](THEORY.md) for proofs, the exact source/cut, cited dependencies,
inheritance ledger, and the remaining realization obligation.

## Reproduce

Python standard library only; all paths are file-relative and correctness
checks use explicit exceptions, never optimization-sensitive assertions.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 run_checks.py
```

The harness reruns each checker normally and with `-O`, compares mode outputs
and retained receipts, and checks that the two complete small-gap pairsets
agree. It produces `verification.json`; SHA256SUMS and MANIFEST.json pin files.

* `verify.py`: all33,867 labelled tournaments of orders1..6; all202,013
  one-vertex deletions and every induced subset; 182,750 positive insertion
  bound checks and2,027,318 strict induced-pair kernel checks; literal paths and odd-cycle deletion formulas through order5.
* `independent_check.py`: a different complete strong presentation cover via
  fixed Hamilton cycles; 549 presentations through order6 and3,251 literal
  deletion checks; no producer import or inherited census. Also retains the
  seven-vertex hostile17->27 to a stronger false linear guess.
* `family_check.py`: literal five/six-vertex parent/apex state and cover counts;
  independently reconstructed tail adjacencies at k6..12; mutated-cut and
  zero-tail-length hostile controls.

Normal and optimized receipts agree byte-for-byte. The infinite conclusions
are the proofs, not extrapolation from these finite controls.

## Scope and ownership

Research base: `04cdf182d9e4543800b421a8893232ccf69b2399` in the read-only
Session133 source worktree. `SOURCES.json` pins exact inherited dependencies.
No inherited frozen file was edited. The parent Session134 task owns repository
integration, publication and remote verification; this child did not push.
`PUBLIC_FILES.txt` is the exact publication list and excludes exploratory probes.
