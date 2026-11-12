# Exact growth syntax and endlessly arriving interfaces

[Proofs](PROOF.md) · [standalone verifier](verify.py) · [finite receipt](receipt.json)

The complete **marked batch-extension** state is:

1. old pair-free core and its ordered transitive chains;
2. newborn tournament;
3. a proper word of newborn-adjacency signatures on every old chain;
4. positive lengths for those signature runs.

This reconstructs every child exactly. For a fixed signature template, a
finite cycle index gives a rational generating function for every future
size. The rules do not require a new case-by-case census at size50 or200.
But the templates are unbounded, and forgetting the old subset needs the
specified symmetry correction, not naive division.

Three useful consequences:

- **Observable evolution:** after `k` uniform labelled arrivals, the number
  of surviving old homogeneous edges is exactly `Bin(n-b,2^-k)`. The child
  core rank differs from `n+k` minus that binomial by a controlled rare-new-edge
  correction. Old transitive compression is lost on a logarithmic scale.
- **Structural births:** complete disjoint `k`-boundary interfaces first
  appear at exactly `k+2^k` vertices:3,6,11,20,37,70,135,264,... . This
  external notion differs from ordinary neighborhood VC dimension, whose
  threshold `2^k` is inherited from THM-1450.
- **Hamiltonian meaning:** the exact skeleton and weights retain the full
  path-cover/walk-content substitution data. There is no claim that a scalar
  path count alone is a complete future state, or that the grammar solves
  the inverse target problem in Chan–Pak.

Run from this directory:

```sh
python3 verify.py
python3 -O verify.py
```

Both saved receipts agree. Finite controls verify69,706 labelled ensemble
instances across the listed `(n,k)` tests,180 random normal forms,32 independent
Hamiltonian count comparisons, marked-family orbit coefficients, and explicit
birth witnesses. The same tournament may occur in different listed ensembles;
this is not69,706 distinct isomorphism classes. Independent review is recorded
in the session's audit directory.
