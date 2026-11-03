# Clock depth, one-vertex readers, and the boundary of observation

2026-11-03. **PROVED** results are in [PROOFS](PROOFS.md); independent finite
algorithms and scope checks are in [AUDIT](AUDIT.md). This packet extends
the newly merged S168 ring simulation, rather than deriving a substitute
rule for the user's quoted period.

The inherited rule is R(a,b)=(b,a+b+Sb) over F3. Its literal tournament
realization has order12m and constant full path-cover profile. Our extensions:

- At m=3^a, every nonzero source of surviving u-adic depth d has exact
  period8*3^ceil(log_3 d), for **every physical lift**. The impulse attains
  8m, whereas spatially constant nonzero states have period8. Complete
  cycle multiplicities follow.
- One new vertex, attached at a fixed named port, reads all six literal
  switch states by the six counts109,115,101,77,73,81. This includes the
  deck bit forgotten by the unmarked isomorphism class.
- At every width m, one boundary site's first2m trit observations recover
  the entire logical state. The other2m-1 physical deck bits stay invisible
  forever under that clock and that fixed probe. This is an iff statement
  about future observations, not just a counting lower bound.
- The observed deck bit carries partial information about the next trit.
  An exact three-trit carry formula counts every shorter physical-history
  class and proves2m is still the sharp worst-case consecutive-history
  length. This is not a bound for adaptively spaced probes.
- The cyclic polynomial `(T^2-T-1)^m-T^m` gives logical periods at every
  width. At the prime-to-three widths divisible by8, its repeated quadratic
  introduces an extra factor3. The physical lift's possible doubling is
  measured by two boundary winding parities. Width11 has logical
  order968=8*11^2 and physical impulse period1936, so odd width alone
  does not prevent doubling.
- The compatible infinite formal-series limit has no nonzero periodic
  states, despite finite periodic truncations at every level. A scheduled
  single-gate clock also differs from a round clock: the impulse has
  8m rounds but32m^2 gate attempts when idle attempts are counted.

[CONNECTIONS](CONNECTIONS.md) records the Anchor/Niche/Wildcard board and
maps to older golden clocks, joint cubic tournament observers, signed-unit
transport, LCA/Zeckendorf carries, July odometer obstructions and fixed-twist
extension laws. It ends with four precise open helper questions.

The open Hamilton-path spectrum classification and Collatz convergence
remain open. In particular, profile-preserving computation does not prove
universal return to a source that realizes3H+2 or3H+4. The new probe shows
exactly one way to retain boundary data that the full profile loses.

## Reproduce from the repository root

```text
python3 05-knowledge/results/clock_boundary_simulation_20261103/check.py
python3 05-knowledge/results/clock_boundary_simulation_20261103/audit.py
python3 -O 05-knowledge/results/clock_boundary_simulation_20261103/check.py
python3 -O 05-knowledge/results/clock_boundary_simulation_20261103/audit.py
```

Both scripts use only the Python standard library and explicit exceptions.
The audit imports no research implementation. `source_pins.json` records the
reviewed inherited files at the integration revision. Optimized reruns must
leave both JSON outputs byte-identical. The audit is by this same agent,
with separate algorithms; no subagent or external referee was used.
