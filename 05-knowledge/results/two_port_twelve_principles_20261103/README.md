# Two-port twelve-vertex principles

**PROVED scoped algebra; FINITE-EXACT controls.** Universal reconstruction,
the general twelve-vertex module residual, and exact even endpoint imbalance
remain **OPEN**. No exact full-response collision is claimed.

[T1-T4](PROOF.md) give:

* Full response modulo four descends through every ambient quotient for
  arbitrary same-size module pairs. The multiplier is the odd H(Q).
* The modular converse contrast has a complete cubic formula using paired
  paths. Its shared-label cubic coefficient at twelve is the start-parity
  sum over supports of sizes 2,4,6,8 omitting the other named vertex.
* The old C-skew modulo-four test is equivalently a parity count of odd-support
  paths between each pair. This supplies a separate implementation.
* Exact core/converse comparison needs at most 1+2n+n^3 queries: 1753 at
  twelve. General tournament pairs can still have a quartic difference.

The formula keeps the full response and its paired endpoints. None of these
results equates it to H,S,E,C or promotes a necessary congruence to a collision.
The inherited exclusion through module size eleven remains credited to
[the complete sieve](../two_port_seven_core_20261103/residual/RESULT.md).

Reproduce the controls with `python check.py --output check.json` and
`python -O check.py --output check.optimized.json` from this directory.
The [board](BOARD.md) names the inheritance, hostile probes, maps and losses.
