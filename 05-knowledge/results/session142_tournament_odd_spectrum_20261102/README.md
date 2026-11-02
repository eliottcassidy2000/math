# Session142: actual tournament counts and the missing legal lift

**The full spectrum conjecture remains OPEN.** For actual finite tournaments,
the target is every positive odd Hamiltonian-path count except7 and21.
This packet recovers the proven exclusions and the exact realization
obligations; rational weights, signed evaluations and modular images do not
count as tournament witnesses.

- [Frontier and proofs](frontier/proof.md): a272-case structural cover for
the small strong spectra needed by the all-order7/21 exclusion; exact
strong carry witnesses49,63,343; target-monoid atoms and maximum factor length.
- [Jacobian bridge audit](jacobian_bridge/PROOFS.md): coefficientwise parity,
continuous-relaxation interval, integral grading theorem, and explicit
unimodular one-arc responses that still fail the legal source guard.
- [Synthesis](SYNTHESIS.md), [connections](CONNECTIONS.md),
[concept board](BOARD.md), [source status](SOURCE_STATUS.md),
[verification](VERIFICATION.md).

The complete realization problem reduces to the two strong prime lanes
`p` for odd primes `p!=7`, and `7p` for odd primes `p!=3`. Those are still
unproved existence statements. The old two-lane mechanism is recovered,
not claimed as a new proof of completeness. All strong carry witnesses
are literal0/1 tournaments and are checked by two exact path counters.

The inherited three-variable polynomial has its exact determinant`-2`
and three-point collision independently rechecked. This does not construct
any prescribed tournament count. Its external formalization is **CITED,
not locally built**; attribution and arithmetic validity are distinct.

Run `python3 -B replay_checkpoint1.py` from any working directory, and
`sha256sum -c SHA256SUMS` from this packet. Frozen source snapshots retain
their original bytes, including historical claims that this packet scopes
or corrects. The narrow elliptic-Jacobian canon repair is independently
reviewed and byte-pinned; no central mistakes ledger was read or edited.

## Checkpoint2: asymmetry and legal constructor language

[Prime-lane proof](minimal_gap/PROOF.md): every prime-count witness H=p>3
is asymmetric; for prime p>=11, H=7p permits group orders1 or7. The first
symmetric7p prime quotient is sharply593, attained by an explicit strong
order10 witness with H4151 and group order7.

[Constructor proof](constructions/BINARY_NO_GO.md): finite actual seed banks
of count at least5 under singleton append/re-mark binary growth yield only
O(log X) outputs, even allowing the specified neutral S1/reset interleavings.
This does not exclude richer gadgets. Two actual positive rays are recovered
and credited to their inherited sources. See [checkpoint2](CHECKPOINT2.md).
Run `python3 -B replay_checkpoint2.py` for the combined portable checks.

## Checkpoint3: module gain and recovered actual witnesses

[Path-cover proof](path_cover_amplification/PROOF.md):
`k! p_k >= H binom(n-1,k-1)` gives `H(Q[T]) >= (n+2)H(T)` at cyclic
homogeneous ports, with an extra1 for even n, even when the context varies.
This does not cover nonhomogeneous or multi-port operations.

[Recovered witnesses](prime80557/README.md): actual strong tournaments of
orders19 and29 have prime counts80557 and100000007. Full subset DP checks
the first; an independent direct adjacency partition formula checks the
second. The constructor is inherited, not a new interval or global proof.
Its22source pins use a published commit. See [checkpoint3](CHECKPOINT3.md).

For the complete packet, run `python3 -B replay_checkpoint3.py` (Python3
standard library plus C++17 compiler required), then `sha256sum -c SHA256SUMS`.
