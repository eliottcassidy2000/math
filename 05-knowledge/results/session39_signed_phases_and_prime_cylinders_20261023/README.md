# Session 39: signed phases and prime cylinders

2026-10-23. **Completed research packet, independently audited.**

Start with [the final synthesis](SYNTHESIS.md) and [research map](RESEARCH_MAP.md).
The two publication checkpoints preserve completed proofs, independently checked
finite controls, and exact reproduction artifacts. Each claim retains its scope.

## Scope and routes

- [Signed valuation cylinders and primes](negative-collatz/SIGNED_CYLINDERS_AND_PRIMES.md): exact finite-word residue classes; prime/composite indistinguishability for finite valuation-only observations. No global Collatz convergence claim.
- [Flux conventions and tournament algebra](flux/FLUX_TWISTS_AND_TOURNAMENTS.md): conditional K-matrix interpretation, with physical conventions and algebraic obstructions separated.
- [Native characteristic refinement](flux-characteristic/NATIVE_CHARACTERISTIC_ADDENDUM.md): graph-lattice characteristic data and exact native parity checks.
- [Critical linking controls](linking/README.md): independent native terminal witness, scaling obstruction and small same-tower quadratic controls.

The native numerator 153/256 and incoming 105/256 are generator-dependent
presentations of the same cyclic pairing: 9^2*153 = 105 modulo 256. Both
normalize to 1/256. This is not a new invariant or a contradiction; see the
[concurrent frontier](../metagraph_five_port_phase_frontier_20261023.md).
No physical phase realization, unrestricted invariant completeness, LRC result,
or global Collatz classification is claimed.

- [Independent audit and producer pins](audit/AUDIT.md). Standalone marked-cell audit controls accompany the later producer packet.

- [Marked cyclic cells and affine transport](marked-cells/README.md): holonomy, phase-blindness, and the scoped static-rooted-observer boundary.
- [Closed-target cylinders](negative-collatz/CLOSED_TARGET_CYLINDERS.md) and [coverage-clock boundary](negative-collatz/COVERAGE_CLOCK_BOUNDARY.md).
- [Converse is not time reversal](linking/CONVERSE_IS_NOT_TIME_REVERSAL.md).
- [Extension audit](audit-extensions/AUDIT_EXTENSIONS.md) and [converse/rooted audit](audit-final/AUDIT.md).

- [All-order rooted boundary](marked-cells/ALL_ORDER_ROOTED_BOUNDARY.md), [independent lift audit](audit-lift/AUDIT.md), and [reconstruction precedents](RECONSTRUCTION_PRECEDENTS.md).

## Portable reproduction

Run from this packet directory, with Python 3 (SymPy for the flux checker).
The previously published Session38 wiring is a required, pinned input and is
not duplicated here. Repeat each Python command with `-O` before its script.

```sh
python3 -B negative-collatz/verify_signed_cylinders.py
python3 -B flux/check_flux_twists.py
python3 -B flux-characteristic/check_characteristic.py ../session38_incidence_ports_electrical_20261023/wiring
python3 -B linking/check.py --wiring ../session38_incidence_ports_electrical_20261023/wiring
python3 -B audit/audit_exact.py
python3 -B audit/audit_graph.py
python3 -B marked-cells/check_marked_connections.py
python3 -B marked-cells/check_rooted_cutoff.py
python3 -B marked-cells/check_lifted_cutoff.py
python3 -B audit-lift/check_lift.py
python3 -B negative-collatz/verify_closed_targets.py
python3 -B negative-collatz/verify_coverage_clock.py
python3 -B linking/check_converse_obstruction.py
python3 -B audit-extensions/check_closed_targets.py
python3 -B audit/audit_terminal.py --edges ../session38_incidence_ports_electrical_20261023/wiring/edges.tsv.gz --witness linking/n8_terminal_lift.json
```

The publisher replayed the listed fifteen scripts normally and optimized (30 executions). Standard
output agrees byte-for-byte with the corresponding frozen receipt or log;
the linking checker also compares its computed JSON with its saved receipt.
Normal-mode assertions, where present, supply validation; optimized execution
agreement alone does not establish that assertions remained active.
Source mathematical files are unchanged by publication. Historical scratch
commands in producer notes retain provenance; the commands above are portable.
Only explicit, approved generated/public-source artifacts are included. Prior
private attachment critiques, copied papers and withheld audits remain excluded.
