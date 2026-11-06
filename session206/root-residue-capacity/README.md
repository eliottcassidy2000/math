# Session 206: root residue capacity and a two-source obstruction

Status: independently audited all-order result. For every fixed odd modulus p, the rooted-edge O/Z parameter count in each residue class is

A_a(X) = (log(2)/p)X + o_p(X).

The exact root converse involution then bounds the number of distinct outputs in each odd class modulo 2p by (log(2)/(2p))X+o_p(X). Consequently the root image misses a positive proportion in every such class. Since the physical JJ source (5,6,2,3) never produces a value congruent to 1 modulo 3, the root-plus-JJ bank is not cofinite among odd integers. The earlier unweighted capacity bound does not exclude this bank.

This does not prove a distinct-output density limit, uniform convergence for growing moduli, or impossibility of every finite physical bank. The general-seed weighted-capacity extension is not part of this checkpoint.

## Contents and reproduction

- ROOT_RESIDUE_CAPACITY.md: unchanged independently reviewed proof, using a finite primitive digit chain and a rigorous endpoint-cylinder cutoff argument.
- verify.py and RESULT.json: finite stationarity/connectivity controls and exploratory parameter residue counts.
- audit/REVIEW.md and audit/REPLAY.json: frozen independent review and optimized author-script replay.
- audit/verify_loops.py and audit/LOOP_RESULT.json: independent exact modular-loop controls, 980 checks across 49 odd moduli through 99.
- SHA256SUMS: frozen payload integrity manifest.

From this directory run:

    python verify.py > /tmp/root-residue-replay.json
    cmp RESULT.json /tmp/root-residue-replay.json
    python audit/verify_loops.py > /tmp/root-residue-loop-replay.json
    cmp audit/LOOP_RESULT.json /tmp/root-residue-loop-replay.json
    sha256sum -c SHA256SUMS

The finite sample for modulus 7 is still visibly far from its limiting coefficient. Convergence can be slow and is not asserted uniformly in the modulus. These computations are regression evidence, not the proof of the limiting statement.

## Established companion results

- [Exact parameter capacity](https://github.com/eliottcassidy2000/math/blob/ea15da8016b2623c3f10365202301f492a7e9a84/session206/parameter-capacity/PARAMETER_CAPACITY.md).
- [Root converse collision and d=1 collapse](https://github.com/eliottcassidy2000/math/tree/e5686802c8998e01db36cf714d98793ef2b13bdd/session206/root-oz-noncofiniteness).
