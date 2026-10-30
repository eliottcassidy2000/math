# Session 109: Pisano clocks, ramified symmetry, and exact attachment

2026-10-30. **Checkpoint 1: ramified golden attachment proved and independently checked.**
The session remains in progress: affine-carry and observation-clock branches are
under independent review. This checkpoint does not claim those branches complete.

## Verified core

In the quadratic ring R_k = Z[phi]/(5^k), phi^2=phi+1, let g=phi^4.
There are exactly 2^(5^k-1) labelled Cayley tournaments invariant under g.
Their independent orientation bits obey

    B_0=0,  B_k=5 B_(k-1)+4.

Each nonzero orbit below has five orbits above it. The reduction kernel adds
four antipodal orbit pairs. A fixed lower tournament has exactly 16 coherent
lifts, where coherence concerns only differences with nonzero reduction.
Arbitrary invariant tournaments need not be coherent. Multiplication by phi^2
reverses every tournament in this bank.

These are specified-symmetry labelled counts on 25^k vertices, **not A000568**
and not unrestricted tournament-isomorphism counts.

- [Independent proof, hostile examples and audit](audit/RAMIFIED_ATTACHMENT_AUDIT.md)
- [Executable audit](audit/check_ramified_attachment.py)
- [Exact audit output](audit/ramified_attachment_audit.json)
- [Historical recovery and general ramified candidate](history/RECOVERED_PISANO_ADJACENT_WORK.md)

The independent audit enumerates all points at k=1,2,3, including 15,625 vertices
at the last level. Normal and optimized Python executions agree, with validation
checks retained under optimization.

The prior repository result `05-knowledge/results/unit_lift_clock_20261015a.md`
already handles all depths in unramified split and inert rings. This work extends
its explicitly excluded ramified boundary. Classical order lifting, Pisano
periods, and orbit counting are credited rather than claimed newly discovered.
No Collatz convergence, Wall–Sun–Sun classification, or LRC result follows.
