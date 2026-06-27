# Observer-Cut Orbit Ledger

S218 gave the right abstraction: every quotient has a next observer, and the
forgotten coordinate is often the cut payload for that observer.  The extra
precision is that the payload is an orbit, not a raw label.

```text
C_q(x,o) = orbit_Aut_q(x)(boundary slice, incidence word, extended shadow)
```

That formula is the common object behind the tournament perspective defect,
AP/GW boundary topology, q=23 Haar squares, K33 state lifts, automaton shadows,
analytic blindness, diagonal-layer rectangles, and matrix sidecar columns.

The controlled-forgetting ladder should therefore be run as a ledger.  Each
row names the quotient, the next operation, the visible automorphism group,
the cut-payload orbit, the predicate it can change, and its discharge:
reconstruct, exactify, dual-annihilate, descend, boundary-stop, or name debt.

This sharpens the question for any stuck LRC quotient.  Do not ask first
whether the scalar looks strong.  Ask which observer it will be forced to
meet next, and what orbit of the old fiber that observer sees.

Tournament Analysis also becomes cleaner.  Vertices are not runners by
default; they are payload columns and discharge modes.  A directed cycle among
payloads is a useful signal: the discharges may not commute, and the proof may
need a bicomplex/fiber-product carrier instead of a linear ranking.

The next concrete job is an HYP-2963 observer-cut ledger over coarse fibers.
The row schema should be:

```text
base_quotient
fiber_id
observer_kind
visible_automorphism_group
cut_payload_orbit_id
changed_lrc_predicate
separating_sidecar
discharge_mode
residual_debt_name
```

That would turn "controlled forgetting" from a slogan into an auditable proof
interface.
