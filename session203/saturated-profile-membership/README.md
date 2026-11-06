# Saturated-state bounded-profile membership

PROVED, independently audited: an exact terminating membership algorithm from the audited
finite weighted profile presentation. For fixed positive start budget K,
only reachable active-control/forest-count states are searched. Counts are
capped at max(X,K)+1, which preserves every requested B<=X and s<=K because
all transitions and the physical B,s readouts are nonnegative.

This decides finite membership questions. It does not establish F(8,137),
cofinite numerical support, or a uniform finite verification of all targets.
The theoretical state bound is enormous. Polynomial dependence on numerical
X for fixed K is not polynomial dependence on the binary input length.

Run the reproducible complete control:

    python3 check_membership.py --K 1 --X 137 --max-states 100000
    python3 -O check_membership.py --K 1 --X 137 --max-states 100000

The outputs match NORMAL137.json and OPTIMIZED137.json byte-for-byte. The
search exhausts 26,760 states and 62,470 transitions and returns the same
41 unit-start bases as the independently published source compiler. Every
returned instruction trace is replayed to literal adjacency rows, and each
physical profile is directly recounted. In particular B=137 is absent only
at start 1; no global Hamilton-count or staggered-profile exclusion follows.

Resource-bound hostile controls:

    python3 check_membership.py --K 2 --X 13 --max-states 1000
    python3 check_membership.py --K 1 --X 1 --max-states 1

These match GUARD_K2.json and GUARD_EMPTY.json. Both return INCOMPLETE,
including the case where some verified witnesses have already been found.
No absence claim is made before the full reachable queue is exhausted.
Unexpected process interruption or an exception yields no certified result.

The proof's core search bounds concern compact membership/rectangle output
and a shared predecessor DAG. The demonstrator also expands and directly
recounts its small control witnesses; that optional validation/output work
is not included in those core complexity claims.

The independent audit, frozen original candidate, literal controls, and
portable reproduction script are under audit/. The release preserves the
algorithm and receipts byte-for-byte; only status and the reviewed scope
clarifications change the proof.
