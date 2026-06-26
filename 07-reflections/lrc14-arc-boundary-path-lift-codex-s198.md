# LRC14 Arc-Boundary Path Lift - Codex S198

The useful sideways pull this session was from the old path-homology and
deletion-persistence work.  I did not put path homology on the runner
tournament.  That would be the wrong quotient for the current LRC predicate.
The live topology object is the closed danger-arc Cech nerve from HYP-3025 and
HYP-3030.

The S198 script adds a GF(2) boundary-lift layer to that nerve.  It scans all
`21913` HYP-2963 default rows for exact boundary/open status, finds the same
two residue-terminal status collisions as HYP-3030, and computes path-boundary
representatives only on those collision rows plus named controls.

The result is clean:

```text
path_lift_rows=41
closed_h1_rows=2
AP and GW only
```

In the AP collision fiber:

```text
size=30
status={'boundary':1,'open':29}
h1={1:1,0:29}
```

In the GW collision fiber:

```text
size=11
status={'boundary':1,'open':10}
h1={1:1,0:10}
```

The boundary rows have explicit representatives:

```text
AP:        C0/C1=1/1, d1=90,  d2=84, rep_edges=58
GW 12->24: C0/C1=1/1, d1=102, d2=90, rep_edges=58
```

For both, deleting any speed in the representative owner support kills the
closed-H1 signal.  That is the old deletion-persistence idea, but applied to
danger arcs rather than tournaments.

The proof implication is modest but useful.  HYP-3030's topology gate should
not stop at the scalar phrase "closed beta1 equals one."  The theorem target is
an owner-essential path-boundary cycle: every zero-open packet should either
carry an AP/GW-type closed danger-arc H1 representative, or route to named
F7/THM-572/harmonic residual debt.

Tournament Analysis used proof carriers as vertices:

```text
arc_boundary_path_lift >
arc_cech_beta >
owner_deletion_persistence >
coarse_et_unit_status >
residue_terminal_word
```

This is not a creativity ranking.  It is proof order.  The arc-boundary path
lift is more expensive than the beta scalar, but it carries the actual cycle
and owner-deletion data needed to make the equality atom local.
