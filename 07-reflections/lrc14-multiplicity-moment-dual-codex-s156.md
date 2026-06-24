# LRC14 Multiplicity-Moment Dual - Codex S156

This pass tries a genuinely different proof lens from the apex-aperture and
few-apex lift machinery.  I ignored where the safe intervals are and looked
only at the integer random variable

```text
X_S(t) = number of runners dangerous at t.
```

The useful dual fact is simple and rigorous: if `P(0)<0` and `P(k)>=0` for
`k=1..13`, then `E[P(X_S)]<0` proves positive lonely mass.  So a strict
counterexample has to satisfy every one of these moment inequalities.

The default audit over qdiv-hard AP-neighborhood rows found:

```text
rows = 17104
zero_safe_rows = 2
positive_safe_rows = 17102
certified_by_barrier_family = 17102
```

The only zero-safe rows are AP and GW.  Every positive row got a negative
admissible barrier using odd binomial barriers or root barriers of degree at
most 7, except for 50 rows that used the universal `B13` separator.

What feels new here is that covering rows blocked at every denominator-14 apex
still separate globally by multiplicity moments.  The row `12->84` has the
`B13` certificate, and `12->168` has a degree-7 root certificate.  That means
the local geometry can fail while the global integer distribution already knows
there is lonely mass.

Assumption challenge: I considered runners, residues, arcs, endpoints,
components, Fourier modes, multiplicity values, and proof obligations.  The
chosen vertices are multiplicity values and moment barriers.  This preserves
the LRC predicate `Pr[X=0]>0`; it destroys endpoint owners and witness
location.  That destruction is intentional, because the goal is a dual
necessary condition for zero-open packets before applying labels.

Next proof target: prove moment rigidity.  A non-AP/GW row in the qdiv-hard
source core should have a negative admissible barrier, or else it is genuinely
moment-indistinguishable from AP/GW and should be handed to NORK/F7 or
HYP-2908/THM-572.
