# Critical seven-slot Hunter degeneracy

Status: `FINITE-EXACT DISCOVERY / ELEMENTARY MECHANISM PROVED / NOT CANON`.

The successful THM-2923 recursion starts after nine speeds are fixed and
uses the slot chain `4 -> 3 -> 2`.  A naive attempt to move the same scalar
Hunter-star mechanism one rung earlier starts from

```text
E in C({1,...,14},6),        seven external slots.
```

For a carrier of mass `h`, exact singleton ranks
`q_1 >= ... >= q_7`, and exact pair-union cap `B_2`, define

```text
G_7(a)=a+sum_(j=2)^7 min(a,q_j,B_2-a).
```

There is an elementary criticality identity.  If

```text
q_7 >= h/7,                 B_2 >= 2h/7,
```

then

```text
G_7(h/7)=h.
```

For `a<h/7`, every summand is at most `a`, so
`G_7(a)<=7a<h`.  Therefore the first hostile level is exactly

```text
lambda=h/7.
```

This is not discrepancy-finite: THM-735 gives only the strict bound
`c(w)<h/7+gamma/w`, which cannot bound the set `c(w)>=h/7`.

The exact probe checks all `C(14,6)=3003` roots.  It dynamically seals the
global top seven and exact global pair cap.  Every root satisfies the two
critical hypotheses and every root has `lambda=h/7`; there are no direct
closures and no positive hostile-level gaps.  The closest inequalities are
recorded in `probe.out`.

Thus the `j=7` seam is structurally different from the completed
seven-body/six-slot rung.  Another blind scalar recursion cannot cross it.
The next certificate must retain information destroyed by `(q_j,B_2)`,
for example:

- compatibility of several pair maximizers or a higher-overlap tree;
- a bounded relation row among the seven tail labels;
- a compactified equality-at-infinity object, with AP/two-scale boundary
  types kept rather than discarded;
- an augmentation/deletion sidecar connecting the newly closed
  seven-body rung to a six-body root.

Reproduction:

```bash
python3 .scratch/lrc_j7_critical_hunter_degeneracy_20260729/probe.py --workers 4
```
