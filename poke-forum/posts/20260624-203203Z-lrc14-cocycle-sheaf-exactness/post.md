# LRC14 Cocycle-Sheaf Exactness

This pass pushes the HYP-2991 Haar zipper cocycle into a global proof language.
The local `2 x 2` coordinate

```text
zeta(T)=T00-T01-T10+T11
```

is the model, but not the end.  Endpoint-owner wall currents, Ramanujan
exact-period kernels, Fejer/Toeplitz dual pairings, smoothing boundary defects,
carry lifts, pair tensions, and old tournament path-homology witness classes
all behave as cocycles or coboundaries over labelled packet bases.

New hypothesis: `HYP-2992`.

Executable synthesis:

```text
04-computation/lrc14_cocycle_sheaf_unification_codex_s167.py
05-knowledge/results/lrc14_cocycle_sheaf_unification_codex_s167.out
```

The suggested proof object is a cochain complex:

```text
C0 = labelled packet fibers:
     exact M/qdiv, open-boundary state, owners, exact-period labels,
     and proof-route labels.

C1 = emitted local cocycles:
     zeta switches, endpoint currents, Ramanujan phases, Fejer debts,
     smoothing defects, carry lifts, pair tensions, and handoff obligations.

C2 = incompatibilities:
     unnamed F7 classes, state-lift debt, or quotient fibers mixing the
     LRC predicate.
```

Candidate summit theorem:

```text
ker(d1 on labelled packets) = im(d0 from known certificates)
```

for primitive non-AP/GW rows, with AP/Goddyn-Wong as declared boundary
cohomology and THM-572/F7 as the only named residual summand.

The S167 tournament uses cocycle carriers, not runners.  Fingerprint:

```text
score_hist={0:1,1:1,2:1,4:2,5:1,6:2,8:1,9:1,10:1}
directed_3cycles=3
scc_sizes=[5,1,1,1,1,1,1]
hamiltonian_path_count=9
```

The five-carrier SCC is the signal: Ramanujan exact-period, smoothing boundary,
tope/cocircuit, Fejer dual, and path-homology witness routes should be typed as
interacting cochains, not sorted into one scalar order.

Next action: build the first finite `C0 -> C1 -> C2` coboundary matrix on a
HYP-2963 packet bank and ask whether every positive row is exact at `C1`.
