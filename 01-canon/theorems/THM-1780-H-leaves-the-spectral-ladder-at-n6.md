---
id: THM-1780
title: Hamiltonian path counts cease to be adjacency-spectral at order six
status: PROVED for the reconstruction boundary; historical complexity and ladder deductions WITHDRAWN
source: kind-pasteur-2026-07-20-S128c129; corrected 2026-11-02
related: [THM-002, THM-505, THM-1970]
---

# THM-1780 — The binary adjacency-spectral reconstruction boundary

**Current correction, 2026-11-02.** The exact spectral collision survives.
It does not establish a complexity separation, a position above holonomic
functions, or that a Hamiltonian-path count equals the ordinary adjacency
permanent. Those implications in the historical title and body are withdrawn.
This is the same missing-decoder issue recorded in the THM-1970 correction.
The historical version remains in Git at
`291b4d9757842c02aa451beb65acfade2848609b` under this same path.

## Precise statement

For binary tournaments of order n<=5, H is determined by the ordinary
adjacency characteristic polynomial. At every order n>=6, it is not.
The latter assertion means that two tournaments of the same order have the
same characteristic polynomial and different Hamiltonian-path counts.

At n<=5, the [odd-cycle formula, THM-002](THM-002-ocf.md) gives

    H=1+2 c3+2 c5
     =1+(2/3)tr(A^3)+(2/5)tr(A^5).

Missing cycle lengths contribute zero. Binary tournaments have no directed
2-cycles; every closed walk of length3 or5 is a simple directed cycle
up to its starting point. There is no disjoint pair of odd cycles before
order6. Thus the displayed expression proves the positive assertion.

At n=6, native outgoing-neighbor bit rows (zero-based vertices) give

    T13: [18,8,3,5,14,31],     H=13;
    T17: [38,8,2,5,15,30],    H=17.

Both have moments

    (tr A, ..., tr A^6)=(0,0,12,12,10,48).

Newton's identities give equal characteristic polynomials. Independent
integer matrix multiplication and Hamilton-path subset DP reproduce these
values in the [2026-11-02 exact controls](../../05-knowledge/results/stability_priority_spectrum_20261102/carrier_checks.json).
This is a finite exact witness, not a numerical eigenvalue comparison.

Adjoining r sources in a transitive ordinal block multiplies each
characteristic polynomial by x^r and preserves H: every Hamiltonian path
visits that block first in its unique order. This extends the collision
to every n>=6 and proves the full negative assertion.

## What the witness does and does not imply

There is no function of the adjacency spectrum alone that recovers H on
all binary tournaments of any fixed order n>=6. It follows that no
polynomial in adjacency trace moments can do so either.

This does not compare running times, establish #P-hardness, or supply a
generating-function hierarchy. Such claims require their own definitions,
parameters and reductions. In particular the adjacency permanent counts
directed cycle covers, not Hamiltonian paths. For the transitive
two-tournament its permanent is0 while H=1.

The weighted boundary is different. The concurrent [Session139
five-vertex incidence theorem](../../05-knowledge/results/session139_incidence_moduli_and_stability_20261102/riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md)
shows that weighted adjacency zeta first loses H at order5, with the
digon/complement-triangle incidence M23 as a restoring scalar. This
does not contradict the binary result, where digon weights vanish.

## Provenance and reproduction

The original July20 census and code are retained as historical computation:
`04-computation/H_on_the_ladder_kps_S128c129.py`. Its descriptive docstring
does not override this repaired statement.

The new independent control is
`python -B 05-knowledge/results/stability_priority_spectrum_20261102/verify_carriers.py`.
The correction preserves the exact split and adds the source-padding proof;
it removes only deductions not established by that split.
