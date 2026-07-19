---
id: HYP-7870
title: The seven-wall strict spectrum crosses the four-comb obstruction; the remaining object is a positioned weighted tree plus deletion circuit
status: RESERVED/ACTIVE. THM-1221 proves the global `15/154` Hunter safe-mass floor and therefore a uniform `1/462` incidence margin over the sharp THM-1203 BAD ceiling. The finite protected-needle/beat transport and six-killer coherence are not proved.
source: codex-2026-07-19-S82
depends_on: [THM-1166, THM-1203, THM-1218, THM-1221, HYP-7678]
---

# HYP-7870 -- from global incidence to a positioned witness

This namespace is claimed for the proof interface exposed by THM-1221.  The
proved global inputs are

```text
mu(Safe_7) >= 15/154,
mu(BAD_4)  <= 2/21,
mu(Safe_7 minus BAD_4) >= 1/462.                         (1)
```

For a non-arithmetic-progression quartet, THM-1218 sharpens the last margin to

```text
15/154-60/637=45/14014.                                 (2)
```

Thus the continuum incidence problem is no longer at equality.  The remaining
claim is not another global measure inequality: it is a transport theorem that
places one of the positive-measure phases from (1) inside the particular
eroded slow gap or consecutive beat block supplied by a hypothetical LRC(14)
cover, while retaining compatibility with the other two killer obligations.

The proof-bearing object must keep at least four layers:

```text
weighted pair graph:       rho(si,sj) and the 1/63 strict/weak colors;
position sidecar:          gcd period, tooth address, and interval phase;
deletion circuit:          which four killer labels define BAD_4;
blocking incidence:        the two complementary killers on the surviving gap.
```

A switched runner tournament preserves the threshold cut and a Hamiltonian
path, but destroys the weights, strict/equality color, gcd positioning, and
deletion-circuit membership.  Alternative vertices to test are threshold
components, reduced ratio channels, wall-crossing events, active-pair
handoffs, eroded safe components, and proof obligations.

The immediate quantitative target is one of the following equivalent-style
suppliers, with every error term explicit.

1. A localized Hunter tree on the protected needle whose loss is strictly
   below `1/462` (or `45/14014` off the AP branch).
2. A beat-grid sampling theorem showing that the relevant consecutive block
   meets `Safe_7 minus BAD_4`, with component/boundary error below the same
   budget.
3. A six-killer coherence lemma showing that a phase outside the chosen BAD
   quartet cannot be jointly consumed by the two complementary killers.

This file deliberately does not claim any of those suppliers.  It records the
new fact that their allowable loss is positive and exact, and prevents future
work from reverting to the obsolete zero-margin formulation.
