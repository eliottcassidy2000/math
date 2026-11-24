# Moore-57: parity, seam descent, and the size of a trade

**Status (2026-11-24): OPEN.** Neither the graph nor an unrestricted sharp row in the surviving Gewirtz model is established. This session builds on the [transport and flag chart](../moore57_switching_flags_20261124/README.md) and [two-seed Gewirtz exact cover](../moore57_gewirtz_exact_cover_20261124/README.md).

| Result | Scope |
|---|---|
| [Parity and fusion frames](moore-tournament-holonomy.md) | A rectangle counts ordinary permutation parity. Pairwise sharpness forces every closed holonomy even. Anchored Gewirtz transports are confined to a subgroup of order 40,320. |
| [One sharp-cover seed](EVEN-SEARCH.md) | Two earlier searches with 21,172 and 20,292 choices reduce to one search with 11,212. |
| [Exact best seam flip](SEAMS.md) | Evaluate and optimize the defect over a whole matching by assignment. Degree seven has positive-energy local minima even with this larger move. |
| [Large anchored trades](ANCHORED-TRADES.md) | Exact anchored completions differ on at least 68 seams. A trade at that bound would remove and replace at least 2940 physical matching edges. |
| [Global energy](moore57-multiedit-energy.md), [cut obstruction](moore-switch-trades.md) | Tangent and curvature identities, at least 167 ordinary edge edits between distinct labeled solutions, and no nonisomorphic trade from one balanced bipartite cut complement. |

A matching edit plays the role of a marked tile flip, but the state stores its full address in two fibers. Parity of transport around a loop survives fiber relabeling; the permutations and the uniqueness of each routed label are lost if only parity is retained. Single-seam descent and simultaneous exact trades are therefore different problems, even on the same metagraph. There is no implication here for tournament Hamilton-path counts.

Run moore-gewirtz-sign.py, moore-gw-mod3-target.py, even-seed.py, and moore_switch_small.py with Python 3 for the finite certificates. The seam-descent.py replay and optional bounded MILP in even-seed.py use NumPy/SciPy.

The normalized rooted matching search itself is also developed in [Less is Moore](https://arxiv.org/abs/2610.04815); the parity, score and trade statements here are about this particular extension, not a solution of the remaining existence question.
