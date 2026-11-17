# Preserve the original endpoints with staggered single-arc routes

This improves the [batch bypass](../chanpak_increasing_square_20261117/README.md) for the original labeled square. Keep both strong order-eight modules, the surrounding nine exterior vertices, and the contacts between them from [the first square](../chanpak_commuting_six_20261117/README.md). No batch reversal or vertex insertion is needed. Choose *different orders* for the five exterior arc reversals in the two columns:

    lower: (1,7), (1,8), (3,4), (6,8), (4,6)
    upper: (1,7), (6,8), (3,4), (1,8), (4,6).

The shared spanning cycle is kept in every intermediate graph. The lower module has H=l, the upper H=u>l, and both have P₂=P. Each stage reads alpha*P+beta*H. The physical coefficients are determined by replacing the module with TT1 and TT2 and counting paths independently; the fixed single exit limits the module to two runs.

| step | lower read | upper read | lower actual | upper actual | arc distance across the row |
|---|---|---|---:|---:|---:|
| 0 | 128P+l | 128P+u | 6103423 | 6103429 | 1 |
| 1 | 256P+5l | 256P+5u | 12207995 | 12208025 | 1 |
| 2 | 276P−9l | 292P−7u | 13156233 | 13919837 | 3 |
| 3 | 475P+109l | 599P+103u | 22689747 | 28600387 | 3 |
| 4 | 823P+13l | 823P+13u | 39245619 | 39245697 | 1 |
| 5 | 896P+7l | 896P+7u | 42723961 | 42724003 | 1 |

The comparison actually holds for every common P≥2u and l<u. In both columns alpha increases and 2Δalpha+Δbeta>0, hence vertical counts increase. At the two staggered rows the horizontal differences are 16P−7u+9l and 124P+103u−109l, both positive by P≥2u; all other comparisons use positive beta with the same exterior. Thus the original four labeled corners are kept and numerical crossing is avoided.

The trade is explicit: the horizontal endpoints at the two middle rows are no longer connected by the original single core reversal. They differ by that reversal and two exterior reversals. On the concrete pair, the independent enumeration of all 120×120 orders finds 322 valid ordered pairs of columns, all differing at two or more intermediate rows. Thus two staggered rows are optimal within these five edits. This does not contradict the earlier proof that a common exterior cannot move both tracks monotonically without a crossing.

Even allowing longer horizontal paths does not remove this trade within the stated setting. The exact search tests all 322 numerical schedules, keeping the fixed module contacts and allowing every exterior arc and the one core arc in its horizontal detours. A rung must stay strong and have strictly increasing counts between its two endpoints. Shared-exterior rungs already work as single arcs. Over the candidate schedules the search needs only five distinct rung checks, since each schedule can be discarded at its first impossible rung. Each negative bidirectional search exhausts a frontier before either guard bound; none of the 322 schedules has every rung. This permits horizontal detours outside the five-edit cube, but not changes to the fixed module contacts or alternate vertical moves.

Replay with `python replay.py`; for the last check run `g++ -std=c++20 -O3 -Wall -Wextra -Werror rungs.cpp -o rungs && ./rungs`. The result is scoped to these original labeled endpoints and allowed vertical schedules, not a global obstruction to other sevenfold constructions.
