# Changing the source: thirteen-count windows at every height

The fixed-port obstruction does not prevent a bank of **different** chronological
sources from sharing a bulk term. From the [fan law](FAN_CHANNELS.md), if
$S_q-K_q=2^e\lambda>0$ and the frozen defect is $\eta$, attach $j-e$
leaves at $q$. For every $j\ge e+1$ this realizes
$H=\lambda 2^j-\eta$. The seed is selected by the desired absolute offset,
rather than deduplicated by its old Hamiltonian count.

## Two physical windows

| Common odd slope | Offsets (inclusive, by twos) | Every exponent |
|---:|---|---:|
| 1 | 1, 3, …, 25 | $j\ge10$ |
| 17 | 23, 25, …, 43 | $j\ge6$ |

Each arm comes from a literal one-parent history on at most nine vertices.
The complete certificates are in [the replay](replay_parallel_fans.py).
It reconstructs each tournament from its parent list, computes its exact
$S,d,K$ by subset DP, and checks a literal first leaf as well as the
slope and offset. The two windows respectively have thirteen and eleven
consecutive odd totals, at every listed height. For example, the first
window contains all odds from 1025 to 1049 at $j=10$.
Transitive prefix padding puts all arms of a given window at common order
$j+9$ without changing their counts; these padded objects need not be
strong or chronological.

These are exact absolute windows, not merely residue covers. But the inverse
is still restricted: to select a source one must first solve
$N=\lambda 2^j+c$ with an offset actually in the bank. A finite bank with
finite offset sets yields only $O(\log X)$ distinct totals below $X$.
The structural search question is whether a source-changing rule can
grow the width proportionally to the gap between centers, or supply enough
different centers with overlapping windows. More heights at fixed width
cannot do so.
