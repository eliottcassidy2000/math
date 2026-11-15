# A coprime $5,9$ return that retains both interfaces

The [count-five switch](CRT_WINDOWS.md#a-literal-unit-switch-at-the-next-multiplier-five)
has two-sided profile $(h,Q,s,e,t,f)=(5,1024,209,5,303,5)$ and local
path counts $67,69$. Here is a second two-sided switch with count-nine
spines and local counts $207,211$:

| Outgoing binary row masks of spine | Mark mask | Local count |
|---|---:|---:|
| (254,220,248,240,224,194,128,0) | 54 | 207 |
| (254,252,184,240,224,192,132,0) | 74 | 211 |

Both splice the same count-nine core with outgoing rows
$(14,28,24,16,1)$ into a three-singleton transitive buffer. In the first,
the buffer and core marks read $0\to(27)\to00$; in the second, $01\to(18)\to0$.
The marked spines have exactly the same full profile
$(9,672,236,0,545,9)$. Thus exchanging the core location and the cut
changes the count by four with no loss of the incoming or outgoing
interface.

Put a count-five spine first, a count-nine spine last, with one common
apex. The base plus the two switches generates exactly

$$
115543+18 a+20b,\qquad a,b\in\{0,1\},
$$

with common profile $(45,688128,141628,0,560807,45)$. In particular,
the paint pairs $(a,b)=(1,0),(0,1)$ give **115561 and 115563**. The
differences are literal: $5\cdot 4-9\cdot2=2$. They are not inferred
from a gcd alone. Both switches pertain to the full exterior signature,
so any common transitive left **and** right tails preserve this exact
two-count window. Prepending a marked singleton and appending an unmarked
singleton makes both whole tournaments strong and yields the same pair
shifted to **1506169 and 1506171**. Two common cyclic contexts are also
allowed, but then the difference scales by both ordinary path counts.

The [bounded modular reset](PRIMITIVE_LANES.md) gives a limit to repeating
this particular coprime pair indefinitely: $3$ divides both the count-nine
spine's ordinary count and its two-cover count, so three occurrences kill
the primitive target lanes even if they sit at different positions and
their marks and outer baselines are changed. Thus the example demonstrates
one exact, fully matched carry. It does not supply an indefinitely
reusable primitive return.

## A scale-stable width criterion

In the [general carried-window decoder](CRT_WINDOWS.md), take pairwise
coprime spine counts $h_i$, $P=\prod h_i$, $w_i=P/h_i$, and contiguous
palettes with half-widths $M_i$ (including $M_i+1$ digits). Even when the
palette sizes are not multiples of their moduli, the same argument gives
the sharp central half-offset interval
$$
[R-P+1,\ M-R+P-1],\qquad R=\sum_i w_i(h_i-1),\quad M=\sum_i w_iM_i,
$$
provided $M_i\ge h_i-1$: the lowest lift in each residue is at most $R$
and the highest at least $M-R$. The two exterior integers fail by taking
all residues $h_i-1$ or by reflection. In particular, if $M_i\ge2h_i-2$
for all $i$, the interval has half-width
$$
2P-2+\sum_i w_i(M_i-2h_i+2).
$$
The *same type of threshold* is thus preserved by gluing: two sufficiently
wide, coprime matched palette modules yield one more sufficiently wide
matched palette module, with surplus amplified. The practical proof
obligation is to build source families that reach this scale or supply
specific quotient fibers for the target. This exact interval controls
carry, not absolute placement at an as yet unknown baseline.

Replay: compile and run `replay_coprime_switch.cpp` as C++17. It computes
each local profile and the whole unguarded and guarded Hamilton counts by
direct subset DP, without using the seam update.
