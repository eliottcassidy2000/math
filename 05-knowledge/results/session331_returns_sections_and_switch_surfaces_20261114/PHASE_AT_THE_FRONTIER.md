# Reusable 99-unit outputs at one four-vertex budget

All quotient masks use the last vertex as the position of the source
module. Every quotient listed here is strong, so its uniform
substitution is strong for any nonempty source U. All outputs add
exactly four outside vertices. Let H=H(U), q=Q_U(2), and P_k be
the ordered nonempty path-cover counts.

## Two-way unconditional phase regeneration

The strong quotients with outgoing masks

| quotient masks | coefficient vector in (P1,P2,P3) | output |
|---|---:|---:|
| (6,28,24,17,1) | (11,6,0) | 6q−H |
| (22,12,24,17,2) | (13,6,0) | 6q+H |

have the displayed responses by literal visit-word counting. The
identity P2=q−2H explains the last column. If H is coprime to 99,
both outputs are 3-units. They differ by 2H, so they cannot both
vanish modulo 11. At least one output remains coprime to 99. The
source can therefore be regrown indefinitely, four vertices per
step, with no need to retain a root or a marked phase.

This is unconditional regeneration of the unit property, not
arbitrary target coverage. Since q≥(n+1)H where n=|U|, every step
multiplies H by at least 6n+5. From any finite seed bank the
number of itineraries up to a count X is
exp(O(log X / log log X)): the k-th step costs at least a fixed
constant times k, hence k! up to fixed exponential factors. This
sparse fixed grammar alone cannot yield a cofinite cutoff.

## Open the third slot without losing the unit phase

At the minimal exterior budget four for reading P3, select one
of the five genuinely three-visit quotients below. All are strong.

| name | quotient masks | (c1,c2,c3) |
|---|---|---|
| A | (14,28,8,16,5) | (9,10,2) |
| B | (14,12,24,16,3) | (11,17,4) |
| C | (6,28,8,17,5) | (13,10,2) |
| D | (14,4,24,18,3) | (13,12,2) |
| E | (22,12,24,1,10) | (15,11,2) |

Since P3 is divisible by six and P2=q−2H has mod-3 phase
((-1)^n−2)H, the output modulo 3 is c1 H for odd n and
(c1−c2)H for even n. Thus (B,C,D) are the valid rows for
odd n, and (A,D,E) are the valid rows for even n, when H is
a 3-unit. Their respective determinants are 60 and −16, both
nonzero modulo 11. Consequently, as long as the triple
(H,P2,P3) is not identically zero modulo 11, at least one
row in the relevant triple gives an 11-unit. In particular a
99-unit source always has a 99-unit successor at exactly this
budget **with P3 active**. This construction still does not
prescribe the absolute output.

The regular cyclic five-vertex quotient C5 gives a fourth even
possibility, (15,13,2). Its two complementary readings with B
are

    Y_C5 = 2Q(3)+7Q(2)−5H,
    Y_B  = 4Q(3)+5Q(2)−11H,
    2Y_C5−Y_B = 9Q(2)+H.

In particular one may isolate a lower coordinate combination
without calculating the third coordinate. The C5 quotient is
also exactly the one-missing-crossing context in the sharp
third-slot collision and the even-cycle source for the
physical lattice theorem. These statements identify a shared
five-vertex core, not an identification of the different roles
played by numbers in external geometry.
