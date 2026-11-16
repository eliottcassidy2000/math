# Independent width-four windows with weights 2 and 6

The exact answer depends on what "factor 3" means. If it means a
second choice whose Hamilton-count step is *exactly three times*
the unit step while the first four choices retain their steps
independently, the following 24-vertex example realizes it. If
it means a fully synchronized new seed of local count 3 in the
existing seed-product architecture, that architecture cannot
retain the unit step. Both statements below are exact.

## Literal positive test: a 4×4 rectangle

Take a 15-vertex transitive upper tape (U), an eight-vertex
block (C), and one head (x). All tape vertices beat every
block vertex. In each of four settings (i=0,1,2,3), the block
outgoing masks and mask of head-to-block arcs are as follows:

| (i) | block outgoing row masks | vertices beaten by head (mask) |
|--:|:--|--:|
|0| 0, 5, 17, 7, 43, 15, 63, 127 |73|
|1| 0, 1, 11, 35, 15, 87, 31, 127 |89|
|2| 2, 8, 3, 21, 7, 31, 63, 127 |77|
|3| 0, 1, 11, 19, 39, 79, 31, 127 |105|

Independently choose head-to-tape contact mask (A_j), read
from the tape's dominant end as a 15-bit binary integer:

| (j) | (A_j) | increment in half-count |
|--:|--:|--:|
|0|28085|3481052|
|1|27854|3481055|
|2|28109|3481058|
|3|28157|3481061|

The Hamilton counts are exactly

    H(i,j) = 6,962,193 + 2i + 6j.

Thus the sixteen settings give thirteen consecutive odd
values from 6,962,193 to 6,962,217; duplicate endpoint
values are the expected carries. The low seam changes only
the block and head-to-block arcs. The high choices change
only head-to-tape arcs. These are genuinely independent
edits on disjoint variable arc sets, but generally not
single-arc edits. In the smaller 16-vertex instance, the
same block bank uses just a seven-vertex tape of contact
mask 108 or 111, giving four counts 27,323–27,329 or four
counts 27,329–27,335. The smaller high switch flips only the
last two tape-to-head arcs.

To derive the displayed formula, use the
[one-head seed-monoid identity](https://github.com/eliottcassidy2000/math/blob/7c1a3ae91cd9ea0f7b847521e42a3a23c73f71cd/05-knowledge/results/chanpak_seed_monoid_line_banks_20261115/README.md).
All four blocks share the unilateral signature (h=9),
(a_o=127), (a_e=118), and their cut entries are
(40,41,42,43). The tape contribution to the half-count is
(9N(A)+245 floor(A/2)+127[A odd]).
Here (N(A)) is the isolated one-apex half-count: on
appending a zero to address (B), increase (N) by
((B+(B mod 2))/2); on a one leave it fixed.
Consequently (i) only affects the cut, regardless of
the tape, and each of the selected high moves raises the
tape contribution by three. The accompanying
`check_two_dials.py` verifies all 24 settings by counting
paths split at the head independently of this algebra.

## Why an actual local count-3 block cannot do that by itself

For fully synchronized seeds in the known one-head monoid,
each new block has fixed diagonal count (h_i); the settings
agree on all entries except the weighted cut. Matrix
multiplication makes the change in (H) for digit (i)
equal to (2 Δc_i ∏_{j≠i} h_j). Thus an
ordinal count-three shield rescales an isolated +2 to +6,
never preserves it. For example an ordinal triangle carries
the known (9,11,13,15) arc bank to (27,33,39,45).

A unit low step beside a second active fully synchronized
bank forces that second bank to have (h=1). Such a block
is transitive; its marked outgoing signature at parity
character (+1) is the binary contact address, which
uniquely determines the word up to relabeling and inert
padding. The cut entry is then already fixed. So there is
no such second active bank. In the positive construction
the upper tape has an effective half-step 3 by cancellation,
not by a multiplication by local count (h=3); the
synchronization is unilateral.

## Why expanding the tape into consecutive free blocks does not iterate it

Write (A(u)) for the usual address of one contact word,
(R_0(v)) for the reversed address of the zeros in the
other, and (b_last(u),b_first(v)) for its
end bits. The cross contribution to the one-apex half-count
of (uv) is

    ( A(u) R_0(v) + b_last(u)·(1−b_first(v)) )/2.

Indeed the long even-path kernel between positions at
distance (d) is (2^{d-2}); at distance 1 the unique
adjacent path contributes the extra (1/2). The outer
block’s own tape term is linear in the concatenated mask,
so the full mixed difference when changing both settings
is (h(Δ A,Δ R_0 +
Δ b_last,Δ(1-b_first))).
The second summand has absolute value at most one. At a
fixed width both addresses are injective. Hence any
neighboring active tape bank fails to interact additively
with a four-setting tape bank: some pair of settings in
the latter differ by at least two in address (at least
three if all four distinct), so their mixed difference
is nonzero. The exceptional two-setting cancellation
requires address steps of magnitude exactly one and
opposite boundary correction, which is a literal local
carry. This is a specific obstruction to the next
pure-tape stacking step, not to nonlocal edits or a new
jointly compensated cut.

## A second bank using 68 as a boundary count

The same construction also works with unilateral signature (15,68,53)
and consecutive cut entries 6–9:

| Block rows | Head contact |
|:--|--:|
| 0,5,17,7,43,79,31,127 | 15 |
| 2,8,3,21,39,15,63,127 | 31 |
| 0,33,11,19,7,29,63,127 | 47 |
| 2,8,3,21,39,15,63,127 | 55 |

The tape offsets are 15N(A)+121 floor(A/2)+68[A odd].
The masks 129,135 at width 8 have offsets 8292,8295, giving
at order 17 all seven odds from 16,611 to 16,623. The masks
258408,284439,265290,262716 at width 19 give offsets
17,400,954+3j. This gives the full rectangle at order 28:
H(i,j)=34,801,935+2i+6j. The 68 is an actual parity-class
boundary count, not an inferred significance from 68=4·17.
The verifier checks these counts by direct split enumeration as well.
