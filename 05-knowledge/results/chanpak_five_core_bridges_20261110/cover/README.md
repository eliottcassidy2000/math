# Five-core image, selector obstruction, and actual switches

2026-11-10. **Chan–Pak and all-height five-core coverage remain OPEN.**
This lane starts from Session251's exact fixed-core compiler and inverse,
then separates an enlarged finite image calculation from all-parameter
proofs about the permitted constructions. No tournament-order census or
finite-state experiment is promoted to a coverage theorem.

## Exact image and the remaining obligation

The bank is K1, C3, S4, S5a, S5b, with outgoing-row representatives

    (0,), (2,4,1), (2,4,9,3), (2,4,9,19,7), (2,4,25,3,11).

Its Hamilton counts are 1,3,5,9,9. Every core cut and arbitrary transitive
left/right flank addresses A,B are allowed. Using the inherited positive
seven-coordinate form, each cut y gives

    H_y(A,B)=h_y+d_y(g(A)-1+g(B)-1)+L_y AB
              +i_y A+o_y B+s_y(A mod2)+e_y(B mod2),

where g(0)=g(1)=1, g(2q+1)=g(q), and g(2q)=g(q)+q+(q mod2).
The sufficient all-height obligation is exactly: for every positive odd
N except7,21, find a bank cut y and nonnegative A,B with H_y(A,B)=N.
Its finite truth does not yet supply a descent, a source return, or a
uniform all-height existence argument.

**FINITE-EXACT:** the target-bounded inverse covers all such N through
16,777,215=2^24-1. The run evaluates419,808,372 parameter pairs across44
distinct cut signatures after the valid A/B-swap quotient. This is complete
for this source grammar at the stated target bound, including arbitrarily
long flanks after the inherited normalization. It is not a search through
all tournaments. The predeclared first test and adaptive stopping rule are
in [PLAN.md](PLAN.md); [support_20.json](support_20.json) and
[support_24.json](support_24.json) retain the results.

The first three cores alone leave322 extra holes through2^24-1. Removing
only S5a from the full bank leaves exactly

    86641, 313075, 1733389, 7391107;

removing only S5b leaves exactly

    21901, 715183, 1633423, 3274267, 8362471, 12880735.

Thus the two H9 charts remain individually necessary within this particular
bank and finite target universe. Their failure sets are disjoint there.
No conclusion about disjointness at all heights follows.

[verify_support.py](verify_support.py) independently runs the original
Python inverse over all90 raw core cuts through4095 and agrees with all32
per-core support-mask frequencies of the C++ implementation. Both larger
runs also agree with the inherited complete2^18 prefix. This is an
independent small-prefix check, not a second complete2^24 implementation;
[support_checks.json](support_checks.json) states that limitation explicitly.

## All-parameter results

* **PROVED, inherited mechanism extended:** [SELECTOR.md](SELECTOR.md)
  rules out a positive-density rational selector retaining both ordinary
  MSB-first binary flank addresses, for any fixed finite nonempty core
  bank. Pumping a mixed input loop creates positive linear drift in g;
  the bilinear AB term supplies only exponentials and cannot cancel it.
  Accepted inputs consequently have bounded runs, giving only polynomially
  many logarithmic-height outputs after the uniform-axis normalization.
  This restricts a proof architecture; it does not refute coverage.
* **PROVED, scoped actual surgery:** [SWITCH.md](SWITCH.md) changes just
  two apex-core arcs in S5a and gives the exact count difference
  2(B-A+1). It yields adjacent odd pairs on A=B and a count-neutral source
  change on A=B+1. The lower endpoint is not prescribed. Exhausting every
  next-apex cut for the sixteen sources at t=0,...,7 on the neutral line
  shows different response supports. A separate rooted two-bin calculation
  finds1-ready roots at t=0 only and no2-ready roots in this finite family.
* **PROVED specialization of inherited divisibility:** [GAP.md](GAP.md)
  gives the exact seven-coordinate observer for C3 -> TT_r -> C3, retaining
  the internal binary word. For every r>=1 all its apex children are
  divisible by3, with arbitrary exterior flanks. More generally any
  nonempty separator D between the triangles has this property. The
  gapless count41 is the necessary hostile boundary. The mechanism is
  already implicit in Session193's ordered-block capacity formula and
  Session228's triangle-bank divisibility; this packet diagnoses the
  proposed omitted-gap escape. Dividing the formula by3 supplies no
  automatic realization of the quotient.

These compare three distinct losses of information: a finite core label
loses the unbounded digit drift; equal current H loses next-apex response;
and a scalar divisor loses an actual quotient tournament. Any proposed
global construction must retain or reconstruct the missing information.

[AUDIT.md](AUDIT.md) independently checks the other lanes' five-module
classification, polynomial B95 central switch, all-additions readiness
obstruction, and all-order H141-source/425 obstruction. Its separate
standalone frontier implementation agrees on all35,840 labelled cuts.

## Reproduction and finite controls

From this directory, run `python3 -I -B` on each of `check_selector.py`,
`check_switch.py`, `check_reset.py`, `check_gap.py`, and `verify_support.py`.
The corresponding `*_checks.json` files record exact scopes. Checks use
explicit failures, so optimized Python does not disable them. In particular:

|Check|Exact finite universe|
|---|---|
|Selector drift|126 nonempty binary loop words;1134 prefix/suffix cases;2816 two-input signature cases|
|Two-arc switch|All8 retained core contexts;128 literal children with A,B in0,...,7|
|Neutral-source response|Every next-apex cut and every rooted two-bin defect of both sources at each t=0,...,7;256 literal next-apex controls at t=0|
|Internal gap|480 middle-word/cut states for r=1,...,4;192 literal children;1184 arbitrary labelled separator/cut states of orders1,...,3;gapless41 hostile|
|Independent support|All90 raw cuts through4095;all32 support-mask counts;inherited2^18 prefix comparison|

`audit_frontier.py` reproduces the independent cross-lane frontier receipt.
All six checkers have byte-identical normal and optimized outputs, with
source/output hashes retained in [lane_checks.json](lane_checks.json).

To reproduce the large exact image, first run `make_signatures.py`, compile
`support.cpp` using C++17, then run the executable with arguments
`signatures.txt 16777215 support_24.json`. The inherited finite inverse
supplies completeness; the retained code uses signed64-bit arithmetic
with an enforced target bound at most2^30-1. The large image test was stopped
after2^24 in favor of the structural investigations above.
