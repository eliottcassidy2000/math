# Independent review of the induced ternary repair

**ACCEPTED:** [repairs/TERNARY.md](../repairs/TERNARY.md), with its explicit
scope: one physical source family has both induced ternary children; the
children are not proved to rearm. No Chan--Pak closure follows.

The generic two-birth lemma is correct. The sink a can only be followed by b.
A Hamilton path avoiding b->u therefore ends a,b and contributes H(D).
If it uses b->u, either a precedes b or a is last. Each case corresponds
bijectively to an ordered two-path cover with the designated path starting
at u; the complementary path may be empty. Thus the full count is H+2G,
with no omitted endpoint case.

Both displayed base embeddings preserve every source arc, including the
module star. Uniform substitution therefore preserves the entire old
inflated source, not just its scalar count. The three frames are strong
and pass the two-run matching-star guard. Exact evaluation on TT1 and TT2
then independently recovers the unmarked vectors (25,12), (77,36), (79,36).

The marked response has a particularly direct two-run justification. Source
module port2 has only one exterior out-neighbour, the named root1. That
vertex is already the first vertex of the designated path. On either path
in the cover, once M is entered it cannot be left for an exterior vertex.
Thus there are at most two module runs in total. The checked marked values
on TT1 and TT2 determine G=27K1+12K2, agreeing with the complete finite
compiler. No higher-run response is silently discarded.

The independent [referee_check.py](referee_check.py) uses a separate path
enumeration for G and subset DP for H. It verifies the generic lemma on all
285 rooted labelled tournaments of orders1..4. Literal module controls are:

| M | h | p | Source / +2 / +4 counts | G at source root1 |
|---|---:|---:|---|---:|
| TT1 |1|2|25 / 77 / 79|27|
| TT2 |1|4|49 / 149 / 151|51|
| TT3 |1|8|97 / 293 / 295|99|
| C3 |3|12|147 / 447 / 453|153|
| Q4 |5|32|389 / 1177 / 1187|399|

The C3 and Q4 cases are hostile controls for treating the scaled digits as
fixed digits. They satisfy increments 2h and 4h, as claimed.

The chosen child frames have no root with G-H=2; the alternate finite H77
and H79 heads do have one. These claims were independently recomputed.
Consequently the proof correctly separates scalar attainment, same-count
alternative heads, and a reusable full-response operation. No change to
the current statement is requested.

Reproduction:

    python -B 05-knowledge/results/prime_terminal_20261113_1610/absorption/referee_check.py

Only the reviewer's own receipt is written by that command. No changes were
made in the repairs lane or to the repository's truth/navigation files.

## The H25 bridge obstruction and its three-port chart

**ACCEPTED:** [repairs/OBSTRUCTION.md](../repairs/OBSTRUCTION.md). The physical
inequalities and congruences are valid: for nonsingleton modules,
K2>=2H and K2=2 modulo4; K3 is even and K_j for j>=4 is divisible by4.
The complement pairing argument for P2 modulo4 needs n>=2 and that guard
is explicitly present.

The reviewer independently checked every sign, lower-bound substitution,
and residue claim against the complete saved coefficient matrices: all21
port pairs partition as claimed into11 uniform-sign cases,6 lower-bound
cases, and4 residue cases. The case distinctions correctly handle a
singleton on either side. The finite exterior-word recurrence has no cycle
at a fixed exterior mask because each two-port quotient is acyclic; it is
an exact run-word enumeration. Thus the finite certificate supports an
all-module exclusion, rather than only the literal sample list.

The three-port chart's exact instance was independently recomputed from
the literal adjacency matrices: (h,p)=(25,304), b=2 gives2725 and8179,
with8179=3*2725+4 and prime by trial division through90. Additional b=1
and b=3 controls agree with the universal formula. This operation is
correctly distinguished from a neutral inflation and an induced return;
the two-port exclusion has not been contradicted.

## The root frame packet

**ACCEPTED F1--F4:** [frames/PROOF.md](../frames/PROOF.md). The fixed-comparison
order bound follows by dividing by h and using p/h>=n+1; when the p
coefficient vanishes, oddness forces h=1 for a fixed digit2 or4.
The Bezout example and shifted family were independently recomputed using
the reviewer's subset DP, yielding199/599/601 and the shifted pairs
199->599,415->1247,847->2543.

The all-order unshifted exclusion correctly retains its primitive hypothesis.
The highest run coefficient supplies a positive leading exponential j^m,
so an affine2^m tail forces at most two runs. Every SCC outside the root's
SCC multiplies all run coefficients by its H. Coprime coefficients therefore
force those factors to have H=1, permitting deletion of singleton SCCs.
The remaining strong frame has H113 or115 and order at most9 by the
inherited strong minimum. The stated finite catalogue plus the verified
order-nine extension universe then closes precisely those two profiles.
This does not exclude the shifted324-coefficient chart.

One draft wording correction concerned F5's filtered universe: after
retaining only parents with H<=601, the exhaustive claim is that every
order-nine tournament with H<=601 is represented, not every order-nine
tournament without qualification. The root author corrected the sentence
and the reviewer verified the repaired text. The intended target counts
satisfy that guard, so F4's argument is unaffected.
