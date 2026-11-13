# A negative terminal digit can select a different smaller parent

**PROVED scoped forest bounds and an induced-parent obstruction;
INHERITED actual negative-digit and positive-digit families;
FINITE-EXACT controls. Uniform physical source selection remains OPEN.**

Keep the eight seeds and multiplication from [PROOF.md](PROOF.md), and
allow a certified terminal operation3h-2 as well as3h+2 and3h+4. Do not
require its output to re-enter the same source class. This is inspired by
choice of branches in mixed Collatz, but uses no new Collatz descent theorem.

For each unseeded critical target N, the exact residual obligation becomes

    N=2 mod3: D2((N-2)/3),
    N=1 mod3: D4((N-4)/3) OR Dminus2((N+2)/3).          (1)

A domain membership means an actual applicable source and a verified gate,
not scalar attainability alone. Source counts7,21 are forbidden; source1
does not occur at an unseeded gate. The existing eight seeds are unchanged,
including67. The iff proof for the declared grammar is the same atom-last-
operation proof as Theorem B: an unseeded atom cannot first appear as a
product, so its last gate must be one of exactly the choices in (1).

## The total proof budget survives the extra choice

**Theorem.** For any nonunit target N with g such gates and nonunit seed
leaves b_i, counted with multiplicity,

    3^g product_i(b_i-1) <= N-1.

In particular3^g*2^l<=N-1, where l is the number of leaves, and
g<=floor(log_3(N-1)).

**Proof.** Use w(n)=n-1 on nonunits. At products,
w(uv)>=w(u)w(v), since uv-1-(u-1)(v-1)=u+v-2>0. At a gate with
d in{-2,2,4}, w(3h+d)=3w(h)+d+2>=3w(h). Induct over the entire forest.
All leaves are at least3, and no unit source gate is used. QED.

This is a bound on total symbolic gate obligations, not physical graph
size, intermediate counts or selector runtime. The positive-only stronger
unshifted product bound is not silently applied to a negative digit.

The divisor recursion can choose either parent in (1). Through100001 it
reduces the minimum formal gate count on14120 allowed targets compared with
the positive-only model; the largest saving in that universe is six at37987.
For example199 has the one-gate formal option3*67-2 rather than the two-gate
positive-only route. No physical67->199 gate is claimed by that calculation.

## One actual target with two different physical ancestries

Put B=2^m, m>=1. The incoming
[positive construction](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md)
and [Session304 negative chart](../../session304_collatz_peeling_20261113/README.md)
have the *same literal target*

    T_m=W(1001 0^m 110),        H(T_m)=45B+1.

The specified positive source S[TT_m] has count15B-1. Its specified +4
construction changes an old core arc. The negative source

    V_m=W(101 0^m 10),          H(V_m)=15B+1

has an induced two-birth transition to T_m, of count3H(V_m)-2. Thus choosing
the alternative parent count h+2 can change the actual operation resource.
These two families are inherited; their role as alternative terminal parents
is explicit here. No optimality among all possible +4 constructors is asserted.

There is no induced shortcut from a count15B-1 source into V_m. The exact
counts after deleting its vertices are as follows:

| Deleted vertex | Remaining H |
|---|---:|
| first1 |3B+3|
| next0 |9B+1|
| third1 |6B+3|
| any of the m middle zeros |15B/2+1 when m>=2;15 when m=1|
| trailing1 |10B-1|
| last0 |5B-1|
| apex |1|

These follow from the inherited streaming inversion formula: appending a
zero run of length m>=1 to prefix value X and last bit b adds(B-1)X+b
to H. The largest deletion count is exactly10B-1, strictly less than15B-1.
Every proper induced subtournament is contained in a one-vertex deletion,
and induced tournament extension cannot decrease H by the standard insertion
injection. Consequently **no** proper induced subtournament of V_m has
H=15B-1, regardless of its presentation or chosen root.

The scalar difference2 between the two parent counts is therefore not an
induced preparation. The successful move is selecting the other actual
parent directly. This is precisely the latitude allowed by (1), while a
renewable same-witness claim would require more.

Run `python -B 05-knowledge/results/prime_terminal_20261113_1610/arithmetic/alternative.py`
and repeat with `-O`. The arithmetic universe has49998 nonunit allowed targets.
Literal one-apex path-assignment DP verifies the two-ancestry target and every
deletion type for m=1,...,8; the old-vertex embedding is checked directly.
The all-m deletion statement is the displayed formula, not extrapolation.
