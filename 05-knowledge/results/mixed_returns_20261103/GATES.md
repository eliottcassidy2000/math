# Compatible gate clocks can still return signed predecessors

Status: **PROVED exact periodic classifications**, with complete finite
invertible-clock certificates. They concern bulk integrality, not sources.

Use the earlier level-raising blocks and their normalized quotient matrices:

    Q_s78^2=3R49,  R49=[[3395,6688],[1232,2427]],
    Q_D=3R100,     R100=[[238350,336547],[89350,126161]],
    D=L2 b L2.

Their row lattices are respectively a-21b=0 mod49 and a+50b=0 mod100.
These are inherited from [the s78 proof](../adaptive_state_induction_20261103/INCOMING.md)
and [the D proof](../return_channels_20261103/INCOMING.md).

Appending k triangles rotates the exact test: a necessary and sufficient
condition for an INTEGRAL BULK predecessor of 3 beta(v) under s78^2 b^k
or D b^k is

    beta(v) Q_b^(-k) in Z^2 R49 or Z^2 R100,
    Q_b=[[5,8],[2,3]], Q_b^(-1)=[[-3,8],[2,-5]].

These conditions remain insufficient for full-state integrality and physical
positivity. Prefixing additional integer return matrices cannot enlarge the
row lattice of the last factor: Z^2 AB is contained in Z^2 B. Operation order
therefore matters when trying to escape a gate.

For each supplied seed the pair period is112 modulo49 and60 modulo100:

| Input H | Accepted k modulo112 for gate49 | Accepted k modulo60 for gate100 |
|---:|---|---|
|31,95,97|1,57|none|
|287|9,65|none|
|289|25,81|none|
|293 (the supplied repair)|33,89|none|
|295|9,65|9,39|

The full pair clocks are enumerated until their first return; invertibility
ensures there is no transient. This proves the classifications for every k.
For H295 both gates pass exactly at k=9 mod840. The full pair modulo4900
has period1680, with accepted indices9 and849. The acceptance clock loses
one binary phase of the actual pair clock and must not replace it.

Simultaneous acceptance for these TWO candidate words is not integrality
for their concatenation. For genuine mixtures over {b,q2,s78^2,D}, the
bulk valuation increases by exactly the number of s78^2 and D blocks:
their normalized matrices are 3-units, as are b and q2. A return to
3 beta(v) therefore requires that number to be at most1+v3(beta(v)).
In particular a word containing BOTH raisers already fails integrality
for the original H31 and all four supplied repaired endpoints, whose
bulk levels are zero. H95/H97 have level one; the projective theorem
still excludes their exact physical returns.

The two normalized raiser matrices generate the full48-element GL2(F3),
certified by complete finite matrix closure. They do not reduce to a single
commutative multiplicative F9 clock. The actual words U s78^2 D and U D s78^2
have the same added vertex cost31, error gain-4 and bulk-level increase2,
but different normalized phases and H values. Literal graph counts check
both. A known normalized bulk vector still has only eight nonzero states;
48 is the action-group size, not a claim that every observer needs48 states.

At k=9 the two exact bulk predecessors are

    R49 branch:  (259559108622, -715260467245),
    R100 branch: (-8088161056805, 21575973029296).

Both are integral, and both violate the nonnegative physical cone. For digit
j=1 the full state predecessors of 3(295,216,260)+j(2,-1,0) are also integral:

    R49:  (55987722323, -255844540473, 259559108622),
    R100: (-6646159140355, 10066985556423, -8088161056805).

The returned words have vertex costs39 and46, including the nine triangles.
These are formal inverse witnesses, NEVER claimed as actual sources. The
all-depth projective theorem excludes every accepted clock phase at once.

The smaller original-H31 example already passes gate49 at k=1, with bulk
predecessor(-1101,3034). It is the clean failed implication: integrality of
a mixed return is not physical realization. The stronger surviving tool is
the exact gate plus the physical cone plus an actual-source inverse.
