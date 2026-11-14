# A nonterrace source through the signed cofactor switch

Combine the [arbitrary-root preparation](ARBITRARY_ROOTED_BRIDGE.md)
with the [two-context corridor switch](SOURCE_SWITCH.md). The latter
only requires a ready two-state signature; it does not need the source
itself to have a terrace drawing. Thus swapping the two final small
corridors gives an exact cofactor difference of \(2a\) (up to sign)
for the source \(a=H(T-v)\), including sources physically outside
the old cone's terminal dyadic bands.

For a small literal example, start with the four-vertex tournament \(C\)
with outgoing rows \((8,1,3,6)\). It has \(H=5,F_2=32\).
Prepend a dominating root and use a width-one bridge. This produces
a ready **seven-vertex** tournament with \((a,b)=(5,64)\);
the complete literal outgoing-mask graph is
\[
                         (62,112,98,102,108,64,1).
\]
It is impossible to obtain those coordinates on the terrace:
\(b/a=64/5\) would need width power 8, but
\(64/(5\cdot8)=8/5>4/3\).

Use \(q=14,r=1\) followed by the \(q=1,2\) and \(q=2,1\)
switch, then the \(q=7,r=2\) divisor. The two cofactors are
\[
                         3539307127,\quad 3539307117.
\]
The first is prime by trial division and the second is divisible
by 3. The actual ordinary Hamiltonian counts are seven times
these numbers, so the repair changes the realized total by 70
without changing the physical starting graph or the multipliers of
its second coordinate. The starting tournament is neither the old
fixed same-reset pair nor a change of the width of its initial gap.

The same code can use any source from the target compiler. To get
a fixed unit correction in the cofactor, this particular repair
requires \(a=1\); arbitrary cores deliberately give a varying
correction size. This example does not choose an arbitrary target
baseline.
