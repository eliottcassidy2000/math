# Three signed branches on the rows of odd numbers

## Two honest Collatz charts

A direct numerical identity is $3q+2=3(q+1)-1$ and $3q+4=3(q+1)+1$. But $q+1$ is even, so the ordinary fixed-sign Collatz rule would **halve** there; this is not a Collatz orbit segment.

The useful chart includes all three tournament digits. Put $x=(H+1)/2$ and write $\Phi_b(H)=3H+2+2b$, $b\in\{-1,0,1\}$. Then
$$
(\Phi_b(H)+1)/2=3x+b .
$$
So the outer digits 0 and 4 are the two signed odd-branch numerators, and the middle +2 is the neutral $3x$ branch. At odd $x$, exactly one of $3x\pm1$ has valuation one and the other has valuation at least two. At even $x$ the two are odd. In particular, on the previous source $H=84\cdot2^s-5$, $x=42\cdot2^s-2$ is even: simultaneous availability of the numerators there is a tournament fact, not enforcement of a Collatz parity guard.

If both signs and dyadic stripping were always accessible, the arithmetic is easy: for any odd $x>1$ choose the sign making $3x+b$ divisible by 4. Its odd part is at most $(3x+1)/4<x$, and its distance from 1 shrinks by at least a factor $3/4$. This does not say a fixed sign converges, nor does it give physical dyadic stripping with the rooted data intact.

## Rows, inverse addresses, and the missing carry

Every positive integer is uniquely $2^a m$ with $m$ odd. In the $H$ coordinate each row is therefore
$$
H_{m,a}=2^{a+1}m-1,\qquad H_{m,a+1}=2H_{m,a}+1.
$$
The dyadic row with $m=1$ contains $1,3,7,15,31,63$; the 21-hole is at the head of the row $m=11$. The **ternary** parent relation is different: the raw $H$ children of 7 are $21,23,25$, those of 21 are $63,65,67$. These are arithmetic parent addresses, not inherited witnesses. In particular 63 belongs to the same dyadic row as 7 but lies two ternary zero steps after it through 21.

On a fixed row $y=2^a m$, the balanced branch label $b$ is the unique representative of $y\pmod3$ in $\{-1,0,1\}$. The address $(x,b)=((2^a m-b)/3,b)$ obeys
$$
(x,b)\mapsto(2x+b,-b)
$$
on going one cell right (except the endpoint where the formal predecessor $x=0$ is not a positive source). This includes the central sheet $b=0$ when $3\mid m$. In count coordinates this is the exact identity
$$
\Phi_{-b}(2H+1+2b)=2\Phi_b(H)+1.
$$
Thus to double a single old branch along its output row, one needs a source $2H-1$, $2H+1$, or $2H+3$ depending on the old sign. One common source lift does not suffice. If instead the input lift is fixed at $2H+1$ and the same branch labels are retained, the three doubled outputs have to be corrected by $+2,0,-2$ to make a new consecutive window. The correction after $a$ common lifts is exactly $2b(1-2^a)$. This gives a signed *target* for a physical repair rather than identifying existing +2 gadgets with an unproved inverse.

## What is physically doubled

In the first ready slot family, $H+5=84\cdot2^s$ and $2E=H+5$. In the reopened family, $H+245=369\cdot2^s$ and $2E=H+245$. They are two different dyadic rows after translation. At their matched handoff $H=493$ the old child has $(E,O,g)=(252,267,15)$, hence $2E-H=11$, while the new parent has $(369,374,5)$, hence $2E-H=245$. The change of continuation is literally a change of dyadic origin. It is not a longer run of the same doubling move.

For comparison, the older fixed-root near-sink block transfer (all nonroots beat the block, and the block beats the root) has $I'=h I$ and $g'=hg$, where $I=H-2O+2g$ and $h$ is the count of the block. Put $\Delta=I+g$. It is zero on the original ready family and positive $du$ on the two nonzero outputs. A near-sink step preserves its sign. For an ordinary ordinal sink block with $q=P(2)\ge2h$, the transfer is $H'=hH,g'=hg,E'=O'=q(E+O)/2$, giving
$$
\Delta'-h\Delta=-qE-(q-2h)O\le-2hE.
$$
Both types of block never decrease $E$ (old even rooted paths can be retained with the block in the complement). For a first nonzero digit $0<\Delta<2E$. Before the first ordinary block, $\Delta$ grows by the product of the $h$'s and $E$ grows by at least the same product, again by retaining old rooted paths. The first ordinary block therefore makes $\Delta<0$, and later ordinary or near-sink blocks keep it negative. Ordinal prefixes only multiply it. So a word of these fixed-root moves cannot reset the old ready equation. This does not rule out a new ready equation, such as the actual 493 re-encoding.

## The charge law and the actual induction obligation

For the sector $O=E+g$, put $k=(H-2E+g)/2$ and $c=g-2k$, so $H=2E-c$. An exact synchronized fresh digit numbered $d=0,1,2$ (scalar increment $2d$) satisfies
$$
(E,g,k,c)\longmapsto (3E,3g,3k+d,3c-2d).
$$
Increasing the length of a marked transitive slot instead doubles $E$ and keeps its family's $g,k,c$ fixed. After $m$ consecutive exact tripled-continuation steps from the old $c=5$ source, the smallest possible $c$ is $5\cdot3^m-2(3^m-1)=3^{m+1}+2$. This follows because the radix-three address ranges from 0 to $3^m-1$. A finite bank of fixed centers cannot accommodate indefinite such steps, even interspersed with their center-preserving slot doublings. The bound stops applying once we change the continuation as at 493.

Conversely, retaining the center alone is not enough. In the [earlier +2 splice](../triple_splice/README.md), fixing its first inflated slot at size one gives $H=9\cdot2^t-1$, $g=3,k=1,c=1$. The +2 setting preserves $c=1$ but sends $(g,k)$ to $(9,4)$. The scalar dyadic row has survived while the readiness can still fail.

There is no unknown Collatz termination issue in the pure forward ternary approach. Each odd $N\ge69$ uniquely reduces to the smaller positive odd $(N-2d)/3$ for the unique $d\in\{0,1,2\}$ making the numerator divisible by 3. That predecessor is at least 23. A **ready** class closed under all three forward digits and seeded at every odd from 23 to 67 would therefore prove coverage from 23 on, with the smaller values handled separately. Starting from an arbitrary scalar witness instead would incorrectly demand $1\to7$. The arithmetic address proves neither physical readiness nor that the signed descent can be run backward through an arbitrary tournament.

## Archive threads used, and what they do not imply

* The October 14 [signed Collatz inverse-fibre note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/signed_collatz_pair_lifts_20261014d.md) gives the alternating-sign lift for the two noncentral sheets. The present three-digit chart adds the central sheet and makes the mismatch of the three physical lifts explicit.
* The [four-vertex signed reading](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md) identifies converse with reversing a local forward/backward reading, not with swapping signs. That is why a converse move is not used as a signed repair here.
* The October 4 [paid carry interfaces](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_carry_interfaces_20261004.md) and [two-clock source addresses](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/source_address_two_clocks_20261004.md) distinguish a surviving address from its native guard. Here the sidecar is the changed continuation; a correct scalar is not a certificate of readiness.
* The October 14 [Mahler reflection and phase restoration](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/recursion_carry_bridges_20261014i.md) likewise restores a bit lost under projection. Its map is not Collatz. The analogy is that the quotient can preserve a current value but lose the information needed for the next move.
* The old [first doubling seam](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-371-first-doubling-unit-pair-seam.md) distinguishes odd-row entry from generic subsequent doubling. The [variable reader](https://github.com/eliottcassidy2000/math/blob/ef469e64524d6aadec6ef90c1458e92d0dd881a6/05-knowledge/results/variable_reader_coverage_20261116/README.md) holds the rooted response through an exact +2 but also shows that a common nontransitive block can scale the unit. The old [triangle-prefix context recurrence](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity.md) keeps signed endpoint defects, not just a multiplicative scalar. The [archival rearm analysis](ARCHIVAL_REARM.md) also examines rerooting and the moving-root tower. The fixed-root blockade above uses this distinction, then the new positive example escapes by switching the rooted state itself.
