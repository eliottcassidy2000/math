# The first exposed port can distinguish identical full path-cover profiles

Write $P=A\to TT_1$ and $B$ for the two eight-vertex blocks of [Session 310](../session310_positive_clocks_20261113/NEUTRAL_TAILS.md), so $Q_P(t)=Q_B(t)$. Uniform exterior tests cannot distinguish these blocks. Here a *marked* test is allowed to orient its edges separately toward one named vertex $v$, but has uniform incidences toward all of $T-v$.

## 1. A two-vertex mixed-difference detector

Add exterior vertices $x,y$ with $x\to y$ and $x\to T-v\to y$. Let $h_{ab}$ be the Hamilton path count when $a=1$ selects $v\to x$ (instead of $x\to v$), and $b=1$ selects $y\to v$ (instead of $v\to y$). Write $F=V(T)-v$, and let $H(\varnothing)=1$. Then
```
h11 − h10 − h01 + h00 = 2 Σ_{S⊆F, |S| even} H(T[S]) H(T[F−S]).
```
In particular if $|T|$ is even, the right side is exactly $Q_{T-v}(2)$.

**Proof.** Use [the odd-cycle collection formula](../../../01-canon/theorems/THM-002-ocf.md), assigning weight 2 to every packed odd cycle. In every setting, a cycle through $x$ must use $v\to x$ (the only possible incoming arc at $x$), and a cycle through $y$ must use $y\to v$ (its only possible outgoing arc). Consequently an odd-cycle packing using both arcs must put them in the *same* cycle, since the two arcs touch $v$. Packings using neither or just one vanish in the mixed difference; the opposite orientations of the two arcs are never in a cycle in the states where they occur. The surviving cycle has the unique form $v\to x\to$(directed spanning path of $T[S]$)$\to y\to v$, with $S$ even because the cycle has $3+|S|$ vertices. The direct edge $x\to y$ handles $S$ empty. Every remaining packing lies in $F-S$, giving its Hamilton count by the same formula. Finally for odd $|F|$ complementation exchanges the even and odd sets with equal summands. Thus doubling the even sum recovers the two-slot partition sum. □

For the Session 310 pair the rooted deletion values in label order are:
```
P: 752,416,608,608,416,752,672,1276
B: 868,128,868,868,612,868,612,868
```
These sets are disjoint: no assignment of one named vertex in each core can preserve all two-exterior marked tests. This is not a statement about general full rooted tomography. For example, the above formula gives only one parity-restricted sum when $|T|$ is odd.

## 2. For this pair, one exterior vertex is enough

Use the same single test everywhere: a fresh vertex $x$ receives from all of $T-v$, and sends only to $v$. Direct subset dynamic programming gives
```
P: 297,325,405,405,409,741,445,95
B: 145,191,215,299,609,577,673,413.
```
Again the sets are disjoint. This is a sharp exterior-size boundary for this pair because the test of exterior size zero is just $H(P)=H(B)=95$. Not every selective test works: using the opposite cut, a fresh vertex that loses only to the mark produces the same 239 in $P$ at mark 2 or 3 and in $B$ at mark 5. So agreement at a favored marked count is not a proof of a reusable selective interface.

## 3. Where marked composability *does* survive recursively

Let $C$ be any quotient with a distinguished slot $u$; mark any number of vertices outside the inflated slot. Replace the uniform module $P$ by $B$. In any further context (allowing arbitrary arcs to the marked vertices) that still sees the whole module at $u$ uniformly, every path-cover count remains equal: collapse consecutive visits to the module in each labelled path slot; the ordered list of all visits is filled by the same ordered cover number for $P$ and $B$. Thus marks can be kept through arbitrarily many layers of uniform wrappers provided no future test uses a mark inside the altered module. This is an actual recursively stable *external* port, not selective substitutability at an arbitrary internal port. For instance repeatedly wrap in a directed triangle of blocks with two singleton blocks, mark any of the wrapper singleton vertices, and retain all exterior labels; all wrapped tournaments are strong at every positive depth. The single-probe separation above is the warning against moving such a port into the hidden core without a new certificate.

All finite lists and the general mixed-difference identity on representative cases have a standard-library replay in `verify_ports.py`.
