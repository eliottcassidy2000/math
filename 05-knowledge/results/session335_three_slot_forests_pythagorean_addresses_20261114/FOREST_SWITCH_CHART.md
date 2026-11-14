# Two hubs turn one path tree into a three-slot forest

Fix two named hubs $x,y$ in a tournament, with $x\to y$, and an ascending transitive spine. Delete the hubs from a Hamiltonian path. What remains is an increasing linear forest; there are **three ordered slots, not always three nonempty components**. If a hub is an endpoint or the hubs are consecutive, a slot is empty. A fixed order of the hubs and endpoint-labelled slots gives the frame polynomial: every free spine vertex independently chooses one eligible slot, contributing 1, 2 or 3. This explains exactly where the $2^\alpha3^\beta$ terms live.

Here we instead forget the slot labels and fix the bare forest. This changes what the inverse problem sees.

## Glue the bare forest by a tiny permanent

For each nonempty component let $p$ and $q$ be the letters at its minimum and maximum, encoded as $2[x\to v]+[y\to v]$. Write $p_x$ for the corresponding bit, and analogously for $q$. For the hub ordering $x,y$, the column of possible positions of this component is
$$
        (1-q_x,\ p_x(1-q_y),\ p_y)^{\mathsf T}.
$$
The first entry tests its outgoing joint into $x$; the second tests both joints; the third tests its incoming joint from $y$. For $y,x$ exchange hub names. An empty slot is supplied by a dummy column $(1,[x\to y],1)$ using the applicable hub order. If there are $t$ nonempty components, add $3-t$ such columns to make a square matrix. The contribution of this bare forest to the ordinary Hamiltonian count is
$$
        J(F)=\frac{\operatorname{per}(M_{xy})+\operatorname{per}(M_{yx})}{(3-t)!}.
$$
The divisor only removes the fictitious ordering of identical dummy columns. Each path is counted once: the increasing internal order is forced and the designated hubs and slot assignment specify the full path. This is not a matching of arbitrary pairings of endpoints; the same component keeps both of its endpoints. Summing $J(F)$ over partitions of the spine into at most three increasing pieces recovers the ordinary count exactly.

In particular the **number associated to one bare forest** need not be smooth: the sum over gluing orders can equal 5 or 7 even though a locked-frame weight is always $2^\alpha3^\beta$. For example, the word 1030 with the two pairs (1,2) and (3,4) on the 1-based spine has precisely seven gluings. Its **total** is 23, not seven.

### Sharp seven-gluing footprint

A bare forest with exactly seven gluings requires two nonsingleton components, with endpoint types $\{(1,0),(3,0)\}$ or $\{(3,0),(3,1)\}$ in digit notation. In *every* two-hub tournament admitting such a forest, the ordinary path count is at least 23, sharply. Indeed the endpoint-permanent test has only four one-point types and sixteen two-endpoint types; exhausting the $t\le3$ lists with the displayed matrices proves the classification. In particular three nonempty components never contribute 7. Delete all nonendpoints of the two components, retaining their four ends and the hubs. This preserves both endpoint types, hence the seven gluings. Hamiltonian-path count cannot decrease on passing from an induced tournament to its supertournament: insert each new vertex into each old path at a legal break, and deleting that vertex recovers the original path, an injection. Consequently it suffices to check the three pair partitions of four ordered sites and their endpoint labels. The replay checks these literal six-vertex tournaments by independent permutation count; the minimum is 23.

Thus there is an exact example of a local 7 sitting inside the first permitted upper part of the spectrum, and a compulsory **at least sixteen other paths** around it. This is a scoped witness, not another proof that 7 or 21 is absent for general tournaments. It explains a failure mode of trying to promote a frame or a matching witness to an absolute count.

## The gauge conservation law behind the three pieces

Fix a switching class of a labelled tournament on the spine and $h$ named hubs. Prescribe a forward internal path in each of $t$ nonempty disjoint spine pieces covering the spine. Then the forest constraints are satisfiable in every switching class. In a compatible representative, the remaining switches may independently toggle the $t$ pieces and the $h$ isolated hubs, modulo global complementation: exactly $2^{t+h-1}$ charts. This count is independent of the lengths of the pieces. In each ordering of the hubs, inject the $t$ named pieces into the $h+1$ ordered gaps. The seams of this candidate completion form a tree on the $t+h$ contracted units, so the required seam orientations fix exactly one of the remaining switch charts. This counts the full binary slice, even if a switch destroys the unused ascending edges between pieces.

Therefore the sum of the gluing counts over these $2^{t+h-1}$ charts is exactly
$$
      h!\,\frac{(h+1)!}{(h+1-t)!}.
$$
For two hubs and three nonempty pieces it is **12 across 16**, an average of $3/4$. Neither local balance nor all slots being available in principle means the charts have equal weight. In fact no forest with three nonempty pieces can use more than **10 of the 16** charts. At the other extreme all 12 gluings can coincide at one chart (each piece has first letter 3 and last letter 0). The sharp support set is $\{1,2,4,5,6,7,8,9,10\}$ when endpoints are allowed to vary. To certify this, introduce a chart bit at each piece and each hub. Traverse each of the 12 permutations, assigning the next bit to agree with or differ from the previous one according to the required seam; normalize the $x$ bit to zero. The complete $16^3$ endpoint assay is finite, because neither vertex counts nor other arcs affect a seam. It gives, by support size, literal numbers of endpoint assignments
$$
1:16,\quad 2:48,\quad 4:288,\quad 5:384,\quad 6:1680,\quad
7:576,\quad 8:816,\quad 9:192,\quad 10:96.
$$
Those numbers are **endpoint assignments**, not numbers of isomorphism classes or tournaments. A construction with three two-vertex pieces realizes any required endpoint assignment, so the support bound 10 is sharp. Singleton components merely restrict some endpoint assignments and cannot raise the bound. The bit-table does not assert that those 10 charts would all preserve the transitive spine.

The gap between the two charts is the transferable point. Frame weight varies with the number of double- and triple-eligible *sites*. Switch-chart rank varies with the number of *components* and can remain fixed while the spine grows. An inverse proof must transport actual seam words as well as their period or rank. 
