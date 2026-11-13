# A sharp obstruction to reusable +2 moves at a single disagreeing contact

Consider tournaments which differ only in the orientation of $pq$.
Assume the two vertices $p,q$ have the same orientation to every other
vertex except $a$, where $p\to a\to q$. Let $R$ be the tournament
with $p,q,a$ deleted, and write $P,Q\subseteq V(R)$ for the
out-neighborhoods in $R$ of $a$ and of either $p,q$, respectively.
Use $H(\varnothing)=1$.

**Theorem.** Replacing $p\to q$ by $q\to p$ increases the Hamilton
path count by at least $2H(R)$. Equality at increment two is possible
at exactly the absolute ordered count pairs
$$
 (1,3),\quad (3,5),\quad (9,11),\quad (13,15).
$$
Consequently the minimum positive count increment in this near-twin
grammar is not a repeatable +2 repair at large targets, regardless of
the size or internal complexity of the reservoir $R$.

## Proof: the contact is absent in one contraction and doubled in the other

A path that does not use $pq$ as consecutive vertices survives
unchanged. Contract a path that does use it to a ghost $z$: incoming
arcs to the ghost go to the first vertex and outgoing arcs leave the
second. The ghost has the same arcs to $R$ in both contractions.
With $p\to q$ neither ghost contact to $a$ is present; with $q\to p$
both are present. The added paths use exactly one of these contacts.
Contracting $z,a$ once more, in each order, gives
$$
 \Delta=D_R(P,Q)+D_R(Q,P),
$$
where $D_R(P,Q)$ counts paths on $R$ and a ghost whose allowed
out-neighbors are $P$ and allowed in-neighbors are $V(R)\setminus Q$.
The two flags are independent; this intermediate is not necessarily a
tournament.

For a directed Hamilton path $r_1,\dots,r_m$ of $R$, put
$p_i=[r_i\in P]$ and $q_i=[r_i\in Q]$.
Across both ghost conventions, the number of legal insertions is
$$
  p_1+q_1+2-p_m-q_m+
  \sum_{i=1}^{m-1}\big((1-q_i)p_{i+1}+(1-p_i)q_{i+1}\big)
 =2+\sum_{i=1}^{m-1}\big(q_i(1-p_{i+1})+p_i(1-q_{i+1})\big).
$$
This is a telescoping lower bound of two **per path**, not per
isomorphism class of endpoints. Deleting the ghost recovers the
underlying path and insertion slot. Other contributions are positive:
they are the ordered two-path covers of $R$ whose seam would not be
an arc if the ghost were deleted, but whose named ends fit one of its
two interface conventions. Thus the square-sum distinction between a
missing seam and a valid detour has an exact count analogue. If $R$
is empty both ghost terms are 1.

## Classifying the sharp case

If $\Delta=2$, the lower bound forces $H(R)=1$. Then $R$ is
transitive: any directed triangle would give three distinct paths
that extend by insertion. Write $R=(1,\dots,m)$ in its transitive
order. Remove any leading positions with $(p_i,q_i)=(0,0)$
and trailing positions with $(p_i,q_i)=(1,1)$. These are universal
sources or sinks in both original tournaments, so deletion preserves
the ordered counts.

In the reduced nonempty sequence, the first pair has at least one 1,
and the last has at least one 0. The start/end contributions to the two
ghost counts already sum to at least 2. They must equal 2, so both
endpoint pairs are mixed, and all interior ghost contributions vanish.
The coefficient for a left path ending at $i$ and a right path
starting at $j$, ignoring the ghost interface, is
$$
 w_{ij}=
 \begin{cases}
  1&j=i+1,\\
  2^{i-j-1}&j<i,\\
  0&\text{otherwise}.
 \end{cases}
$$
Indeed when $j<i$, each label strictly between the endpoints can be
assigned independently to either path, and all others are forced.
For $m\ge2$ the mixed endpoint pairs cannot coincide, since
$w_{m1}>0$. Up to swapping the names of the bit strings they are
$10$ first and $01$ last. The first 1 in the $P$ string requires
$q_i=1$ for every $i>1$, lest $w_{i1}>0$ contribute. The last
0 in $P$ requires $q_j=0$ for every $j<m$, lest $w_{mj}>0$
contribute. An interior index cannot satisfy both; hence $m\le2$.
For $m=2$ the endpoint sequence is $10,01$ or $01,10$;
for $m=1$ the sole pair is either mixed pattern. The empty case is
also allowed. Their exact count pairs, respectively, are $(9,11)$
or $(13,15)$, $(3,5)$, and $(1,3)$.
These tiny evaluations are reproduced by the companion script;
the same script verifies the exact mixed seam formula on all
transitive reservoirs through order five.

The obstruction suggests a sharp design rule for the next construction:
**a large-target constant repair cannot be a one-sided contact gain.**
Use contacts of both signs, or change the reservoir, so the two
ordered seam channels can actually cancel. The theorem makes no
assertion that such a signed repair always exists.
