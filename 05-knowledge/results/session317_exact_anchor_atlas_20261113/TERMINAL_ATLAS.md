# Exact terminal atlas for arbitrary two-apex endings

Let $a\to b$ be two apexes and let the other vertices form an increasing
chain. Each letter records whether $a,b$ beat that vertex, in bits 0,1.
Take any prefix $P$ in letters 0,1 and **any** suffix $u$ in letters
0,1,2,3. Put $W=T(b,u)$, the tournament on $b$ and the suffix alone.
All prefix vertices beat all of $W$. The out-cut $U$ of $a$ in $W$
contains $b$. Use these four physical terminal parameters:
- $h=H(W)$;
- $s$, the full paths of $W$ starting in $U$;
- $E$, the ordered two-slot covers with one distinguished nonempty
  path ending outside $U$ and the other slot unrestricted;
- $c=H(a+W)$ with the actual cut.

If $B$ is the binary integer encoded by $P$, $\ell=B\bmod2$, and
$K(0)=0,\ K(2B)=K(B),\ K(2B+1)=K(B)+B+(B\bmod2)$, then
$$
\boxed{H(Pu)=c+(h+E)B-hK(B)-s\ell.}
$$
Leading zeroes are universal sources. For $B>0$ use the expansion
beginning in 1. The result is strong precisely when every vertex of $W$
can reach $a$ in $a+W$: from $a$ we reach the first prefix vertex,
hence the chain and $W$; a return to the prefix must go through $a$.

**Proof.** Project the
[common-apex block law](../session313_zeckendorf_target_selection_20261113/ADDITIVE_PORTS.md)
to two ordinal blocks $I,W$. The transitive prefix $I$ has full count
1, start response in two slots $B$, and end response in one slot
$1-\ell$. If $s_I$ is its one-slot start count, its two nonempty
split count around $a$ is $B-s_I-K(B)$: from all two-slot starts
subtract the empty other slot and the splits whose other endpoint
is also in the out-cut. The latter count is $K(B)$. When a 0 is
appended, it can only go in the slot with marked start, hence the
count is unchanged. When a 1 is appended, placing it at the end of
the other slot contributes $B$ and placing it as the entire
start-marked slot contributes $B\bmod2$; the other arrangements
give the old count. Thus the empty-left and same-$I$ terms sum
to $h(B-K)$. The empty-right and same-$W$ terms, with the forward
adjacent contribution, give $c-s\ell$. The backward adjacent term
is $EB$. The empty prefix is immediate.

In particular the curvature coefficient $-h$ can **never** vanish.
Equivalently, flip the incidences at $a$ of the first two chain
vertices independently. The mixed difference comes only from paths
traversing $a$ between those two vertices. The two contractions are
a universal source followed by $W$, so the mixed difference is
exactly $-2h$. No straight affine decoder is possible in this
binary-prefix family, even with arbitrary suffix length.

## An unbounded absolute triple corridor

The following suffixes have identical $(h,E,s)=(9,245,0)$:

| suffix | $c$ |
|---|---:|
| 021102 | 91 |
| 0102330 | 93 |
| 0102101 | 95 |

The three strong outputs for each $B\ge1$ are exactly
$$
 254B-9K(B)+91,\quad 254B-9K(B)+93,\quad 254B-9K(B)+95.
$$
They are consecutive odd targets at unbounded scale, with one
unchanged prefix. The standalone replay computes physical terminal
parameters by subset DP and checks strongness with a prefix of one
letter; the reachability proof then applies to every nonzero prefix.
Another short pair, suffixes 120 and 312, has counts
$13+15B-3K(B)-3\ell$ and $15+15B-3K(B)-3\ell$.
This does not assert overlapping triples across different prefixes.

## Why a fixed-count reservoir cannot provide a wide offset table

The exact terminal anchor obeys
$$
 h\le c\le E+s\le E+h.
$$
The lower bound is by insertion of $a$. If $e$ counts full paths
ending outside $U$ and $D$ the splits around $a$, then $c=s+e+D$
and $e+D\le E$ by forgetting a constraint on the second slot.
So the offset is a smaller exact search at scale $E$, not an
arbitrary residue.

**For fixed $h$ and $E\ge1$, the number of possible anchors $c$ is
polylogarithmic in $E$**, uniformly over order, even allowing arbitrary
$W$. Order the strong components $W_i$. The nontrivial ones have at
most $O(\log(h+1))$ total vertices: their counts multiply to $h$
and [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
bounds their size logarithmically in those counts. Thus at fixed
$h$ there are finitely many lists of such components and cuts on them.

Take the last component containing an in-cut vertex $d$. Any subset
of vertices of earlier components can be assigned to the first slot
and followed by $d$, with the rest assigned to the second. Hence at
most $\log_2 E$ vertices lie in those earlier components. Singleton
components later than the last in-cut component can be deleted:
they are ordinal singletons in $a+W$, and must be in the second slot
in an end-marked cover. Deletion preserves $h,E,s,c$. After trimming,
the possible singleton-run lengths are bounded by a power of
$\log(E+2)$ at fixed $h$.

For each choice of run lengths and nontrivial components/cuts,
there is at most one placement of the cut on singleton vertices at
the specified $E$. Indeed
$$
 E=\sum_i E_i\prod_{j<i}Q_{W_j}(2)\prod_{j>i}H(W_j).
$$
The local $E_i$ marks the end in component $i$; earlier components
can use both slots and later ones only the other slot. A singleton
contributes its in-cut bit times the displayed weight. Consecutive
singleton weights grow by at least a factor two, because the earlier
singleton gives two and intervening factors have $Q(2)/H\ge2$.
Each weight exceeds the sum of its predecessors, giving unique
binary decoding. This rules out linearly wide anchor tables at fixed
reservoir count by stretching neutral runs. It leaves changing $h$
and changing interface types open.
