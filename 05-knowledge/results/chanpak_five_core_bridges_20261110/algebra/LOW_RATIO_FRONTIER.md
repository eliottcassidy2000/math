# A complete finite source-preserving test below the two-ear threshold

**PROVED elementary consequence of the strong-ear argument in this packet.**
This is an exact test for extensions of one specified source. It neither
selects a source for every integer nor proves Chan–Pak. The H141 instance
and its independent computation are in
[BIRTHS](../alternatives/BIRTHS.md) and [the frontier audit](../cover/AUDIT.md).
No priority claim is made for the classical relative strong-extension lemma.

## 1. General theorem

Fix an actual strong tournament S with n>=3 vertices and h=H(S). Let N be
an integer satisfying

    h <= N < (n+2)h.                                    (1)

Start at S. Repeatedly adjoin one vertex with every possible incidence
pattern, retaining precisely the children which are strong and have H<=N.
Keep the old named source as an induced subtournament. This finite search
produces, up to naming the added vertices, **every strong supertournament
of S with H<=N**. It terminates after at most N-n added vertices. In
particular an empty frontier is a proof of exhaustion at every graph order,
not a choice to stop an order-bounded search.

Proof: suppose R is a proper strong induced subtournament of a strong U.
If some outside vertex has both an incoming and an outgoing arc to R,
adjoining it is a strong one-vertex enlargement. If there is no such vertex,
outside vertices split into those dominated by all R and those dominating
all R. Strongness forces an arc from the former class to the latter. Those
two vertices and R form the cyclic two-ear R->x->y->R.

That two-ear has H(R)+Q_R(2) paths, where Q_R(2) counts two ordered path
bins, including empty bins. Cutting each Hamilton path at all |R|+1
positions injects (path,cut) into these covers, so the two-ear has at least

    (|R|+2) H(R) >= (n+2)h > N.                         (2)

Hamilton insertion is monotone under induced inclusion: an old Hamilton
word always has an insertion position for a new vertex, and deleting that
vertex recovers the word. Thus a two-ear inside a target U of count<=N is
impossible. There must be a strong one-vertex enlargement inside U. Repeat
until U is reached; every intermediate is one of the retained children.

For termination, a strong tournament of order m has at least m Hamilton
paths: from each prescribed first vertex, repeatedly insert an outside
vertex after the last current vertex beating it. Strongness guarantees an
outgoing edge from the current set until it spans. This constructs a path
with that first vertex, and the m starting vertices distinguish the paths.
Consequently H(U)<=N implies |U|<=N. The algorithm's stated bound follows.
Sharper Moon–Busch order bounds may improve runtime but are not necessary
for correctness. No assumption about the conjectured spectrum was used.

The strict threshold is sharp. For S=C3, n=h=3 and Q_S(2)=12. Its cyclic
two-ear has H15=(n+2)h, but each of the two added vertices alone is a
uniform source or sink relative to S. No strong one-vertex chain reaches
this particular target. Replacing the strict inequality in(1) by a
nonstrict one would therefore make the stated completeness theorem false.

## 2. The ternary implication and the nonstrong boundary

For every h>=3 and n>=3, both targets 3h+2 and3h+4 satisfy(1). Therefore
the actual-source question

    Can any number of births, retaining S, reach a strong target3h+2 or3h+4?

is exactly settled by this finite one-vertex strong-frontier search.
This does not require each birth in the originally proposed program to
have strong intermediates: the theorem supplies another ordering of the
target's extra vertices which does. It also does not replace a multi-birth
extension by a single birth, as the explicit141->215->425 trace illustrates.

For an arbitrary, possibly nonstrong target of count N containing S, the
strong component C containing S has H(C)>=h and H(C) divides N. Hence if

    every proper positive divisor of N is smaller than h,                 (3)

then H(C)=N, and arbitrary-order target existence is equivalent to the
strong target test above. Any other components are neutral singletons.
This includes prime targets, and the present425/427 instances at h141.
Without(3), a containing component of smaller divisor count may be joined
to other components; that factorization and its physical witnesses must
be retained. Dropping(3) would be an incorrect all-target inference.

## 3. What the new source obstruction establishes

The two C4 sources at k2 both have n8,h141. The threshold is1410, so
N427 is well inside(1). Their exact strong frontiers have layer sizes
25,4,0 and28,4,0. All actual intermediate states and all cuts are retained
in the certificates; two independent counting methods agree. Neither
frontier contains425. The first also excludes427, while the second has
a one-birth427 witness. The proper-divisor condition upgrades both
negatives to arbitrary tournaments containing the respective source.

A different, explicit five-core-grammar source with the same H141 has
the two-birth425 trace. Thus source replacement is a concrete mathematical
obligation even within the present grammar. A neutral scalar identity,
an unpointed profile equality, or an unbounded supply of unused vertices
does not perform that replacement.

The general frontier theorem provides a reusable cheap hostile test for
future proposals: test the actual source before spending time on a larger
birth budget or a more elaborate clock. A positive target must still be
followed through its actual marked response; finite membership at one
source is not a reusable induction on all target counts.
