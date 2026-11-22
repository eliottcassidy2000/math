# A Paley obstruction to universal physical binary and quaternary digits

**Proved.** Let \(P\) be the Paley tournament on seven vertices, with \(H(P)=189\). It has no induced supertournament with \(H\in\{2H(P)\pm1,\;4H(P)\pm1,\;4H(P)\pm3\}\), irrespective of the number of added vertices or choice of root. If old orientations may be changed, the minimum numbers of old arc reversals to reach the two binary digits are exactly **one for 379** and **three for 377**, still allowing arbitrarily many new vertices. This is a retained-source obstruction, not an assertion about rebuilding a different tournament of the same count.

## Preliminary facts

The odd-cycle packing formula expresses \(H(T[S])\) as a sum of nonnegative weights of packings each contained in a vertex subset of \(S\), including weight one for the empty packing. Thus \(S\mapsto H(T[S])\) is monotone and supermodular:
\[
 H(T[S\cup S'])+H(T[S\cap S'])\ge H(T[S])+H(T[S']).
\]
This applies to groups of added vertices, not just single apices.

If every new vertex sees an old tournament \(T\) uniformly, contract \(T\) to a marked quotient vertex. If its quotient SCC is a singleton, \(H\) is an odd multiple of \(H(T)\). If not, cutting a directed Hamilton cycle of that SCC gives at least three spanning quotient paths, and following the cycle from the marked vertex back to itself gives a two-visit skeleton. Fill the latter with two nonempty ordered paths covering \(T\); there are at least \((|T|-1)H(T)\) such covers by cutting every full path at each nontrivial position. Together the output is at least \((|T|+2)H(T)\). The same holds if only a subset of external vertices makes such a quotient SCC, by monotonicity.

For a single apex with tag set \(S\), the exact cut formula is a directed cut with symmetric weights \(W(u,v)\) equal to the sum of the two ordered cover counts with swapped seam endpoints. In a regular arc-transitive \(n\)-vertex tournament all weights are equal and the directed cut has \(s(n-s)/2\) arcs, \(s=|S|\). If \(Q\) is the two-slot cover count including empty slots, the weight is \((Q-2H)/\binom n2\). In particular any nonuniform one-apex extension has gain at least \((Q-2H)/n\). Whenever this exceeds \(H+1\), no induced supertournament gives either binary digit: a mixed apex alone overshoots, and the homogeneous case above excludes the other possibilities for \(H>1\) and \(n\ge3\).

## The exact window for the Paley source

Use the orientation \(i\to j\) iff \(j-i\bmod7\in\{1,2,4\}\), with outgoing masks \((22,44,88,49,98,69,11)\). Exact subset DP gives \(H=189\), \(Q=2016\), and \(W=78\). The one-apex counts are \(189+39s(7-s)\), so for mixed \(S\) they are \(423,579,657\). Vertices labelled mixed here and below see the old tournament, not a later partial state, in both directions.

There cannot be two mixed new vertices at any height at most 759. By supermodularity, if one has \(s=2,3,4,5\), the two already contribute at least \(189+390+234=813\). Otherwise each beats just one old vertex or loses to just one old vertex. Move the first exceptional vertex to 0 by translation, and move the second to \(0,1,3\) by the square multipliers, covering equality, outneighbor and inneighbor. The two exceptional contact types each have two possibilities, and their mutual arc has two orientations. The 24 exact induced count cases give minimum 855 (included in replay). Any bigger set with two such mixed vertices overshoots too.

With no mixed vertex the homogeneous dichotomy leaves only \(189,567\) under 759. With exactly one mixed vertex \(x\), write the other vertices as \(B,A\), before or after the old \(P\) respectively. Every \(B\)-vertex must beat every \(A\)-vertex: otherwise after deleting \(x\), the quotient has a cycle through \(P\). If \(A\) or \(B\) contains a directed triangle, it gains \(2H(P)=378\) in the restriction to \(P\) and that triangle, so supermodularity with \(P+x\) yields at least \(423+378>759\). Thus \(B,A\) are transitive.

After stripping ordinal extremes, the two transitive flanks of \(x\) have binary addresses \(u,v\ge0\). The established star formula has a nonnegative term \(H[p(u)+p(v)]\), a mixed term \(Q uv\), and one unsigned/signed read on each flank of the form
\[
 c w-z[w\text{ odd}],\qquad c=e+o,\quad z=e-o,\quad e,o\ge0 .
\]
Here \(e,o\) are the corresponding weighted even/odd rooted path reads. Both flanks nonempty overshoot by \(Q\). Any positive even address gives \(p(w)\ge2\), so also overshoots. For an odd address one, the linear read is \(2o\). For odd \(w\ge3\) the read is \((2o)+(w-1)c\ge 2(2o)\). Vertex transitivity supplies the same per-tag odd end/start mass 144, verified by subset DP and the converse isomorphism. Thus in either nonempty flank \(2o\ge288\). A nonempty flank must have address one and can yield only 711 below 760: with \(s=1\) or 6 the smaller return is 288 on base 423, while all other returns have larger minimum and base. Altogether the exact induced window up to 759 is
\[
 \{189,423,567,579,657,711\}.
\]
All six entries have literal realizations: 567 is the ordinal product with a directed triangle, and one-address 711 follows a singleton mixed apex.

## How many old arc changes are needed?

Arc transitivity fixes the first changed old edge to \((0,1)\). After this change the old count is 159. The mixed one-apex counts below 380 are \(315,353,373,375,379\), and the smallest relevant positive one-address flank gain is 222. Two mixed apices overshoot 379 by supermodularity. If there is only one mixed apex, the same before/after split applies; a directed triangle or either active flank overshoots. Thus no 377 appears after one old flip, while tag mask 8 gives 379 exactly.

To exclude 377 after two old flips, fix the first as before and test each of the other 20 old arcs for the second. These old tournaments remain strong since each original directed cut has at least three arcs in each direction. The exact finite replay gives the following simultaneous bounds, all across these 20 labelled hosts:

* old count \(h\ge111\), and none of the old counts divides 377;
* every mixed one-apex count is at least 189, and its gain over its own old host is at least 78;
* every positive one-address flank gain is at least 144, and each \(Q\ge1428\);
* any two mixed vertices have count at least 333;
* no count 377 occurs by adding zero, one or two vertices to any of these hosts.

These checks imply a global exclusion, not only a short search. No mixed vertex is impossible by the homogeneous dichotomy and divisibility. With at least three mixed vertices, supermodularity gives \(333+78>377\). With one or two mixed vertices, the uniform remaining vertices split into before and after groups as above, and either such group containing a directed triangle would add at least \(2h\) by supermodularity, overshooting. They must both be transitive. With exactly two mixed vertices, a uniform vertex that fails to be ordinal relative to one of them would add at least 144 by the one-address gain and supermodularity against the two-mixed host. Thus the uniform vertices are all truly ordinal and can be removed. With exactly one mixed apex, the two-flank star formula leaves at most one effective uniform apex: both active already add \(Q\), any positive even address adds at least \(2h\), and any odd address at least 3 adds at least \(288\). All remaining candidates are within the exhaustive two-addition replay.

Finally 377 is attained by reversing old arcs \((0,1),(0,2),(0,3)\) and adding one apex with tag mask 59. Thus the lower bound of three old changes is sharp; changing \((0,1)\) and adding tag 8 is the sharp plus construction. These repairs concern ordinary outputs, not one common physical seam between the two repaired endpoints.

The four quaternary tests further show the obstruction is one of retained structure, not scalar arithmetic. Target 755, 757 or 759 is attained after old flip \((0,1)\) with two consecutive apex tags \((1,247),(4,191),(4,239)\), respectively. Target 753 uses the two old flips \((0,1),(0,2)\) and tags \((32,177)\). In each pair of tags the second mask includes its bit for the first apex.

## Replays

binary-audit-check.py is standard-library Python for profiles, exact reduced Paley controls, and all repair witnesses. binary-old-two.cpp uses subset endpoint/cover weights to check every one- and two-apex output after at most two old reversals (up to the stated first-edge normalization); it also independently recounts the displayed two-apex witnesses by ordinary Hamilton subset DP.
