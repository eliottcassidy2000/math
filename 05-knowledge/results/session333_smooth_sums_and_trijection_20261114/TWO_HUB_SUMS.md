# Smooth sums are exact; an independent smooth-frame compiler is impossible

The problem is the proposed converse in [Wave 6](../chanpak_waves_20261113/README.md): the two-hub path count splits into a sum of \(2^a3^b\) terms. Can an arithmetic divisibility antichain be installed as the list of literal frames of one tournament? The answer is **no in that literal form**, even with one baseline path separated out. A successful proof has to regroup interacting frames, or use another physical interpretation of the arithmetic digits.

## The arithmetic lemma, with an exact decreasing inverse

Every positive integer is a sum of powers \(2^{a_i}3^{b_i}\) such that the \(a_i\) strictly increase and the \(b_i\) strictly decrease. In particular no summand divides another. This is the classical primitive 3-smooth representation setting ([Blecksmith, McCallum and Selfridge](https://www.tandfonline.com/doi/abs/10.1080/00029890.1998.12004921)); here is an elementary constructive proof of existence. For even \(N\), represent \(N/2\) and double every term. For odd \(N\), let \(3^b\le N<3^{b+1}\); either equality gives one term, or write
\[
N=3^b+2M,\qquad 0<M<3^b.
\]
Represent \(M\) inductively, then double its terms. Their 2-exponents are positive and their 3-exponents are all below \(b\) because each term of \(M\) is less than \(3^b\). So adding \(3^b\) at the beginning preserves the entire antichain. The inverse is positive and strictly decreases the absolute target, not just a residue.

This gives \(34=18+16\), \(45=27+18\), and \(69=27+18+24=45+24\). These are actual arithmetic sums, with no missing hidden integer carry. Requiring *consecutive* 3-exponents instead is an extra condition and false: 17 has \(9+8\), but a consecutive three-exponent form needs at least \(9+6+4=19\), and a consecutive two-exponent form cannot use \(3+2^a=17\) or \(9+3\cdot2^a=17\).

## Positive frame formula, including empty runs

Let deleting two hubs \(x,y\) leave a transitive spine on ordered labels. For hub order \(x,y\), a path is \(P,x,M,y,S\), where every spine run is increasing. A frame records \(a=\max P\), \(b=\min M\), \(c=\max M\), \(d=\min S\); an empty run uses a sentinel rather than a vertex. If \(b=c\) that one middle anchor is used only once. Test the hubs' boundary arcs at the anchors, and test \(x\to y\) if the middle run is empty. Across nonempty runs the anchors must be distinct. Each other spine vertex independently chooses one of the runs whose interval contains it: the left open ray below \(a\), the open interval \((b,c)\), or the right open ray above \(d\), omitting empty runs. Within each run the order is then forced. Configurations with no choice at some vertex contribute zero. The two hub orders and all frames partition the full path set, so
\[
H=\sum_{\text{admissible frames}}2^\alpha3^\beta ,
\]
where \(\alpha\) counts exactly-two coverage and \(\beta\) exactly-three coverage. This is the Wave 6 formula with the endpoint conventions made explicit. The ternary sites form the common interval
\[
(\max(b,d),\,\min(a,c))
\]
on the spine, excluding anchors, not a freely scattered register of digits.

### A compulsory weight-one obstruction, even after one baseline is removed

Start from the unique increasing spine path and insert hubs one by one. For any ordered directed path, writing a 1 when the new vertex beats a vertex of the path, the possible insertion gaps are the ascents \(0\to1\) after adjoining a leading 0 and a trailing 1. There is always at least one, and the number is one plus the number of descents. At the first insertion stage producing a nontransitive induced tournament, the preceding tournament was transitive, and its ordered path was unique. The new bit word then has a descent: otherwise the new vertex could be put at its unique sorted-order cut, giving a transitive orientation. Thus there are at least two paths preserving the spine order at that stage. Every subsequent insertion preserves their distinctness on restriction.

Any spine-order-preserving completed path is a unit frame: the occupied runs have disjoint ranges as the spine labels increase across the path, so no free vertex has multiple placements. It follows that **every nontransitive two-hub tournament has at least two weight-one frames**. The literal list can never be a divisibility antichain; even removing one designated unit frame leaves a unit among at least two remaining contributions (the Hamilton count is odd and at least three). This is independent of resource limits or finite enumeration. It does not rule out mathematically meaningful carrying or grouping of these frames.

### Ternary sites impose a paired contribution

If \(b<c\) and \(d<a\), take the frame \((a,b,c,d)\) for order \(x,y\) to \((c,d,a,b)\) for order \(y,x\). Exactly the same forced hub arcs are used. For any nonanchor \(i\) the old number of choices is
\[
[i<a]+[i>b]+[i<c]+[i>d]-1,
\]
which is the same in the reversed hub order (use \(b<c\) or \(d<a\) to expand the middle interval). Thus the two frames carry identical weights. Every frame with \(\beta>0\) satisfies both inequalities, so its exact contribution occurs in a pair. An independent unpaired ternary monomial cannot just be installed as one physical weight.

There is an all-order size cost: a frame with \(\beta>0\) forces
\[
H\ \ge\ 2\cdot3^\beta+4\cdot2^\beta+1.
\]
Restrict to the four anchors, the \(\beta\) triple sites and the two hubs; insertion afterwards only increases the path count. The paired frames give \(2\cdot3^\beta\) paths, each with a forced backwards spine bridge through each hub, \(a,x,b\) and \(c,y,d\). Neither hub can be deleted from those paths while retaining a valid path. Now delete \(y\) from the graph, use only the forced backwards bridge \(a,x,b\), and let each intermediate spine vertex choose either side of \(x\). This gives \(2^{\beta+e_x}\) paths, where \(e_x\) counts the other two anchors inside \((b,a)\). Insert \(y\) into each, keeping the paths distinct. Do the reverse construction with \(x\) deleted to get \(2^{\beta+e_y}\) outside the paired frames. These two lifted families are disjoint: if the hubs are not adjacent, deleting \(x\) from a path whose deletion of \(y\) uses \(a,x,b\) leaves the illegal backwards edge \(a,b\); if the hubs are adjacent, the external predecessor and successor would have to be both pairs \((a,b)\) and \((c,d)\), impossible. Each of the two groups of anchors (left of every triple site, right of every triple site) contributes one interior anchor to exactly one interval, so \(e_x+e_y=2\) and \(2^{e_x}+2^{e_y}\ge4\). An increasing-spine path is in none of these families. This proves the bound, which already excludes triple depth three at **either 45 or 69**.

## A stronger absolute obstruction for fixed runs and fixed frame budgets

Here **spine colour runs** mean maximal consecutive strings with the same two-bit orientation to the hubs. They are not the three path runs \(P,M,S\). This distinction matters for growth.

For every fixed maximum \(r\) of spine colour runs, the number of possible Hamilton counts at most \(X\) is \(O_r((\log X)^r)\), even allowing tournaments of arbitrarily large order. More quantitatively, fewer than \(O(X^{0.9})\) counts at most \(X\) have any witness in this class with at most \((\log X)/5\) colour runs, for all sufficiently large \(X\).

Proof: strip every singleton SCC. The count is unchanged by ordinal product. Every surviving strong component contains a hub, since the spine is transitive, so there are at most two. A strong tournament of order \(n\) has \(H\ge 5^{(n-1)/3}\): for the near-transitive case use \(1+2^{n-2}\ge5^{(n-1)/3}\) starting at \(n=3,4\); otherwise use Busch's strong one-vertex overlap and multiply the factor bounds, since \((n_1-1)+(n_2-1)=n-1\) ([earlier overlap argument](../session327_clock_orbits_and_activation_20261114/OPTIONAL_SLOT_FLOOR.md)). By the SCC product the remaining order is at most \(2+3\log_5 X\). If a hub was removed, a dummy neutral vertex may be reserved as the second mark without increasing the number of spine runs. The surviving spine word is a subsequence of the original, hence never gains colour runs. For \(L=O(\log X)\) positions with number \(L\le 3\log_5 X+4\), the total number of words with \(k\) runs and length at most \(L\) is \(4\cdot3^{k-1}\binom Lk\), by summing positive run lengths. A factor two accounts for the hub edge. This is \(O_r((\log X)^r)\) for \(k\le r\). For \(r\le(\log X)/5\), the elementary bound \(\binom Lk\le(eL/k)^k\) gives exponent at most \(0.2\log(15e/\log5)<0.87\) to leading order in \(\log X\). Polynomial logarithmic factors fit into \(O(X^{0.9})\).

Independently, even an unconstrained multiset of positive smooth terms needs length at least \((1-o(1))\log X/\log\log X\) for almost every target through \(X\): there are \(O((\log X)^2)\) available terms through \(X\) and at most \(\binom{C(\log X)^2+\ell}{\ell}\) multisets of length at most \(\ell\). At \(\ell=(1-\epsilon)\log X/\log\log X\) this is \(X^{1-\epsilon+o(1)}\). In particular a bounded number of frames with unbounded sizes cannot prove a cutoff. In the stronger spine-run test, a successful dense absolute construction would have to implement **linearly many switches in log target size**, not a fixed collection of growing constant runs.

This rules out two simple versions of the proposed inverse. It does not prove that the unrestricted two-hub language is sparse: its changing frame configurations are precisely the unresolved part.
