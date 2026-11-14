# What the three 1/3 questions can and cannot select

All tournament statements below are independent of the conjectures. We use
[a recent primary paper on 1/3--2/3](https://arxiv.org/abs/2609.30888)
for background, and
[Dorbec--Henning, The 1/3-Conjectures for Domination in Cubic Graphs](https://arxiv.org/abs/2401.17820)
for the Verstraëte and Kostochka formulations and special cases. The latter two concern the ordinary domination number of undirected graphs, not independent domination.

## 1. Balance is exact inside a locked frame, not a reduction of all paths to a poset

Write a path with fixed hub order as \(R_0,x,R_1,y,R_2\); the spine is in ascending order internally in each run. Once the four endpoints or empty-run markers are fixed, each free site independently chooses from its eligible set \(E_i\). Thus the uniform law on the frame is the product of uniform distributions on the \(E_i\).

This really is a linear-extension model: for each site take a chain on \(\lvert E_i\rvert-1\) auxiliary symbols and an independent distinguished symbol. Its possible positions, uniform over the extensions, list the eligible slots in increasing order. Take the ordered sum over sites, making all symbols from an earlier site less than all symbols from a later site. The extension count is exactly the original product. Comparing a free site to a hub separating its eligible runs gives precisely \(1/2\) when it has two choices, and \(1/3\) or \(2/3\) when it has three choices. For two distinct fully eligible sites \(i<j\), their relative order in the path is \(i\) before \(j\) in six of the nine equally likely run assignments: the chosen run of \(i\) is no later than that of \(j\). By contrast the corresponding probability is \(3/4\) for two sites with the same two eligible runs.

The unmarked tournament path set is not the linear extensions of any poset on **those tournament vertices** if it has more than one path. In any linear extension of a nonchain poset some adjacent pair is incomparable, as otherwise transitivity makes the whole extension a chain. Swapping those adjacent symbols produces another extension. Two Hamiltonian paths of a tournament cannot differ by one adjacent swap, since their common intervening unordered edge would need opposite directions. The auxiliary one-site posets above encode one frame; they do not eliminate the global union over frames or its unknown weight.

## 2. Exactly what a genuine triple overlap certifies

Fix hubs \(x,y\) and distinct spine anchors \((a,b,c,d)\) in slot order (last of \(R_0\), first and last of \(R_1\), first of \(R_2\)). When there is a triple site \(i\), necessarily \(b,d<i<a,c\). The enforced directed arcs give two vertex-disjoint triangles
\[
                 x\to b\to a\to x,\qquad y\to d\to c\to y.
\]
Conversely let these two disjoint hub triangles be given and suppose the open spine intervals \((b,a)\) and \((d,c)\) share at least one site. They produce exactly **two** positive-triple frames for this pair: the indicated \(x,y\) frame and the reverse hub order frame with anchor tuple \((c,d,a,b)\). Every nonanchor site has an available run because the three eligibility intervals together cover the spine outside the anchors. Let \(u<v<w<z\) be the sorted four endpoints. Outside \((u,z)\) a free site has one run; in \((u,v)\) and \((w,z)\) it has two, and in \((v,w)\) it has three. Thus both orders give the same
\[
    \beta=\min(a,c)-\max(b,d)-1,\quad
    \alpha=|b-d|+|a-c|-2,\quad \text{pair contribution }2^{1+\alpha}3^\beta.
\]
These are non-overlapping frame pairs: the two triangles recover the anchors and the hub order recovers the frame. This describes exactly the ternary portion of the ordinary count, not just an existence test for one frame. In particular, no-disjoint-triangle tournaments such as the 385 and 573 examples have \(\beta=0\) in **every** two-hub frame. They may have three nonempty runs; they cannot have a free site common to all three. A generic pair of disjoint triangles need not yield triple overlap: sharing an **interior spine site** is essential.

## 3. A real cubic graph with domination <= 1/3 and seven invisible overlap patterns

The naive frame incidence graph, with free sites on one side and the three actual slots on the other, cannot directly import either cubic-domination conjecture. A cubic instance would be \(K_{3,3}\) with girth 4, and even two fully eligible sites already create a 4-cycle.

One can preserve all three labels by splitting slot vertices into copies. This yields a sharper test. For \(r=3\) or \(r\ge5\), use points \((i,j)\) with \(i\in\{0,1,2\}\) and \(j\in\mathbb Z_r\), and lines \((d,t)\) with \(d\in\{0,1,2\}\) and \(t\in\mathbb Z_r\). The incidence is \(j=t+di\). This is a connected cubic bipartite graph on \(6r\) vertices. For connectedness, the line with slope 0 links points with the same \(j\); moving along slopes 0 and 1 connects \(j\) to \(j-1\). Parallel lines are disjoint, and if two distinct points share two slopes then \((d-d')(i-i')\equiv0\pmod r\), impossible since its nonzero absolute value is \(1,2\), or \(4\). Thus its girth is at least 6. All the points with \(i=0\), together with all the lines of slope 0, dominate the graph, giving a dominating set of size \(2r\).

Now project line copies to one of the three physical slot names. Assign role 0 to slope 0, role 1 to slope 1, and role 2 to exactly \(q\) of the slope-2 lines, \(0\le q\le r\); assign role 0 to the remaining slope-2 lines. The graph and dominating set do not change. But the lists after projection are \(\{0,1\}\) at \(3r-3q\) vertices and \(\{0,1,2\}\) at \(3q\) vertices. This gives distinct local multiplicities
\[
                            2^{3r-3q}3^{3q}.
\]
The case \(r=6\) gives **36** vertices, an explicit dominating set of size **12**, and seven possible overlap weights on the same undirected graph. This is nonidentifiability of the undirected **rank/degree/girth/size** data in a concrete model, not a disproof of either domination conjecture or a claim about all labels required by any possible reduction. The equality only concerns named-slot choices, not matching the cubic graph's adjacency to the tournament.

Those lists have genuine physical locked frames: for \(\alpha=3r-3q,\beta=3q\), arrange the double sites before the triple sites and use the spine word below. Copies of a slot do not consume capacity in a Hamiltonian frame; enforcing a one-to-one matching to copies would define a different counting problem.

## 4. A constructive exact inverse on a new nonsmooth absolute family

Use the convention \(x\to y\), spine digit \(2[x\to i]+[y\to i]\). The strong word
\[
                        2\,0^\alpha\,1\,0^\beta\,00
\]
has a locked frame with anchor **positions** \(b,d,a,c\) in that increasing order: the leading 2, the 1, and the final two zeros. The \(\alpha\) sites between \(b,d\) have roles 0 and 1, while the \(\beta\) sites between \(d,a\) have all three. For strongness, the hubs and the site of digit 1 form a triangle; the leading 2 is reached from \(x\) and reaches that site; every 0 is reached from the leading 2 along the spine and beats \(x\).

For \(\alpha,\beta\ge1\) its **ordinary** path count is
\[
 N=1+2^{\beta+3}+2^\alpha(2^{\beta+3}+18\cdot3^\beta-1).
\]
For a homogeneous nonempty block, fix exactly which runs it occupies. Its completions inside the block are surjections onto those runs, yielding the basis \(1,2^e,3^e\) in each positive block length \(e\). Ordinary counts on the two blocks therefore have nine coefficients, uniquely certified by the \(3\times3\) grid at 1,2,3 in each coordinate. The independent subset-DP grid is \((155,419,1163)\), \((293,805,2261)\), \((569,1577,4457)\). This supplies the claimed whole-count formula, not a fixed-frame substitution. The frame construction itself works also with zero exponents; the displayed ordinary-count formula is only asserted for positive exponents.

**Absolute input, overlap output.** For prescribed \(N>1\), take \(v=v_2(N-1)\) and write \(F(\beta)=2^{\beta+3}+18\cdot3^\beta-1\), which is odd. Exactly three branches are possible:

- If \(\alpha<\beta+3\), then \(\alpha=v\). Integer bisection on the strictly increasing \(H(v,\beta)\), starting at \(\beta\ge\max(1,v-2)\), decides the branch.
- If \(\alpha>\beta+3\), then \(\beta=v-3\). The quotient \((N-1-2^v)/F(\beta)\) must be an integer power of two with exponent greater than \(v\).
- If \(\alpha=\beta+3\), then \(N-1=2^{\beta+4}(2^{\beta+2}+9\cdot3^\beta)\) with odd bracket, so \(\beta=v-4\), \(\alpha=v-1\), and direct substitution decides the branch.

Only positive exponents are legal, and every output gives both the absolute count and the actual overlap anchors. For example 805 gives \((\alpha,\beta)=(2,2)\), although \(805-1\) itself is not 2,3-smooth. It can also return multiple patterns: the total 8849 has both $(\alpha,\beta)=(4,3)$ and $(7,1)$ in this single family. The inverse supplies alternatives rather than a unique prescription. This does not cover arbitrary \(N\), but it is a target-to-overlap compiler on an explicit two-parameter family.

## 5. The numerical echoes, with their roles kept separate

[Session 334](https://github.com/eliottcassidy2000/math/blob/2fbf5edce71ae1278acbf5e07d902d7f72683b68/05-knowledge/results/session334_switching_compass_and_three_strands_20261114/README.md) used a period of **6** in the phase differences of a fixed block skeleton, and a different **36** as coefficient in the scalar rail \(1+36\cdot2^k3^\ell\). [Session 335](https://github.com/eliottcassidy2000/math/blob/056218f6a12c81850a22149e2c32c44473a29d57/05-knowledge/results/session335_three_slot_forests_pythagorean_addresses_20261114/FOREST_SWITCH_CHART.md) had **12** gluing orders for three nonempty runs and two hub orders, distributed over 16 binary cut charts. In the cubic construction, 6 is the girth lower bound, 12 the exhibited domination budget when \(r=6\), and 36 the graph's actual order. The exact common mechanism is that labels/weights are lost under a projection, not an identification of the numerical invariants.

The archival [three-slot factorial-moment theorem](https://github.com/eliottcassidy2000/math/blob/1b91742ce/01-canon/theorems/THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary.md) is another instructive analogy: a third statistic detects distinctions left after a projection to a two-dimensional mean-zero plane. The preceding labelled cubic construction specifies the actual lost coordinate here, so it avoids treating that analogy as a deduction.
