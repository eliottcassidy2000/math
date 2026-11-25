# The 7–13–21 Singer alphabet and a degree-57 orbit filter

This is about the conditional **free order-13 action** on the missing Moore graph, not a proof that the graph or such an action exists. It complements the [odd alphabet](THIRTEEN-ODD-ALPHABET.md) and the [cyclic crossing constraint](SQUARE-REFLECTION-CROSSING.md). It does not identify any of these phase parameters with a universally missing tournament height.

## Four shifts force a projective plane and nine witness orbits

Between two moving orbit labels \(A,B\), let \(D_{AB}\subset\mathbb Z/13\) give the edges \(A_x B_{x+d}\). A repeated ordered nonzero difference of two shifts gives two common neighbours of the same pair of vertices, hence a forbidden square. Therefore \(|D_{AB}|\le4\). At equality its differences cover all twelve nonzero residues, so every pair within either endpoint orbit already has a common neighbour in the other. Both orbit interiors are independent, and neither endpoint can have another multiple contact.

The exact finite check finds 52 sets in the affine orbit of \(\{0,1,3,9\}\), including 16 containing zero, and these exhaust the four-shift possibilities. The induced 26-site block is the point–line incidence graph of \(\mathrm{PG}(2,3)\). This is the middle example in the actual [Singer difference-set construction](https://personal.ntu.edu.sg/bernhard/Publications/pub/Singer.pdf):

| Phase size \(q^2+q+1\) | Shifts | Order of projective plane |
|---:|---|---:|
| 7 | \(\{0,1,3\}\) | 2 |
| 13 | \(\{0,1,3,9\}\) | 3 |
| 21 | \(\{0,1,4,14,16\}\) | 4 |
| 57 | \(\{1,6,7,9,19,38,42,49\}\) | 7 |

The former 7 and 21 discussion had no common mechanism at the level of Hamilton counts; here they have a genuine common mechanism at the level of cyclic phases. The latter two rows are design comparisons, not proposed cyclic automorphisms of the Moore graph. For example, cyclic order 21 is excluded by the [group restrictions](https://www.sciencedirect.com/science/article/pii/S0024379509003735). A saturated order-7 block, if it occurs in the order-7 branch, similarly gives the Heawood graph. A fixed vertex cannot touch either orbit in any saturated block: it would be adjacent to the entire orbit, duplicating within-orbit common neighbours.

In the free order-13 branch, each endpoint has 53 external singleton neighbours at the orbit level. For a common outside orbit \(R\), the two singleton shifts determine one relative phase \(\delta\) for a path \(A_x-R_y-B_{x+\delta}\). It must be a nonedge phase in \(\mathbb Z/13\setminus D_{AB}\); every such phase needs exactly one common outside label, since no internal edge can complete the route. Thus the 53-label lists intersect in **exactly nine** labels, with all nine phases different and none in \(D_{AB}\). Of the 248 outside labels, 44 touch just \(A\), 44 just \(B\), nine touch both and 151 touch neither. At the vertex level the 117 vertices over those nine labels are in bijection with the 117 nonincident point–line pairs. This is a direct SAT/CSP constraint on the marked quotient, not just a bound on its degrees. In the cyclic 125-label crossing model the crossing block is binary, so a saturated pair must be same-side.

## Triple contacts form alternating paths and long even cycles

Write a gap as a nonzero shift difference modulo sign, represented by 1 through 6. Each independent orbit must spend each of its six gaps exactly once. A cyclic orbit of step \(s\) has already spent \([s]\) on an edge and \([2s]\) on its internal two-step path. Hence it spends the other four via external contacts. A contact of \(m\) shifts spends exactly \(\binom m2\) different gaps at both endpoints. If there are \(t\) internally cyclic orbits of each step class, the total contact usage of each individual gap, counting each unordered contact once, is exactly \(125-t\).

There are eight possible triple gap sets:
\[
123,\ 134,\ 145,\ 156,\ 235,\ 246,\ 256,\ 346.
\]
If an orbit has two triple contacts, their gap sets partition all six gaps. The **only complementary pair is 134 / 256**. The triple-contact graph therefore has maximum degree two and alternates these two sets at every degree-two vertex. All cycles are even.

A 4-cycle would carry two shift triples of each type. After choosing directions, translate each set so that it contains zero; translation does not affect whether its fourfold sumset is the whole field. The exact check of \(6^4=1296\) combinations shows the sumset is always all of \(\mathbb Z/13\). One choice of shifts then closes a physical four-cycle through the four distinct orbits, impossible. Thus triple-contact components are paths or even cycles of length **at least six**.

For an internally cyclic orbit of step 1, the only triple type avoiding gaps 1 and 2 is 346, and its remaining double gap is 5. Scaling by \(s\) gives all cyclic cases. None of their types is 134 or 256. A triple cannot join two cyclic orbits: for unequal steps \(s,t\) it would have to avoid the gap classes \(s,2s,t,2t,s+t,s-t\), leaving at most two gaps (checked for all 15 pairs). The last two exclusions correspond to four-cycles using one internal edge in each orbit. Equal steps forbid even a single cross shift, by the parallel-step square. In particular any cyclic orbit incident to a triple is joined to a **degree-one independent** vertex in the triple-contact graph. Its triple edge forms an isolated two-label component.

## Two tournament pulses can read the full cyclic census

Enumerating all 64 translation-invariant 13-site tournaments gives this complete marked height spectrum (counts of sign masks in parentheses):
\[
3669497\ (12),\quad 3683797\ (12),\quad 3703011\ (4),\quad
3704857\ (12),\quad 3707483\ (12),\quad 3711175\ (12).
\]
Start with forward differences \(-1,\dots,-6\); read state zero without changing it, or read state \(s\in\{1,\dots,6\}\) by reversing the entire step class \(s\). Repeat with a second baseline obtained by first reversing classes 1, 2, 3. The two height readings, in order of states 0 through 6, are
\[
(3711175,3669497,3683797,3669497,3703011,3704857,3707483),
\]
\[
(3711175,3683797,3703011,3707483,3704857,3669497,3669497).
\]
The pairs distinguish all seven marked states. The first reading has \(v_5=2\) only at state zero, and \(v_{13}=2\) at states 1 and 3, with \(v_{13}=1\) otherwise. Thus an ordered join of first readings has \(v_5\) equal to twice its independent-orbit count, separately for each half if desired. Across the full census \(v_{13}=250+2t\), giving 260, 290 or 320, with only 260 left in the cyclic crossing subbranch. This counts a constructed census tournament, not an adjacency decoder for all external shifts, and a pulse reverses 13 arcs, not one.

## Phase-channel crossing rank

In the cyclic crossing model on 125 labels, \(f\) is a transversal modulo 25. Over \(\mathbb F_{13}\), \(\Phi_{125}\) is irreducible since the order of 13 modulo 125 is 100. It cannot divide \(f\): divisibility would make all five coefficients in each residue class modulo 25 equal, whereas exactly one is 1. The order-5 and order-25 packets still vanish, and the constant mode \(f(1)=25\) is nonzero modulo 13. Thus the crossing block has rank exactly 101 over \(\mathbb F_{13}\). Every nontrivial Fourier crossing block over \(\mathbb Q(\zeta_{13})\) reduces to this block at the prime \(1-\zeta_{13}\), so its rank \(h\) is at least 101. Its rank is at most 117, because it is proportional to the off-diagonal crossing block of the rank-117 negative projection in the surviving \(t=5\) stratum. Galois conjugacy makes the twelve nontrivial ranks equal.

Therefore the full 1625-site crossing block has rank **\(101+12h\), \(101\le h\le117\)**, not just an interval with endpoints 1313 and 1505. The rank floor and the Singer/triple rules are independent ways to reject proposed phase wiring. The replay checks the two finite modules, all height entries, and the modular packet arithmetic.
