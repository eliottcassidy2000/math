# A rank-three condition for three disconnected color contacts

## General theorem

Let F be a fixed rainbow forest with r+1 components and r free colors. Suppose at most three free-color contact graphs on the quotient by F are disconnected. If there are three, suppose one of them, called A, has graphic rank at least3. Then the rainbow spanning-tree completion fiber is connected under at-most-two-color moves, with diameter at most r+4.

The earlier sharper bounds apply when there are at most two bad colors. Graphic rank means the number of quotient vertices minus the number of connected components of the specified color-contact graph, including its isolated vertices. Actual edge realizations inside F components are retained throughout.

## 1. Preparing the other two bad colors

Call the other potentially disconnected colors B,C. Any acyclic selection of one B edge and one C edge extends to an acyclic three-bad-edge selection: the two-edge forest has rank2, while all A contacts together have rank at least3, so some A edge crosses two of its components.

As in the two-bad theorem, reconfigure the B,C type pair to the target type pair using at most three one-coordinate changes. A type is an unordered pair of F components, and the pair is acyclic exactly when its two nonloop types differ. The exceptional swap when both available type sets consist of the same two types preserves the uncolored quotient forest and can be made directly. Use the actual target edge realizations whenever the corresponding target type is installed. Type-preserving realization changes are also one-coordinate changes; the total partial-phase bound remains three.

The remaining issue is lifting each one-coordinate partial change without disturbing the full tree. This is addressed next.

## 2. Lifting one B/C change in at most two full moves

Suppose B_old changes to B_new while C remains fixed. Contract F and C. Both B edges are nonloops. The rank of the A contact graph after this contraction is at least2, since contracting one edge lowers graphic rank by at most1.

If B_old and B_new have the same unordered endpoint pair in this further quotient, replace B directly: the full quotient tree is unchanged. Otherwise they define different parallel classes X,Y.

The current A edge is not a loop or in X, because the current three-bad-edge forest is acyclic.

- If the current A edge is also outside Y, the B change leaves the bad forest acyclic. Lift it as a single-bad-coordinate change in one width-two full-tree move.
- Otherwise the current A edge lies in Y. If there is a nonloop A edge outside both X and Y, first replace A by such an edge, then replace B. Each proposed bad forest is acyclic, so each single-coordinate change lifts in one full-tree width-two move.
- If no nonloop A edge lies outside X∪Y, the A rank≥2 ensures that it has an edge in X as well as its current edge in Y. Replace A by an edge in X and B by B_new simultaneously. Before and after, the two uncolored quotient edge types are exactly X,Y. Hence the full quotient tree is unchanged; this is one valid width-two move with all good edges retained.

Here “single-bad-coordinate change lifts” uses the connected-contact repair: insert its proposed edge. Either the old same-color edge is on the fundamental cycle, or the cycle contains a good-color edge, since the proposed entire bad forest is acyclic. Remove the old same-color edge and that good-color edge, then reconnect the remaining cut with the good color. The other two bad edges stay fixed.

The exceptional B,C type swap from the first section also preserves the full tree with A fixed, because its uncolored quotient forest over F is unchanged. Thus at most six full moves install the exact target B and C edges. Then install the target A edge in at most one further full move: the target bad triple is acyclic. Once all three bad target edges are installed, the r−3 good target edges can be installed monotonically, one width-two move each. Total≤6+1+(r−3)=r+4.

## 3. Every five-free-color integer-difference fiber

On labels0,…,m, freeze any acyclic selection of colors6,…,m, leaving1,…,5 free. For every m≥5 the completion fiber is connected under width-two moves, with diameter at most9.

For m≥7, Q1 is connected because all distance-one edges form a spanning path. Q2 is connected because all distance-two edges form the two parity chains, joined by the frozen odd color-seven edge. Thus only Q3,Q4,Q5 can be disconnected.

There are six F components. The full distance-three graph has three residue-chain components. Adding F and then contracting F cannot increase that component count, so Q3 has at most three connected components. Therefore rank(Q3)≥6−3=3. If all three of Q3,Q4,Q5 are bad, use A=3 in the general theorem; if fewer are bad, the previous criterion suffices. This gives the bound r+4=9.

At m=5 the upper tail is empty; at m=6 it consists only of the compulsory edge{0,6}. These finite boundaries are checked exhaustively rather than inferred from the parity-chain argument.

## 4. Why a rank assumption matters

A generic three-bad-color obstruction exists on star components H,a,b,c. Let A offer Ha or Hb, B offer Hb or Hc, and C offer Hc or Ha. The only acyclic selections of all three colors are (Ha,Hb,Hc) and(Hb,Hc,Ha). They differ in all three colors, so width-two moves cannot connect them. Every bad color has rank2.

Additional vertices with connected good-color contact graphs do not remove this partial-forest obstruction: every full completion must still choose one of the two bad triples, and no width-two move can change between them. In the actual five-free-color integer-difference problem, the rank bound for Q3 rules out this pure cyclic-star realization. This is an arithmetic restriction on the obstruction, not a generic rainbow-tree connectivity claim.

## Limits

All moves preserve the frozen forest and one edge per required color, but can change the abstract tree shape. Global one-color GT1 and the Graceful Tree Conjecture are not proved. No theorem for arbitrary numbers of free colors is asserted. Elementary fundamental-cycle repair and finite quotient rank are used as tools; historical novelty is not claimed without further literature work.

## Exact checks and reproduction

Run `python verify_three_bad_rank.py` and `python -O verify_three_bad_rank.py`. The identical receipts contain86115 checks and8409 individually validated moves. The checker covers all752 source states of the two actual three-bad tails at m7, using three fixed targets per fiber (2256 source–target transformations). It checks both finite boundaries exhaustively:40 states at m5 and164 at m6, each with exact width-two diameter3. The generic rank-two cyclic-star hostile has exactly two bad selections. A separate abstract parallel-contact control verifies the rank-three swap's unchanged component partition and full-tree connectivity. Parallel contacts refer to the quotient by F, where different original edges can have different colors between the same components.

The finite controls support the implementation and boundary claims. The all-m theorem follows from the quotient rank proof. The [two-bad-color theorem](../graceful_two_bad_color_contacts_20261019/TWO_BAD_COLOR_CONTACTS.md) supplies the two-coordinate partial-state recipe and connected-color repair. No dependency on mutable solver output is required.
