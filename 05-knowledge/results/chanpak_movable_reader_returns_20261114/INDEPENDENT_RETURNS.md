# Independent switch blocks: a residue decoder with an exact carry

A local change that preserves both exposed boundary responses is not
automatically a repeatable absolute-target return. In a common-apex
composition, it behaves like a CRT digit with a separate carry condition.
This generalizes the two-sided unit switch in [bilateral.md](bilateral.md)
from one switch center to arbitrarily many independently changed centers.

## Exact gluing

Take nonempty tournaments B_i with all cross-block arcs forward in the
listed order and a common apex v. Inside each block v beats U_i and loses
to its complement. Write h_i=H(B_i), f_i=F_2(B_i), and K_i for the count
after adding the apex to just B_i. The exposed profile has four more
coordinates:

* P_i is the weighted sum over subset paths that start in U_i; Q_i is
  its dual sum over subset paths ending outside U_i. A subset path on S
  has weight H(B_i−S); sums include all nonempty S.
* S_i and T_i denote the corresponding weighted sums with an additional
  sign (−1)^|S|. No endpoint-identity convention is needed below.

Let P=∏ h_i. Then the full ordinary Hamiltonian count is
\[
 N=\sum_i(P/h_i)K_i-(m-1)P+
 \sum_{i<j} P_iQ_j
     \prod_{i<t<j} f_t \prod_{t<i\text{ or }t>j} h_t+
 \sum_{j=i+1}S_i T_j\prod_{t<i\text{ or }t>j}h_t .
\]
In the final sum i ranges over adjacent block pairs. (The symbols P and
P_i distinguish the baseline product from an exposed coordinate.)

For the proof use the [Grinberg–Stanley odd-cycle packing expansion](https://arxiv.org/abs/2307.05569).
A packing not using v contributes P. A cycle using v within a single
block contributes (K_i−h_i)P/h_i. For a cycle visiting block indices
between i and j, the remainder of the cycle after deleting v is an
ascending sequence of block paths. At its two ends the contact
requirements give the P_i,Q_j channels; middle blocks use an arbitrary,
possibly empty selected subset and hence F_2. An odd cycle requires an
even total size of selected subsets and contributes two. The parity
projector therefore gives exactly the sum of unsigned and signed
channels. In any nonempty middle block, the signed contribution
vanishes: expand both complementary Hamiltonian counts by odd-cycle
packings and color each odd cycle or singleton selected/unselected.
Each has opposite signs for its two color choices. This leaves only
the adjacent correction. The same proof works if the graph is not
strong; strong witnesses are obtained, for instance, if the apex beats
the whole first singleton block and loses to the last singleton block.

## Target first, not a gap test

If a block is replaced by any alternative with the same six coordinates
(h_i,f_i,P_i,Q_i,S_i,T_i), the base count changes through K_i alone. Fix
one permitted reference option in every block, yielding N_0. Define
the digit set D_i to consist of the actual values (K_i−K_i^0)/2.
The displayed formula now gives
\[
              N=N_0+2\sum_i(P/h_i)d_i,\qquad d_i\in D_i.
\]
Every such independent choice is a literal tournament; no sign
subtraction between path counts is being used.

Suppose the switch-bearing h_i are pairwise coprime; other blocks
can be absorbed into N_0 and a common multiplier on the gains.
Remove that multiplier first (reject if it does not divide the
desired gain), and now let P be the product over switch-bearing
blocks. For target N, set k=(N−N_0)/2 after this common scaling.
Reject a noninteger. Set r_i to the least nonnegative representative of
\[
               k(P/h_i)^{-1}\pmod{h_i}.
\]
The switch can succeed **if and only if** there are choices
\[
\begin{split}
        j_i&\in\{(d_i-r_i)/h_i:d_i\in D_i,\ d_i\equiv r_i\pmod{h_i}\},\\
        \sum_i j_i&=\frac{k-\sum_i(P/h_i)r_i}{P}.
\end{split}
\]
The right side is an integer by CRT. This is an exact finite carry test;
saving one attaining option at each carry sum recovers the actual
chosen blocks. Even perfect residue coverage is not an absolute
compiler without coverage of the carry on the right. In particular
a unique representative for every residue gives just one possible
carry per residue. The h_i=1 case has no nontrivial switch at fixed
profile: by the packing expansion H=1 implies transitivity, and on a
transitive block P_i is a binary expansion determining all contacts.

Without coprimality all gains lie in the single lattice
\[
              2P/\\operatorname{lcm}(h_1,\ldots,h_m)\;\\mathbf Z .
\]
This follows from the gcd of the coefficients and is not by itself a
claim that all points of the lattice are attained. For repeating
one unit switch center with h_i=15, the exact gain set of m copies
is 2·15^{m−1}{0,1,…,m}, whatever fixed transitive wings are used.
No parameter hidden in the other contact words removes that
granularity while those words and the baseline N_0 are held fixed.
For coprime h_i≥3 with two options per block, at most 2^m of the P
residues are accessible. Independent binary choices at one same
baseline are therefore not a way around target-dependent selection.

The useful design target is now more precise: change the reader
and the cores jointly so that mod h_i choices are available *and*
the remaining quotient is supplied by a physical carry. It is
not enough to exhibit infinitely many independent fixed-gap pairs.

## An unbounded family of literal neutral returns

There is nevertheless an exact repeatable return once the baseline
has been chosen by the carry test. Use m copies of the 15-path switch
center, and put a transitive 2-vertex wing before and after all the
centers. The apex beats both vertices of the first wing and loses to
both of the last. Orient all cross-block arcs forward. Every choice is
strong, and deleting any non-apex vertex leaves a strong tournament:
one vertex remains in each wing, so every remaining interior vertex
is reachable from the first wing, reaches the last wing, and both
wings lie in the cycle through the apex. Deleting the apex leaves
an ordinal sum. Thus the apex is intrinsically identifiable, and
its puncture has its m strong centers in recoverable order. The two
local choices have different K, so have no isomorphism fixing their
apex. As a result, the `m choose j` words of weight j really give
that many pairwise nonisomorphic strong tournaments with the same N.

Changing a low block to high and a high block to low is an exact
neutral return, regardless of m, j or the positions. Such moves
connect every weight-j binary word, by swapping a 1 to a position
where the other word needs it. All six local boundary coordinates
are preserved throughout. But if one deletes a designated block
whose option differed between two words, the resulting literal
counts differ by 2·15^{m−2} (for m≥2); deleting a block tests the
selected hidden digit. This separates repeatability and readout
from absolute-baseline control: the large connected family still
sits on the explicit coarse lattice above.
