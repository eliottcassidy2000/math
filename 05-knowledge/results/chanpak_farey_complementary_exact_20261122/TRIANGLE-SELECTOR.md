# A triangle selector and its higher-slot state

Let q_m(T) count covers by m labelled directed paths with empty paths allowed;
write h=q_1, Q=q_2, and q_0=0 on nonempty tournaments. Take two nonempty
modules A,B and a singleton c, with c→A→B→c. Its exact count is
$$
H=2h_Ah_B+(Q_A-h_A)(Q_B-h_B).
$$
Indeed, a spanning path with c at either end uses intact paths in A and B.
Otherwise, cut the path at c. Each side consists of an A part followed
by a B part; the left B part and right A part must be nonempty. Within
each module all ordered splits with one named part nonempty are allowed,
yielding the two factors. Every case is reversible.

Set P_A=Q_A-h_A, positive since every cut of a Hamilton path gives a
distinct ordered two-cover. For fixed A an exact target N reduces to
choosing a feasible smaller h_B with
$$
2h_A h_B\equiv N\pmod {P_A},\qquad
Q_B=h_B+(N-2h_A h_B)/P_A.
$$
Any literal profile match is sufficient; a residue hit alone is not.
Q_A is even because complementation acts freely on the corresponding
vertex splits for a nonempty tournament, whereas h_A is odd by Rédei.
Thus P_A is odd, and if gcd(h_A,Q_A)=1 there is precisely one residue
class for h_B. Also Q_B≥(|B|+1)h_B implies
$h_B\le N/(2h_A+|B|P_A)$, an optional pruning bound.

## State propagation

Let W(A) be the same construction with B a singleton. Reversing B→c
produces an ordinal sum and changes the literal height from h+Q to h.
The full update is
$$
q_m(W(A))=m^2q_m(A)+m(q_{m+1}(A)-q_{m-1}(A)),\quad m\ge1.
$$
If the two anchors have different labels, their positions relative to
the A vertices are forced, contributing m(m−1)q_m(A). With a fixed
common label, the order B,c permits two optional ordered A parts
around the anchors, contributing q_{m+1}; the order c,B takes a
nonempty A part between the anchors, contributing q_m−q_{m−1}.
These cases prove the update, including m=1. An ordinal singleton
multiplies q_m by m; and q_m(A→B)=q_m(A)q_m(B) by splitting each
labelled path at the module boundary.

The coefficient of q_{k+1}(A) in H(W^k(A)) is k!: the only sequence of
k shifts reaching index k+1 is the all-upward sequence, whose successive
coefficients are 1,2,…,k. All other terms use smaller slots. Thus the
sequence of wrapped heights recovers every q_m, and equality of the
entire profile suffices for all those future heights. This is a concrete
operational minimality statement for repeated W observations, not a
claim every arbitrary numeric profile has a tournament realization.

There are two literal strong order-nine tournaments with the same
(h,Q)=(2231,36088), but q_3=302037 and 302277. See the row-generation
data in the replay. After one wrap both counts are 38319; after two
they become 782283 and 782763. In fact H(W²A)=−H(A)+5Q(A)+2q_3(A).
The replay verifies each original profile by two independent methods
and recounts both actual wrapped graphs by subset Hamilton DP.
The profile is a principal specialization of the
[Rédei–Berge symmetric function](https://arxiv.org/abs/2307.05569).
This selector does not prove the required smaller-profile existence or
Chan–Pak completeness.
