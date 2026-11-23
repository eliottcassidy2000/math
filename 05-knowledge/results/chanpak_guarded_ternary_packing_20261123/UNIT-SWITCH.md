# An exact renewable unit switch in the five-read guarded apex state

Work with the positive/signed transfers for the deletion word of a marked apex $z$. The first block is the singleton O, beaten by $z$, and the final block is the singleton I, which beats $z$. This makes every output tournament strong. For a guarded word use the five coordinates $(D,A,B,Q,H)$: deletion mass, scaled positive start and end reads, scaled positive unrestricted read, and the actual Hamilton count. The two signed outside entries are zero. The slot counts and probe definitions agree with the five-read construction.

## A physical switch fixing four coordinates

Between O and I insert one six-vertex ordinal block $C$ with outgoing row masks $(58,4,9,18,38,14)$. The apex beats mask 28 in that block, using zero-based labels. Reverse the actual internal edge $(2,5)$, preserving every contact and both guards. The two complete states are exactly
$$
(29,898,964,1280,751),\qquad(29,898,964,1280,753).
$$
Thus the switch acts by a pure translation on $H$ and does not perturb any of the four continuation reads.

For verification let $h=H(C)$ and $q=Q_2(C)$, the allowed-empty ordered two-slot count. Within $C$ weight a nonempty path by the Hamilton count of its complement. Write $U_e,U_o$ for paths starting at an apex out-neighbor and $V_e,V_o$ for paths ending at an apex in-neighbor, splitting by parity of the number of vertices. Let $W_e$ be the even read with both restrictions. Each orientation has $(h,q,U_e+U_o,V_e+V_o)=(29,320,129,162)$, while $W_e+U_o+V_o$ changes from 201 to 202. Splitting each apex-containing odd cycle into its path in the deletion word, its contribution to the height is twice that complement-weighted even-path count. The possibilities are both endpoints in $C$, the first at O, the last at I, or both guards. The last contributes $q$: the complement-weighted even unrestricted slot read including the empty path is $q/2$ by signed slot cancellation. Consequently the full state of O-C-I is
$$
(h, 2(q+U_e+U_o), 2(q+V_e+V_o), 4q, h+q+2(W_e+U_o+V_o)).
$$
The companion replay computes these counts from literal path permutations and the Hamilton counts of complements, and separately evaluates the complete tournaments by subset DP.

## Guard renewal and both ternary signs

Appending a singleton O sends $(A,Q)$ to $(2A+D,2Q)$ and leaves $(B,H)$ unchanged. Appending I sends $(A,B,Q)$ to $(2A,B+Q,2Q)$ and adds $A+D[p=O]$ to $H$, where $p$ is the previously last singleton. These are the positive transfer and its nearest-neighbor signed correction. Keep the old first guard and require the new final letter to be I; temporary letters O need not themselves be final guards. The recurrences now preserve the unit difference and all shared coordinates for every subsequent common right singleton word, however often the mark is switched. Equally, a common arbitrary left/right ordinal singleton context preserves the unit gap since the old word has zero signed outside reads and all new terms depend only on its shared four coordinates.

Two specific finite guard restorations make the lower or upper reading, respectively, divisible by three at every further width:
- Append OIOI to both cores, then $m$ further I's. The lower height is $14658\cdot2^m-4695=3(4886\cdot2^m-1565)$; the upper height is greater by two.
- Append IOI to both cores, then $m$ further I's. The upper height is $7242\cdot2^m-1941=3(2414\cdot2^m-647)$; the lower height is less by two.

These are exact strong one-arc neighbors at every $m\ge0$. In particular the middle height in each line is really $3H$ of a separately attainable input at $m=0$: the two input heights are $3321=3^4\cdot41$ and $1767=3\cdot19\cdot31$. The factors 19, 31, 41 have respective explicit rows $(48,1,3,39,46,6)$, $(52,33,34,7,46,8)$, $(48,33,35,7,14,24)$ and the factor 3 is a directed triangle. Taking ordinal products proves the claimed roots, not merely their divisibility. Two other directly certified input roots are $251$ (the original upper state 753) and $2085=3\cdot5\cdot139$ (append OII to the lower state to get 6255). The replay includes their masks. This note does not certify all members of either infinite quotient family as attainable.

## Why a mere mass-three graft cannot perform the full uniform ternary step

The guarded states do form a closed composition semigroup. For two guarded words write their five coordinates with indices 1,2. Concatenation is
$$
\begin{aligned}
D'&=D_1D_2,&Q'&=Q_1Q_2,\\
A'&=D_1A_2+A_1Q_2,&B'&=B_1D_2+Q_1B_2,\\
H'&=H_1D_2+D_1H_2-D_1D_2+A_1B_2.
\end{aligned}
$$
These follow from triangular transfer multiplication, with zero signed start/end interactions. The same algebra shows that an isolated switch of height difference two in one guarded word becomes a signed height difference of $2D_\text{ext}$ when surrounded by arbitrary fixed contexts of combined deletion mass $D_\text{ext}$.

There is a stronger no-go for a context that leaves the original guarded word literally unchanged and contiguous and tries to upgrade it by total external deletion mass three. The odd-cycle packings in the old apex tournament survive; the external deletion components independently contribute their factor three, already producing $3H$. If there is any left context, its first guard must be a singleton O. That new guard and the old final I make a new directed triangle with the apex. Its complementary Hamilton read is $3D$, since removing guard singletons changes neither inside nor outside deletion mass. If instead there is only right context, use old initial O and new final I. In either case the new triangle contributes at least $6D$, so the new height is at least $3H+6D$. It cannot be any of $3H-2,3H,3H+2$. This proof uses actual packings, not a modular approximation; it does not cover replacing or retagging part of the old word.

The construction and obstruction thus separate two senses of renewal: the pure physical unit switch can be re-used after a guard-restoring context, but a context graft of mass three cannot itself transport an arbitrary guarded input to its correctly centered ternary child.

## Both signed exits from the same attainable guarded input height

There is also a common height \(h=12693\), with an explicit guarded
parent. Use \(O\,P\,I\) with
\[
 P=(214,68,72,51,38,71,152,62),\qquad \text{apex contact mask }71.
\]
Its literal Hamilton count is \(12693\). From the original six-block
\(C\) and the same marked pair \(\{2,5\}\), use two guard-preserving
contexts:
\[
\begin{array}{c|cc|rr}
 &\text{left of }C&\text{right of }C&\text{before flip}&\text{after flip}\\\hline
\text{lower}&OIOIO&IOI&38077&38079\\
\text{upper}&OIIOOIO&I&38079&38081
\end{array}
\]
The shared center is \(38079=3\cdot12693\). Both pairs preserve
their four accompanying measurements and both guards. The two
center hosts are not the same guarded state (their \(A\) coordinates
are 54850 and 49538 respectively). The independent parent, likewise,
is not asserted to sit as an induced tournament in these hosts.
This is a genuine two-context ternary fork from a common attainable
height, not a uniform rewiring of an arbitrary source or a lift
from one common center. The companion replay recounts all five
literal tournaments.
