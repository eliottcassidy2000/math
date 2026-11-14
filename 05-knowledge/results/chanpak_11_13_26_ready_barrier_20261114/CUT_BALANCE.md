# Balanced cut flips: an ordinary-target quadratic alternative

This construction is separate from the affine source lines for terminal
\(7p\). A changed middle contact, with two balanced outer binary
contexts, produces arbitrarily many exact **ordinary** two-count
repairs. It does not prove that a pair starts at every target.

## Exact binary readout for a transitive spine

A bit word \(W\) gives the arcs of a tournament with one apex over an
otherwise transitive spine. Bits are read in the spine's forward order;
a 1 means the apex beats that vertex. Write \(W\) also for its binary
integer (most significant bit first). Number bit positions from the
right by \(0,1,\ldots\). Let
\[
 E(W)=\sum_{\substack{i<j\\W_i=W_j=1}}\omega(j-i),\quad
 \omega(1)=1,\quad \omega(d)=2^{d-2}\ (d\ge2).
\]
This is the even-endpoint count in the two-contact formula, generalized
to any contact set: an ascending path whose two endpoints are contacted
has either no interior slots or a freely chosen even subset of its
interior slots. A lone rooted start over a transitive spine has
\(R_+=2^i\) at bit position \(i\), and \(R_-=0\) except for the last
spine vertex where \(R_-=-1\). The odd-cycle packing identity (or
counting the two ascending runs split by the apex directly) gives
\[
 H(W)=1+W-(W\bmod2)-2E(W).
\]
The zero-follower formula is important: use both raw signs, not the
half-of-all-subsets argument for positive follower width. In that
positive-width variant the first three terms are instead \(1+2^r W\).

If \(\overleftarrow Y\) denotes the reverse of \(Y\) in its declared
length, then
\[
 2E(XY)=2E(X)+2E(Y)+X\overleftarrow Y+X_{\rm last}Y_{\rm first}.
\]
The final correction is exactly the extra half per adjacent marked pair,
multiplied by two. The cross-distance factor for every other pair splits
as the product of its high and reflected low bit weights.

## A rank-two balance and a genuine single-arc unit switch

Take \(W=A\,101\,B\), where \(b\) is the declared length of \(B\).
Change **one arc**, that of the central vertex, by replacing 101
with 111. The three-bit integer and its reversal both rise by
2, the internal \(E\) rises by 2, and the first and last central bits
do not change. Write \(U=2^b-\overleftarrow B\). Direct substitution gives
\[
                        \Delta H=2(U-A-2).
\]
To get any prescribed even difference \(2d\), take
\[
 U=A+d+2,\qquad B=\overleftarrow{(2^b-A-d-2)}.
\]
Every \(A\in[2^{b-1},2^b-d-2]\) is legal when this interval is nonempty
(for integer \(d\ge1\)); the first spine bit is 1 and the last is 0,
so both full tournaments are strong. Every spine vertex belongs to the
apex's strong component via the first and last spine vertices.
In particular for \(d=1\) a single contact change gives exact \(+2\)
in exponentially many contexts.

This repairs a real ordinary count, not just a cofactor or residue.
Examples, verified by a direct two-run Hamiltonian path DP and trial division:

| \(b\) | \(A\) | \(B\) | lower \(H\) | upper \(H+2\) |
|---:|---:|---:|---:|---:|
| 9 | 434 | 420 | 1,517,507 prime | \(1,517,509=7\cdot216,787\) |
| 10 | 646 | 954 | \(3,354,407=7\cdot479,201\) | 3,354,409 prime |

The respective cofactors 216,787 and 479,201 are prime. This is not
an infinitude claim for these prime-pair types.

## Absolute target inverse and genuinely distinct baselines

For fixed \(b,d\), substitute the balance into the lower count:
\[
 H(A)=(8A+5)(A+d+2)-5A+
 H(B)-2E(A)-2-A_{\rm last}-B_{\rm first}.
\]
For fixed \(d\) the leading term is quadratic in the **integer**
context address \(A\); its growing derivative is the cost of the
balanced unit, not an artifact of using two source coordinates.

For each high marked bit at position \(j>0\), the sum of its pair
weights with *all* lower bits is \(2^{j-1}\). Hence
\(0\le 2E(A)\le A\), and \(1\le H(B)\le2^b+1\).
The lower-order term
\(H(B)-2E(A)-2-A_{\rm last}-B_{\rm first}\)
lies in \([-(2^b+3),2^b]\), a common interval of width
\(2^{b+1}+3\). Consecutive leading terms on the strong range
increase by \(16A+8+8(d+2)\ge8\cdot2^b+8\).
Thus the count is **strictly increasing** in \(A\). Given an absolute
target \(N\), binary-search this exact expression at each allowed
length \(b\), then accept only literal equality (or equality to
\(N-2d\) when looking for an upper-member repair).
There are only \(O(\log N)\) lengths since \(H\) is bounded below
by \(2^{2b}\), and each search uses \(O(b)\) exact readout evaluations.
No enumeration of tournaments is needed. The decoder does not claim
that any width succeeds for every target.

For each fixed \(d\) and all large \(b\) there are \(\Theta(2^b)\)
distinct lower counts of size \(O_d(2^{2b})\) and hence
\(\Omega_d(\sqrt X)\) distinct ordinary pairs with both endpoints
at most \(X\). There are also only \(O_d(\sqrt X)\), since at each
allowed width there are at most \(2^b\) contexts and one sums a
geometric series to \(2^b\le\sqrt X\). Thus the true order of
growth is \(\Theta_d(\sqrt X)\), not positive density.

## Why the named 11/13 reflection cannot balance

Use central 1011 vs 1101 (11 vs 13), with common ends.
These reverses have the same internal \(E=4\), but the ordinary
integer rises by 2 while the reversed integer **falls** by 2.
The same derivation gives \(\Delta H=2(U+A)\), with both signs
now positive. There is no free cancellation: \(+2\) occurs
only when \(A=0,U=1\), i.e. all right bits are ones,
yielding counts 3 and 5 up to neutral padding. With one forward
follower the two raw local terms are \(2\cdot11=22\) and
\(2\cdot13=26\); that follower only increases the monotone gap.

The underlying signed test extends to arbitrary tournament reservoirs
around the central transitive 4-block. With the relative ordering
\(L\to TT_4\to B\), retain any fixed contacts on \(L,B\). Write
\(h_L,h_B\) for their path counts, \(q_B=F_2(B)\),
\(R_L\) for the sum of positively weighted rooted *start* responses
on the contacted vertices of \(L\), and \(Q_B\) for the analogous
rooted *end* response of the contacts in \(B\). Empty reservoirs
have \(h=F_2=1\) and empty contact response zero. Then the
ordinary 11-to-13 shift with \(r\ge0\) actual followers is
\[
 \Delta H=2\{h_L(2^r q_B-Q_B)+h_B R_L\}.
\]
To prove it, in the even-endpoint sum \(E\) the cross contribution
from \(L\) into the central block changes by \(-h_B R_L\):
the center's reflected integer drops 2, and the even projector
divides the unsigned change by 2; the signed contribution is
unchanged since the first central bit stays fixed. The contribution
from the central block into \(B\) changes by \(h_L Q_B\):
the center's normal integer rises 2 and the last central bit
stays fixed. Other contact pairs are unchanged, including those
entirely inside the center since both have \(E=4\). The positive
rooted start response rises by \(2h_Lq_B\), and the negative
rooted start response is unchanged. At \(r=0\) use the latter
in addition to the positive response; at \(r\ge1\) it vanishes
under summing the signs of the follower subsets. This proves
the formula including both boundaries.

Here \(q_B-Q_B=h_B+\sum_{v\notin U_B}Q_+(B,v)\), where the
\(h_B\) term is the empty selected subset. Thus the lower bound
is \(2h_Lh_B\) and the unit case again has only counts 3 and 5.
More generally a fixed gap \(2d\) forces \(h_Lh_B\le d\),
\(R_L\le d\) and \(q_B-Q_B\le d\). A contacted vertex in
\(L\) cannot have more than \(\log_2 d\) vertices in subsequent
strong components (collectively they multiply its start response
by their \(F_2\), which is at least \(2^{|V|}\)).
The analogous restriction holds for missed contacts in \(B\)
looking backwards; for \(r\ge1\) all of \(B\) already has
size at most \(\log_2 d\) from the \(q_B\) term, and
\(r\le\log_2(d+1)\). At \(r=0\), all remote leading
singleton components of \(L\) and remote trailing singleton
components of \(B\) are neutral and can be deleted. Remaining
nonsingleton components have bounded size from Busch's exponential
strong tournament path lower bound and their \(H\)-product is at
most \(d\). Therefore this named reflection with fixed gap cannot
yield unboundedly many absolute pairs, even with arbitrary reservoirs.
The useful change is to make both radial coefficients have the
**same** sign, as in 101-to-111, rather than optimize a positively
signed reflection.

The two literal five-vertex cores below provide a separate exact
11/13/26 equality, not a substitute for the cut calculation:
\[
 (\mathbf a,v;H,R_+,R_-)=
 ((24,1,3,6,14),0;11,13,-1),\quad
 ((18,8,3,5,14),4;13,26,-2).
\]
Positive one-contact bridges of widths \(r+1\) and \(r\)
respectively then yield \(11+26\cdot2^r\) and
\(13+26\cdot2^r\). These are physically correct, but just
one-width rails, whereas the cut balance allows a growing
choice of context and an absolute inverse in the ordinary spectrum.

## Provenance and replay

The change-of-contact count is a specialization of the first
[two-contact extension](https://github.com/eliottcassidy2000/math/blob/3c4ef16d160baa9aca5ccbb7bbb86e8d52b37905/05-knowledge/results/chanpak_literal_prefix_compiler_20261114/TWO_CONTACT_REPAIR.md);
the starting one-apex inversion count is also recorded in the
[nonuniform-apex transfer](https://github.com/eliottcassidy2000/math/blob/9717032b8bf0afc5371763d0f26b3810d0332f0d/05-knowledge/results/chanpak_terminal_seven_nonuniform_contacts_20261114/NONUNIFORM_APEX.md).
The new ingredient here is solving the *signed interface difference*
by exact bit reversal at the ordinary-target level, together with
proving its strict quadratic ordering. No novelty is claimed for the
binary path counter itself. The cited odd-cycle identity is that of
[Grinberg and Stanley](https://arxiv.org/abs/2307.05569); the exponential
bound used only for the 11/13 bounded-gap claim is
[Busch's](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).

Run replay_cut_balance.py.
