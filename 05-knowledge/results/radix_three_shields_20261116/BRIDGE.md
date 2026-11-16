# A four-state bank with exact factor-three sensitivity in strong tournaments

For a root x of U, put
\(P_\epsilon(U,x)=\sum_{S
i x, |S|\equiv\epsilon(2)}
 H_{\rm start}(U[S],x)H(U-S)\), with \(H(\varnothing)=1\).
Write \(E=P_0,O=P_1\). For a right module V use paths terminating
at its root y, equivalently take its converse before calculating.
Join U to V with all cross arcs forward except for y→x.

By the odd-cycle-pack formula for tournament Hamiltonian counts,
\[
 H=h_Uh_V+2(E_UO_V+O_UE_V).
\]
Indeed, no crossing cycle gives the independent-cycle-pack product.
Any crossing cycle must use the lone backward arc y→x exactly once.
After that arc its path consists of a root-start path in U, one
forward cross arc, and a root-end path in V. Such a cycle is odd
exactly when the two path orders have opposite parity. The remaining
vertices have independent internal cycle packs because both roots
have been used. If U and V are strong and both have at least two
vertices, the full join is also strong: all paths from right to left
can route through y→x, and an unreversed forward arc provides paths
from left to right. The formula needs neither module to be strong.

Four strong left modules in row-bit adjacency format (0-based)
realize \((h,E,O)=(109,131,141+j)\), for j=0,1,2,3:

| j | Root | Outgoing row masks |
|---|---:|---|
| 0 | 0 | 94,60,56,112,96,65,6 |
| 1 | 3 | 96,65,3,39,15,22,60 |
| 2 | 4 | 72,97,3,6,15,29,60 |
| 3 | 2 | 30,60,120,80,96,73,3 |

Use a root of the triangle on the right for \((h,E,O)=(3,1,2)\):
the four outcomes are \(1133+2j\). Use instead the converse of
(\(8,1,3,6\)), root 3, to get \((h,E,O)=(5,3,4)\):
the outcomes are \(2439+6j\). Each of the eight is a strong
tournament. Thus a factor-three shield for the *difference* exists
even with nontrivial cross-cycle contributions. It changes the
baseline, and its ordinary Hamilton count is 5, not 3. The adjacent
four-state bank is attenuated neither by disconnection nor by uniform
multiplication of every Hamilton count.

Run `python3 check_bridge.py` to enumerate the local ports by
literal permutations and independently verify the complete joins by
terminal-vertex subset DP. This does not yet give a second independent
four-state bank, a radix-3 cascade, or a proof of absolute placement.


## Retained-root transfer

Keep x as the root for another reverse bridge and write (D=H(U-x)).
A path starting at x cannot use y→x, and its complement has no
backward cross arc because x was used. Set
(K_i=\sum_{|S|\equiv i(2)}H(V[S])H(V-S)), allowing empty S.
If every forward arc were available, root path counts would give
((K_0 E+K_1 O,K_1 E+K_0 O)).
But x→y is missing: the forbidden concatenations consist of the
singleton left path at x followed by a path starting at y. Their
complement has left count D. If right-start profile is ((e_s,o_s)),
the exact return is
\[
 (E',O',D')=(K_0E+K_1O-D o_s,\ K_1E+K_0O-D e_s,\ h_V D).
\]
Here (K_0=K_1=:K) for nonempty V. Expand both H factors by odd-cycle
packs. Treat unused vertices as singleton cycles; for each full pack,
S can be the union of an arbitrary subset of its odd components. The
signed sum over these choices is ((1-1)^c=0) for c≥1 components.
Thus the continuation depends only on ((E+O,D)).

With a right triangle, (K=6) and start and end profiles both have
((e,o)=(1,2)). Repeated right triangles multiply D by 3 and send
(E+O\mapsto12(E+O)-3D). Two initial sources with common S=E+O and
D but different outputs after the first gluing retain identical E,O,D
at all later steps. Their Hamilton-count difference multiplies by 3
at every later step, since their cross-cycle contributions are equal.

### Three fully reusable shifted counts

The next three strong order-eight sources share \(H=657\), \(D=153\),
and \(E+O=1157\); their \(E\) coordinates are \(540+j\).

| j | Root | Outgoing row masks |
|---|---:|---|
| 0 | 7 | 28,37,56,194,106,201,135,23 |
| 1 | 6 | 26,44,97,148,102,201,139,23 |
| 2 | 7 | 28,41,82,164,106,197,139,23 |

One triangle gives \(H=5365+2j\) and the same retained data
\((E',O',D')=(6636,6789,459)\). A second gives \(56217+6j\).
For every \(k\ge1\), the same formula and induction give three
strong outcomes at step \(k\) whose successive differences are
exactly \(2\cdot 3^{k-1}\). The numerical witnesses at steps 1 and
2 are independently checked in `check_bridge.py`.

This realizes a single ternary digit at arbitrary height. It does
not yet give independent digits at several heights in one tournament.
The needed operation must deposit a fresh small change locally without
changing E, O, D or interfering with any previously deposited digit.
