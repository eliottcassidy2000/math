# An exact graded crossing bridge for arbitrary module counts

This construction uses an auxiliary count without identifying it
with the ordinary Hamiltonian count. It gives positive exact
families and a target inverse, but not cofinite coverage.

For a tournament \(T\), put \(h_T=H(T)\) and let
\[
 Q_T=\sum_{U\subseteq V(T)} H(T[U])H(T-U)
\]
be its ordered weak two-path-cover count; the two slots may be
empty, with \(H(\varnothing)=1\). If \(T\) is nonempty, exactly
half of its weak two-slot covers have an even first-slot size:
use the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569)
for both slot counts and color the resulting odd-cycle-or-singleton
atoms by their slots. For every nonempty underlying set of atoms,
the signed color sum vanishes, since each atom has odd size.
Thus the even and odd weighted subset sums are both \(Q_T/2\).

## A literal crossing construction and its exact formula

In a tournament \(A\), mark a nonempty set \(X\). Let \(c_A\) count
weighted selected paths starting in \(X\), with the unused set
weighted by its Hamiltonian count. Let \(E_A\) count the same
weighted paths with both endpoints in \(X\) and even cardinality,
and define the signed grade
\[
                            d_A=h_A-2E_A.
\]
For any nonempty tournament \(S\), start with the ordinal join
\(A\to S\to b\), where \(b\) is a single vertex, and reverse
exactly the arcs \(x\to b\) for \(x\in X\). The new count is
\[
                         d_A h_S+c_A Q_S.                 \tag{1}
\]

**Proof.** All odd cycles that leave either component must contain
\(b\), since the only backwards edges are \(b\to X\).
The unchanged odd-cycle packings have weight \(h_Ah_S\).
A new odd cycle is \(b\), a nonempty path in \(A\) starting in \(X\),
and then an optional path in \(S\). With \(S\)-part nonempty its
parity must agree with the selected \(A\)-path's parity. Parity
balance gives \(Q_S/2\) choices by weight, minus \(h_S\) when the
\(A\)-path is even to exclude the empty slot. The empty \(S\)-part
reinstates that \(h_S\) term exactly if the even \(A\)-path ends
outside \(X\). This leaves one subtraction of \(E_Ah_S\).
The crossing cycle carries weight two, and all untouched cycles
are counted on the two complements; hence (1). If \(X\) meets the
first strongly connected component of \(A\), the resulting
tournament is strong: \(b\) reaches all of \(A\) and \(S\), while
all vertices reach \(b\) via \(S\). ∎

An untagged ordinal block \(D\) inserted after \(A\) changes
\((d_A,c_A)\) to \((d_Ah_D,c_AQ_D)\). The even path whose
endpoints are both tagged cannot enter \(D\), while the selected
tag-start path and the complementary path split \(D\) freely.
Therefore any marked seed with \(d_A=1,c_A=\lambda\) gives, for
**arbitrary independently chosen nonempty modules** \(D,S\),
\[
                         h_Dh_S+\lambda Q_DQ_S.           \tag{2}
\]
Positive fixed seeds include:

| \(\lambda\) | Outgoing rows | Tag mask | \(h_A\) | \(E_A\) |
|---:|---|---:|---:|---:|
| 1 | \(0\) | 1 | 1 | 0 |
| 6 | \(2,4,1\) | 3 | 3 | 1 |
| 11 | \(2,4,9,3\) | 3 | 5 | 2 |
| 12 | \(2,4,9,3\) | 6 | 5 | 2 |
| 29 | \(18,4,9,19,6\) | 18 | 11 | 5 |

Every displayed seed is strong or a singleton, so each resulting
host is strong. Padding the seed on the right by untagged
singletons permits any further dyadic multiple of \(\lambda\).
The ordinary range (2) is therefore a genuinely graded pairing:
its second coefficient need not be one of the dyadic choices of
a transitive marked singleton.

There is also a symmetric \(\lambda=3\) without using an
independent fixed seed. Wrap \(D\) with singleton ends
\(x\to D\to y\to x\), marking both \(x,y\). Its ordinary count
is \(h_D+Q_D\): all cycles use the backward edge and an odd-size
selected middle path, with total weight \(2(Q_D/2)\). Its even
double-end count is \(Q_D/2\), since the direct path \(y\to x\)
supplies the empty-middle contribution, and all other admissible
paths go \(x\to U\to y\) with nonempty even-size \(U\).
Its unsigned start mass is \(3Q_D\): starts at \(x\) contribute
\(2Q_D-h_D\), and starts at \(y\) contribute \(Q_D+h_D\);
pick the middle subset and decide which endpoint follows, in
accord with the fixed orientations. Thus (1) yields
\[
                         h_Dh_S+3Q_DQ_S
\]
in a host with only three extra vertices total.

## An explicit two-parameter absolute selector

For \(D\) transitive of width \(a\ge1\), take \(S\) transitive
of width \(b+2\), but reverse its last-to-first extreme arc,
with \(b\ge1\). Then \(h_D=1,\ Q_D=2^a\), and
\[
                  h_S=1+2^b,\qquad Q_S=2^{b+2}+2\cdot3^b-2.
\]
For the latter read, distribute the two extremes into distinct
slots (weight \(2^{b+1}\)) or the same slot (weight
\(2(2^b+3^b-1)\) after choosing the middle subset).
Thus the symmetric bridge gives the exact strong count
\[
 N=1+2^b + 3\cdot2^a(2^{b+2}+2\cdot3^b-2).
\]
A prescribed odd target has a complete \(O(\log N)\) candidate
test in this family: try each \(b\) with \(2^b<N-1\) and divide
\(N-1-2^b\) by \(3(2^{b+2}+2\cdot3^b-2)\), accepting exactly
a positive integral power of two at least two. This is
equality at the absolute target, not a modular selector.
More generally for fixed input grades \((h_D,Q_D)\) and
\((h_S,Q_S)\), the auxiliary multiplier in (2) is determined
by one exact division and then checked against the physical
seed catalog. None of these inverse statements entails
coverage of arbitrary odd targets.

## Why a fixed-height apex cannot freely supply the auxiliary grade

Fix an ordinary apex height \(K=H(A+v)\), with \(v\) beating
the tagged set of \(A\). Factor \(A+v\) into its ordered SCCs.
Let \(C\) contain \(v\), set \(D=C-v\), and write
\(A=L\to D\to R\). All vertices of \(L\) are untagged; all
of \(R\) are tagged. The weighted start profile is therefore
\[
                     c_A=h_L[(c_D+h_D)Q_R-h_Dh_R].         \tag{3}
\]
A selected tagged-start path omits \(L\), leaving it in the
other slot. If the selected path starts in \(D\), distribute
\(R\) over the two slots freely. If it starts in \(R\), the
whole \(D\) lies in the other slot and the first \(R\)-slot
is nonempty, of weight \(Q_R-h_R\). These are all possibilities.

After deleting all singleton SCCs except \(v\), only finitely
many normalized marked graphs remain at fixed \(K\).
A nontrivial strong SCC of order \(n\) has at least \(n\)
full paths by rotating a directed Hamilton cycle; all
component counts multiply to \(K\). Thus their orders and
number are bounded purely by \(K\). The discarded left
singletons leave (3) unchanged. The discarded right
singletons only double \(Q_R\) once each. Thus **for fixed
\(K\), every possible auxiliary \(c_A\) lies in a finite
union of exact affine dyadic rays**
\[
                               c_A=\alpha 2^m-h_A .
\]
Here \(\alpha,h_A\) depend on the normalized marked graph,
and \(m\ge0\) is the number of restored right singletons.
Every point of each ray is attained by appending those
singletons; their apex contacts follow from the ordinal
order. Therefore a proposed compiler keeping its ordinary
apex height fixed cannot transmit arbitrary auxiliary
integers through this mark. The arbitrary-module positive
construction above instead keeps the *signed grade* fixed
and allows ordinary heights to vary.

## Slot-neutral choices and a larger exact family

The strong four-vertex module \((2,4,9,3)\) has grades
\((h,Q)=(5,32)\); these share the auxiliary grade of a
transitive five-vertex module, whose ordinary grade is 1.
The strong six-vertex module \((38,4,8,19,39,14)\) similarly
has grades \((19,256)\), matching the auxiliary grade of
a transitive eight-vertex module. Hence in an ordinal
deletion one may choose independently, in each prescribed
slot, to carry a nontrivial ordinary grade or carry 1
with exactly the same auxiliary grade, not just its phase.
Both weak path-cover counts and ordinary path counts
multiply under ordinal joins.

For every fixed seed multiplier \(\lambda\) in (2), this
gives the exact strong family
\[
                            5^i19^j+\lambda 2^k,
 \qquad i,j\ge0,\quad k\ge5i+8j,\quad k\ge1.
\]
Use \(i\) copies of the four-vertex module, \(j\) copies of
the six-vertex module, and \(k-5i-8j\) singleton components
on the untagged side of the bridge. The host stays strong
even though these internal inputs are ordinal sums.
An exact target test tries the finitely many pairs
\(5^i19^j<N\) and checks whether
\((N-5^i19^j)/\lambda\) is an integral power of two of
the required exponent. There are \(O((\log N)^2)\) such
pairs for fixed \(\lambda\). This supplies a physical
alternative to guessing an auxiliary second count against
an already fixed ordinary count. It does not claim that
the tested target language is cofinite.

## A universal grade congruence and a sharp neutral edge

Let \(Q_T(j)\) count ordered weak covers by \(j\) path slots
(\(j\ge2\)). Expanding each slot by the odd-cycle packing
formula colors every odd atom independently by its assigned slot.
A packing with \(a\) atoms, of which \(r\) are nontrivial cycles,
thus contributes \(2^r j^a\). As every atom has odd size,
\(a\equiv |T|\pmod2\), and consequently
\[
                      Q_T(j)\equiv j^{|T|}H(T)\pmod{j^2-1}.
\]
In particular equal \(Q=Q(2)\) at fixed order implies equal
ordinary count modulo three. Since ordinary counts are odd,
any nonzero same-order \(Q\)-neutral change of \(H\) is at least
six in absolute value. Thus a \(Q\)-neutral module switch
cannot realize a \(+2\) jump in the grade-one bridge.

The bound is sharp even for one literal arc. The lower
rows \((244,81,210,39,8,150,56,90)\) and their reversal
at \((0,5)\) have grades \((383,5960)\) and \((389,5960)\).
Use either module in a fixed grade-one bridge with multiplier
\(\lambda\), appending \(k\ge0\) singletons to the untagged
module if desired. The resulting strong hosts differ only
on that same arc, and their exact counts are
\[
             383+5960\lambda 2^k,\qquad389+5960\lambda2^k.
\]
The literal replay checks both source grades and small
uncompressed grafts. This is an exact \(+6\) family preserving
its auxiliary contribution, not a \(+2\) construction.
