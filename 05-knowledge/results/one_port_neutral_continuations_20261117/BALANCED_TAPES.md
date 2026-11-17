# A genuinely growing unit switch with one live port

Keep the triangle \(D=(b,c,y)\), oriented \(b→c→y→b\), and
the two ports x,a from the original switch: x→b,y; c→x;
a→b,c; y→a. Minus has x→a, plus has a→x.
Let L and R be transitive tapes in an ordinal chain \(L→D→R\).
**Only a has nonuniform contacts.** The inert port x loses to
all L and beats all R. A left gate is a vertex that a beats,
and a right gate is one that beats a. Nongates respect the
ordinal orientation around a. Remove any initial left
nongates and terminal right nongates: they contribute only
global sources or sinks and do not change either count.

Let \(P_L=\sum_{\text{left gates }i}2^{|L|-i}\) and
\(P_R=\sum_{\text{right gates }j}2^{j-1}\), indexing inside the
tapes from one in forward order. A nonempty left tape thus
starts at its highest gate; a nonempty right tape ends at
its highest gate; their low bits describe the vertices
nearest to the triangle. The exact switch difference is
\[
 H_+-H_-=2+9P_R-6P_L+(P_R\bmod2)-2(P_L\bmod2).
\]
So **every pair \((P_L,P_R)=(6t,4t)\), \(t\ge0\), gives gap two**,
and these are exactly the trimmed address pairs giving gap two.
Indeed modulo 3 both low bits must vanish, and then
\(3P_R=2P_L\) gives the stated divisibility and parameterization.

No tape can be visited by a Hamilton path in more than two
runs: removing a leaves \(L→(D\cup\{x\})→R\), an ordinal list.
The one live port supplies one backward crossing per tape.
In contrast to the common-contact architecture, there is
no third free run.

## Proof of the gap formula, in a general fixed tile

For any nonempty D attached the same way to the inert port
relative to external tapes, let \(I_x,I_a\subseteq D\)
be their arrival sets, i.e. the vertices beating each port.
Contract consecutive ports in the affected paths to q.
Minus reads arrival set I_x and departure complement I_a on D;
plus reads I_a and departure complement I_x.
Let \(K_\pm\) be their Hamilton counts on D+q; let \(Z_-\) be
the minus count ending at q and \(Y_+\) the plus count starting
at q. Count paths with both whole tapes on their ordinal sides
of q: the only lost paths are minus ending at q when the first
R vertex is a gate, and plus starting at q when the last L
vertex is a gate. These give the signed difference
\(K_+-K_-+rZ_--lY_+\), where r,l are the low bits.

Minus can also place some L vertices after q but cannot place
any R vertex before q: no R vertex has an exit from the right
tape. The first L gate after q determines the mandatory left
vertex, and any subset of later L vertices can follow it.
All other L vertices precede q. The core D vertices before q
must make an optional path ending in I_x; the remainder makes
an arbitrary path after the returning L vertices. The choices
therefore factor as \(C_L P_L\), where
\[
 C_L=H(D)+\sum_{\emptyset\ne S\subseteq D}
 H(D-S)\,\#\{\text{Hamilton paths of D[S] ending in }I_x\}.
\]
Plus can put some R vertices before q but cannot put any
L vertex after q. Exactly symmetrically the last R gate
before q is mandatory; earlier R vertices can be chosen
independently. The core D after q must start outside I_x
if nonempty, so these choices factor as \(C_R P_R\), where
\[
 C_R=H(D)+\sum_{\emptyset\ne S\subseteq D}
 H(D-S)\,\#\{\text{Hamilton paths of D[S] starting outside }I_x\}.
\]
This exhausts the contraction patterns. In the original triangle
\(K_+-K_-=2, Z_-=1,Y_+=2,C_L=6,C_R=9\);
the five values are independently checked in verify_balanced.py.
With one tape absent, the general formula visibly restricts a
fixed gap to at most two addresses: for the right tape,
\(C_R\ge H(D)\ge Z_-\), so \(C_R P+Z_-(P\bmod2)\) increases
strictly except possibly equal adjacent pairs. Likewise for L.

## Exact absolute counts and their scale

Let p be the binary carry checksum
\[
p(0)=0,\quad p(2v+1)=p(v),\quad
p(2v)=p(v)+v+(v\bmod2)\quad(v>0).
\]
Then at \((P_L,P_R)=(6t,4t)\),
\[
 H_-(t)=11+768t^2+138t+5[p(6t)+p(4t)],\qquad H_+(t)=H_-(t)+2.
\]
Thus this route supplies unbounded exact neutral switches,
not density one. The main quadratic coefficient outweighs
all checksum jumps: p is even by its recurrence, and induction gives \(0\le p(v)\le v\) (if v is odd, the bound improves to v−1). So
the jump is at least \(768(2t+1)+138-50t>2\).
Both branches together have \(\Theta(\sqrt X)\) distinct
counts through X, and the exact monotone center can be
searched by integer binary search.

For a derivation delete a, leaving the ordinal chain
\(L→C→R\) with four-vertex \(C=D\cup\{x\}\). Its outgoing masks
are (2,12,1,5) in order b,c,y,x. It has H=5, optional ordered
two-cover count 32, weighted arrivals (paths ending at y or x)
15, and weighted departures (paths starting at b or c) 12.
Inserting a by the positive odd-cycle packing identity yields
the bare minus count 11 plus contributions from nonempty
even paths in the deletion starting at an outcontact of a,
ending at an incontact of a, and weighted by H of the complement.
Each weighted path contributes twice because it closes an odd cycle.

* A cross path from a left gate to a right gate contributes
  the optional core factor 32 times the product of binary
  gate weights. The parity split is exact: the signed optional
  two-cover count of a nonempty C vanishes, since one can
  allocate the odd blocks of a positive cycle packing
  to either slot. Sum to get \(32P_LP_R\).
* Paths from L into C contribute \(15P_L\), and from C into R
  contribute \(12P_R\). In each case a gate has at least one
  vertex available to choose freely in its tape between that
  gate and C, because neither low bit is a gate. That free
  choice balances parity exactly.
* Within one tape, one chooses endpoints of opposite contact
  type in the indicated order. At distance d an endpoint pair
  contributes \(2^{d-1}+[d=1]\) before weighting by H(C).
  With low bit nearest the core, the low endpoint is a
  nongate and the high endpoint is a gate; summing over both
  endpoint positions is the stated recurrence for p. This gives
  \(5[p(P_L)+p(P_R)]\).

Equivalently \(p(6t)+p(4t)=p(t)+p(3t)+6t+2(t\bmod2)\),
making the effect of a binary carry under multiplication
by three explicit. None of these counts are asserted to
place an arbitrary external target.

## A general entropy barrier on this route

For any fixed D in the general gap law, a fixed gap constrains
\(C_RP_R-C_LP_L\) to a bounded set of integers (the parity bits
are bounded). The two positive coefficients then force
\(P_L\) and \(P_R\) comparable above fixed constants, and given
one address the other has only O(1) possibilities.

For both growing tapes, split a Hamilton path of the fixed
core \(D\cup\{x\}\) into two nonempty paths C0,C1.
The physical trace \(L_0 C_0 R_0\,a\,L_1 C_1 R_1\) has
at least \(P_LP_R\) realizations: choose the last R0 gate
and first L1 gate, and independently choose arbitrary allowed
vertices in their assigned tape portions. The two gates
and the subsets can be recovered from the path. It works
on either switch state since a is adjacent only to the tapes.
Thus even allowing all addresses, the exact fixed-gap subset
contributes O(sqrt(X)) distinct counts up to X for each fixed
tile D (and finite union of such tiles). The explicit triangle
case above attains that exponent. This is not an obstruction
to variable tiles or more active state.
