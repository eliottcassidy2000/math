# A bridge of smaller tournaments needs its third read on the second pass

Use the following fixed five-vertex quotient, with two vertices replaced homogeneously by arbitrary nonempty tournaments \(D,S\). The other three vertices are \(x,y,b\). Its arcs, as out-neighbor lists, are
\[
x:\ D,S;\quad D:\ y,S,b;\quad y:\ x,S;\quad S:\ b;\quad b:\ x,y.
\]
The quotient is strong, so the resulting tournament \(\mathcal B(D,S)\) is strong even when either input is not. In a strong quotient, a cycle through a substituted vertex supplies a route between any two of its replacement vertices.

Let \(Z_T(q)\) count ordered \(q\)-slot covers by directed paths, allowing empty slots, and let \(W_T(i)\) be the corresponding count with all slots nonempty. Write \((h,Q,R)=(Z_D(1),Z_D(2),Z_D(3))\) and \((h',Q',R')\) for \(S\).

## Exact compression of the two smaller parts

Compress every maximal run of \(D\) and \(S\) to its symbol, separately in each slot. In \(q\) slots, \(D\) can occur at most \(q+1\) times: a run either starts a slot or is entered from the singleton \(x\). Symmetrically an \(S\) run either ends a slot or exits to the singleton \(b\). The skeleton records the order and slot of every run. Filling the \(D\)-runs independently of the \(S\)-runs gives
\[
  Z_{\mathcal B(D,S)}(q)
   =\sum_{1\le i,j\le q+1}C_q(i,j)\,W_D(i)W_S(j).
\]
A fixed skeleton is filled by the ordered nonempty path covers, not by arbitrary products of scalar counts. The coefficient matrices at \(q=1,2\) are
\[
 C_1=\begin{pmatrix}13&6\\6&3\end{pmatrix},\qquad
 C_2=\begin{pmatrix}104&100&28\\100&92&26\\28&26&8\end{pmatrix}.
\]
For a finite exact computation, give the \(i\) occurrences of \(D\) and \(j\) occurrences of \(S\) distinct names. Use the quotient arcs between differing symbols and no arcs between identical symbols, count the ordered covers by subset/endpoint DP, and divide by \(i!j!\). Every unlabeled skeleton has exactly that many namings. The included replay performs this count.

Converting \(W_1=h,\ W_2=Q-2h,\ W_3=R-3Q+3h\), the ordinary count and next read become
\[
\begin{aligned}
 H(\mathcal B(D,S))&=hh'+3QQ',\\
 Z_{\mathcal B(D,S)}(2)&=-6(hQ'+Qh')+8QQ'+2(QR'+RQ')+8RR'.
\end{aligned}
\]
The negative signs in the second read arise only from the change of basis: the skeleton kernel itself is nonnegative. At arbitrary \(q\) the top coefficient \(C_q(q+1,q+1)\) is positive: use \(D,S,b,x,D,y,S\) as one slot and \(D,S\) as each of the other slots. Thus the uncompressed hierarchy uses the next slot count at each level, not simply repeat applications of the scalar law.

## An actual collision proves that two reads do not close under recursion

The two eight-vertex inputs with outgoing row masks
\[
(2,104,147,213,67,29,165,51),\qquad
(110,184,210,212,1,156,178,17)
\]
both have \(h=345,\ Q=5640\), but \(R=45225\) and \(44505\). Use a singleton for \(S\), so \((h',Q',R')=(1,2,3)\). The resulting strong tournaments both have \(H=34185\), yet their two-slot counts are \(1352400\) and \(1332240\). Applying the same bridge once more with a singleton gives distinct ordinary counts \(8148585\) and \(8027625\). The replay checks the input profiles, both exact literal 12-vertex tournament outputs, and the final formula.

This is a precise obstruction to any iteration of this bridge that retains only the two input reads \((H,Z(2))\). It does not rule out a marked invariant with more information or a restricted class of inputs where the third read is a function of the first two.
