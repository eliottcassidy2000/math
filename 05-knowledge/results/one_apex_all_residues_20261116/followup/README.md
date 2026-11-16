# Logarithmic residue controls and the one-gap fork obstruction

Follow-through to the [one-apex construction](../README.md). The underlying congruence density is inherited from Session 25. This note improves the control size and puts a precise limit on a proposed exact child rule.

## 1. Three pairs make one binary place

Use the same notation \(M=2^am\), \(m\) odd, \(D=\max(2,a+1)\), \(L=\operatorname{ord}_m(2)\) (set \(L=1\) for \(m=1\)). In a pair of well-separated losing marks, swapping phases \((0,d)\) to \((d,0)\) adds
\[
\delta_d=\tfrac12(2^d-2^{-d})\pmod m.
\]
The external contributions depend only on phase sums and do not change. All affected gaps are long, so there is no dyadic change. Give a switch three pairs with labels \(d+1,d,d\). In its off setting swap the two \(d\)-pairs; in its on setting swap only the \(d+1\)-pair. The increment is
\[
\delta_{d+1}-2\delta_d=3\,2^{-d-2}\pmod m.
\]
Let \(g=\gcd(m,3)\), \(q=m/g\), and \(k=\lceil\log_2 q\rceil\). For \(d=0,\ldots,k-1\), these are the binary weights \(2^{k-1-d}\) times the common generator \(3\,2^{-k-1}\) of the multiples of \(g\). Since \(2^k\ge q\), any target in the prescribed coset is attained by one ordinary binary expansion. The old short binary pairs and ternary coset gadgets are unchanged. Every spine length at least \(B(D+L)\) works, where
\[
B=2\max(a-1,0)+6[3\mid m]+6k+1.
\]
The final \(+1\) is the terminal losing mark, which makes the tournament strong. An arbitrary fixed suffix after a long gap is still allowed; external cross terms do not change, so the leading-bit \(O_M(\sqrt N)\) placement argument also survives. This gives logarithmically many marks in the odd modulus, not logarithmic total order unless the period is also controlled.

Example: modulo 8191 the explicit uniform construction now uses 79 marks and a spine of 1185 positions; the earlier square-root-tail bound used 183 marks and 2745 positions. These are constructive bounds, not claims of minimality. The verifier constructs every admissible class for \(M\le140\), checks a second fixed order, and includes larger mixed and odd-power moduli.

## 2. The arithmetic of a shared gap

Write a strong one-apex spine word as \(UV\), with bit 1 for an apex win. Set
\[
X=X(U),\quad Z=Z(V),\quad K=H(U)+H(V)-1+[U_{\rm last}V_{\rm first}=10].
\]
Here \(X\) is the ordinary binary value and \(Z\) the value of the reversed complemented word. The concatenation rule in the parent note gives \(H=K+XZ\). Adding a zero or one at the gap gives deviations from \(2H\)
\[
s_0=X+\epsilon_0-K,\qquad s_1=Z+\epsilon_1-K,
\quad\epsilon_0=[11],\quad\epsilon_1=[00].
\]
A boundary gap of a strong nontrivial word cannot provide both children, since at least one insertion would merely pad with a source or sink. At an interior gap giving both signs, \(|X-Z|\le3\) and \(K\le\min(X,Z)+2\). If the cut reads \(00\), the alternatives are \((X,Z)=(K-1,K)\) or \((K+1,K-2)\). If it reads \(11\), use reverse-complement. Otherwise \(X,Z\) are \(K-1,K+1\) in some order. Necessarily
\[
H\in\{K^2,\ K^2-2,\ K^2+K-1\}.
\]
So there are only \(O(\sqrt N)\) possible scalar parents up to \(N\) for this particular two-child rule. No converse claim follows from the arithmetic.

## 3. Even allowing different next gaps, the fork cannot renew immediately

**Claim.** For a shared-gap fork as above with \(H\ge63\), each physical child has at most one legal \(2H'-1,2H'+1\) child under any one-vertex induced extension staying in the one-apex class. The second step may use any gap.

First \(X,Z\ge7\): if their smaller value were at most 6 then \(XZ+K\le6\cdot9+8=62\). For any word, the zero-insertion deviation is nondecreasing as the cut moves right, constant when crossing a 0. For a crossed 1 with adjacent old exterior bits \(p,q\) and remaining right suffix \(V'\), its jump is \(X_{\rm before}+1+q-p+Z(V')\); this is at least \(X_{\rm before}+Z(V')\). To see it, the respective forward values are \(X_{\rm before}\) and \(2X_{\rm before}+1\), the \(K\) coordinate decreases by \(Z(V')\), and the \(11\) correction changes from \(p\) to \(q\). Take \(p=0\), \(q=1\) at empty exteriors. Reverse-complement implies one-insertion deviations are nonincreasing, constant across 1s.

In \(U0V\), the one-insertion deviations immediately before and after the new zero are respectively at least \(Z-2\) and at most \(3-X\). Before it the suffix value is \(2Z+1\) and the seam coordinate is \(K+\epsilon_0\); after it the suffix value is \(Z\) and the seam coordinate is \(K+X+\epsilon_0\). The bounds use \(|Z+\epsilon_1-K|=1\). These consecutive cuts straddle both legal signs strictly, so no gap gives a legal child by adding a 1. Another 0 in the same run gives \(s_0-2\epsilon_0\). Leaving this run crosses a 1, if there is one: to the left the right suffix still contains \(V\) as a suffix, giving \(Z(V')\ge Z\); to the right the left prefix contains \(U\) as a prefix, giving \(X_{\rm before}\ge X\). Each such jump is therefore large enough to miss both signs, as the deviation on the run is between \(-3\) and \(1\). For \(U1V\) apply reverse-complement; the surviving deviation is \(s_1-2\epsilon_1\).

This holds for the unpointed tournament too. By the archived [two-apex classification](../../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md), a strong tournament of order at least four with two eligible apices is the cyclic blow-up with two singleton parts and a transitive part. Its two pointed words are \(10^{n-1}\) and \(1^{n-1}0\) (where \(n\) is spine length). For \(10^{n-1}\) the gap after its leading 1 gives a fork only at \(H=5\); an interior zero-run gap after \(j\ge2\) bits needs \(2^{n-j}-2^{j-1}=2\), which gives \(H=9\). The other word is its reverse-complement. If one of the large fork children had two apices, deleting the inserted bit in its pointed word would give a parent of this same form or a transitive parent, impossible here. Each child thus has a unique old apex, which any next one-apex extension must retain: deleting another old vertex would already give a transitive deletion in the child, and the new vertex cannot be the apex because the child is nontransitive. This completes the claim.

The one-apex restriction matters. This does not forbid distinct gaps for the *initial* fork, editing old arcs, multi-vertex operations, or moving to a larger feedback core. The published 121 fork is now diagnostic: 241 and 243 both have only their minus continuation at the next induced one-apex step.

## Replay

Run python3 verify.py. The distinct exact bitwise count formula checks every strong spine through length 13. It returns eight large shared-gap fork presentations, with 14 next children having exactly one legal continuation and two having none. The finite census checks the algebra and the construction; it is not the proof of the all-order bound.
