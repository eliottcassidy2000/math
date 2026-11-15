# A reverse-prefix bound for carrier rows

Let \(G(u)=cu+kp(u)+\phi[u\ {\rm odd}]\), with fixed \(c,k>0\), \(\phi\in\mathbb R\), and \(p\) as in the carrier packet. The physical examples have integer coefficients. Write \(V=2^\ell\). An ordered collision whose largest differing binary position is \(\ell-1\) can be counted up to a factor 2 by writing
\[
 u=Vw+b,\qquad u'=Vw+b',\qquad b\ge V/2>b'.
\]
The full [row collision lemma HYP-9752](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/hypotheses/HYP-9752-row-collision-lemma-a-carrier-row-cu-plus-kp-u-has-linear-energy-which-would-give-the-chan-pak-spectrum-positive-lower-density.md) asks for \(O(V/(1+w))\) collisions here. We prove a stronger restriction than the initial suffix-pair uniqueness in the large-\(w\) regime, but not that full conjecture.

**Result.** For fixed \((\ell,w)\), let \(N_{\le B}(\ell,w)\) count collisions with \(0<|d|\le B\), where \(d=\operatorname{rev}_\ell(b)-\operatorname{rev}_\ell(b')\). Put
\[
 D=\min\bigl(B,V,CV/(1+w)\bigr),\qquad C=3+c/k+|\phi|/k.
\]
Then \(N_{\le B}=0\) if \(D<1\), and otherwise
\[
 N_{\le B}(\ell,w)\ll D^2(1+\log(V/D)).
\]
There are two consequences at dyadic \(N=2^L\), counting inputs in \([N/2,N)\):
- All collisions with \(\ell\le(2/3-\epsilon)L\) contribute \(O_{c,k,\phi}(N^{1-3\epsilon}\log N)=o(N)\).
- All collisions with \(|d|\le B\), at any \(\ell\), contribute \(O_{c,k,\phi}(\sqrt N B^{3/2}\log N)\); in particular reverse displacement \(\le N^{1/3-\epsilon}\) gives \(o(N)\).

Additionally, anywhere in \(0\le u,u'<N\), collisions with \(|u-u'|\le H\) contribute \(O(H^2\log N)\), so ordinary displacements \(\le N^{1/2-\epsilon}\) are negligible on the linear-energy scale. These results leave large ordinary split levels with genuinely large reverse displacement. They do not imply linear total energy or positive density.

## Proof of the reverse-prefix count

The published split gives
\[
 k d w+k(w\bmod 2)=c(b-b')+k(p(b)-p(b'))+\phi([b\ {\rm odd}]-[b'\ {\rm odd}]).
\]
The left reverse displacement \(d\) is odd (the most significant differing bit becomes the least significant bit) and nonzero. Since \(0\le p(v)\le v\), the equation gives \(|d|(1+w)\le CV\).

Now take \(Q\) the smallest power of two at least \(\lfloor D\rfloor\), with \(Q\le V\). Work in the *reversed* arguments \(x=\operatorname{rev}_\ell b,x'=\operatorname{rev}_\ell b'\), which differ by at most \(Q\). Their quotients by \(Q\) agree or differ by exactly one. In the first type they share a prefix with residue width \(A=Q\). In the latter types, expose the entire single carry or borrow: if the smaller quotient ends in \(j\) ones, their last \(j+1\) bits are respectively \(2^j-1,2^j\), their common prefix is unrestricted, and \(A=Q2^{j+1}\le V\). Direction and \(j\) determine all extra exposed bits of both words. Thus for each of \(O(1+\log(V/Q))\) types and each of at most \(Q^2\) base suffix pairs,
\[
 x=A z+B_0,\qquad x'=A z+B_1 ,
\]
where \(B_0\ne B_1\) are fixed and \(z\) ranges over common prefixes.

Let \(A=2^a\), \(s=\ell-a\). Undoing the reversal,
\[
 b=2^s\operatorname{rev}_a(B_0)+t,\quad
 b'=2^s\operatorname{rev}_a(B_1)+t,\quad
 t=\operatorname{rev}_s(z).
\]
Equivalently the full inputs are \(2^sH_0+t,2^sH_1+t\) with fixed distinct
\(H_i=2^a w+\operatorname{rev}_a(B_i)\). If \(s=0\) there is only one \(t\). Otherwise use the direct binary split at their **common low suffix** \(t\). The terms \(p(t)\) and the explicit parity of the full inputs cancel. The difference \(G(u)-G(u')\) is
\[
 c 2^s(H_0-H_1)+k\{p(H_0)-p(H_1)
 +(2^s-1-\operatorname{rev}_s(t))(H_0-H_1)
 +[t<2^{s-1}]((H_0\bmod 2)-(H_1\bmod2))\}.
\]
On each of the two halves for \(t\) its dependence on \(\operatorname{rev}_s(t)\) is strictly affine with nonzero slope \(-k(H_0-H_1)\). Hence each of the previous types and suffix pairs contributes at most two actual coincidences. This proves the bound, with no assumption on \(w\) beyond being fixed.

## Summing the bound

In a dyadic input annulus both inputs have leading bit 1, so \(\ell\le L-1\), \(N/(2V)\le w<N/V\), and there are \(O(N/V)\) eligible prefixes at one level. For levels \(V\le N^{2/3-\epsilon}\), the contribution per level is
\[
 O_{c,k,\phi}((N/V)(V^2/N)^2\log N)
 =O_{c,k,\phi}(V^3\log N/N).
\]
A geometric sum gives the first assertion. For reverse displacement at most \(B\), use the per-level bound
\[
 O_{c,k,\phi}((N/V)\min(B,V^2/N)^2\log N).
\]
Sum upward on the side \(V\le\sqrt{NB}\), and downward for \(V>\sqrt{NB}\). Both geometric sums are \(O_{c,k,\phi}(\sqrt N B^{3/2}\log N)\).

Finally, in the ordinary displacement condition \(b-b'\le H\), both \(b,b'\) lie within \(H\) of the midpoint \(V/2\). There are \(O(\min(H,V)^2)\) such suffix pairs. Each gives at most two \(w\) by the previously known top equation, one per parity of \(w\). Summing over levels proves the third bound.

**What remains.** A single fixed positive power decay in \(w\) beyond a bound \(O(V)\), for example \(N(\ell,w)\ll V/(1+w)^\eta\) for any \(\eta>0\), would already imply linear energy. The new squared-displacement bound only gives this kind of control far out in the prefix and cannot replace that missing estimate.


## Extremal reverse distance one

If \(c\ge 2k+|\phi|\), as in the three calibration rows, \(d=-1\) is impossible. Indeed \(\operatorname{rev}_\ell b\) is odd, and adding 1 to it makes \(\operatorname{rev}_\ell b'\). If the odd reversed word ends in \(j\ge1\) ones, the two ordinary words have top blocks \(1^j0\) and \(0^j1\) and a common low suffix. Their difference is \(h=(2^{j+1}-3)2^{\ell-j-1}\), and \(b'<2h\). The right side of the top equation is then at least \((c-|\phi|)h-kb'>0\), but its left side is \(k(-w+[w\ {\rm odd}])\le0\). At \(\ell=1\) the negative direction is already excluded by parity and order.

For \(d=1\) and \(\ell\ge2\), the reversed pair is even to odd with no carry. Hence \(b=b'+V/2\); put \(S=V/2\), \(t=b'\), \(z=\operatorname{rev}_{\ell-1}(t)\). The exact condition is
\[
 w+[w\ {\rm odd}]=\frac c k S+S-1-z+[z\ {\rm even}].
\]
This follows from the shared suffix of length \(\ell-1\), high parts \(2w+1,2w\), and \(p(2w+1)-p(2w)=-w-[w\ {\rm odd}]\). The map \(z\mapsto z-[z\ {\rm even}]\) on \([0,S)\) has multiplicity at most two, so the entire distance-one sector for any \((\ell,w)\) has multiplicity at most two in this orientation.

The [replay](verify_rows.py) checks every literal collision for the three calibration rows through \(u<2^{15}\), including derivative parity, top equations, and uniqueness on each half of every reverse carry chart.


The reverse-displacement bound extends to the entire initial interval \(0\le u,u'<N\), at the cost of a second logarithm. For collisions within equal dyadic input annuli, sum \(O(\sqrt{2^j}B^{3/2}\log(2^j))\) geometrically over their maximal dyadic cutoff \(2^j\). Otherwise the largest differing bit is the leading bit of the larger input, so the top prefix is \(w=0\). Each of these \(O(\log N)\) levels contributes \(O(B^2\log N)\) by the reverse-prefix lemma (and less if \(B>V\)). Since \(B^2\le\sqrt N B^{3/2}\) for \(1\le B\le N\), the whole initial interval has at most \(O(\sqrt N B^{3/2}\log^2 N)\) such collisions.

