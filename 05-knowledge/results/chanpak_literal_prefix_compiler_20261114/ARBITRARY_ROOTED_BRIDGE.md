# A physical bridge that accepts arbitrary rooted tournament seeds

The [previous square-root inverse](https://github.com/eliottcassidy2000/math/blob/9717032b8bf0afc5371763d0f26b3810d0332f0d/05-knowledge/results/chanpak_terminal_seven_nonuniform_contacts_20261114/OUTER_DEPTH.md)
gave two integral coordinates, but not a graph. The new
interface below no longer requires the seed to be on the
minimum-triangle terrace. The same corridor matrices really act
on arbitrary rooted cores **after** a physical preparation.

## 1. Exact gate and three rooted readouts

For a tournament \(T\) with marked vertex \(v\), let
\(s_v(S)\) be the number of directed Hamiltonian paths of \(T[S]\)
starting at \(v\), and set
\[
 R_+(T,v)=\sum_{S\ni v}s_v(S)H(T-S),\qquad
 R_-(T,v)=\sum_{S\ni v}(-1)^{|S|}s_v(S)H(T-S),
\]
with \(H(\varnothing)=1\).
Append \(r\ge0\) forward transitive followers and a last vertex
\(w\). All old vertices and followers beat \(w\) except
\(w\to v\). The followers beat no old vertex. Mark \(w\).
Call this a bridge with width \(r\).

Every odd-cycle packing using \(w\) has one cycle that begins
\(w,v,\ldots\); its old part is a \(v\)-start path on some
\(S\ni v\) and its remaining vertices are an ascending subset
of the \(r\) followers. Other cycles lie wholly in \(T-S\).
For \(r\ge1\), exactly half of the follower subsets make the
cycle odd. The potentially missing return \(v\to w\) affects
only the one-vertex old part with no followers, an even
two-cycle already excluded. For \(r=0\), the old set must
instead have even size. By the
[odd-cycle packing identity](https://arxiv.org/abs/2307.05569),
which contributes a factor 2 for the new cycle, the bridge's
literal total is therefore
\[
 H'=\begin{cases}
 H+2^r R_+,&r\ge1;\\
 H+R_++R_-,&r=0.
 \end{cases}
\]

For the root-profile sums at \(w\), the singleton subset
\(\{w\}\) contributes \(H\), and every other path starting at
\(w\) immediately uses \(v\). Each follower independently
belongs to the selected slot or the other slot; sign-weighted
summing multiplies by \((1-1)^r\). Hence the simultaneous
updates are
\[
 \begin{array}{c|cc}
 &R_+'&R_-'\\\hline
 r\ge1&H+2^rR_+&-H\\
 r=0&H+R_+&-H-R_- .
 \end{array}
\]
No hypothesis is made about triangles inside the old tournament.

## 2. Exactly the old corridor matrices, on physically larger sources

A positive bridge always gives the ready signature
\[
 (H',R_+',R_-')=(a+b,a+b,-a),
 \qquad a=H(T),\quad b=2^rR_+(T,v).
\]
Here \(a=H(T'-w)\) and \(b=H(T')-H(T'-w)\) are literal deletion
coordinates. From such a signature take \(q\) zero-width bridges
and then a positive bridge of width \(r\).
The zero update on \((H,R_+,R_-)\) has characteristic polynomial
\(x^3-x^2-x-1\), determinant 1.
Just before the final positive bridge, the \((H,R_+)\)
readouts for \(q=0,1,2\) are respectively
\[
(a+b,a+b),\quad(a+2b,2a+2b),\quad(3a+3b,3a+4b).
\]
Consequently for every \(q\ge0\), the new deletion coordinates
are exactly \(D_{2^r}C_q(a,b)^t\) with the previously used
terrace kernel. This **does not** assert that a bare arbitrary
root with only its two deletion counts known is already ready.
The preparation supplies the missing rooted information.

For a broad explicit seed family, let \(C\) be any tournament
and let the initial \(v\) dominate \(C\). Then \(H=H(C)\),
\(R_+(T,v)=F_2(C)\), and the first positive bridge gives
\[
                           (a,b)=(H(C),\,2^rF_2(C)).
\]
The bridge tournament is strong: every vertex of \(C\) and
every follower lies on a triangle through \(v,w\).
When an already strong prefix is extended by another bridge,
strong connectivity persists, since all new vertices reach
the old root through the new last vertex and are reachable
from the old graph. Thus the construction outputs a strong
ordinary witness at the terminal seven seam.

These bridge operations may also be expressed using neutral
followers and a one-outneighbor root extension. The new point
here is the exact interface matching to the later corridor
decoder, and the use of nonterrace sources, not a claim that
those elementary operations have no precedent.

## 3. Choose a real source, or reject the whole bridge-source class

For a prescribed pair \((a,b)\) test a source \(C\) with
\[
            H(C)=a,\qquad b/F_2(C)\in\{2,4,8,\ldots\}.
\]
This is complete for the dominating-root bridge class above.
A finite exact search on the strong cores suffices:
\(H\) and \(F_2\) both multiply under ordinal sums; remove any
singleton components of \(C\) and absorb their powers of 2 in
the bridge width. If \(a>1\), the remaining strong component
counts are factors of \(a\) at least 3. For each possible factor
one only needs the strong classes up to the logarithmic order
bound from [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
Retain the pairs \((H,F_2)\) whose second coordinate divides
\(b\), and use exact divisor recursion on the two coordinates.
The class includes an empty \(C\) when \(a=1\); then
\(F_2(C)=1\) and \(b\) must be a power of two at least 2.
For a bank of named cores, the same test needs only the
saved pair list and literal arcs. Failure of that bank is
not general impossibility.

The old greedy square-root inverse can be made to recover any
one of these actual source cores (not to recover every target).
Let \((u_k,v_k)=(47,56)M^k\), \(A_k=u_k+2v_k\),
\(M=D_2C_1\). For any \(C\), at all sufficiently large \(k\)
we have \(a=H(C)\le A_k/2\) and \(F_2(C)<3A_k\).
Choose the least \(r\ge1\) with \(b=2^r F_2(C)\ge3A_k\);
then \(b<6A_k\). Set \(p=u_ka+v_kb\).
Since \(A_k/3\le v_k\le A_k/2\) and \(A_{k+1}\ge5A_k\),
\[
 A_k^2\le p<\tfrac72 A_k^2<A_{k+1}^2,\qquad 0<a<2v_k.
\]
So the target chooses the same \(k\), and the least positive
residue of \(a\) modulo \(2v_k\) is exactly the constructed
one. This works with unboundedly varied actual cores. It is
not a theorem that a core exists for every target-selected pair.

### A literal nonterrace inverse hit

Take the six-vertex core with outgoing bit masks
\[
                           (2,24,3,53,37,7).
\]
Its \(H=23\) and \(F_2=296\). Use the dominating root and four
followers, then reverse the last arc to the root. The **actual**
12-vertex prefix has outgoing masks
\[
(2046,3972,4016,3974,4074,4042,3982,3840,3584,3072,2048,1),
\]
where bit \(j\) specifies an arc to vertex \(j\).
Its deletion coordinates (last vertex marked) are
\((23,4736)\). Literal subset path DP gives its total 4759
and verifies \((R_+,R_-)=(4759,-23)\).

No terrace prefix can have those coordinates. For any
terrace positive-ending state the ratio \(b/a\) is
\(2^r\) times a number in \([1,4/3]\): use
\(t\ge2\) in \((f_q+g_qt)/(e_q+f_qt)\), or the seed ratio 1.
The only possible exponent for \(4736/23\) would be 7,
but its residual ratio \(37/23>4/3\).
Thus this is a genuine new physical source, not a disguised
terrace representative. One \(q=1,r=1\) stage produces
\((9495,19036)\); the \(q=7,r=2\) seam then reads
\[
                    H=7\cdot1512281=10585967.
\]
The cofactor is prime by trial division, and feeding 1512281 to
the previous greedy inverse gives the actual pair \((23,4736)\).

Run `python bridge.py` for the literal small-root interface
checks and the witness.
