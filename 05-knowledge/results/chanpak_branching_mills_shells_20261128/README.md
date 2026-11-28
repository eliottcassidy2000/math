# Literal +2 pairs in cubic shells, and a branching quadratic Mills theorem

2026-11-28. The [sparse-digit obstruction](../chanpak_mills_sparse_digit_rigidity_20261128) showed why fixed contact complexity cannot sustain a Mills chain. Here the words change with the target. The same symmetric contact split yields two positive constructions:

1. Every integer \(N\ge100\) has a pair of strong one-apex counts \(L,L+2\) in \((N,N+16\sqrt N)\), with one literal arc reversal between the witnesses. This gives genuine branching cubic Mills chains.
2. **Every odd \(h\ge5\) has two distinct strong one-apex counts \(k\) in its strict square shell \(h^2<k,\ k+1<(h+1)^2\).** Thus there are uncountably many \(A\in(5,6)\) for which all floors \(\lfloor A^{2^j}\rfloor\) are strong one-apex counts. The two square successors need not differ by two or be arc neighbors.

There is no prime claim and no assertion of exact coverage of arbitrary integers.

## 1. Exact binary contact laws

A contact word \(w\) means that the spine edges point forward and the apex beats a spine site iff its bit is 1. If the word starts at 1 and ends at 0, it is strong: every spine site can reach the apex through the last site and be reached from the apex through the first site. Let \(H(w)\) be the path count. Equivalently,
\[
H(w)=1+\sum_{\substack{i<j\\w_i=1,\ w_j=0}}\bigl(2^{j-i-1}+[j=i+1]\bigr).
\]
To count, classify paths with both an initial and terminal spine segment by the sites immediately to the left and right of the apex. Each inversion supplies \(2^{j-i-1}\) choices, distributing the intermediate sites between the two increasing segments. For active words there are also two empty-side paths and one path for each adjacent 01, whose total is one plus the number of adjacent 10; initial zeros and final ones are harmless ordinal padding. The all-padding case has count one. This is the [one-apex recurrence](../chanpak_contact_run_complexity_20261128) with the full inversion weights retained. Starting from \(H=1\), appending a bit \(d\) adds \((1-d)(x+b)\), with \(x\) the old binary prefix and \(b\) its last bit (zero for an empty prefix).

Write \(F(x)=H(\operatorname{bits}(x))\) for \(x\ge1\), and set \(F(0)=1\). Reversing and complementing a word gives the converse tournament up to relabeling, so it preserves \(H\). For \(u=\operatorname{bits}(x)\), \(v=\operatorname{rc}(\operatorname{bits}(y))\) and \(x\) even, the exact product-plus law is
\[
H(uv)=xy+F(x)+F(y)-1 .
\]
In fact prefix \(u\) contributes an extra \(x2^i\) at each suffix zero at position \(i\) from the left starting at zero, with no first-bit boundary correction since \(x\) is even. These suffix zeros encode \(y\). Inserting a switch yields the exact flip gap
\[
H(u1v)-H(u0v)=y+[y\text{ odd}]-x-[x\text{ odd}].
\]
Changing the switch from zero to one loses \(x+[x\text{ odd}]\), but each later zero at suffix position \(i\) gains \(2^i\) and the first adds one extra if it is zero.

The basic estimates used below are
\[
F(x)\le x+1-[x\text{ odd}],\qquad
F(2m)\ge m+1\quad(m\ge1).
\]
Indeed \(F(2m)=F(m)+m+[m\text{ odd}]\) and \(F(2m+1)=F(m)\) (including \(m=0\)). The upper inequality starts with equality at 1; doubling preserves its slack and appending a 1 strictly increases its slack once \(m\ge1\). In particular for positive even \(x\), equality \(F(x)=x+1\) occurs exactly at powers of two.

## 2. A literal +2 pair above any target, then a cubic tree

For \(N\ge100\), take the smallest positive even \(x\) with \(2x^2\ge N\), then \(y=x+2\). Both words \(u0v,u1v\) start at 1 and end at 0, and differ only at the middle apex contact. Their heights are \(L\) and \(L+2\), where
\[
L=2x^2+5x+F(x)+F(x+2)-1.
\]
Since \(F(x)\le x+1\), one has the convenient loose bound \(2x^2+5x+1\le L\le2x^2+9x+3\). With \(x<\sqrt{N/2}+2\), this gives \(N<L<L+2<N+17\sqrt{N/2}+31<N+16\sqrt N\). The number of vertices is \(O(\log N)\).

For any odd \(h\ge63\), set \(N=h^3+1\). The inequality \(16\sqrt{h^3+1}<3h^2\) puts both children and their upper floor endpoints strictly inside the cubic shell. The word \(1001100\) gives the seed \(h=63\). Its first pair is 252,965 and 252,967 on 20 vertices, and the lower child leads to 16,187,557,770,191,969 and 16,187,557,770,191,971 on 56 vertices. More generally for any fixed real \(c>2\), take \(N=\lfloor h^c\rfloor+1\); for sufficiently large \(h\), the pair and its upper endpoints fit strictly, because the shell width grows as \(c h^{c-1}\) and the error is \(O(h^{c/2})\).

For either construction the rooted intervals \([h_j^{1/c^j},(h_j+1)^{1/c^j}]\) are strictly nested and shrink to a point. The two children have disjoint closed intervals since their counts are distinct and odd. Choosing an infinite binary path therefore fixes a distinct \(A\) with every required exact floor, not merely a shared excluded endpoint. For integer \(c\) a fixed computable branch gives a computable \(A\) by rational enclosure of the roots. This is analogous to [Mills's construction](https://doi.org/10.1090/S0002-9904-1947-08849-2) with a different admissibility predicate.

## 3. Why square shells can also branch

A generic square-root interval bound does not fit cleanly into a square shell. Instead center the two contact flanks at the parent height. For \(t\) odd with \(0<t<h\) put \(x=h-t,\ y=h+t\), both even. Their literal count is
\[
K_t(h)=h^2-t^2+F(h-t)+F(h+t)-1.
\]
For \(t=1\), \(h\ge5\), this is strictly above \(h^2\) by the lower estimate on \(F\). Also the neighboring evens \(h-1,h+1\) cannot both be powers of two: the only powers of two differing by 2 are 2 and 4. Their total \(F\) is therefore at most \(2h\), using oddness of each \(F\) to save at least two from \(2h+2\). Thus \(K_1(h)-h^2\le2h-2\). This proves strict square inclusion. For \(t=3\), \(h\ge9\), the excess is between \(h-8\) and \(2h-8\); for \(h=5,7\) the explicit children are 27 and 53, which also qualify.

**Distinctness except at 37.** Let \(D(n)=F(n+1)-F(n)\) and \(a=(h-1)/2\). Since the linear terms and parity corrections for four adjacent even arguments cancel,
\[
K_3(h)-K_1(h)=-8+D(a+1)-D(a-1).
\]
The binary recurrence gives
\[
D(2m)=-m-[m\text{ odd}],\qquad
D(2m+1)=D(m)+m+1+[m+1\text{ odd}].
\]
In particular \(D(n)\ge-n/2-1\) by induction. If \(a\) is odd, \(D(a+1)-D(a-1)\) is 0 or \(-2\). If \(a=2m\) with \(m=2r\), \(r\ge1\), the difference is \(-2r-2[r\text{ odd}]-D(r-1)+2<8\) by the same bound. If \(a=2m\) with \(m=2r+1\), the difference is \(D(r)+2r+2\). Equality to 8 forces \(D(r)+2r=6\), hence \(r\le4\); checking \(r=0,1,2,3,4\) gives only \(r=4\), or \(h=37\). There \(K_1=K_3=1425\); use \(t=5\) instead, giving \(K_5=1407\), strictly in the square shell by direct calculation. In particular \(5\) branches to \(27\) and \(33\); the same strict nesting argument as above gives uncountably many genuine **square** floor constants with strong tournaments at every step.

## 4. Separate absolute-size bound

Even without the square-shell symmetry, every integer \(N\ge9\) has a strong one-apex height in \((N,N+(3/\sqrt2)\sqrt N)\). Choose \(x\) from the two bounding powers of two around \(\sqrt N\) to minimize \(x+N/x\) (if \(\sqrt N\) is itself a power of two, choose it). Both choices are even and \(y=\lceil N/x\rceil-1\ge1\). The product-plus law reads \(H(uv)=(y+1)x+F(y)>N\). Its excess is at most \(x-1+y+1=x+y<x+N/x\). If \(s\in[1,2]\) is the ratio of \(\sqrt N\) to the lower power, the two choices give \(s+1/s\) and \(2/s+s/2\) in units of \(\sqrt N\). Their minimum has maximum \(3/\sqrt2\) at \(s=\sqrt2\). This is an absolute-size selector, not an exact integer selector.

## Replay

Run node replay.js. It checks the literal +2 construction through target 100,000, checks both strict square children for every odd parent below 100,000, and verifies four cubic steps with exact integers. It cross-checks the prefix count against the explicit inversion sum.
