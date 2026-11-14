# Free address windows and affine rigidity

## A physical four-count window with an unrestricted parameter

Write g(0)=1, g(2u)=g(u)+u+p(u), g(2u+1)=g(u), where p(u)=u mod 2.
For each positive integer u, there exist strong tournaments at each of the four consecutive odd counts
$$
908u+9g(u)+{752,754,756,758}.
$$
Unlike a window indexed by powers of two, the entire positive integer parameter is permitted here. The first window is 1669,1671,1673,1675. The construction is not an all-target theorem: for example 2027 is not in this bank, and its upper relative density in the odd numbers is at most 8/908, since all four counts at index u are greater than 908u.

For a literal marked ordinary core (rows, U), use the bilateral apex formula of the movable-triples note:
$$
K+fxy+Qx+Py+h(g(x)+g(y)-2)+Tp(x)+Sp(y).
$$
The four choices of rows, outneighbor mask U, multiplier M, appendage w, and right address b (use x=Mu+w,y=b) are:
$$
egin{array}{c|r|r|r|r|r}
\text{rows}& U&M&w&b&\delta\\
(2,8,3,21,7)&4&4&3&2&752\\
(2,8,3,21,7,31,63)&77&2&1&1&754\\
(2,8,3,21,7)&22&2&1&5&756\\
(2,8,3,21,7,31)&19&4&3&1&758 .
end{array}
$$
Use the true binary width in each outer wing. Both addresses are positive, so the apex beats the first ordinary vertex and is beaten by the last, and the output is strong.

If M=2^r and 0≤w<M, repeated use of the defining recurrence gives, with r-bit reversal,
$$
g(Mu+w)=g(u)+g(w)-1+u(M-1-\operatorname{rev}_r(w))
                    +p(u)(1-\operatorname{msb}_r(w)).
$$
Every witness has h=9, most significant appended bit 1, and slope M(fb+Q)+h(M-1-rev_r(w))=908. The displayed constants are obtained by substitution in the bilateral formula; linear_window_replay.py evaluates them from independent path DP, then recounts 40 whole tournaments directly. It also checks a complete inverse within this chart. The inverse uses dyadic residue branches and binary search on each once the linear rise dominates the maximal possible drop of g on that branch.

## There is no non-dyadic fixed affine reanchoring of free linear charts

**Theorem.** Let F_i(u)=D_i u+h_i g(u)+B_i(u) for real D_i, positive h_i and bounded B_i. For an integer a≥1 and b∈ℤ, the difference F_1(u)−F_2(a u+b) is bounded for all sufficiently large integer u if and only if
$$
a=M=2^r,quad 0≤b<M,quad
h_1=h_2,quad D_1=M D_2+h_2(M-1-\operatorname{rev}_r(b)).
$$
At r=0 use b=0 and an empty reversal.

**Proof.** Write φ(w) for the binary radical inverse of a nonnegative integer, obtained by reflecting the binary digits after the binary point, and ψ(w)=1−φ(w). For any fixed A≥1 and w≥0 the suffix formula gives
$$
\lim_{k\to\infty} \frac{g(A2^k+w)}{2^k}=A\psi(w).
$$
For example select a fixed bit width R covering w, apply the suffix identity with quotient A2^(k-R), and use g(A2^j)=g(A)+(2^j−1)A+p(A) when j≥1. Set u=2^k+w and compare the two putatively bounded-distance charts after division by 2^k. For every sufficiently large fixed w:
$$
D_1+h_1\psi(w)=aD_2+a h_2\psi(aw+b).
$$
Write a=Mq where M is dyadic and q is odd, and b=Md+c where 0≤c<M. Binary concatenation gives
$$
\psi(aw+b)=(M-1-\operatorname{rev}_r(c))/M+\psi(qw+d)/M.
$$
Both ψ(w) and ψ(qw+d), over all sufficiently large integer w, are dense in [0,1]: the odd multiplier q is invertible modulo every power of two, so every finite initial fractional digit pattern occurs. Equality of the two closed range intervals, using positivity of the h_i, forces h_1=q h_2 and D_1=aD_2+q h_2(M-1-rev_r(c)). The original relation then reduces to ψ(w)=ψ(qw+d). Radical inverse is injective on nonnegative integers, so q=1, d=0. This gives the necessity. The suffix recurrence gives sufficiency, including any bounded parity sidecar. ∎

This is a rigidity of uniform alignment. It does not forbid isolated exact intersections, or an address chosen non-affinely from the target. As a control, the apparent related banks at slopes 1816 and 1825 are exactly the old bank with u restricted to odd and even addresses respectively, because
$$
908(2u+1)+9g(2u+1)=1816u+908+9g(u),qquad
908(2u)+9g(2u)=1825u+9g(u)+9p(u).
$$
They do not supply new absolute targets. The discarded formal phase is analogous to the older MRB branch-jump example: asymptotic and local charts can create multiple presentations, but only data carrying the full observable certifies an independent exact shift. The silent square-birth/Pell note offers a complementary lesson: an exact cancellation relation can classify exceptional coincidences without being a universal existence mechanism.

## A finite menu of non-dyadic affine maps still cannot align uniformly

Extend ψ to the 2-adic integers by reflecting the infinite binary expansion after the point. This is continuous; fixing the first r digits changes the value by at most 2^(-r). Ordinary nonnegative integers remain a dense subset, on which the radical inverse is injective.

Fix one reference F and a finite collection F_j as in the theorem, and suppose every positive integer multiplier a_j has an odd factor greater than 1. For all sufficiently large u, one cannot keep |F(u)−F_j(a_j u+b_j)| bounded by adaptively choosing j. In fact a common arithmetic progression of dyadic scales u=2^k+w has discrepancies of order u for all j, for some fixed nonnegative w.

For the proof use the continuous limiting discrepancies
`L_j(w)=D+hψ(w)−a_j(D_j+h_jψ(a_jw+b_j))`. If one zero set contained a cylinder w=Lz+c with L dyadic, strip the low bits on both sides and compare the remaining reflected ranges and values as in the proof above. After stripping, the first address is z and the second is q_j z+d for q_j the odd part of a_j. Their reflected ranges are both dense in [0,1], and equality after normalization implies z=q_j z+d at every sufficiently large ordinary z, contradicting q_j>1. Thus all the zero sets are closed with empty interior, so there exists an open cylinder disjoint from their finite union. Choose one ordinary w in it and large enough for nonnegative arguments. Along u=2^k+w the discrepancies divided by 2^k tend to fixed nonzero L_j(w), as asserted.

This excludes a finite menu solely of non-dyadic uniform maps. It does not exclude a dyadic bank acting on an exact partition into residue cylinders, a true variable-core recurrence, or target-dependent affine parameters of unbounded complexity.

## Independent-wing transitive rigidity

For the old all-transitive spine, originally split in any place into two contact words, consider any dyadic suffix schedule x=2^r u+a,y=2^s v+b. The two short suffixes concatenate into one inner transitive old core of order n=r+s, with word the fixed-width binary a followed by the reverse complement of fixed-width b. This leaves the independently free positive outer addresses u and v. The coefficient of uv fixes f=2^n. After removing the universal g(u)+g(v), the independent linear coefficients determine P and Q, including if u,v each have restricted parity. In a transitive core P is exactly the binary value of the marked word: for each marked possible minimum, independently select an arbitrary subset of later vertices. Consequently P fixes the entire mark and therefore the constant K, S and T too. So no two such pure suffix schedules with common full two-variable response can have a nonzero constant shift; only actual nontransitive marked interfaces or coupling the two addresses can avoid the obstruction.
