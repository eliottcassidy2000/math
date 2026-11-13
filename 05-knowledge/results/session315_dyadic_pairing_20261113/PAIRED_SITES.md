# The paired-site model for self-converse tournaments and blue lines

A self-converse tournament on $2m+\epsilon$ vertices, with $\epsilon\in\{0,1\}$, can be encoded as **an ordered pair consisting of a tournament and an undirected graph on $m$ sites**, with one additional binary site color when $\epsilon=1$. This carries actual blue path-reversal lines, as well as a scalar bivariate class polynomial.

## Canonical directed toothpicks

Choose an involutory anti-automorphism $\rho$. Existence is inherited from [THM-024](../../../01-canon/theorems/THM-024-involution-anti-aut.md) or [THM-647](../../../01-canon/theorems/THM-647-anti-redei-theorem.md). Every automorphism has odd order: an involution would swap and reverse some arc. The extended group of automorphisms and anti-automorphisms has twice this odd order, so all anti-involutions are Sylow-2 generators and are conjugate by automorphisms (an odd-coset conjugator can be multiplied by the involution it is conjugating). An anti-involution fixes at most one vertex, so it has $m$ transposed pairs and a fixed center $c$ iff $\epsilon=1$.

In each pair name the vertices $a_i,b_i$ intrinsically by $a_i\to b_i$. Encode a tournament $A$ by $A_{ij}=[a_i\to a_j]$. Encode an undirected graph $G$ by $G_{ij}=[a_i\to b_j]=[a_j\to b_i]$ for $i\ne j$. If a center exists let $u_i=[a_i\to c]=[c\to b_i]$. The equalities are forced by $\rho$; the orientations on the $b_i$ form the converse of $A$. Conversely all of these bits can be chosen independently and specify the tournament, with the intended anti-involution. Relabeling preserving $\rho$ must send each $a_i$ to an $a_j$ because of the directed twin arc. Thus marked objects modulo isomorphism are exactly the simultaneous $S_m$-orbits of $(A,G,u)$, and the forgetful map from these orbits to ordinary self-converse classes is bijective by conjugacy of anti-involutions. The stabilizer of $(A,G,u)$ is exactly $C_{\operatorname{Aut}(T)}(\rho)$, not in general the entire automorphism group.

This is a genuine two-vertex growth recursion. Add a new site, choose its two bits per old site (tournament direction and graph edge), and its center color if applicable. Canonical augmentation can retain the old-site deletion witness, then take simultaneous site-permutation classes. It does not require trying arbitrary new old-new arcs and testing for an anti-involution afterward. In labelled coordinates the inventory multiplier is $2^{2m+\epsilon}$.

## Decorated cycle-index generating function

Let $z_\lambda=\prod_d d^{m_d}m_d!$ for a partition $\lambda$ of $m$; restrict to odd parts since a fixed tournament $A$ has no even-order automorphism. Let $\mathcal O(\lambda)$ be the multiset of lengths of orbits on unordered site pairs. For a cycle of length $d$ there are $(d-1)/2$ internal orbits of length $d$. Between two cycles of lengths $d,e$ there are $\gcd(d,e)$ orbits of length $\operatorname{lcm}(d,e)$. Denote the total number of orbits by $o(\lambda)$. The exact joint class polynomial in undirected edge count and center weight is
$$
 F_{m,\epsilon}(x,v)=
 \sum_{\lambda\vdash m,\ \lambda\ {\rm odd}}
 \frac{2^{o(\lambda)}}{z_\lambda}
 \prod_{\ell\in\mathcal O(\lambda)}(1+x^\ell)
 \prod_{d\in\lambda}(1+v^d)^\epsilon .
$$
Indeed on each unordered-pair orbit one freely chooses the tournament orientation (a return cannot reverse a pair since the action has odd order), and independently chooses all or none of its undirected graph edges. Center colors are independently constant on each site cycle. Burnside with invariant weights proves the polynomial. Replacing $G$ by its complement leaves this formula unchanged. At $x=v=1$ the result is $\sum_\lambda 2^{2o(\lambda)+\epsilon\ell(\lambda)}/z_\lambda$, giving the known self-converse counts via a new decorated-site derivation; no literature-priority claim.

## Exact folded paths and their blue partners

An anti-palindromic Hamilton path with respect to $\rho$ is specified by a directed word $W=(w_1,\ldots,w_m)$ picking exactly one vertex $w_i$ from each pair (in arbitrary order), followed by the center if present, then $\rho(w_m),\ldots,\rho(w_1)$. The middle junction is allowed iff the last vertex is a canonical tail $a_i$ at even order, or $w_m\to c$ at odd order. The full path is closed iff the first vertex is a canonical head $b_i$, because its terminal-to-initial arc is $\rho(w_1)\to w_1$.

An exact subset transfer $D(S,w)$ on sets of selected sites starts at $D(\{i\},a_i)=D(\{i\},b_i)=1$, and sums $D(S-\{i\},v)$ over all signed vertices $v$ in earlier sites with $v\to w$, where $w$ lies at site $i$. Retain first sign to count closed ports, or the actual first vertex for marked start data; then filter the final signed vertex by the middle junction. Thus the exponential dimension for blue path counts is in $m$, not $2m+\epsilon$. Quotient these $\rho$-fixed paths by $C_{\operatorname{Aut}(T)}(\rho)$ to obtain exactly the blue tiling darts: every blue path has a unique positional reversing involution, and two $\rho$-fixed paths isomorphic in $T$ can only be related by an automorphism commuting with that involution. The action is free on directed paths. In particular divide the fixed-$\rho$ count by the **centralizer** order, not by the full automorphism order. Counting paths fixed by all choices of involution, with each path assigned to its unique involution, gives the alternate full-automorphism normalization.

Use the path-reversal line convention of [THM-4729](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md). For consecutive vertices in the half-word, toggle the $A$-edge of their sites if the signs agree (both tails or both heads) and the $G$-edge if signs differ. The second half of the path reverses precisely the mirror arcs of those toggles. At odd order flip the center color of the last site. At even order the middle twin arc reverses, so switch tail and head labels at that last site to restore the canonical gauge. For each other site $j$, with $i$ the switched site and $x=A_{ij}$, $y=G_{ij}$ measured *after* the half-word toggles, the regauging is
$$
 (A_{ij},G_{ij})\longmapsto(1-y,1-x).
$$
Indeed the new $a_i\to a_j$ was $b_i\to a_j$, and the new $a_i\to b_j$ was $b_i\to b_j$. Thus the full physical blue partner, not just its count, is computed inside the paired-site space. Recover all blue edges with parallel multiplicities by pairing the resulting darts as usual. The fold does not handle black edges, which need the general disagreement recurrence or full tiling data.

## A signed recurrence that isolates strong nodes

Write any self-converse tournament uniquely as an ordinal sum $L\to C\to L^{\rm op}$, with $C$ the middle strong SCC or empty and $L$ an arbitrary tournament (possibly empty). The anti-involution reverses the SCC sequence and pairs left with right. In the paired-site model, all $|L|$ corresponding tails lie in $L$, so every $G$-edge touching an exterior site is present, and the remaining graph is exactly that of $C$. Thus the number of missing $G$-edges is unchanged. At odd order the exterior sites each have center color 1, so center weight is increased by $|L|$.

The Hamiltonian paths in an ordinal sum concatenate, so $H(T)=H(L)^2H(C)$. A fixed-$\rho$ path has arbitrary left path, the anti-palindromic central path, then the forced paired right path. The centralizer likewise has an arbitrary automorphism on $L$, its forced mirror, and the centralizer in $C$. Hence blue multiplicities satisfy $b(T)=t(L)b(C)$, where $t(L)=H(L)/|\operatorname{Aut}L|$ and the empty factors are 1. This gives an exact joint factorization: in the power series by pair-sites $z$, missing graph edges $x$, center weight $v$, with basis $[h,b]$ multiplying by componentwise multiplication, the all-self-converse series is the product of the central strong series and
$$
 \sum_{[L]}z^{|L|}v^{\epsilon|L|}[H(L)^2,t(L)].
$$
The central series includes an empty term only at even parity (and the singleton at odd parity). This is a triangular deconvolution for the strong sector. It does not imply that the full unlabelled graph is determined by the resulting coefficients.

At odd order every blue edge changes center weight by exactly one. It follows that the odd blue graph is graded by this intrinsic weight, with no loops. For $m>0$, setting $v=-1$ in Burnside gives zero, because every site cycle in an odd partition has odd length. Thus the two parity shores have equal size, whether or not the graph is connected. More sharply, forget $h,b,x$ in the factorization and keep the sign $(-1)^{\sum_i u_i}$. The signed odd strong self-converse inventory is exactly
$$
 \sum_{m\ge0}z^m\!\!\sum_{\substack{[T]\ {\rm strong,SC}\\ |T|=2m+1}}
     (-1)^{\sum_i u_i}
 =\frac1{\mathcal T(-z)},
 \qquad
 \mathcal T(z)=\sum_{j\ge0} A000568(j)z^j
$$
with $A000568(0)=1$ and singleton considered strong. The numerator is 1 because every positive pair-site coefficient in the total signed series vanishes. This scalar series is for one specific signed observable, not the full graph.

## Reproduction and scope

The companion script independently constructs all paired-site states through order 7, compares every fixed-twist folded partner against literal reversal and regauging, checks the folded counts and closed subcounts, and checks each weighted orbit coefficient against a separate Burnside sum. Orders 2 through 7 give class counts 1, 2, 2, 8, 12, 88. The proofs do not depend on these finite checks. The model gives exact self-converse growth and all blue edges, but not a small scalar model for black wiring or an eightfold period.

## A two-state joint edge generating function

There is an exact coarse edge polynomial that does not need class canonicalization. Mark the old and new undirected-edge counts by $x,y$, the old and new center weights by $v,w$ at odd order, and closure of the old dart by $t$. Fix the half-word site order to $1,\ldots,m$: every actual blue dart has a unique such coordinate presentation because ordering the sites kills the remaining site relabeling. Choose only the sign word. Its first sign contributes the row $(1,t)$ (tail is open, head is closed).

For consecutive sites not at a regauged boundary, the sign-pair matrix is
$$
M(x,y)=\begin{pmatrix}1+xy&2x\\2y&1+xy\end{pmatrix}.
$$
Equal signs force the $A$-direction and leave the $G$-bit free and unchanged. Tail-to-head forces old $G=1$ and new $G=0$ (factor $2x$, with two choices of $A$). Head-to-tail forces the reverse change (factor $2y$). A nonadjacent pair neither toggled nor regauged contributes $2(1+xy)$.

At odd order $n=2m+1$ with $m\ge1$, all center bits except the last are free and unchanged. A final tail forces the old center bit 1 and flips it to 0; a final head forces 0 and flips it to 1. Thus the exact blue dart polynomial is
$$
[2(1+xy)]^{\binom m2-(m-1)}(1+vw)^{m-1}
(1,t)M(x,y)^{m-1}\binom{v}{w}.
$$
In particular from a closed endpoint, an increase in center weight forces an even change in the undirected edge count, and a decrease forces an odd change. Indeed the number of $G$-toggles modulo 2 is the XOR of first and last signs. The first is a head iff closed, and the last is a head iff the center weight increases. Therefore the two node grades determine the orientation of every odd blue line joining them: the lower-center-weight endpoint is closed iff the two undirected edge counts have the same parity. There cannot be opposite closure orientations among parallel odd blue lines with the same endpoint classes.

At even order $n=2m$ with $m\ge2$, every pair meeting the final site is regauged, and the final sign must be tail. For a nonadjacent pair meeting that site, the old $G$-bit and the new $G$-bit (the post-toggle $A$-direction toward the final site) are independent, giving $(1+x)(1+y)$. For the final adjacent pair, an entering tail forces old $A=1$ toward the final site, toggled to zero, giving $1+x$. An entering head forces old $G=0$ and leaves $A$ free, giving $1+y$. The exact polynomial is consequently
$$
[2(1+xy)]^{\binom{m-1}2-(m-2)}[(1+x)(1+y)]^{m-2}
(1,t)M(x,y)^{m-2}\binom{1+x}{1+y}.
$$
In both polynomials the coefficient of $t$ counts each actual blue line once, since each line has exactly one closed dart. At unit weights, the total dart counts recover the ordinary blue tiling census. The replay checks every coefficient by literal signed-word enumeration through order 8. The polynomials retain real adjacency between grades; they do not identify individual class pairs within a grade.

The transfer is genuinely order two: $M=(1+xy)I+N$ with $N^2=4xy I$, so after removing the prefactors its scalar matrix contractions satisfy the recurrence with coefficients $2(1+xy)$ and $-(1-xy)^2$. The two formal eigenvalues are $(1+\sqrt{xy})^2$ and $(1-\sqrt{xy})^2$. In particular at odd order with $m\ge2$, forgetting the graph weights gives exactly
$$
2^{m^2-m-1}(1+t)(1+vw)^{m-1}(v+w).
$$
Thus the projected blue lines form a binomial ladder between center grades $k,k+1$: there are $2^{m^2-m-1}\binom{m-1}{k}$ closed-source lines in each direction, once graph-edge parity is forgotten. This is a count of actual line instances, not merely of possible transitions. The direct replay also checks the first four signed strong coefficients $1,1,0,1$ against the reciprocal law.
