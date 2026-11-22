# Selecting both endpoint counts from an absolute target

**PROVED.** Fix the tagged six-vertex middle and its designated edit from the [exponent-15 certificate](E15-FAREY.md). Given a positive odd lower target \(N\), the following finite inverse returns a pair of tournaments with counts \(N,N+2\) differing in that edge if one occurs in the class of arbitrary nonempty left and right tagged tournaments and two arbitrary exterior words. The denominators are bounded by the absolute height, and they determine both endpoint numerators before either flank is built. The theorem gives a decision procedure for this grammar. Its search size at exponent 61 is large.

Let \(a=H(A), b=H(V)\), let \(u\) count the full paths of \(A\) ending at a nontag, and let \(r\) count the full paths of \(V\) starting at a nontag. In the indicated direction the physical half-change is \(ar-ub\). A unit edit requires
\[
       29ab\le N,\qquad \gcd(a,b)=1,\qquad ar-ub=\epsilon\in\{-1,1\}.
\]
The first inequality holds because deleting the apex and exterior strings leaves an induced ordinal sum of height \(29ab\). Every Hamilton path of an induced subtournament extends to the full tournament by inserting the missing vertices one at a time. At each step insert at an end when possible, otherwise at a consecutive transition from a vertex beating the new vertex to one beaten by it. Deletion recovers the original path, giving an injection.

For each pair of odd denominators satisfying the bound, and each sign, if \(a,b>1\) set
\[
       u=(-\epsilon b^{-1})\bmod a,\qquad r=(ub+\epsilon)/a .
\]
Here the first quantity is the least nonnegative residue. This is the unique possibility in \([0,a]\times[0,b]\), since extremal numerators would force the other denominator to be one. If either denominator is one, take all solutions in the closed rectangle directly. Every accepted numerator is coprime to its own denominator.

### Why two strong components suffice for the boundary test

The strongly connected components of \(A\) form an ordered sum. If \(P\) is the product of the heights of all components before its last, then \(P\mid a\) and \(P\mid u\): all full paths terminate in the last component and each path there has exactly \(P\) preceding choices. Consequently \(P=1\). A nonsingleton strong component has at least three paths, so all preceding components are singletons. Dually all components of \(V\) following its first are singletons since their height product divides both \(b\) and \(r\). Move these chains into the respective exterior words, keeping their apex contacts in order. If the entire side is transitive keep its boundary singleton. Exterior leading zero padding is harmless for the height and can be dropped in this existence test.

For \(a>1\), enumerate strong cores of order at most \(1+3\log_5 a\) with height \(a\), using Busch's bound. Select nontag subsets of their terminal-path vectors summing to \(u\). Perform the dual construction with starting-path vectors and \(r\) for \(b\). A tag-preserving automorphism group acts freely on full paths because fixing an ordered spanning path fixes every vertex. It also acts freely on the subset counted by \(u\) or \(r\). Its order therefore divides two coprime integers and is one on each side. Ordinary automorphisms may still identify *different* taggings. [Busch, The number of Hamiltonian paths in strong tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)

For each pair of retained tagged cores compute the joined lower chart \((h,q,c,t,e,r',k)\) in the orientation selected by the sign. All coefficients are nonnegative. The middle has positive two-slot readings \(c_D,t_D\) in both orientations, so the join formula and positive heights therefore give \(c,t>0\). For \(k\le N\), enumerate the finite rectangle
\[
 0\le x\le (N-k)/t,\qquad 0\le y\le (N-k)/c.
\]
Accept at \(k+t x+c y+qxy+r'[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y))=N\). Necessity follows from component reduction and the exact chart, and sufficiency from the tagged construction with this very middle arc. This generalizes the fixed \(3/5\) Berggren branch by selecting *both* numerators and both strong cores. The residue controller acts later, on the exterior words. [Fixed-branch selector](BERGGREN-RECEIVERS.md) · [Unrestricted inverse on canonical tower levels](FAREY-HYPERBOLA.md) · [Controlling the exterior words](LOCAL-WORDS.md)
