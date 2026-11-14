# Exact overlap on the seven rail, and the state that survives it

**Proved:** disjoint multiplicative coordinates at every scale, positive same-count overlap among strong tournaments, and an exact marked endpoint decoder for these overlaps. The result does not make the whole odd spectrum an interval.

## 1. A reproducible lateral source pass

Starting from the lexicographically listed top-level results Markdown of the inherited local \`origin/main\`, Python \`random.Random(325).sample(paths,18)\` sampled, among others:

- [Modular Fourier and Zeckendorf routes](../modular_fourier_zeckendorf_routes_20261003.md), where a common spectral signal cannot replace the second carry register and actual route pointer.
- [Mixed Smith uniform residues](../third_20260906_smith_density.md), where exact collision depth and the number of colliding sources differ, even after the main valuation is known.
- [Mazur positive-density audit](../mazur_positive_density_20260928_audit.md), whose exact fixed-word affine identity keeps both the 2-adic source cylinder and 3-adic inverse congruence. Neither resonance nor a mean is an integer inverse.
- [The complete planar Jacobian boundary audit](../planar_jc48_sep08_four_three_one_audit.md), where a local degree can exceed the entire global pole budget; existence of each local trace does not supply a lawful global primitive.

The native transfer is deliberately limited. First apply two simultaneous valuation gates, as in the Collatz audit, to rule out accidental scalar overlap between different chart exponents. Then build and count actual overlapping sources, retaining the marked coordinate needed for continuation; this is the Zeckendorf/Smith distinction in a tournament one can inspect. The route itself comes from [the seven rail](../session324_busch_inverse_and_gap_structure_20261113/FERMAT_SEVEN_RAIL.md), not from a Collatz or Jacobian correspondence.

## 2. The absolute exponent lattice has no collisions

For \(a\ge0,\ s\ge2\), all numbers
\[
7\cdot 3^a(2^s+1)
\]
are distinct as pairs \((a,s)\) vary. In particular the seven rail \(7(2^s+1)\) enters a different multiplication-by-three chain for each \(s\); it is not a disguised single ladder.

Suppose an equality uses different exponents. By swapping the two pairs if necessary it gives \(2^s+1=3^d(2^t+1)\) with \(s>t\ge2\) and \(d>0\). If \(t=2\), the right side is 5 or 7 modulo 8, while the left side is 1. For \(t\ge3\), reducing modulo 8 forces \(d\) even. Subtracting yields
\[
2^t(2^{s-t}-3^d)=3^d-1,
\]
whose first factor after \(2^t\) is odd. Therefore \(t=v_2(3^d-1)=2+v_2(d)\). The valuation identity follows by starting at \(3^{2u}-1=(3^u-1)(3^u+1)\) for odd \(u\), with valuations 1 and 2, and adding one at each subsequent doubling.

The original equality forces \(s\) odd. The identity \(v_3(2^s+1)=1+v_3(s)\) follows by reduction modulo 9 when \(3\nmid s\), then factorization \(x^3+1=(x+1)(x^2-x+1)\) when lifting an exponent by three; the second factor has exact valuation one for \(x\equiv-1\bmod3\). Hence \(s\ge3^{d-1}\). But \(3^d(2^t+1)<2^{t+2d+1}\) gives \(s\le t+2d\), and \(t=2+v_2(d)\le2+d/2\) for even \(d\ge2\). For \(d\ge4\), \(3^{d-1}>2+5d/2\) (already at 4 and increasing faster), a contradiction. If \(d=2\), then \(t=3\) and \(2^s+1=81\), again impossible. Equal \(s\) of course forces equal \(a\).

Thus residue closure alone gives no overlap here. In particular the independent strong ray \(63\cdot3^j\) intersects the unmultiplied seven rail only at \(s=3,j=0\): it is the \(s=3\) chain since \(2^3+1=3^2\).

## 3. Strong overlaps within every suitable exponent pair

The literal star of the S324 rail has blocks \(I\to U\to W\) and an apex \(v\). On \(I\) its outgoing rows are \((4,1,10,67,15,31,55)\), with apex out-cut mask 81. On \(W\) its outgoing rows are \((8,1,3,6)\), with apex out-cut mask 13. All of \(U\) beats \(v\). The exact arbitrary-module law there is
\[
H(T)=7\bigl(512 Q_U(2)+H(U)\bigr),
\]
where \(Q(2)\) counts ordered two-path covers with empty slots allowed. The quotient is strong for every nonempty \(U\).

Take \(U\) to be any ordered join of \(a\) directed triangles and \(b\) singleton components, \(a,b\ge0\), \(a+b>0\). A triangle has \((H,Q(2))=(3,12)\); a singleton has \((1,2)\); both coordinates multiply on ordinal joins. Thus every such order gives the literal **strong** count
\[
H(T)=7\cdot3^a(2^{\,9+2a+b}+1).
\]
This proves a strong subcone of the \((a,s)\) lattice. In particular it is not justified to assert this star realizes every 3-multiple of a fixed \(s\): it specifically requires \(s=9+2a+b\).

At fixed \(a,b\) there are \(\binom{a+b}{a}\) choices of the SCC word of \(U\). If \(3a+b>12\), these yield that many nonisomorphic strong tournaments at the same count and order. Here is the precise recognizer. Let \(x\) denote the root in the 13-vertex quotient formed by making \(U\) a singleton. On this quotient no proper nonsingleton module contains \(x\); the supplied finite verifier checks the full set of candidates. In any inflation, a module containing \(U\) and other vertices therefore equals the whole graph. If a module \(M\) intersects \(U\), meets its complement, and omits some vertex of \(U\), then \(M\cup U\) is a module: any vertex outside the union has uniform orientation to each of the two intersecting modules and hence the same orientation to both. It projects to a module containing \(x\) and another vertex, so \(M\cup U\) must be the whole graph. But a vertex in \(U\setminus M\) loses to every vertex of \(I\) and beats every vertex of \(W\); both are subsets of \(M\), a contradiction. Thus every proper module either lies inside \(U\) or is disjoint from it, and the latter has at most 12 vertices. This makes \(U\) the unique largest proper module. Its own SCC word then distinguishes the displayed isomorphism classes.

## 4. The exact lost continuation coordinate

Put \(S=N^+(v)\cap I\), and let \(\alpha_0\) be the number of spanning paths of \(I\) starting in \(S\). From the displayed mask \(\alpha_0=12\). Let \(\beta=14\) be the number of ordered two-piece covers of \(W\) whose first piece is empty or ends at a vertex beating \(v\), as in the inherited star law.

For a marked singleton SCC \(r\) of \(U=D\to r\to E\), let \(S_T(r)\) count spanning paths of the full inflated tournament starting at \(r\). Then
\[
S_T(r)=\alpha_0\beta H(D)Q_E(2)
      =168\cdot3^a\cdot2^{2a_R+b_R},
\]
where \(a_R,b_R\) are the numbers of triangle and singleton components of \(U\) strictly after \(r\). To prove it, split the full path at \(v\). A path starting at \(r\) has no \(I\)-prefix, since \(I\to U\to W\) and the only return is through \(v\). Its entire \(I\)-part is therefore the post-apex first block and has \(\alpha_0\) choices. The \(W\) split contributes \(\beta\). All of \(D\) must lie in the second \(U\) slot (otherwise it precedes \(r\)); it contributes \(H(D)\). The singleton \(r\) itself starts the first slot, and either slot may traverse any allowed part of \(E\), giving \(Q_E(2)\). Conversely these choices construct the whole path and are read off uniquely.

Moving a marked singleton from just before to just after a triangle divides its actual start count by four, while leaving total \(H\) exactly unchanged. Moving it past a singleton divides it by two. Moreover the starts at the singleton SCCs reconstruct their entire placement: list them right-to-left, divide by \(168\cdot3^a\), and read the binary exponents \(e_j=2a_R+j-1\). This recovers \(a_R\) at each singleton, hence the gap allocation. When \(U\) is intrinsic, the unordered list of these values has the same recovery power because the values strictly increase from right to left.

**A complete marked start corridor at fixed count.** If \(b\ge1\), one distinguished singleton may have arbitrary \(a_R=i\in[0,a]\) and \(b_R=j\in[0,b-1]\), realized by placing the chosen numbers of components on each side. The possible start exponents are therefore exactly
\[
\{2i+j:0\le i\le a,\ 0\le j<b\}.
\]
For \(b\ge2\) the intervals \([2i,2i+b-1]\) touch or overlap, giving every integer exponent \(0,\ldots,2a+b-1\) without changing the absolute \(H\). For \(b=1\) the set is only the even progression \(0,2,\ldots,2a\). Multiple marked starts are not independent: a proposed list \(e_1,\ldots,e_b\), in right-to-left singleton order, is realized iff \((e_j-j+1)/2\) are integers, nondecreasing in \(j\), and all in \([0,a]\). Necessity follows from the suffix counts; sufficiency constructs the intervening triangle runs from the successive differences, with the two end runs determined by the extremes.

Thus this chart has two different kinds of overlap: \((a,s)\) has no arithmetic overlap at all, whereas the actual strong inverse fiber at fixed \((a,s)\) has many physically distinct members. The extra state does not buy a new absolute target by itself, but it provides a tunable **marked** start exponent at an exact count. A proof seeking an endpoint-ready continuation can vary these coordinates without paying an \(H\)-drift. This does not assert a readiness profile beyond the proved start count.

Replay: python3 replay_rail_overlap.py.
