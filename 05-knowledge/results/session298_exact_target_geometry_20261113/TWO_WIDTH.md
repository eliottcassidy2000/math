# A two-width neutral repair chart

Let the cyclic quotient have parts C=TT_m, X={x}, B=TT_b and arrows C→x→B→C. Move a named internal root through C while preserving its total-order structure. If it has l predecessors and r successors, put U=2^b and V=3^b. For m,b≥1, the count and root-start count are
```
h=(2^m−1)(U−1)+2,  s=2^r.
```
For m≥2, deletion gives q=(2^(m−1)−1)(U−1)+2. The exact root cover coordinate is
```
F=2(2^r·3^l−2^l)(V−U)+3·2^l U−2·2^l−[r=0](U−1).
```

Split a Hamilton path before and after x. When both pieces are nonempty, the first is B*C+ and the second B+C*; independently choose a nonempty set of C vertices for the first and a nonempty set of B vertices for the second. Each empty-side possibility is unique. If the path starts at the marked root, B is all in the second piece and the first resource run selects an arbitrary subset of the root's successors.

For F use ordered pieces (P,Q) with P ending at root. If x is absent from P, it is an optional B-prefix and a C-run. Choose i predecessors and j B vertices for P. Its complement count is
```
(2^(l+r−i)−1)(2^(b−j)−1)+2−[l+r−i=0 or b−j=0].
```
The correction handles the transitive case when a part is empty. Summing binomially gives (2^r 3^l−2^l)(V−U)+2·2^l U−2^l−[r=0](U−1). If x lies in P, the post-x B run is nonempty and the terminal C run contains the root. For an empty pre-x part there are 2^l(U−1) choices. Otherwise each predecessor has three assignments (two runs or Q), each successor two, and the pre-x C run must be nonempty; this gives 2^r3^l−2^l assignments. The B allocation similarly gives V−U. Q is transitive in this case. Adding yields the claimed formula.

## A nondyadic exact source: H47

At m=4,b=2 the graph has seven vertices and H=47. At l=0 the state is (47,8,23,80), giving δ=18 and t=57·2^k+5, never divisible by 3. Moving the same root two positions gives (47,2,23,180), δ=112 and d=157. Its returned pair is

    t=157·2^k−89,  h'=3t+224.

For every k≡1 (mod 6), t is odd, divisible by 3 and nonzero modulo 7. Since 224 has only primes 2 and 7, the pair is primitive; h' is a unit modulo 3. Therefore this exact neutral preparation supplies infinitely many actual entries into the balanced controller at an initial target not of the form 2^m+1.

The chart still contains only O(log² X) distinct initial targets at most X: each factor in (2^m−1)(2^b−1) is at least one, so each parameter is O(log X). This is an exact source with renewable outgoing branch, not density or arbitrary-target supply.
