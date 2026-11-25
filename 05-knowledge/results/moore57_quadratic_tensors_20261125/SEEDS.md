# A twenty-branch partial construction and a test for unrestricted extensions

**FINITE-EXACT.** The [twenty-branch array](seed20.txt) constructs a connected graph on 1,120 sites, with a complete matching between every two of its twenty branches, and no triangles or quadrangles. On adjoining the root and these twenty branch heads it gives a 1,141-vertex rooted partial Moore chart. Its nonconstant fiber spectrum lies strictly between \(-8\) and \(7\), as required by interlacing. Yet **even one site of a twenty-first branch cannot be attached with one neighbor in each of the original twenty branches** without introducing a square or triangle. The projection onto the [first nineteen branches](seed19.txt) admits at most one further full branch in any girth-five extension, even with arbitrary matchings on every missing branch. This is a certificate of a search boundary, not a full Moore graph.

## The extension test

In a cyclic partial chart on \(k\) branches let the transport \(i\to j\) add \(a_{ij}=-a_{ji}\) in \(\mathbb Z/56\). The supplied text files give successive rows of the strict upper triangle. No short cycle occurs exactly when the \(k-1\) offsets
\[
                  a_{ij},\qquad a_{i\ell}+a_{\ell j}\quad(\ell\ne i,j)
\]
are distinct at each pair. These offsets comprise the set \(F_{ij}\). A prospective exterior site chooses **one physical neighbor** \(x_i\) in each branch. It can be attached precisely when \(x_j-x_i\notin F_{ij}\) for every pair: a direct old edge would make a triangle, and an old two-step route a square. Subtracting \(x_0\) normalizes the possible words to \(w_0=0\). One normalized legal word gives an entire additional branch by taking all 56 simultaneous translates. Distinct translates have disjoint neighborhoods, so there is no new short cycle involving two new sites. Conversely, even an arbitrary extra branch must contain a legal word. The step from one site to a whole branch is special to the symmetry of the *given seed*. No symmetry of a putative final graph is assumed.

For more than one site, let \(m_w\) be the number with normalized word \(w\). At an unused cell \((i,j,d)\), at most 56 sites can have \(w_j-w_i=d\), since there are only 56 physical pairs of those slots. In a full Moore completion each physical pair requires its own midpoint, so every one of these capacities is attained. Thus simultaneous packing, or real-valued nonnegative weights relaxing the \(m_w\), gives stronger tests than checking each missing offset independently. This is the repeated-label distinction in the attached paper's list-distinctness argument in a finite form.

## Exact outcomes

On the 19-branch seed there are only six normalized words; they miss 5,964 required difference cells. Each word hits at least two of \((0,8,12),(2,5,11),(3,12,3)\), with zero-based labels. The three capacities therefore give at most \(3\cdot56/2=84\) exterior sites, fewer than two branches. Every pair of the six words has a repeated coordinate difference \(w_i-v_i=w_j-v_j\), \(i\ne j\). Consequently once the full orbit of one word has been inserted, a new site of any normalized word would agree in two positions with a translate already used. This is a forbidden square. The twentieth branch in the displayed seed is the orbit of the last word in the replay order. Both the nineteen- and twenty-branch graphs are connected: their triangular holonomies generate the full cyclic group (the anchor row is zero and the remaining increments generate \(\mathbb Z/56\)).

On the [separate 18-branch seed](seed18.txt) all remaining 39 offsets at all 153 pairs occur among its 2,479 legal words. Even these local choices cannot be allocated simultaneously. [Exact capacity certificates](SEED-BALANCE.md)

To certify the spectral test without approximate eigenvalues, let \(B\) be the partial matching adjacency and restrict to the subspace \(W\) of vectors with mean zero in each branch. In a full degree-57 chart the compression on those complete fibers has eigenvalues \(7,-8\); the partial chart is a principal compression on the direct sum of the selected fiber subspaces. Thus a necessary condition is that \(C=2B+I\) on \(W\) have all eigenvalues in \([-15,15]\). With \(u_i=56e_{i,0}-\mathbf{1}_i\), translation invariance gives
\[
                 {\rm tr}_{W}(C^{2m})=\frac1{56}\sum_i\|C^m u_i\|^2 .
\]
The positive even trace is strictly smaller than \(225^m\) at an explicitly found \(m\) for each seed, so every eigenvalue satisfies the strict bounds. [Exact integer trace replay](seed_spectrum.py) · [Full fiber spectrum](../moore57_switching_flags_20261124/LOCAL-SMITH.md)

From this directory run

~~~sh
c++ -O3 -std=c++17 seed_transversals.cpp -o /tmp/moore57-transversals
/tmp/moore57-transversals --words < seed19.txt
/tmp/moore57-transversals < seed20.txt
/tmp/moore57-transversals --words < seed18.txt | python3 -B seed_balance.py
python3 -B seed_spectrum.py
~~~

The [enumerator](seed_transversals.cpp) starts the first coordinate at zero, takes the next coordinate with the least available values, and propagates the forbidden differences after each assignment. Each normalized solution has a unique path in this exhaustion. The explicit arrays are positive certificates irrespective of the heuristic used to choose them. They were obtained by repeatedly adding a cyclic branch with randomized backtracking and selecting candidates by the spectrum of their Fourier matrices.
