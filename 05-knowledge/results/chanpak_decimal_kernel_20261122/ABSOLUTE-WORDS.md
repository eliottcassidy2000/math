# From a residue automaton to many distinct strong classes at one height

Fix one of the eight canonical Farey spines \(n=4,8,\ldots,32\), its four-vertex left flank with tags 12 and its physical middle flip. Let \(L(x,y)\) denote the lower count for binary addresses \(x,y\); the upper count is always \(L+2\). Let \(q>0\) be the coefficient of \(xy\). Given a fixed modulus \(m\) and lower residue \(a\), let \(\delta\) be the exact probability of \(L\equiv a\pmod m\) on two independent uniform register states, as in the [local inverse and Ogg-prime calculation](LOCAL-WORDS.md). All lower residues represented by odd integers have \(\delta>0\) in these charts.

**PROVED.** For every sufficiently large order \(N_v\) there exists a single lower height \(H\equiv a\pmod m\) with at least
\[
             \left(\frac{m\delta}{7q}-o(1)\right)N_v
\]
pairwise nonisomorphic **strong lower tournaments of that order** admitting a one-arc reversal to \(H+2\). Their upper tournaments are pairwise nonisomorphic as well. The height can vary with the order. In particular this holds for the exponent-61 residue at the product of the 15 Ogg primes, or at any fixed product of their powers. This is a multiplicity result for one prescribed congruence class. The exponent-61 absolute height still has to satisfy its own inverse.

### Recovering the two addresses from the ordinary class

Use the words \(11\alpha,11\beta\) as the binary representations of the two addresses. On the left, the first two physical vertices are beaten by the apex. On the right, reverse complementation makes the last two physical vertices beat the apex. The apex deletion is the ordinal sum of singletons, the three strong blocks \(A,D,V_n\), and the remaining singletons. Its largest strong component has size \(n+2\) (including the size-six tie when \(n=4\)). Both middle orientations are strong. The entire host is strong: the apex, first and last vertices form a directed triangle and include every intermediate vertex on a directed return.

Delete any vertex other than the apex. There is a surviving left vertex among the first two, and a surviving right vertex among the last two. Every vertex strictly between them is reached from the first and reaches the last. The strongly connected component containing the apex and these two vertices has order greater than \(n+2\), even after a deletion inside the largest block. Thus the unmarked graph itself identifies its apex as the unique vertex whose deletion has all components of size at most \(n+2\). The initial and terminal singleton strings of that deletion have unique orders, and their apex contacts recover both words. This argument works for each orientation.

### Counting inside a fixed order

Take \(|\alpha|=i\), \(|\beta|=j\), \(i+j=K\). The order of the host is \(N_v=K+n+17\). For one split there are \(2^K\) address pairs. The [finite register chain](LOCAL-WORDS.md) mixes to uniform measure exponentially in the number of free bits, from any initial state. Indeed it is doubly stochastic, and there is a common positive power of its transition matrix by the fixed-length control theorem. Restrict to \(i,j\ge\lceil\sqrt K\rceil\). Uniformly over these splits the proportion of pairs producing the prescribed phase is \(\delta+o(1)\). Different splits give different pairs. The combined count is \((\delta+o(1))K2^K\).

Now \(3\,2^i\le x<4\,2^i\) and \(3\,2^j\le y<4\,2^j\). Since \(0\le f(x)\le x\), the terms other than \(qxy\) are uniformly \(o(2^K)\) in this range. The heights lie between \(9q2^K\) and \((16q+o(1))2^K\). There are at most \((7q/m+o(1))2^K\) integers in the prescribed residue in that interval. Pigeonhole and address recovery prove the theorem. The sign, order, height and physical adjacency are retained. The restriction to the prefix 11 is the additional data needed to control unmarked isomorphisms.

### The height-weighted frequency of any local phase

**PROVED.** Count all nonnegative ordered address pairs with \(L(x,y)\le X\) and \(L(x,y)\equiv a\pmod m\), now without the leading 11 restriction. Their number is
\[
                      \frac{\delta}{q}X\log X+O_{m,n}(X).
\]
In particular the computed enrichment at 17 on spine 12 survives ordering the constructions by their actual Hamilton count. This formula counts presentations; the preceding fixed-order argument supplies a separate injective subfamily. [Earlier two-tape multiplicity work](../chanpak_rule30_campaign_evaluation_20261114/README.md)

For a single register state \(s\), its count on \([1,t]\) is \(d_s t+O_m(t^{1-\epsilon})\) for some \(\epsilon>0\), with \(d_s\) its uniform mass. To see this, partition a binary initial interval into at most one cylinder at each free suffix length. On a cylinder of length \(2^k\) the error is \(O_m(2^k\rho^k)\) for some \(\rho<1\) by mixing. Summing the geometric bounds gives the assertion after reducing \(\epsilon\) if necessary. A union of acceptable state pairs is finite. For a pair of states \(s,t\), summing over \(x y\le Z\) gives \(d_s d_t Z\log Z+O_m(Z)\) by applying the one-variable count at \(Z/x\) and partial summation. For the absolute height apply this with \(Z=X/q\). The chart lies between \(qxy\) and \(q(x+C)(y+C)\) for an integer \(C\) depending on the chart. For positive integers the two hyperbolic regions differ in \(O_C(Z)\) lattice points, by translation and the boundary rows and columns. The two axes add only \(O(X)\).

### Finite automata, tape computations and the absolute graph

For one address the update on an appended bit \(b\) is \(x'=2x+b\), \(f'=f+(1-b)(x+(x\bmod2))\). On two parallel inputs, \((xy)'=4xy+2b_y x+2b_x y+b_xb_y\). The integer registers \((1,x,y,xy,f(x),f(y),x\bmod2,y\bmod2)\) therefore have a fixed linear update for each of the four bit pairs. Reducing the registers modulo a fixed integer gives a finite machine. Exact variable-length registers give a deterministic tape verifier: process the bits in order, store the registers in binary and add shifted registers on each step. With \(k\) input columns this uses \(O(k)\) binary cells and \(O(k^2)\) steps on a fixed number of tapes, for a fixed chart. The same update admits a one-dimensional cellular implementation: encode the tracks and a head marker in each cell, perform each addition by a sweep carrying a finite carry state, and return to the input marker for the next column. A tape head sweep moves one cell per tick, so only neighboring cells exchange state. [Linear representations and regular sequences, Allouche and Shallit](https://www.cambridge.org/core/books/abs/automatic-sequences/kregular-sequences/474383BEAD7CF7CFFC5C406CDE05D1FF)

**PROVED.** No finite automaton recognizes the exact graph \((x,y,L(x,y))\) as three synchronous binary tracks read most significant bit first, padded on the left to the shortest common width. Restrict to \(x=y=2^k-1\). Then both digital reads vanish and
\[
 L=q(2^k-1)^2+(t+c)(2^k-1)+K+e+r .
\]
Its bit length is \(2k+d\) for all sufficiently large \(k\) and a fixed integer \(d\): the leading term is \(q4^k\), and the lower-order terms decide the side if \(q\) is a power of two. Intersect a supposed recognizer with the regular condition that the first two tracks agree and consist of zeros followed by ones. Project to the first track. Apart from finitely many small words, the resulting language is \(0^{k+d}1^k\), which fails the finite-state pumping lemma. This explains the distinction between successful local machines and an exact computation that needs growing memory.
