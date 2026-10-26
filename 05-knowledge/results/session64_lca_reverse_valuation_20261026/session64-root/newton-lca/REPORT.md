# LCA depths give the Newton precision losses

Session64. **PROVED:** the following elementary bridge for the specific
evaluation nodes p,p^2,... at a prime p. **FINITE-EXACT:** an independent
standard-library checker, importing no production counting engine.

Session56's least-significant-digit trie gives

    depth(LCA(x,y)) = v_p(x-y).

For distinct positive indices r,i,

    v_p(p^r-p^i)=min(r,i).

Consequently, for q_j(X)=product_(i=1)^j(X-p^i) and r>j,

    v_p(q_j(p^r))
      = sum_(i=1)^j depth(LCA(p^r,p^i))
      = 1+2+...+j = j(j+1)/2.

For r<=j the product vanishes exactly. This is a concrete transport of the
prior LCA observer into the new Newton bank: the precision lost by Newton
coordinate j is the sum of its evaluation node's pairwise ancestor depths
with preceding nodes. It is not an analogy based on the word "valuation".

At absolute precision p^K on nodes p,...,p^S, the monic integer Newton basis
therefore gives the observation image

    direct_sum_(0<=j<S, j(j+1)/2<K) Z/(p^(K-j(j+1)/2)).

Proof of exactness: at row1 only coefficient c_0 contributes. After its
required divisibility is imposed, it vanishes at every row. Inductively row
j+1 forces c_j divisible by p^(K-j(j+1)/2) when that exponent is positive;
the corresponding diagonal after removing its p-power is a unit. Conversely
those divisibilities annihilate every evaluation. For degree<h, stop at
j<min(h,S). This proves both directions and the image cardinality, not just
a sufficient collection of zero vectors. Sparse exponent supports are a
different sublattice and need not attain this full image.

The continuation operation P(X) -> T*p^J*X^delta*P(p^a X) descends whenever
T is an integer, J,delta>=0, and every queried source exponent a+r lies in
1..S. Its image at destination node p^r uses exactly a retained source
observation and integral multiplication. This is the operation-specific
kernel-stability proof. At p=2 these are the nodes and transitions used in
the signed clone-excess implementation; extending that native recurrence
itself to other primes is not claimed here.

The map from labelled trie data to these invariant factors keeps the named
evaluation nodes, their order, and pairwise depths. It loses unit residues,
polynomial coefficients and any native counting semantics. Those must remain
available for actual transition arithmetic. In particular an arbitrary tree
shape alone does not determine a general evaluation lattice, and this proof
does not assert a formula for arbitrary node sets or Hermite data.

## Independent finite controls

    python -B check.py > normal.json
    python -O -B check.py > optimized.json
    cmp normal.json optimized.json

Both runs agree. The checker verifies 588 diagonal/LCA evaluations at primes
2,3,5; 58,968 horizon-compatible kernel-transition evaluations; and 72
exhaustive full coefficient-lattice image cardinalities at primes2,3.
The source uses explicit exceptions so its checks survive optimized Python.
This independently validates the observation-module and LCA bridge, not the
entire native signed clone recurrence or a runtime advantage.
