# From order of 2 to a physical route: quotient graphs with integer return shifts

This note separates two logically different proofs: mixing on a finite
quotient, and moving an actual tournament across a prescribed absolute
count. The analogies below motivate a definite certificate for the
second problem, not reductions from unsolved conjectures.

## A finite certificate for absolute placement

Suppose a finite connected directed graph of **marked physical
interfaces** is given, with reverse edges. A physical realization at
vertex v carries an integer coordinate h, read by the distinguished
Hamilton count as r+Mh (M>0 fixed in this paragraph). An oriented edge
e:v→w has signed voltage d(e) in Z, with d(reverse e)=−d(e).
The necessary physical hypothesis is strong but testable: **from any
realization of the source interface at height h≥B, the edge has an
actual lift** to a realization of the destination at h+d(e), preserving
whatever designated +2 choices are needed. Not all heights need be
realized in advance.

Pick a base vertex u and a spanning tree. Let D(v) be the shift
along its tree route from u. For each non-tree edge form the integer
return shift s(e)=D(v)+d(e)−D(w). These numbers generate exactly the
possible shifts of abstract closed walks at u: subtract the tree
potentials and telescope; conversely the tree-edge-tree walk realizes
each generator abstractly. The subgroup is gZ, where g is the gcd of
the s(e). This is a necessary obstruction before checking any guards.

If g=1, the extended Euclidean algorithm gives a finite signed
combination of these walks of total shift +1. Reverse edges turn
negative coefficients into actual walks. Compute the minimum prefix
shift of the resulting walk, say −D₀ (or use zero if all prefixes
are nonnegative). A single physical realization at u of height
h₀≥B+D₀ makes the entire walk executable, since its source at each
step meets the threshold. It ends at u at h₀+1. Repeat the same walk
to obtain every height at least h₀, with the same interface.
For g>1 the same argument gives a +g walk and a lifted half-ray
in one congruence class, not full support. Intermediate shifts need
not be unit shifts. This certificate complements the existing
first-crossing criterion, which needs small intermediate jumps.
Uniform guard control and a physical anchor are essential; allowing
formal negative operations or losing the marked interface invalidates
the conclusion.

Proof design: look for a **finite period of interface changes** whose
return paths have coprime *absolute* shifts, rather than enlarging a
tree of increasingly many isolated witnesses. No such complete
voltage certificate for Chan–Pak is supplied here.

## An exact Artin and ×2×3 countermodel to quotient optimism

Take the affine maps f(x)=2x+a, g(x)=3x+b, with integral a,b. For any
modulus m coprime to 6 put δ=2a−b and d=gcd(m,δ). The finite
positive-word orbit of x is exactly the inverse image of the
multiplicative ⟨2,3⟩ orbit of x+a modulo d, with the translation by a
understood in this statement. Indeed in the y=x+a coordinate the
maps are 2y and 3y−δ. Reduction modulo d leaves precisely the
multiplications. The two words gf and fg have the same linear
part and differ by δ. Inverses of the finite permutations are
positive repetitions, so composing gf with the inverse of fg
supplies translation by δ modulo m and generates every allowed
lift over a residue modulo d. The argument also covers nonunit
input residues.

For a concrete firewall take a=0,b=1, x>0 fixed. Every residue is
reachable modulo **every** m coprime to 6, yet absolute outputs
are power sparse. For each word the output is at least x times
the product of its letters' slopes. Let s>0 satisfy
2^(−s)+3^(−s)=1. In the abstract word tree stop each branch at the
first product >X; stopping leaves have products in (X,3X], and
their Bernoulli probabilities, the inverse s-th powers of their
products, sum to one. Hence there are Θ(X^s) stopping leaves
and only O(X^s) words of product at most X. Here s<1. This
is not a tournament construction, but is a literal model in
which positive moves, commutator translation on every finite
quotient and exceptionally good local coverage do not produce
integer translation. Using an inverse as a positive permutation
power only works modulo m.

Even without inverses, replacing fg by gf inside a positive word
shifts the count by δ times the product of the slopes of the
*subsequent* letters. This is the actual price of a local
commutation. Fixing the word at the same place modulo m does
not certify this price at the absolute scale.

The [repo's older Artin/Collatz audit](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_artin_corrections_20260927.md)
already distinguishes population statistics from orbit coupling.
[Hasse's original parity paper](https://gdz.sub.uni-goettingen.de/id/PPN235181684_0166?tify=%7B%22view%22%3A%22info%22%2C%22pages%22%3A%5B25%5D%7D) concerns whether ord_p(2) is odd, not
whether 2 generates all of F_p*: with p−1=2^s u, u odd, odd
order is exactly 2^u=1; a primitive root has full even order.
For example ord_17(2)=8 while ord_73(2)=9 despite both primes
being 1 modulo 8. In the cyclic group of order p−1, failing the
Artin test means 2 is an ℓ-th power for some prime ℓ|p−1; the
Hasse test instead asks if 2 is a 2^s-th power, with the exponent
selected by the prime. These are coordinated tests, not a fixed
congruence table.
[Hooley's Artin paper](https://doi.org/10.1515/crll.1967.225.209)
gets prime density conditionally under suitable GRH. These are
different finite orbit tests, neither certifies a uniform physical
lift. [Furstenberg's original disjointness paper](https://math.ucsd.edu/~asalehig/F_Disjointness.pdf)
provides a topological rigidity comparison for ×2×3 on the
circle; the measure-rigidity question is distinct. Neither
topological density nor a finite quotient calculation imposes a
prescribed terminal height.

## What the other analogies demand of the certificate

- **Ramanujan and class numbers: a local genus is not a global
  class.** The particular form x²+y²+10z² misses 7 and 21 even
  though both are locally eligible; 7 has z=0 only, and for 21
  the remaining sums of two squares are 21 and 11, neither
  possible. The exception-list coincidence is not a shared
  mechanism. [Ramanujan's original discussion](https://ramanujan.sirinudi.org/Volumes/published/ram20.html)
  and [Ono–Soundararajan](https://uva.theopenscholar.com/files/ken-ono/files/025_8.pdf)
  distinguish genus eligibility from actual representation by
  this form. Likewise a residue is a label on a physical move
  graph, not its component. The gcd of return shifts is a
  concrete candidate for an additional component invariant.
  The repo's [class-number/Paley ranking identity](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md)
  reads a signed ranking statistic, not the Hamilton-count
  spectrum; it demonstrates why the choice of reader matters.
- **Tunnell and BSD: replace a proxy with a certificate.**
  [Tunnell's original paper](https://doi.org/10.1007/BF01389327)
  relates the square of a specific signed theta coefficient to
  the central L-value of a congruent-number curve. The inference
  from vanishing to positive rank invokes BSD; it is not simply
  equality of local counts. The [repo's exact terminals](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/bsd_tunnell_terminals_20261013c.md)
  illustrate the alternate nonvanishing certificate. The
  battlefield for tournaments is similar in shape: prove a
  signed object detects the wanted physical lift, rather than
  interpreting signed evidence as the witness itself.
- **Friedlander–Iwaniec: prove correlation estimates for one
  constrained family.** Their [x²+y⁴ primes theorem](https://annals.math.princeton.edu/articles/13033)
  proves an asymptotic for a thin family rather than deducing
  primes from congruence admissibility. For an exact target the
  analogue is stronger than average mixing: on a *literal*
  family of counts in an interval shorter than a chosen
  modulus Q, a uniform bound Σ_(k≠0)|μ̂(k)|<1 gives positive
  mass at the target residue by Fourier inversion, and hence
  at that integer (when the interval also contains the target).
  If Q is shorter than the interval, the same bound only sees
  a union of absolute levels. A usable tournament analogue
  requires control past that aliasing scale, or a physical
  return walk instead.
- **Hadamard: preserve all cross constraints on the same
  object.** In a ±1 Hadamard matrix, row orthogonality is an
  exact dot-product requirement. In a one-apex tournament,
  simultaneous disjoint mark changes interact via an *exact
  weighted* pairing; zero mod M and zero as integers are
  distinct. The [repo's Hadamard multitoggle result](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance.md)
  gives a parallel example where higher mixed information
  is essential. Import the exact cancellation test, not
  the unsupported assumption that a Hadamard completion
  realizes the required tournament interfaces.

Taken together these point to the same missing object: physically
admissible paths with a nonzero absolute index, together with guards
and a certified finite return. More residue coverage alone does not
make that object available.

## A small class-number example with every local witness

The discriminant −23 form x²+xy+6y² cannot represent 2 over the
integers: with |y|≥1 it is at least 23/4 by completing the square,
and at y=0 it is a square. But (1/2,1/2) represents 2 over every
Z_p with odd p. At p=2, take x=1 and lift y=1 mod 2 as a root of
6y²+y−1, whose derivative is odd. The other integral form
2x²+xy+3y² of that discriminant represents 2 at (1,0). This is
an explicit local/global class distinction, not a reduction to
tournaments. Compare it with the separately proved affine
[prime-power firewall](CRT_FIREWALL.md).
