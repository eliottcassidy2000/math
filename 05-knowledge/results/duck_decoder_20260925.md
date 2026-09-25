# Four colors, ten divisor letters, and decoders that retain their carries

**Date:** 2026-09-25. **Status:** PROVED scoped constructions; FINITE-EXACT
controls; inherited results identified below. The proposed global Collatz
certificate and square-sum/graceful-tree equivalence remain **OPEN**.
No claim of literature priority is made for the elementary group actions,
Fourier count, Fibonacci charge, or block transform.

## Inheritance and working board

Closest proved mechanism: [THM-3339, Fibonacci three-ray Berggren transplant
and moving-owner obstruction](../../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md).
It already joins Fibonacci parameters, primitive triples, three matching
channels, and a four-state translation sidecar. The
[ternary residue clock](ternary_bridge_20260925.md) separately identifies
decimal repunit residues with sparse inverse-fibre Berggren heights.

Canonical hostile: the three matching channels do not select a four-state
origin in THM-3339. Corrected near miss: canonicalizing every Collatz sibling
to its smallest representative discards available certificates, as the
[241 and483 examples](creative_sibling_20260925.md) show. Least-used sidecars
here are overlap carries, entry/exit ports, and the full family of pairings.
The historical Fibonacci-depth-reset claim was already retracted; the
[Zeckendorf companion](duck_zeckendorf_20260925.md) recovers that boundary.

Anchor: compatible arithmetic decoding after tripling. Niche: ambiguity of
representations and diagonal color carries. Wildcard: construct an actual
C3 action whose cardinality is the decimal sequence, then exploit its
four-color and K5 realizations.

| Live concept | Preserved object | Essential extra coordinate |
|---|---|---|
| Fibonacci representations | value and two-register charge for distinct atoms | overlap carry for addition |
| Ternary residue trees | refinement, valuation of differences | ordinary height and arithmetic guard |
| Four color states | XOR and cyclic permutation of three nonzero states | distinguished basis/owner |
| Ten divisor letters | C3 action and7+3 split | actual prime labels and multiplicities |
| Tournament compression | module quotient or directed reachability | child blocks, ports, or block contrasts |
| Root certificates | an actual finite arithmetic orbit | exact halving clocks and a terminating dependency |

The research moves used are controlled forgetting, recovering the missing
coordinate, and testing a quotient before treating its size as a rank.
No new META-PATTERNS card is promoted from this session alone.

## 1. What the prime exclusions actually produce

The first twelve primes sum to197. Their sum after deleting11 is186;
after deleting2,3,11 it is181. If instead one asks for the first twelve
**retained** primes, the sums are227 and312 respectively.

More strongly,196 is not a prefix sum of any of these three ordered lists:

| Ordered prime list | Consecutive prefix sums straddling196 |
|---|---|
| all primes | 160 < 196 < 197 |
| omit11 | 186 < 196 < 227 |
| omit2,3,11 | 181 < 196 < 222 |

There are genuine nearby square identities. After deleting2,3,11, the first
three retained primes sum to25 and the first eight sum to144. Deleting those
three primes at a fixed cutoff beyond11 subtracts16. If both unfiltered and
filtered sums were squares, their difference16 would force25-9 or16-0;
neither is compatible with that cutoff. This is an exact fixed-cutoff
obstruction, not a claim about arbitrary square prime sums.

The [prime companion](duck_primes_20260925.md) gives the definitions, direct
checks, and distinctions from arbitrary subset sums. The earlier displayed
decimal run has seven initial prime terms31 through33333331; none of these
prime-sum changes alters that statement.

## 2. An explicit four-color / ten-letter / ternary construction

Let V=F2^2={0,R,G,B}, with R+G=B and addition interpreted as XOR. Choose the
order-three linear map rho cycling R,G,B and fixing0. Let

    D = Sym^2_set(V) = {{u,v}:u,v in V},

where unordered pairs **allow repetition**. This is the symmetric square
of a four-element set, not the symmetric tensor square of a vector space.
It has ten elements: (4^2+4)/2=10, by quotienting ordered pairs by exchange.
The four fixed ordered pairs explain the correction to naive halving.

Split it into

    U = {{R,R},{G,G},{B,B}},                       |U|=3,
    S = {{0,0}} union {{u,v}:u != v},              |S|=7.

The C3 action on D has one fixed letter{0,0} and three free three-cycles.
In particular, the three repeated nonzero pairs are exactly the leftover
diagonal channel. This is the part that raw XOR sends to zero.

For N=p^2*q*r with three distinct primes, the inherited divisor-balance
classification gives F=10,U=3,S=7. Its ten proper nontrivial divisors are
equivariantly identified with this D as follows. Assign p,q,r the three
nonzero color labels c_p,c_q,c_r:

| Divisor type | Unordered color pair |
|---|---|
| single prime i | {0,c_i} |
| squarefree product ij | {c_i,c_j} |
| pqr | {0,0} |
| p^2, p^2q, p^2r | {c_p,c_p}, {c_q,c_q}, {c_r,c_r} respectively |

The last row uses the canonical bijection from the three prime labels to
the nonsquarefree remainder. These letters are **not** the original three
prime divisors; this is a typed correspondence between sets of equal size.
The transported rotation is combinatorial, not an automorphism of the
divisor poset. The [prime companion](duck_primes_20260925.md) proves this
construction and supplies explicit arithmetic hostiles.

Now form length-k words over D and remove the seven constant words whose
letters lie in S. Rotate every letter simultaneously by rho. The only
fixed word before deletion was the constant{0,0} word, so the residual
action is free. Consequently, for every integer k>=1,

    number of residual C3 orbits = (10^k-7)/3.          (1)

Equivalently Burnside gives(10^k+2)/3 before deletion; deleting one fixed
orbit and two free orbits subtracts3. These are rotations of **color
labels**, not rotations of word positions. Equation(1) gives1,31,331,3331,...
and explains both10 and7 from the divisor split. It does not explain or
force primality of those counts.

### The fourth triangular number is now a map, not only a coincidence

Adjoin a fifth vertex infinity to V. Map each distinct pair{u,v} to the
K5 edgeuv, and each repeated pair{v,v} to the edge v--infinity. This is an
explicit C3-equivariant bijection

    Sym^2_set(V) <-> E(K5).

The action fixes vertices0 and infinity and cycles the other three. S is
the six K4 edges on V together with0--infinity; U is the three remaining
spokes. Thus the fourth triangular number10, the7+3 divisor split, the
four-state color action, and the decimal formula share one precise object.
The edge incidence relation is an extra structure; the divisor bijection
does not transport multiplication or a square-sum predicate to it.

### A binary observable inside the ternary count

Give pair{u,v} charge u+v and a word its total XOR charge. The ten letters
have charge multiplicities(4,2,2,2). Character sums on V are10 for the
trivial character and2 for each of the other three. Therefore the numbers
of length-k words of charge0 and of each specified nonzero charge are

    Z_k=(10^k+3*2^k)/4,
    W_k=(10^k-2^k)/4.                                 (2)

Among the removed S-constant words, one has charge0 for odd k; all seven
have charge0 for even k. The residual **orbit** counts split as

    neutral: (10^k+3*2^k-4)/12   if k odd,
             (10^k+3*2^k-28)/12  if k even;
    nonzero: (10^k-2^k)/4 - 2    if k odd,
             (10^k-2^k)/4        if k even.             (3)

For the nonzero part, each orbit visits all three nonzero charges, so its
orbit count is the count of any one charge, not that count divided by3.
The first splits are1=1+0,31=7+24,331=85+246,3331=835+2496.
This supplies a genuine binary/ternary interaction with an explicit parity
guard. It still contains no arithmetic Collatz transition.

## 3. Fibonacci ambiguity has an exact carry repair

The phrase three-color Zeckendorf construction did not identify a unique
current theorem during the bounded repository search. The companion states
the following explicit convention instead of silently inventing attribution.

Use Fibonacci weights f_i=F_(i+2)=1,2,3,5,... . Canonical Zeckendorf
representations are unique. If adjacency is allowed but each atom occurs
at most once, representations need not be unique; nevertheless, their
two-register charge is invariant:

    Q(n)=(n-2*b(n),b(n)),  b(n)=floor((n+1)/phi^2).

The atom charges obey v_0=(1,0),v_1=(0,1),v_(i+2)=v_(i+1)+v_i. Reducing
modulo2 gives R,G,B,R,G,B,... . Every merge of adjacent distinct atoms
preserves charge. A complete proof and independent direct-subset/floor
checks are in [the Zeckendorf companion](duck_zeckendorf_20260925.md).

Allowing repetitions changes the answer:2=2 and2=1+1 have raw XOR charges
G and0. The same issue can arise when adding two canonical representations.
For example, on the strict diagonal a+b=5,

    c(1)+c(4)=R+G=B,
    c(2)+c(3)=G+B=R=c(5).

So a diagonal-constant coloring fails, but it has an exact repair. Define

    delta(a,b)=b(a)+b(b)-b(a+b) in {-1,0,1}.

Then

    Q(a)+Q(b)-Q(a+b)=(-2*delta,delta),
    c(a)+c(b)=c(a+b)+(0,delta mod2).                    (4)

The carry obeys the cocycle identity

    delta(a,b)+delta(a+b,c)=delta(b,c)+delta(a,b+c).

Thus a single carry bit repairs XOR color for every pair, and the signed
carry repairs the full charge. This is stronger than merely selecting a
canonical representative. The correction channel(0,1) is distinguished by
the initial weights1,2: a color rotation must also rotate this basis/owner.
Three nonzero colors form part of a four-state XOR group, not Z/3, and
some positive integers have neutral color (for example6).

### How far the ternary correspondence goes

The inherited clocks R_b(j)=(b^j-1)/(b-1), for b=4,10, satisfy
v3(R_b(j)-R_b(i))=v3(j-i). A Collatz inverse fibre selects ray heights
H(j)=2^(k0-1)*R_4(j); R_10 selects the decimal clock. Their finite
residue-refinement trees are isomorphic, while the ordinary heights remain
sparse. The color action above is an order-three **linear symmetry over
F2**, not an identification of Fibonacci addition with3-adic addition.

THM-3339 supplies a second, genuinely geometric connection: primitive
normalization of Fibonacci parameters produces three interleaved Berggren
rays, with words(BA)^r,A(BA)^r,C(BC)^r. These are not the parabolic inverse
Collatz-fibre rays. Its moving-owner obstruction matches the present need
to retain a color origin, but does not identify the two dynamics. The
remaining comparison must retain the value/height chart and carry at every
step; matching the number of colors alone supplies no commuting map.

There is an exact **finite-residue repair**. Represent an index j modulo
3^d by d color digits, least significant first, encoding ternary0,1,2 as
R,G,B. Increment the first digit and carry to the next precisely on B->R.
This color odometer has order3^d. The charts R_10 and
2^(k0-1)*R_4 send it to the two inherited residue clocks; the digit carry
is the missing transition data. This is a genuine commuting finite-level
decoder, and the charts commute with reduction to smaller depths.

At depth2 it cannot be replaced by simultaneous color rotation: that
operation has order3, whereas the clock has order9. After three correct
increments,(R,R) becomes(R,G), not(R,R). The Fibonacci atom color
rho^i(R) records only i mod3; it does not retain the higher index digits.
The Fibonacci addition carry(4) and this ternary odometer carry are two
explicit, different carry laws. A common color alphabet does not equate
them. The new quotient-word count in(1) uses simultaneous rotation and
is therefore also distinct from this odometer.

## 4. Use tournament structure without demanding an unavailable pair

The [tournament companion](duck_tournament_20260925.md) proves two concrete
escapes and audits their limits.

* In Q[H,H,H,1], with regular odd H of orderq>=3, the three q-vertex blocks
  are intrinsically recoverable strong modules. The quotient recovers Q
  and the three child tournaments, even when Q is cyclic. This is a
  zero-reversal inverse of tripling. Numerically N->(N-1)/3 is a guarded
  **reverse construction edge**, not a forward Collatz step; certification
  uses an actual supplied child orbit containing3q+1=N.
* Contract a strongly connected block while retaining both crossing
  directions. The quotient is a semicomplete digraph and may have digons.
  Directed reachability lifts exactly, with entry/exit ports available for
  witness paths. A strong four-core can contract a cyclic triangle to a
  two-vertex digon. Replacing the digon by one arbitrarily oriented arc
  loses a direction. Hamiltonicity additionally requires path-cover data,
  handled by inherited
  [THM-3121, path-cover/walk-content kernel](../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).

Cut switching does not manufacture a pair: the companion proves that every
cut gauge of Q[H,H,H,1] still has no two-vertex module. Nor is a canonical
perfect matching a free choice. The automorphism rotating one C3 block and
fixing everything else has no invariant matching: an internal edge would
force all three triangle edges, and an external edge would force three
edges sharing its outside endpoint. One can retain the full matching
family equivariantly, but cannot select a member equivariantly in general.

### A tested tournament-specific shortcut and its precise failure

Every tournament has an odd number of Hamiltonian paths; see
[THM-002, OCF](../../01-canon/theorems/THM-002-ocf.md) and
[Irving--Omar's primary account](https://arxiv.org/html/2412.10572v3).
For an even tournament, pair successive vertices of each path and sum
these matching incidence vectors over F2. This gives an intrinsic symmetric
zero-diagonal matrix W. Every vertex has odd degree in its support, since
each path contributes one incident matching edge. The construction is
equivariant and avoids selecting a path.

The hope that W is always nonsingular is **REFUTED** at six vertices.
Use vertices0..5, every edge pointing from larger to smaller except0->4.
There are nine Hamiltonian paths. W has rank4, with null vectors e0+e4
and e2+e5. Exhausting all32768 labelled six-vertex tournaments gives ranks
2,4,6 in counts1680,17520,13568. Orders2 and4 always have full rank.
The surviving universal invariant is an odd-degree edge-chain, not a
symplectic halving. This is a different kernel from
[THM-2290's hafnian/Pfaffian obstruction](../../01-canon/theorems/THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete.md),
although both warn against erasing the pairing family at order6.

There is also a precise K3,3 connection, with a useful hostile control.
Any rank-two symmetric zero-diagonal F2 matrix with odd row degrees on an
even number of vertices is the adjacency matrix of K_(a,b) with a,b odd.
To prove this, factor its alternating rank-two form through F2^2. Each
vertex has a nonzero coordinate because it has odd degree. Vertices with
equal coordinates are nonadjacent and all unequal-coordinate pairs are
adjacent. Thus the graph is complete multipartite with at most three
parts. Each nonempty part must have odd size, since its degree is n minus
that size. Even n forces exactly two parts.

In the six-vertex census,960 kernels are K1,5 and720 are K3,3; their first
orientation codes are72 and83 in the script's explicit edge convention.
Thus K3,3 genuinely appears in this decoder obstruction. The planar star
also fails nonsingularity, so nonplanarity is not the obstruction's exact
predicate. The ten-letter K5 realization and this six-vertex K3,3 kernel
are two specified maps, not an identification of all their graph data.

### Exact halving into a richer object

Choose a perfect matching and order the two members of each pair. Keep
the internal arc and, between each two pairs, the full signed block

    [[a,b],[c,d]], a,b,c,d in {+1,-1}.

Its four integer channels are

    s=a+b+c+d, r=a+b-c-d, q=a-b+c-d, h=a-b-c+d.

They reconstruct the block by

    (a,b,c,d)=((s+r+q+h),(s+r-q-h),(s-r+q-h),(s-r-q+h))/4.

Ordinary pair-module contraction is exactly the special case r=q=h=0 for
every crossing block. Carrying these three contrasts gives a reversible
half-sized **block object** even for the raw tripled tournament and cyclic
cores. Erasing them does not: balanced row-striped and checkerboard blocks
have the same mean and the same repair cost but different arcs. The
matching, internal order gauge, and contrasts are the required sidecar.
Keeping all choices gives an equivariant family of such objects.

This halves the number of containers while retaining the original bit
budget. It is not a decreasing complexity rank or a Collatz proof.

## 5. What the common structure suggests next

Across these lanes the useful object is a family of representations with
legal local rewrites and recorded defects of an observable. Fibonacci
addition has an explicitly solved defect: the cocycle(4). Prime
factorizations have unique exponent vectors, whereas factor groupings
carry prime-support overlap. The square-sum/graceful sign transform carries
a vertex-sign gauge plus the permitted vertex interval and edge-label
set. Tournament quotients carry ports or contrasts. Those are concrete
data, not interchangeable metaphors.

The proposed decoder should therefore have four components:

1. an exact arithmetic value map from the representation family;
2. legal local transformations, with clocks/carries retained;
3. transport of actual root witnesses, not only reachability in an
   unrelated graph;
4. a well-founded dependency for supplying those witnesses on every input.

Items1--3 now have explicit local models, including intrinsic triple-block
recovery, reversible block halving, and the previously proved guarded-join
compiler. Item4 is still the unsolved global obligation. The raw block
transform gives a decisive hostile to using container count alone. The
inverse tripling decoder gives a second: it recovers the smaller parent
but does not automatically supply that parent's root certificate or make
the result compatible with E(N/2).

For square-sum/graceful transport, the next admissible object is a
constraint-labeled representation family, retaining vertex values and
edge predicates along every sign change. The existing spanning-star and
endpoint-label hostiles show why neither Hamiltonian-path existence nor
the three-color quotient can certify the other problem. No new implication
between those conjectures is asserted.

## Reproduction and audit

Run from the repository root:

    python3 04-computation/duck_decoder_20260925.py
    python3 -O 04-computation/duck_decoder_20260925.py

The [script](../../04-computation/duck_decoder_20260925.py) uses explicit
checks that remain active under optimization. The
[frozen output](duck_decoder_20260925.out) records all10^k words for1<=k<=4,
all even tournaments of orders2,4,6, all sixteen signed2x2 blocks, and all
3/945 matchings at orders4/10, plus the three compatible clocks at ternary
depths1..5. Independent direct path permutations check
the DP through order4 and the six-vertex hostile. In a second all-order6
route, each of720 Hamiltonian orders is completed to its1024 compatible
tournaments; all32768 kernels agree with the subset DP. A separate agent
independently enumerated the four charge-orbit counts in(3) and reproduced
the full rank census by path completion and column-pivot elimination. The companion
notes state their own universes, controls, and reproduction commands.
