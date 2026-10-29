# Reciprocal composition, Lucas prime clocks, and observable transport

2026-10-28. **PROVED elementary statements P1--P6**, **CITED** classical
inputs where specified, and **FINITE-EXACT** for enumerated controls.
No external priority or independent-agent review is claimed. In particular,
the trace identity, finite-field order laws and group-action descent are
classical mechanisms; the contribution is the assembled interface, exact
controls, and the resulting distinctions between minimality questions.

## P1. Reciprocal traces have one joint branch bit per tree edge

Work in a field of characteristic other than2. For a nonzero u put
t(u)=u+u^(-1). Then t(u)=t(v) iff v is u or u^(-1): multiply the equality
by uv and factor (u-v)(uv-1)=0. Assume u²!=1 when asking for two distinct
branches. For x=t(u), y=t(v), the two possible product traces are

    z_+=t(uv),       z_-=t(u/v),
    z²-xy z+(x²+y²-4)=0,
    (z_+-z_-)²=(x²-4)(y²-4).

These identities follow by expansion. Thus two scalar traces determine an
unordered pair of product traces, not a single product trace. The two values
are distinct precisely when both reciprocal pairs are nondegenerate. If u
or v is self-inverse, the branch bit is unnecessary. Example: traces3 and7
of phi² and phi⁴ permit product traces18 and3. The difference is the relative
choice of conjugates, not a failure of multiplication to be associative.

**Tree reconstruction theorem.** Fix k reciprocal pairs {u_i,u_i^(-1)} with
u_i²!=1. Observe all individual traces and all Laurent-word traces
t(product_i u_i^(a_i)), a_i in Z. Among the2^k possible choices of roots,
exactly simultaneous inversion remains invisible. The joint fiber has
2^(k-1) classes. Pair-product traces on any spanning tree distinguish them.

Proof: choose reference roots. Each tree-edge trace distinguishes equality
versus inequality of its two branch signs by the nonzero discriminant above.
The k-1 edge signs reconstruct all signs after one arbitrary root sign is
chosen. Simultaneously reversing all signs inverts every Laurent word and
preserves its trace. Any two assignments not related in this way disagree
on a tree edge and are separated. This proves both sufficiency and the exact
fixed-length lower bound k-1 bits over the fixed individual traces.

The global sign is a genuinely free bit only for this inversion-invariant
observer menu. A named root value or an oriented comparison destroys that
freedom. On a graph with cycles, edge signs must have product+1 around every
cycle. These are realization guards, not extra independent degrees of
freedom. An LCA representation stores tree-edge signs and multiplies along
the unique path; it need not store all pair signs.

If previously retained relations already specify consistent relative signs
on a graph with c connected components, only c-1 further bits are needed.
This is when an alignment bit is free. For example a retained relation
v=u² rules out the misaligned product in the3,7 example. One must not forget
that relation and then claim the bit was inferred from the two traces alone.
Writing vertex signs as bits, the edge signs are s_i+s_j over F2: this is
the graph's cut code, with cycle parity as its realization guard. This
particular alignment problem is graphic; it does not turn arbitrary
arithmetic square-class matroids into graphic matroids.

For a quadratic unit with conjugate c/u, the analogous single trace is
u+c/u and powers obey T_(r+1)=T_1 T_r-c T_(r-1). Keep c: the golden unit has
c=-1, whereas phi² has norm+1. Erasing norm parity confuses ordinary Lucas
traces with reciprocal traces at odd indices. P1 is stated for norm+1 to
avoid that change of type.

The field hypothesis also matters. The concurrent
[trace-clock theorem](../session93_reversal_lucas_clocks_20261029/arithmetic/TRACE_CLOCK_CARRIERS.md)
has four norm-one preimages of trace 3 modulo 33, from independent local
inversions. Modulo 5 an entire trace history can hide a five-state nilpotent
fiber. Neither is a counterexample to P1 over a field; both prevent extending
its two-branch statement to arbitrary quadratic rings without CRT and
ramification guards.

Finite control: both path and heap tree charts, all2^k lifts for k=1..8,
plus49 rational pair-polynomial checks. The proof is not restricted to these
examples or to characteristic zero.

## P2. Lucas divisibility is an exact projective half-turn condition

Let L_0=2,L_1=1, L_(n+2)=L_(n+1)+L_n, and

    M=[[1,1],[1,0]],       L_n=tr(M^n),
    A_l=F_l[X]/(X²-X-1),  phi=X, psi=1-X, alpha=phi/psi=-phi².

For a prime l other than2,5 this is an etale quadratic algebra. It is split
when epsilon=(5/l)=+1 and a quadratic field when epsilon=-1. Conjugation
interchanges phi,psi and sends alpha to alpha^(-1). In the inert case that
conjugation is Frobenius. Therefore d_l=ord(alpha) divides l-epsilon.

For every n>=0,

    l | L_n  iff alpha^n=-1
             iff d_l is even and n = d_l/2 modulo d_l.             (P2)

Indeed L_n=psi^n(alpha^n+1), and psi is a unit. In the cyclic group generated
by alpha, -1 occurs exactly halfway around an even cycle. In the split
algebra its two coordinates are inverses; its unique order-two element is
(-1,-1), so the same argument applies there.

The projective action x->1+1/x of M on P1(F_l) has exact order d_l. It has
two fixed eigenlines in the split case and none in the inert case; every
other cycle has length d_l. For the latter assertion, a nonidentity power
fixing a rational point would supply an eigenline for the same nonsplit
torus, which is impossible. Thus P2 says precisely that M^n is a nontrivial
projective involution. This is a native arithmetic clock, not an arbitrary
tournament orientation imposed on residues.

Exceptions are retained: L_n is even iff3 divides n, and L_n never vanishes
modulo5 (its period-four values are2,1,3,4).

### Two different eleven-point-related involutions

The biplane frame map R:x->6/x acts on twelve points P1(F11). It has no
fixed points because6 is nonsquare. Its determinant class is -6=5, a square,
so it lies in PSL2(F11). By contrast

    M^5 mod11 = [[8,5],[5,3]]

fixes the two golden roots4,8 and has nonsquare determinant -1. It is outside
PSL2(F11). These involutions are not conjugate even in PGL2(F11). The compact
frame symmetry is not literally the Lucas-divisibility half-turn.

For general primes l other than2,3, x->6/x has1+(6/l) fixed points and
(l+2+(6/l))/2 orbits on P1(F_l). Membership in PSL2 is instead controlled
by(-6/l). Those are honest quadratic prime-splitting predicates, but neither
is the exact golden order d_l. This separates a useful shared projective
language from an invalid identification of arithmetic actions.

### Prime-index and power-of-two lanes

If an odd prime l divides L_n for odd n, then l splits in Q(sqrt5). This
follows from L_n²-5F_n²=-4: modulo l, 5F_n²=4, and F_n is nonzero. If n=p
is an odd prime, P2 gives d_l=2 or2p. The first case would force phi²=1 and
then phi=0 from its quadratic relation, impossible. Hence

    l | L_p, l odd  ==>  d_l=2p and l=1 modulo2p.

The prime p=3 has L_3=4 and no odd prime factors. The criterion on l is
necessary, not sufficient: an allowed residue class does not fix d_l.

For n=2^k, k>=1, phi^(2n)=-1 whenever l divides L_n, so phi has exact order
4n. In the split case this gives l=1 modulo4n. In the inert case
phi^(l+1)=-1, hence l=2n-1 modulo4n. These two factor classes are kept
separate. There is no claim that all such factors occur or that every
power-of-two index gives a prime value.

The identity L_a | L_(am) for odd m follows from factoring x^m+y^m.
Consequently if L_n is prime and n>=2, then n is a prime or a power of two.
For an odd composite n choose a proper divisor a>1; for an even n with an
odd factor use its power-of-two part. In either case L_a>1 is a proper
factor. This is another necessary gate, not a prime-value classification.

Let r=d_l/2 when d_l is even. For every positive odd m the lifting law is

    v_l(L_(rm)) = v_l(L_r)+v_l(m).

In the unramified quadratic local algebra, apply the usual binomial/LTE
argument to (-alpha^r)^m-1, where -alpha^r=1 modulo a prime over l. The
remaining unit psi^(rm) changes no valuation. Split conjugate places give
the same valuation for the rational integer L_(rm). The initial valuation
v_l(L_r) must be retained; no unproved absence of exceptional lifting primes
is assumed.

**CITED input, scoped consequence.** Bilu--Hanrot--Voutier, Theorem1.4,
gives a primitive divisor of F_(2n) for n>15. Since F_(2n)=F_n L_n, that
divisor occurs in L_n and in no earlier L_j. This proves an endless supply
of prime divisors, not an endless supply of prime Lucas values. Their term
"Lucas sequence" means the U-sequence, here F_n; this distinction matters.
See [the primary paper](https://www.researchgate.net/profile/Paul_Voutier/publication/2397582_Existence_of_Primitive_Divisors_of_Lucas_and_Lehmer_Numbers/links/559999c308ae99aa62cc6983/Existence-of-Primitive-Divisors-of-Lucas-and-Lehmer-Numbers.pdf?origin=publication_detail).

Finite control: all166 primes l<=1000 other than2,5 and n=0..300;
460 divisibility incidences; full trial factorizations of L_2..L_60;
declared odd-multiple valuation probes. L_23=139*461 is a hostile to a wheel
of all primes<=100 even after the prime-index gate. P2c below gives a separate
all-size statement about what a fixed bank can observe.

### P2b. The exact minimum for a named finite divisor bank

Observe, for each prime l in a fixed finite bank, the Boolean l|L_n, and
permit n->n+1 (or any specified integer affine index maps). For l!=2,5,
the individual minimum has d_l states when d_l is even; otherwise the
output is constantly false and its minimum is one state. For2 the minimum
is3; for5 it is one. These are behavioral states, not merely distinct
present-time Boolean outputs.

The exact joint minimum is M=lcm of the active periods. Sufficiency is
n modulo M. Necessity: an individual nonconstant test has one hit position
per period, so any two distinct phases are separated by a future increment.
If n,m have the same entire vector of future outputs, their difference is
divisible by every individual period, hence by M. All M phases occur.
Affine index maps also descend to this carrier. This proves minimality for
the named test vector, not for the weaker OR of its coordinates.

The independent product of clocks is usually larger because all clocks
read the same n. Exact data on overlaps supplies the missing correlation:
the reachable tuples are the diagonal cyclic subgroup, equivalently the
compatible generalized-CRT tuples. Examples verified for their full periods:

| Divisor bank | Individual periods | Joint minimum | Independent product |
|---|---|---:|---:|
|2,3,11|3,4,10|60|120|
|2,5,7|3,1,8|24|24|
|4,2,5,3,7,13|6,3,1,4,8,1|24|576|
|all primes<=29|3,4,1,8,10,1,1,18,24,14|2520|5806080|

The composite test4|L_n is exactly n=3 modulo6, by the recurrence modulo4.
The tests for5,13,17 are constantly false. Thus the familiar two prime
banks differ for a precise reason, and the larger mixed bank has24-fold
overcounting if its coordinates are treated independently.

This is a robust analogue of minimizing a metagraph observer. To decide
primality rather than named divisibility, retain the value/height and the
self-divisor exception L_n=l. A certified trial sieve must reach sqrt(L_n).
Neither this finite clock theorem nor the exceptional PSL2(11) frame gives
the frequency or infinitude of Lucas prime values.

### P2c. Divisor composition and eventual silence of every fixed sieve bank

Put A = {prime indices} union {2^k : k >= 1}, with repeated indices removed.
These are all possible indices of prime values L_n for n >= 2, although many
values at these indices are composite. Distinct values indexed by A are
pairwise coprime.

Here is a proof directly from P2. For an odd prime l dividing L_n, the first
hit r=d_l/2 divides n and n/r is odd. The value r=1 is impossible because
L_1=1. If n is prime, r=n. If n is a power of two, the odd-quotient condition
again forces r=n. A fixed l cannot therefore divide two different candidate
values. The prime 2 divides exactly the candidate L_3=4. This proves the
pairwise-coprimality statement without a statistical independence hypothesis.

Consequently every fixed finite bank of divisors greater than one has an
all-size failure of *coverage*: each divisor hits at most one candidate,
and all its named outputs are false at every sufficiently large index in A.
For the bank of all primes <=29, the only hit indices are {2,3,4,5,7};
L_16=2207=31*71 is already a composite survivor. The bank <=100 also misses
L_23=139*461. These examples are witnessed false positives, while eventual
silence is proved for every fixed bank. We do not deduce from this alone
that infinitely many candidate values are composite or prime.

The restriction to A is not closed under n->n+1. Thus this theorem does not
replace P2b's full-index minimal observer by a one-state increment automaton.
It explains exactly why a perfect finite modular observer for divisibility
can become uninformative on the prime-value target after the index filter.
The growing data must include new divisor clocks, higher lifting data when
needed, and the height/self-divisor guard.

For completeness, P2 and its valuation law give the full positive-index gcd
composition rule. Let g=gcd(m,n). Then

    gcd(L_m,L_n) = L_g, if m/g and n/g are both odd;
                   2, if one quotient is even and 3 divides g;
                   1, otherwise.

If both quotients are odd, every common odd prime has first hit dividing g
with odd quotient, and LTE identifies the minimum valuation with v_l(L_g).
If the quotients differ in parity, the 2-adic orders of m,n differ, so no
odd prime clock can hit both. The remaining 2-adic calculation is

    v_2(L_n)=0 if 3 does not divide n,
             2 if 3 divides n and n is odd,
             1 if 6 divides n.

For odd n this follows from the recurrence modulo 8 (period 12). For even
n iterate L_(2n)=L_n^2-2(-1)^n; once an even index is reached the valuation
is exactly one. This proves the displayed gcd rule, including its exceptional
factor 2. It is a classical type of divisibility identity, not a claim of
priority. Exact controls cover all 7,260 unordered pairs in 1..120 and all
distinct admissible indices through 300.

## P3. Three eigenvalues can support exponentially large return rank

For the constraint star K_(1,m), use the inherited native Gram
D=PB(1-P)BP, where B sums all single-bit flips and P keeps independent sets.
The independence complex is a leaf simplex disjoint from the centre point,
hence has homotopy type S0 for every m>=1. Its legal state count is2^m+1.

The empty state has zero return. On singleton states, D is the signless
Laplacian of the star: centre diagonal m, leaf diagonals1, centre-leaf
entries1. Its spectrum is0 once,1 with multiplicity m-1, and m+1 once.
Each leaf subset of size at least2 has a unique forbidden insertion of the
centre, with no other legal predecessor. Those2^m-m-1 states contribute
independent diagonal1 entries. Thus, on the full legal space,

    spec(D) = {0^(2), 1^(2^m-2), (m+1)^(1)},
    rank(D)=2^m-1,
    D(D-I)(D-(m+1)I)=0.                                         (P3)

For m>=2 the minimal polynomial has degree3; for m=1 it has degree2.
The homotopy type, number of distinct eigenvalues, rank, and dimension of a
specific observed recurrence are four different quantities. A compact
formula for D^r is compatible with an exponentially large rank needed by
any exact real Gram factorization for arbitrary source amplitudes.

The native flip construction is checked for every m=1..9, including the
matrix polynomial on every coordinate vector. Independent modular
elimination supplies rank lower bounds through m=6. This strengthens the
inherited rank example; it does not claim a general circuit lower bound.

## P4. The Lucas graph's cancellation extends to a six-class rule

Let Lambda_n, n>=3, be the toggle graph of independent sets of C_n, with
each edge oriented by adding one occupied position. The empty set is the
unique degree-n vertex: a nonempty independent set S has degree
n-|N(S)|<n. Graph distance from empty is exactly |S|. Every automorphism
therefore preserves weight and the bit-addition orientation. All Lambda_n
are Royle sign-even, regardless of the parity of their automorphism orders.

Counting insertions at one fixed position leaves a path of n-3 positions,
which has F_(n-1) independent sets. Thus

    |V(Lambda_n)|=L_n,   |E(Lambda_n)|=n F_(n-1),
    |E(Lambda_n)| odd iff n=3 or5 modulo6.                         (P4)

For a tournament T large enough to contain copies, sum the products of arc
signs over all embedded unlabelled Lambda_n copies. The reference orientation
is respected by graph automorphisms, so this is well defined. Converse
multiplies every copy sign by(-1)^(n F_(n-1)). An anti-automorphism is a
bijection of copies. Hence the signed sum vanishes on every self-converse T
in the odd-edge cases, even if a particular embedding is positive.

At n=5 this recovers the inherited eleven-vertex cancellation, with
Aut(Lambda_5)=D5 of order10 and15 edges. At n=7 there are56 edges and no
forced cancellation from converse. At odd prime index p>3 the forced-zero
lane is p=5 modulo6. Sign-evenness of the source graph and converse parity
of its signed invariant are distinct predicates. Finite controls cover
all toggle graphs at n=3..13.

## P5. The nine-vertex square-root rule needs a trit for moved ports

Use the inherited W=C3[C3,C3,C3], with vertices3i+j and anti-involution
t0(i,j)=(-i,-j). Its full automorphism group has order81, and the centralizer
C of t0 has order3. For reversal markings s,t, the inherited natural map is
g(t,s)=sqrt(ts), using the unique odd-order square root.

Let H be the subgroup generated by all closed-path transports based at t0.
A vertex observable f is independent of the chosen transport path iff it
is constant on H-orbits. Necessity follows by traversing a loop; sufficiency
follows because any two paths differ by a loop. This is the exact descent
criterion for this observer, not a statement that all transports are flat.

One actual triangular holonomy is

    h=(0,1,2,4,5,3,8,6,7).

It has order3 and lies in C, so H=C. Accordingly the finest vertex observer
that needs no history has the five cells

    {0}, {1}, {2}, {3,4,5}, {6,7,8}.                              (P5)

If a retained port lies in either three-orbit and its exact image is to be
observable, the residual transport requires three distinguishable values.
One Boolean bit cannot encode those three values; a trit, or a coherent
reference frame removing path ambiguity, suffices. This is different from
P1's two-branch reciprocal ambiguity and from the Bruck loop's27 markings.

The checker directly builds81 actual automorphisms,27 anti-involutions and
all729 based triangles. Their three holonomy frequencies are297,216,216.
The all-size flatness and Bruck-loop construction remain credited to the
[prior transport proof](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md).
The twelve-vertex81/27 pair remains a distinct hostile, not this W.

**Current all-size refinement, credited.** The incoming
[recursive centralizer proof](../session93_reversal_lucas_clocks_20261029/holonomy/FULL_CENTRALIZER_GENERATION.md)
and [independent full-holonomy packet](../lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md)
now establish full centralizer holonomy throughout W_k. The
[joint folded-LCA theorem](../session93_reversal_lucas_clocks_20261029/ports/PORT_ORBITS_AND_LCA.md)
gives 2k+1 vertex orbits, 8k²+1 ordered two-port orbits, and the exact minimum
(3^(k-1)-1)/2 named anchors. P5 is the k=2 vertex specialization. Thus it
does not leave recursive-family generation or tuple classification open;
the remaining attachment question concerns other atoms and operation menus,
particularly native complement which changes the underlying adjacency.

## P6. A precise common theorem, with a free-action boundary

Let G be a finite group. Suppose k locally chosen frames are represented by
g_1,...,g_k in G, and changing the common frame acts by simultaneous left
multiplication. Assume the target can distinguish every relative frame
g_i^(-1)g_j. Then the exact joint carrier, over the given local data, has
|G|^(k-1) states. Relative frames on any spanning tree are sufficient and
minimal. If consistent known relative frames connect the sites into c
components, the remaining carrier has |G|^(c-1) states.

Proof: root the tree and set its frame to the identity. Multiplying the
oriented edge labels along its paths reconstructs every other frame.
Changing the initial root frame changes only the common left multiplier.
Every tree-edge assignment is realizable, giving |G|^(k-1) possibilities;
by hypothesis the target separates them. Known edges reconstruct each
component internally, leaving one free frame per component. For a full
edge-labelled graph, realization by vertex frames is equivalent to ordered
product equal to the identity around every cycle. This argument works for
noncommutative G, provided the order of products and inverse orientations
are retained.

This is a standard finite gauge/transport mechanism, not an external
priority claim. P1 is its C2 case. It explains free-frame Frobenius/LCA
compression when the action really is free, and the compatible-CRT example
shows how a shared index can force additional relations before minimization.

There are two essential boundaries. First, an atom with stabilizers has
transporters rather than freely chosen group coordinates: retain joint
incidence and stabilizer actions instead of replacing it by a torsor.
Second, the natural square-root connection in P5 is not flat. Its triangle
product is actual observable holonomy, not an illegal input to be rejected.
One may choose a different coherent reference-frame rule, or quotient the
observer by its holonomy action, but may not enforce flatness and claim to
have preserved the original transports. Finally, no theorem here says all
arithmetic realization constraints are graph-cycle constraints; the repo's
Fano square-class example is precisely a boundary to that extension.
