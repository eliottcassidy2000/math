# Two triangular prime pictures and a priority-preserving PRIMEGAME chart

**PROVED M1–M4; independently audited PASS.**
[Independent audit](../audit_machines/REPORT.md). Finite checks are
**FINITE-EXACT**. The classical prime-generation theorem, Euclid parameters,
and Conway trial costs are inherited, not claimed as new. No Goldbach or
twin-prime infinitude result is asserted.

## Inheritance and the active board

The closest mechanism is Conway's exact division phase, recovered in
[THM-4853 / PRIMEGAME step counts](../../../../01-canon/theorems/THM-4853-primegame-step-counts-in-closed-form-oeis-a007546-is-the-15-14-list-a007547-the-conway-1987-list-and-their-difference-is-a-sum-of-largest-proper-divisors.md).
The canonical hostile is executing a legal fraction after an earlier rule
already applies. The corrected near miss is
[THM-4852's November 2 denominator-precision repair](../../../../01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md):
prime support alone does not justify a one-digit selector. The least-used
sidecar is **which ordered odd number encodes the trial divisor**.

Anchor: a faithful triangular representation of the prime tests. Niche:
the exact local sieve kernel and its source weights. Wildcard: the odd
register pair as a rational linear chart of every literal PRIMEGAME step.
Board: pair orientation / sum row / divisor ray / priority / source height /
physical clock. A point of the same triangle can carry different predicates;
the inverse map must decide which one is being used.

[Conway 1987, printed pp.4 and9–10](https://gwern.net/doc/cs/computable/1987-conway.pdf)
was read directly for the ordered program and division phase. The existing
[TM/CA dictionary](../../machines_automata_fractran_tournaments_dictionary_20261101/README.md),
[halting programs](../../fractran_turing_circulant_goldbach_20261101.md), and
[exact translations](../../lenses_tm_ca_fractran_tournaments_20261101.md)
are inherited. The new claims here are elementary coordinate and kernel
calculations, not new universality or simulation-size theorems.

## M1. The odd-prime pair triangle, including its boundaries

For unordered odd integers 3<=p<=q set

    h=(p+q)/2,   d=(q-p)/2.

This is a bijection onto integers h>=3, 0<=d<=h-3 with h+d odd;
its inverse is p=h-d, q=h+d. Thus h,d have opposite parity. Under this
map:

- odd Goldbach pairs of 2h are the prime points in row h;
- odd twin-prime pairs are the prime points on d=1;
- prime pairs whose sum is a square lie on h=2t², t>=2.

The prime2 is a separate boundary: an even target has a prime pair using2
only at 4=2+2. A prime sum p+q=r necessarily has one input2, and is exactly
a twin-prime relation 2+p=p+2. This differs from two odd primes summing to
an even target. No orientation into a tournament is supplied by this
symmetric relation.

For distinct odd primes p<q, h>d>=1, gcd(h,d)=1 and opposite parity give
the primitive Pythagorean triple

    (A,B,C)=(pq, (q²-p²)/2, (p²+q²)/2).

Conversely, the standard Euclid factors A=(h-d)(h+d), with h-d and h+d
both prime, recover that pair. The nonunit-factor guard is indispensable:
A=15 has both (15,8,17), from factors3,5, and (15,112,113), from1,15.
This reuses the full odd-leg fibre distinction in the
[prior twin/Pythagorean packet](../../twin_pythagorean_124_20261026/tree/RESULT.md).
The diagonal p=q gives d=0 and is outside positive primitive triples.

### The exact row290 comparison

For h=145, there are56 unordered positive odd pairs p+q=290 with
p<=q and gcd(pq,145)=1. The involution p↔290-p has no fixed unit,
so this is phi(290)/2. One pair is (1,289); requiring p>=3 leaves55.
Exactly ten are prime pairs:

    (7,283), (13,277), (19,271), (61,229), (67,223),
    (79,211), (97,193), (109,181), (127,163), (139,151).

The count56 therefore has an actual pair-folding meaning, but coprimality
is not primality. The first prime1451 congruent to±1 modulo290 belongs to
the distinct increasing-prime arrival clock in
[PRIME-MRB](../../jacobian_eth_prime_wait_20261102/PRIME-MRB.md).
Neither these ten pairs nor their order computes that arrival time.

## M2. The local sieve is an exact weighted tensor kernel

Let M>1 be odd and squarefree. On h,d modulo M define

    S_M(h,d)=1[gcd(h²-d²,M)=1].

The change of cell coordinates (h,d)↔(p,q)=(h-d,h+d) is invertible modulo
M, since2 is a unit. It identifies this predicate with p and q both units.
For each prime ell dividing M, its forbidden set in a fixed h row is

    d=h or d=-h mod ell.

There is one forbidden residue when ell|h and two otherwise. Hence the
exact row size is

    sum_d S_M(h,d)= product_(ell|M) (ell-1 if ell|h, else ell-2),
    sum_(h,d) S_M(h,d)=phi(M)².                              (1)

These are counts over a complete residue rectangle, not densities for
prime pairs at a growing height.

**M2.** Over Q the matrix S_M has rank

    product_(ell|M) (ell+1)/2.                               (2)

Its distinct rows are classified exactly by h² modulo M. For ell prime,
negation orbits of residues have r=(ell+1)/2 members. On those orbits the
kernel is J_r-I_r: two orbits are incompatible precisely when equal.
Its eigenvalues r-1 and -1 show invertibility over Q. Duplicating rows and
columns on their original orbits proves rank r. CRT gives the tensor
product over ell, proving (2) and the distinct-row classification.

For literal sums, the quotient must retain the orbit weight:1 at0 and2
at every nonzero orbit. If C is the quotient kernel and w the column-orbit
sizes, its action on orbit-constant functions is C diag(w), not C. CRT
multiplies these weights. This does not change rational rank; it changes
row sums and the operator's spectrum.

For M=145 the full145-by145 matrix has rank3*15=45. Its quotient is the
unweighted graph K3 tensor K15, with45 vertices and degree28. The original
weighted row sums are112,108,84,81 according as h is divisible by both5,29,
by5 only, by29 only, or by neither. Replacing these sums by28 would erase
source multiplicities. Also the45 local-sign orbits are not the73 orbits
of the single global involution h↔-h on Z/145Z, nor the56 unit-sign
orbits used in the root packet's doubling calculation.

In (p,q) coordinates the unit mask is the outer product u(p)u(q), of
rank1. The rank45 does not contradict invertibility of the change of
coordinates: that change mixes row and column variables and is not a
separate row and column permutation. Matrix rank here belongs to a chosen
separation of ports, not to the bare set of admissible cells.

### The correctly weighted spectrum

For an odd prime ell, the full ell-by-ell symmetric matrix S_ell has
characteristic polynomial

    X^((ell-1)/2) (X+2)^((ell-3)/2)
       [X²-(ell-3)X-(ell-1)].                            (3)

Indeed the functions odd under d↔-d form a zero eigenspace of dimension
(ell-1)/2. On the even functions, use the negation-orbit quotient with
weights(1,2,...,2). The functions zero at0 whose values on the nonzero
orbits sum to0 have eigenvalue-2, with dimension(ell-3)/2. The remaining
two-dimensional space consists of a value x at0 and a common value y at
all nonzero orbits. Its matrix is

    [ 0   ell-1 ]
    [ 1   ell-3 ],

whose characteristic polynomial is the bracket in(3). These spaces exhaust
the carrier. In particular S5 has nonzero eigenvalues1+sqrt5,1-sqrt5,-2;
S29 has13+sqrt197,13-sqrt197 and-2 with multiplicity13. CRT makes S145,
after the same permutation of rows and columns, S5 tensor S29. Its nonzero
eigenvalues are therefore the45 products of those nonzero local values,
with multiplicities; the remaining100 eigenvalues are0.

This is an exact source-weighted spectral calculation. Its zero odd sector
occurs because the sieve predicate deliberately forgets either local sign;
it is not the signed-unit doubling observer, whose odd sector carries
orientation transport. Nor is S145 a tournament adjacency matrix. The
appearance of sqrt5 here follows from the displayed weighted2-by-2 block;
it is not evidence that the native operations of the other lanes coincide.

### Height and self-hit guards

At a fixed row h, primality is decided exactly by testing every odd prime
ell<=floor(sqrt(2h-3)), and retaining a coordinate divisible by ell only
when that coordinate equals ell. These self-hits are necessary: at h=7,
the prime pair (3,11) would be removed by a raw divisibility mask for ell3.
This is the earliest such row. The square-root height bound follows from
the least-prime-divisor bound on a composite positive integer.

No fixed wheel can decide all rows. If Q is the product of its odd test
primes, choose p=(2Q+1)² and q=(4Q+1)². Both are odd composites, are coprime
to Q and to each other, and give a valid primitive Euclid pair (h,d).
Thus even primitivity of the associated right triangle does not repair the
missing prime-height test. This is a concrete obstruction to a finite sieve
being a total prime recognizer.

## M3. Conway's divisor triangle is a different chart on odd pairs

A Conway trial is indexed by n,d with 0<d<n, at integer state

    N_(n,d)=5^n 7^d 13.

Put

    p=2d+1,   q=2(n-d)+1.                                (3)

This is a bijection from these trials to **ordered** odd pairs p,q>=3.
Its inverse is n=(p+q-2)/2 and d=(p-1)/2. The divisor test becomes

    d|n  iff  (p-1)|(q-1).                               (4)

A failed trial sends d→d-1 and therefore (p,q)→(p-2,q+2), staying in
one sum row. Divisor points lie on the integer-slope rays
q-1=k(p-1) from(1,1). This is an exact triangular arithmetic picture of
trial division. Its predicate is not that p and q themselves are primes.

The ordered port matters: (3,5) and (5,3) encode n=3,d=1 and n=3,d=2.
The former divides and the latter fails. Folding to an unordered pair
would make the next branch ambiguous. Further, the same row p+q=290
represents candidate n=144 in this chart, while it represents center145
in M1. A visually identical row is not a shared input clock.

The inherited literal endpoints and clocks are

    d∤n: 5^n 7^d 13 → 5^n 7^(d-1) 13,
           time=6n+2 floor(n/d)+2;
    d|n: 5^n 7^d 13 → 2^n 7^(d-1),
           time=4n+2n/d+2.                               (5)

They are Conway's division mechanism and THM-4853's checked costs. Thus
one triangle move is a macro-step with its own nonconstant integer cost.
The actual prime-generation orbit tests d=n-1 downward and stops at the
largest proper divisor. The full prefix of failed trials is the realization
guard for saying that d=1 certifies a prime.

The guard cannot be deleted: beginning directly at5^4*7*13 produces2^4
after26 legal updates, although4 is composite. This does not contradict
Conway's theorem, whose initial state is2. It is the smallest composite
candidate with this direct-start hostile. No finite quotient of the triangle
may infer prime output from the terminal port alone.

## M4. A literal odd-register-pair chart, with one missing exponent exposed

There is a second, literal encoding that advances every PRIMEGAME update,
not merely the division macro. On the carrier of positive integers supported
on primes2,3,5,7,11,13,17,19,23,29 write uniquely

    N=2^a u v,
    u=3^b 5^c 7^e,      v=11^f 13^g 17^i 19^j 23^k 29^l.

Both u,v are odd. Retain a, and set h=(u+v)/2, d=(v-u)/2 with signed d.
The inverse u=h-d,v=h+d and the support restrictions recover N. This is
a sparse, typed subset of the integer cone h>|d|; it is not the whole
prime-pair triangle. On the actual orbit from2, v is either1 or one control
prime, by the inherited one-hot invariant.

For a reduced instruction A/B split its factors by these same supports:
A/B=2^delta rho sigma, with rho acting on u and sigma on v. Selection
requires a>=max(0,-delta), denominator(rho)|u and denominator(sigma)|v,
**and failure of every earlier instruction**. On that exact cell,

    a' = a+delta,
    [h'] = 1/2 [rho+sigma   sigma-rho] [h],
    [d']       [sigma-rho   rho+sigma] [d].             (6)

The image has integral h',d', and reconstructs precisely A*N/B. This is
immediate by applying rho,sigma to u=h-d,v=h+d. The quadratic form obeys
h'^2-d'^2=rho*sigma*(h²-d²); the power of2 restores the whole fraction.
Consequently (6), with all priority guards and a retained, is conjugate
to literal PRIMEGAME on this carrier, one update for one update.

Two minimal information losses are concrete. States1 and2 have the same
odd pair(u,v)=(1,1), but select55/1 and15/2; a cannot be dropped. The
matrices for11/13 and13/11 are reciprocal, yet the proposed two-step word
need not be realized:39→33 uses11/13, and33 then selects29/33, producing29,
before13/11 is reached. The legal loop13↔11 is a positive control. Products
of fraction matrices alone are not a language of legal words.

This value-plane chart does not produce a radius-bounded CA with one tick
per update. On states N=3^b the first selected rule is55/1. The h-coordinate
moves by2*3^b+5, unbounded in b. A one-token encoding on this literal lattice
with quiescent background cannot make that jump in one fixed-radius CA
step. Digit encodings, exponent counters, extra tracks or additional time
are legitimate alternatives, already treated in the inherited simulations;
their costs and guards must remain attached.

## Connection contract and reproduction

M1 preserves prime-pair incidence and distinguishes rows, twin columns and
square-sum rows. M2 gives an exact local sieve kernel with full orbit weights;
its quotient loses height, self-hits and primality. M3 preserves the division
branch and, with (5), its native clock; losing ordered ports or the incoming
trial prefix destroys that preservation. M4 is a literal guarded conjugacy
on its stated prime-supported carrier, and loses nothing if a,h,d and the
priority list are kept. None identifies a prime pair with a tournament
isomorphism type or proves that a residue trajectory is a prime generator.

Run `python3 -I -B check.py` and `python3 -I -B -O check.py` in this folder.
The receipts declare all finite universes; the source manifest pins the
current correction and inherited mechanisms. Independent audit is required
before status promotion.
