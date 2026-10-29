# Eleven guarded states, a Paley subgraph, and a golden Gram spectrum

**PROVED exact finite identities; inherited golden reader identified explicitly.**
This note does not identify a physical packing of eleven squares with a graph
merely because both have eleven objects. It supplies an actual operator map
between a cyclic no-adjacent-one language and a finite-field tournament.

Dependencies: the inherited eleven-square/Fibonacci torus reader
THM-4592 `the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`,
its corrected consumers, and [the guard-rank theorem](GUARDED_FLIP_LEAKAGE.md).
The incoming `odd_sign_spectra_trijections_20261028/THEORY.md` is the free-flip
spectral comparison. Source pin90e704a8f35614d2085778d1b518f03ff0d750ca.

## 1. The exact eleven-state reader

Let W5 be the binary words (x0,...,x4) with no adjacent ones cyclically.
These are the independent sets of C5: one empty set, five singletons, and
five independent pairs. Define

    rho(x) = sum_(i=0)^4 x_i 4^i mod11.

Because 4²=4+1 mod11 and 4^5=1, this is the specified golden-generator reader. Coordinates here are indexed in increasing exponent order; reversing the indexing gives the inherited parity-word convention without changing feasibility.
The residues for integer bit masks are

    0→0, 1→1, 2→4, 4→5, 5→6, 8→9,
    9→10, 10→2, 16→3, 18→7, 20→8.

They are all distinct, so rho is a bijection onto F11. This is a concrete
state bijection, not a statement that all operations on F11 are legal.

An insertion of bit i, when the bit is0 and both cyclic neighbors are0,
changes the residue by4^i. The set of these increments is

    {1,3,4,5,9},

exactly the nonzero quadratic residues modulo11. Therefore every positive
insertion is an arc of Paley11, with convention a→b iff b−a is a square.
The legal insertion graph is a DAG graded by number of occupied bits, with
15 arcs. Its undirected graph has a central empty vertex joined to five
singletons; each of the five independent pairs subdivides an edge in a
five-cycle on those singleton vertices. This is a wheel with each rim edge
subdivided.

By contrast, allowing all additions by ±4^i gives every nonzero difference
in F11, hence the complete graph K11 with55 edges. Even allowing only the
positive increments gives the full Paley11 tournament with55 arcs, not the
15 guarded insertion arcs. The empty state has legal degree5 while each
independent pair has degree2. These are cheap decisive operation hostiles.

The reader transports each **guarded** operation exactly. It does not turn
legal-bit dynamics into an unrestricted Cayley action or identify their
spectra. The parent session's even-length fibre note records a separate
failure when the reader itself ceases to be injective.

## 2. Exact golden-ratio blocks in the return defect

For the operators P,B,A,L,D of the companion note, order W5 by particle
number0,1,2. Direct forbidden-neighbor incidence, or the exchange/pin proof,
gives

    D = [0] direct-sum (2I + A5) direct-sum (3I + A5),

where A5 is a five-cycle adjacency matrix in each block. In the singleton
sector, each forbidden adjacent pair joins its two legal predecessors.
In the pair sector, each state has two exchange neighbors and one extra pin,
giving diagonal3. No forbidden state has predecessors at different particle
numbers, so the cross blocks vanish.

Let phi=(1+sqrt5)/2. The singleton block has eigenvalues

    4 once; phi² twice; phi^(-2) twice.

The pair block has eigenvalues

    5 once; (5+sqrt5)/2 twice; (5−sqrt5)/2 twice.

Thus D has rank10, trace25, smallest positive eigenvalue phi^(-2), and product
of its nonzero eigenvalues4·125=500. The exact characteristic polynomial is

    z(z−4)(z²−3z+1)²(z−5)(z²−5z+5)².

This is a genuine golden-ratio occurrence forced by the pentagonal exchange
operator. It is distinct from the golden ratio as growth rate of the word
count, and from any geometric equation in square packing. To identify those
objects, an additional map preserving the claimed metric or operation is needed.

For general cyclic length n>=3, the number of feasible states is L_n and
there are exactly n F_(n−1) legal undirected edges. For a chosen occupied
coordinate, its two neighbors are0 and the remaining path has n−3 vertices,
giving F_(n−1) choices. Every edge has a unique occupied endpoint at its flipped
coordinate. Consequently

    trace D = n(L_n − 2F_(n−1)).

This agrees with the exact rank formula L_n−1−[n even], but trace and rank are
not interchangeable statistics.

## 3. Exact integer and characteristic-five behavior

The two positive blocks have Smith normal forms

    SNF(2I+A5) = (1,1,1,1,4),
    SNF(3I+A5) = (1,1,1,5,25).

The checker verifies every determinantal divisor by taking the gcd of all
minors of each size. Therefore the full integral cokernel is

    Z direct-sum Z/4 direct-sum Z/5 direct-sum Z/25
      ≅ Z direct-sum Z/5 direct-sum Z/100.

The free Z is the empty-state direction. In particular the Gram has rank10
over Q, rank9 modulo2, and rank8 modulo5.

The characteristic-five degeneration is even more concrete. Write P5 for a
cyclic shift. Since A5=P5+P5^(-1), modulo5

    3I+A5 = P5^(-1)(P5−I)²,
    (P5−I)^5=P5^5−I=0.

Hence the pair block is nilpotent of exact index3: its square is the all-ones
matrix and its cube is0. Its ranks are3,1,0 at powers1,2,3. Over R the same
integer block is positive definite. The boundary incidence itself still has
rank10 modulo5; multiplying by its transpose creates a rank8 Gram because
isotropic cancellation replaces positivity.

This does not contradict Session90's Fourier lattice theorem: that change
of basis is invertible modulo every odd prime. The exceptional5 here belongs
to a particular **guarded Gram operator**, not to the coordinate transform.
Nor does nilpotency index3 explain the unrelated cancellation depths5,21,25
of a prior tournament without an explicit intertwiner and matching source.

## 4. Equal rank does not identify the even-length reader quotient

For even n>=4, the inherited golden reader q:W_n→R_n has exactly one
three-word fibre: empty, alternating-even, and alternating-odd. Its linear
pushforward q_* on real word amplitudes therefore has kernel

    span(e_alt-even − e_empty, e_alt-odd − e_empty).

Its rank is L_n−2, exactly the rank of the leakage Gram D. But the latter
kernel is

    span(e_empty, sum_i (−1)^i e_singleton-i).

The two kernels intersect only in zero: the alternating words have particle
number n/2>=2, the singleton vector has particle number1, and the empty state
has particle number0. Thus neither equality of rank nor equality of target
cardinality supplies a factorization. In particular D cannot factor through
q_*, since q identifies an alternating word with the empty word while D kills
only the latter.

There is also a useful exact repair for **linear amplitudes**. The only
maximum independent sets of an even cycle are its two alternating words.
Each has n/2 forbidden insertions and no exchange neighbor at the same particle
number. Their rows of D are therefore isolated diagonal entries n/2. Given
q_*f and the two scalar observations (Df)_alt-even and (Df)_alt-odd, recover

    f_alt-even = (2/n)(Df)_alt-even,
    f_alt-odd  = (2/n)(Df)_alt-odd,
    f_empty   = (q_*f)_0 − f_alt-even − f_alt-odd.

Every other coefficient is already in a singleton reader fibre. The joint map
(q_*,D) is injective, and these two extra linear observations are minimal because
q_* has a two-dimensional kernel. This is over Q or R; division by n/2 is not
uniformly valid modulo a prime. It is not the same claim as a minimal discrete
three-way tag for named-toggle dynamics in the parent session's fibre note.

## 5. Reproduction and scope

Run `python3 check_guard_leakage.py` and `python3 -O check_guard_leakage.py`.
The code checks the reader, all15 insertion arcs, the five-cycle block
polynomials, exact Smith minors, modular rank and nilpotency, alongside the
all-small-graph guard-rank census. No floating-point eigenvalue is used.

The positive result is a guarded map plus a complete return operator. The
negative result is that state equinumerosity cannot replace feasibility or
source weights. There is no new square-packing optimum, major conjecture
closure, unit-distance realization, or general tournament classification claim.
