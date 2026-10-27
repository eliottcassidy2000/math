# Clock arrows, precision and tape guards

Status: **PROVED** elementary statements C1--C3 and T1--T3 below; **FINITE-EXACT**
bounded controls in `normal.json`; **CITED / inherited** mechanisms explicitly
linked below. These bridges give no Collatz convergence or native-K universality.

## C1. Three different meanings of a clock drawing

For a sequence a_n and modulus m, distinguish:

1. Index arrows n mod m -> a_n mod m. They form a function exactly when
   equal index residues always give equal output residues on the chosen domain.
2. Successor arrows a_n mod m -> a_(n+1) mod m. They form a function exactly
   when every repeated value has the same next value.
3. An iteration x -> f(x), reduced modulo m. It descends exactly when
   x=y mod m implies f(x)=f(y) mod m. For a guarded rule the same implication
   must also preserve legality, branch, or any other declared observation.

These are immediate fibre-constancy criteria, but forgetting which one is
intended produces false clock dynamics. Fibonacci modulo 5 is a small hostile:
the pair (a_n,a_(n+1)) has period 20, while the scalar value 1 has successors
0,1,2,4. The pair is an exact autonomous state; its first coordinate is not.

Such diagrams are directed relations or functional graphs. Retain loops,
two-way pairs, missing pairs and multiplicity. A tournament requires an
intrinsic skew, tie-free observable on every distinct pair, which a sequence's
successor relation generally does not supply. Its vertices here are residues,
its pairwise observable is the declared arrow relation, its orientation is
time direction, and the lost data are the lifted state and branch. No arbitrary
tie-breaking is licensed by calling the picture "tournament-like".

## C2. Exact affine composition and division between different clocks

The affine arrows A_(a,b)(r)=ar+b mod m compose by

    A_(a,b) after A_(c,d) = A_(ac,ad+b).

They are permutations exactly when a is a unit modulo m. Otherwise their
functional components have cycles with in-trees. This recovers the October
[functional-graph route](../coloring_clock_session_20261020/FUNCTIONAL_GRAPHS_AND_REGULARITY.md)
with an explicit transition monoid, not an inferred tournament.

For positive integers q,m and integers a,b, the map

    D_(a,b,q): Z/(qm) -> Z/m,    [x] -> [floor((ax+b)/q)]

is well-defined. Replacing x by x+qm changes the floor by am. In general it
does not factor through Z/m: for floor(x/2) modulo 5, the representatives
1 and 6 have the same input residue but outputs 0 and 3. The extra precision
is operational input, not a dispensable label.

This distinguishes a cycle on an odd modulus, where 2 has an inverse, from
integer division selected by parity. Multiplication by 2 modulo an odd m
does not itself determine which integer lift and parity branch is realized.
Statements about the number of integer Collatz or 3n-1 cycles require those
guards and a global orbit argument. A first edge in an ordering supplies neither.

## C3. Exact finite-depth Collatz observations

Let T(x)=x/2 for even x and (3x+1)/2 for odd x, the shortcut map. The first k
parity observations are in bijection with residues modulo 2^k.

**Proof.** Work on 2-adic integers, where 3 is a unit. Each output y has one
even inverse branch 2y and one odd inverse branch (2y-1)/3. Both branches
double precision: congruence modulo 2^j becomes congruence modulo 2^(j+1).
Induction assigns each parity word one residue modulo 2^k, and every input
has a unique parity word. This is also a finite modular induction, without
assuming convergence of any integer orbit. QED.

Thus named depth-k parity queries have exactly 2^k classes. For arbitrary m,
knowing x modulo 2^k m determines these k branches and T^k(x) modulo m:
after one equal-parity step, a difference divisible by 2^k m becomes divisible
by 2^(k-1)m, possibly multiplied by 3. Induct on k. A fixed clock can therefore
answer a bounded horizon, while exact repeated division needs a precision
tower or an incoming digit stream.

The repo already has a genuine arithmetic tournament in
[THM-4754 / three-adic Sierpinski and Collatz](../../../01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md).
Its skew sign is built from the leading nonzero 3-adic digit of a difference.
The Collatz pair-transport statement has a synchronized-parity hypothesis.
It is a separate pairwise construction, not an automatic interpretation of
one-outgoing-edge clock graphs.

## T1. A literal two-stack arithmetic simulation

Consider a binary tape with finite support, head h, and control q. Encode
the cells strictly left of h and the head/right cells by nonnegative integers

    L = sum_(i<h) tape(i) 2^(h-i-1),
    R = sum_(i>=h) tape(i) 2^(i-h).

The read bit is R mod 2. After writing w in {0,1}, the exact updates are

    move right:  (L',R',h') = (2L+w, floor(R/2), h+1),
    move left:   (L',R',h') = (floor(L/2),
                              (L mod 2)+2w+4 floor(R/2), h-1).

The program supplies the new control from (q,R mod 2). **Proof:** expand
the positional digits: a right move pushes the new head-cell digit onto L
and pops R; a left move pops L, then prefixes that digit and w to the old
right tail. The formulas agree cell by cell, including an all-zero tail. QED.

This is one native transition, with time preserved exactly. `check.py`
independently edits finite sets of occupied cells and compares 3,584 cases
against these formulas. This is a classical stack representation proved here
for the exact conventions used; no new universality theorem is claimed.

## T2. An explicit growing-memory lower bound

A right scanner that reports the read bit distinguishes all 2^k words in a
k-cell window by its first k reads. A machine which moves right on 0 and
halts on 1 supplies a sharper hostile: R=0 and R=2^k agree modulo 2^k but
the latter halts after k right moves and the former never halts. The differing
read occurs at zero-based time k. Consequently no fixed residue R mod 2^k,
even with the same head and control, gives an exact all-input autonomous
halting observer for this fixed machine.

Across all these finite tapes, different distances to the first 1 give
infinitely many continuation classes. This requires no appeal to general
undecidability: a continuation of the appropriate length separates each pair.
It is the same *delayed guard access* mechanism as the tournament alignment
examples, with a different operation language and observation algebra.

Finite-support realization matters too. The integers 2^k-1 approach the
2-adic value -1. A scanner that moves on 1 and halts on 0 halts on every one
of those finite tapes but never on the infinite all-ones 2-adic tape. Passing
to the completion creates states outside the integer/finite-support domain.

## T3. A finite quotient can exist without being uniformly computable

Fix one deterministic program and one input. Restrict the state space to its
reachable run, allow only STEP, make halting absorbing, and observe only the
Boolean predicate "has halted". This particular continuation quotient is
always finite:

- if it never halts, every reachable state has the all-zero future trace,
  so the quotient has one state;
- if first halting is at time N, the N+1 distinct remaining waiting times
  N,N-1,...,0 are precisely its classes.

Nevertheless there is no algorithm that, for every program and input,
returns a finite exact quotient with an initial state, effective transitions
and the halt observer. Reachability in that returned finite graph would decide
halting. For completeness, a decider H would let a program D(x) loop if H
says that x halts on input x, and halt otherwise; D on its own description
contradicts either answer. The machine model admits effective program coding
and conditional execution, as in the classical argument initiated by
[Turing's original paper](https://www.cs.virginia.edu/~robins/Turing_Paper_1936.pdf),
§§1--2 and 8. The displayed halt-decider contradiction is given directly here;
it is not a verbatim statement of Turing's circularity formulation.

The scope is essential. A quotient for one fixed run differs from one over
all configurations or all inputs, and halt-only observation differs from
observing every tape symbol. Finite per-instance existence, uniform size
bounds, and an effective construction are three separate assertions.

## The inherited tournament tape and the LRC carry connection

Session65 already encodes a binary word w as an ordered join E(w) of TT3
and C3 blocks. It proves

    d_iso(E(u),E(v)) = Hamming(u,v),
    H(E(w)) = 3^(number of ones),    c3(E(w)) = number of ones.

Its [tape report, §3](../session65_cryptid_tournament_interfaces_20261026/session65-cryptids/REPORT.md)
retains the word order, window origin, head, program and read guards. A write
is a controlled local arc reversal; two-stack arithmetic supplies the literal
integer implementation of the same tape update. Scalars H and c3 erase
order: 01 and 10 already differ in future behavior.

The [native-K parity obstruction](../session65_cryptid_tournament_interfaces_20261026/session65-k-write/REPORT.md)
prevents a literal one-arc write from being realized by any fixed-label,
fixed-order sequence of native complement moves when n=1,2 mod 4. Thus this
faithful tape bridge cannot be promoted to autonomous native metagraph
universality. The parity obstruction does not exclude simulations with
different encodings, relabelings, auxiliary vertices or changed orders.

An underused July mechanism is even closer to the current composition target:
[THM-2211 / carry-regime root transducer](../../../01-canon/theorems/THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index.md)
has a finite 39^9-state upper bound when the next nine carry digits are edge
inputs, but an infinite autonomous continuation index for exact root propagation.
It distinguishes named sheet ownership from its signed scalar sum. The
analogy is exact at the level of missing transition input: carry digits for
arithmetic, read digits for tapes, and joint phase/incidence for attachments.
Its own weaker terminal zero-Haar question remains outside that obstruction;
none of these connections settles LRC(14).
