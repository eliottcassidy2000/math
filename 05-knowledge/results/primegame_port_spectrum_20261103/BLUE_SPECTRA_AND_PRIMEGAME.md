# Blue spectra, quartet clocks, and the difference between a product and a move

**PROVED elementary formulas / FINITE-EXACT n5 census and word cells /
CITED recovered mechanisms.** The spectrum conjecture is not deduced from
these comparisons.

## 1. Correct the operator, then explain the radicals

The current n5 census has three pure-blue and five mixed self-converse
classes, plus two further nodes representing non-self-converse pairs.
Thus3+5 describes eight of the ten merged nodes. This is inherited from
[the owner's-gauge census](../groups_of_five_primegame_merged_metagraph_20261031.md),
and independently reconstructed here from all64 fixed-path presentations.

For the user's first five proposed eigenvalues the polynomial is

    x^5+x^4+(3/2)x^3+2x^2-x.

The two roots(-1+/-sqrt(3))/2 are not algebraic integers, so they cannot be
eigenvalues of any integer adjacency matrix, even as part of a larger
spectrum. The second proposed multiset has polynomial

    x^5-x^4-2x^3+x^2-3x+2.

Its trace is1. Both lists fail the trace-zero condition for a loopless
tournament matrix. The second can be realized by an integer companion
matrix, but that does not specify a metagraph operator. The incoming
[PRIMEGAME path-orbit packet](../primegame_path_orbits_20261102/CONNECTIONS.md)
independently records these same corrections; they are not claimed as
exclusive discoveries here.

There is a useful positive repair. Take the **ordinary symmetric adjacency
matrix of the blue-line subgraph** of the merged metagraph. Each blue line
has multiplicity one, and none is a loop. This declared operator has:

| Induced vertex set | Actual graph | Characteristic polynomial |
| --- | --- | --- |
| Three pure-blue nodes | Three isolated vertices | x^3 |
| Five mixed nodes | P4 plus one isolated vertex | x(x^4-3x^2+1) |
| Eight self-converse nodes | C4 with one leaf at each vertex | (x^2-1)^2(x^4-6x^2+1) |
| All ten nodes | Previous graph plus two isolates | x^2(x^2-1)^2(x^4-6x^2+1) |

The mixed-node spectrum is therefore

    0, +/-phi, +/-1/phi,       phi=(1+sqrt(5))/2.

The eight-node spectrum is

    +/-(1+sqrt(2)), +/-(sqrt(2)-1), +1 twice, -1 twice.

The three/five partition is not an invariant block decomposition of the
full adjacency matrix: blue edges cross it. Taking an induced subgraph
changes the operator. These golden and silver values are genuine, but
they are not the supplied two lists or tournament Hamiltonian counts.

### A uniform leaf-attachment recurrence

For any simple graph G on m vertices, attach r>=1 new leaves at every
vertex. Eliminating the leaf coordinates by a Schur complement gives

    chi_new(x)=x^(m(r-1)) det((x^2-r)I_m-x A_G).             (8)

For x!=0 this follows directly from the block matrix; equality then holds
as a polynomial, including at zero. For r=1, every old adjacency eigenvalue
theta gives the equation x^2-theta*x-1=0. The eigenvalues2,0,0,-2 of C4
give exactly the silver spectrum above. Repeating one-leaf attachment gives

    chi_(k+1)(x)=x^N_k chi_k(x-1/x), N_(k+1)=2N_k.          (9)

These are standard block-determinant mechanisms, with the proof written
out for this application. No claim is made that the actual n6 metagraph
is obtained by leaf attachment from the n5 graph. That would be a separate,
false-to-assume identification of vertex growth with metagraph growth.

## 2. Retain the fourteen instructions and ten valuation coordinates

Use [Conway's original1987 list, Theorem1](https://gwern.net/doc/cs/computable/1987-conway.pdf),
in its first-enabled order, with letters A,...,N:

    17/91,78/85,19/51,23/38,29/33,77/29,95/23,
    77/19,1/17,11/13,13/11,15/2,1/7,55/1.

These are fourteen instructions on the ten supported primes
(2,3,5,7,11,13,17,19,23,29), not fourteen graph vertices identified with
ten prime factors. The input integer may have multiple tokens at a prime;
primes outside the list are inert spectators.

The supplied partition is the already audited
[Session129 partition](../session129_primegame_tournament_partition_20261031/partition/README.md):

    X0=JK, X1=AIMN, X2=CEFH, X3=BDGL.

Its unique primitive valuation grading, under the stated block conditions,
has coefficients(-4,0,-4,2,2,2,2,3,-1,3). The increments are respectively
0,-2,+1,0 on X0,X1,X2,X3. We independently recover its one-dimensional
solution space. For an executed return word, #X2=2#X1; this necessary
balance does not specify legal order, nonnegative registers, or priority.

## 3. Classify every contiguous use of a complete quartet

The inherited exact prefix-cell compiler is
[P4 in FRACTRAN_PRECISION](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md).
For a proposed word w and prefix multiplier A/B, rule i has input
divisibility requirement d_i=B*b_i/gcd(A,b_i). Let D be the lcm of the
chosen requirements and retain all requirements of earlier rules as
forbidden divisors. Writing the input as D*t, remove their gcd with D.
The word executes iff none of the resulting divisors divides t. A divisor1
means the cell is empty; otherwise t=1 is its least input witness.

Applying this exact all-positive-integer test to all24 orders of each
quartet and both orders of JK gives the complete table below. Listed
prime sets mean t is coprime to their product.

| Block | Executable ordering | Least input D | Primes forbidden in multiplier t |
| --- | --- | ---: | --- |
| X0 | JK | 13 | 3,7,13,17,19,23,29 |
| X0 | KJ | 11 | 3,7,13,17,19,23,29 |
| X1 | AIMN | 637 | 2,3,5,7,11,13,17,19,23,29 |
| X2 | CEFH | 1683 | 2,3,5,13,17,23,29 |
| X2 | CHEF | 153 | 2,5,11,13,17,23,29 |
| X3 | None | — | Every one of its24 orders has empty source cell |

This is an all-size classification of **these fixed words**, not merely a
bounded input scan or a classification of all program runs. It includes
multi-control-token inputs. The minimum executions include

    AIMN: 637 -> 119 -> 7 -> 1 -> 55,
    CHEF: 153 -> 57 -> 231 -> 203 -> 539,
    CEFH: 1683 -> 627 -> 551 -> 1463 -> 5929.

The impossibility of X3 has a short invariant proof. B,D,G preserve the
total number of control tokens, and L does not change it. B requires a
control token. But L can execute only when there are no control tokens:
each of the six control primes enables an earlier rule, with a bare
control denominator available if the conditional rule fails. Hence no
ordering containing all four instructions B,D,G,L can execute contiguously.

The formal products, recovered from Session129, are

    product X0=1, product X1=55/637,
    product X2=539/153, product X3=585/34.

In particular34*(585/34)=585 is integral, but the X3 quartet has **no
executable ordering from any input**, a stronger boundary than failure
only at34. Unordered cancellation loses the source guard.

## 4. The map that actually transfers

The shared mechanism is operation-sensitive information loss:

* A graph Schur complement preserves a declared boundary response, with
  its leaf ports and rational x-dependence; it does not preserve every
  other graph operation.
* A PRIMEGAME block product preserves the valuation displacement; its
  exact source cell and instruction order restore executability.
* The tournament four-query growth identity preserves an exact edge
  derivative; its four actual source extensions and signs must be retained.
  Their signed difference is not itself a tournament construction.

The ten-prime/ten-node coding remains a declared arithmetic encoding as
in Session129. No natural PRIMEGAME execution/tournament-growth conjugacy
has been produced, and no number-field or spectral coincidence fills the
prime-lane realization gap.
