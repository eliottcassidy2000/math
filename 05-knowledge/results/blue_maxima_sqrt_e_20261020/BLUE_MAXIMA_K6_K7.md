# Exact blue maxima at orders thirteen and fifteen

Status: **FINITE-EXACT; TWO EXHAUSTIVE ALGORITHMS AGREE; INDEPENDENT AUDIT APPROVED**. The general closed form remains open. This packet does not modify the inherited session14 packets.

## Results and inherited meaning

Let `D=(A,b)` be a hollow reflection frame on `k` pair indices, with lexicographically ordered index pairs. Its two masks specify the tournament `A` and the symmetric cross graph `b`. Let `h(D)` count directed signed index permutations ending at an `x` vertex. The inherited regular-subset reduction identifies the maximum blue mass at order `2k+1` with

    M_k = max_D h(D).

The new exact values are

    M_6 = 1185,       M_7 = 8295.

Together with the inherited census this gives

    k:     1  2  3   4    5     6     7
    M_k:   1  3  9  39  195  1185  8295.

There is one maximizing hollow-frame isomorphism type in each of the two new orders. The data are:

| k | order | maximizing `(A,b)` | labelled maximizing frames | frame group | maximum-blue tournament classes |
|---:|---:|---:|---:|---:|---:|
| 6 | 13 | `(314,3618)` | 144 | C5 | 12 |
| 7 | 15 | `(268,39375)` | 1008 | C5 | 24 |

The maximum classes split evenly between the triangle-parity shores: 6+6 at order 13 and 12+12 at order 15.

## Why the reduced census is complete

A frame has `m=binom(k,2)` tournament bits and `m` independent cross bits. Relabel the pair indices to put `A` into one representative of its `S_k` orbit. For that representative, still enumerate **every** one of the `2^m` cross graphs. Every original labelled frame is represented, since relabelling maps the full cross-graph cube bijectively onto itself.

The generator scans every tournament mask, marks the complete permutation orbit of each first-unseen mask, and records its orbit size and automorphism order. An independent Python implementation verifies that the resulting orbits are disjoint, have the recorded sizes, and cover every tournament mask. No external tournament catalogue is assumed.

There are 56 tournament types on 6 indices and 456 on 7. Accordingly the exact searches need

    k=6: 56 * 32768 = 1,835,008 reduced frames;
    k=7: 456 * 2097152 = 956,301,312 reduced frames.

The latter represent all `2^42=4,398,046,511,104` labelled external frames. Weighting each cross graph by the orbit size of its tournament gives the complete labelled histogram and moments.

## Two exact transforms, plus dynamic-programming checks

Fix `A`. A candidate signed prefix consists of an index permutation `p` and polarities, with the final polarity `x`. Reject it if a same-polarity edge contradicts `A`. Otherwise its mixed-polarity edges require a set `S` of distinct cross bits, with required-one subset `R`.

Its indicator is

    product_(e in R) b_e * product_(e in S\R) (1-b_e).

### Walsh computation

Write `chi_T(b)=(-1)^(sum_(e in T)b_e)`. The indicator is

    2^(-|S|) sum_(T subset S) (-1)^|T intersect R| chi_T(b).

All denominators divide `2^(k-1)`. Accumulate the scaled integer Fourier coefficients and apply the Walsh-Hadamard transform. The resulting value at each `b`, divided by `2^(k-1)`, is exactly `h(A,b)`.

### Independent ordinary-monomial computation

Expand the same indicator in the ordinary monomial basis:

    sum_(U subset S\R) (-1)^|U| product_(e in R union U) b_e.

A subset-zeta transform evaluates this polynomial at every binary vector `b`. `replay_zeta.cpp` uses polarities indexed by the vertices, whereas the Walsh producer indexes them by positions in the permutation. Its coefficient arithmetic and transform are different from the Fourier calculation.

The two programs stream one tournament's array at a time through a FIFO and compare **every entry**. No multi-gigabyte result array needs to be retained. At `k=7` all 956,301,312 comparisons agree; the streamed replay took about 19.4seconds in the recorded environment. The main Walsh census took about 11.8seconds. Peak working memory is comfortably below 1GB.

At `k=6`, a third method, the forward pair-subset DP, computed every reduced frame. A direct subcube indicator summation also agrees at every frame. The Walsh program was then checked pointwise against that DP array before the order 15 computation.

Additionally, the order 15 Walsh run independently checks 256 deterministically sampled cross graphs for every tournament type and **all** cross graphs attaining that type's maximum, using the subset DP. These 117,355 checks are supplementary; they are not substituted for the complete zeta replay.

The full labelled moments agree with the inherited all-order formulas:

| k | mean h | variance h | sum h | sum h squared |
|---:|---:|---:|---:|---:|
| 6 | 720 | 89685/2 | 773094113280 | 604777029304320 |
| 7 | 5040 | 1875825 | 22166154415964160 | 119967383853151027200 |

The order 15 square sum is accumulated in 128-bit integer arithmetic. Every value is checked positive and odd.

## Extremizer structures and all maximum classes

For `k=6`, the representative tournament has a regular five-vertex tournament on indices `{0,1,3,4,5}`, all pointing to the sink2. The cross graph is exactly the star joining2 to those five indices.

For `k=7`, indices `{0,1,2,3,4}` induce a regular five-vertex tournament. Vertex 5 dominates those five, and 6 dominates all six other vertices. The cross graph is exactly the clique on `{0,1,2,3,4}`.

The census has just one maximizing cross graph for the one maximizing tournament representative in each order. Therefore it is fixed by the entire automorphism group of that tournament. Direct automorphism enumeration shows that group is `C5`, acting as a five-cycle and fixing the remaining `k-5` indices. The complete frame orbit consequently has `k!/5` elements, exactly the recorded number of labelled maxima.

For binary matching completions, the number fixed by `C5` is `2^(k-4)`. Every other vector has trivial stabilizer and a five-element orbit, since the group has prime order. The maximum-blue classes thus number

    (2^k - 2^(k-4))/5.

There are 4 additional classes of blue mass237 in the order 13 maximizing hollow type, and 8 additional classes of blue mass1659 in the order 15 type. No other hollow type can attain the global maximum, since `B=h/|H_a|<=h<=M_k`.

Toggling a matching coordinate fixed by the group preserves triviality of its stabilizer and reverses triangle parity, so maximum classes pair across shores. This concerns the inherited auxiliary matching-cube operation; no new adjacency claim about the actual blue metagraph is made.

`witness_n13.json` and `witness_n15.json` contain the full adjacency row masks, explicit automorphisms, a trivial-stabilizer matching subset, and completion counts. Independent direct enumeration checks all 1185 and 8295 symmetric Hamiltonian words, respectively. SHA-256 digests of those deterministic word streams are included; the full word lists are not needed to reproduce the count.

## General methods and boundaries

The transform identities above are valid at every `k`; they give a general exact finite algorithm, not a closed form for `M_k`.

All seven known exact values happen to satisfy

    M_k/k! = sum_(j=0,...,floor(k/2)) 1/(2^j j!).

This is **only a finite pattern**. No proof of this formula, a universal `sqrt(e)` upper bound, or an asymptotic equivalence is claimed. The separate `SQRT_E_LOWER_BOUND.md` proves a one-sided asymptotic lower bound by a different construction.

A deterministic, gradient-checked local search supplies independently verified witnesses

    M_8 >= 66123,        M_9 >= 597249,        M_10 >= 5975877.

These are **certified lower bounds only**, not maxima. The tempting finite-pattern values would be 66465, 598185, and 5982795. A 100,000-start order 17 search failing to reach 66465 does not refute or prove that formula. Exact extension searches of the previous maximizing frame are also only restricted-family maxima.

Every exact maximizing frame through `k=7`, as well as the search witnesses, has

    d_A(i)+d_b(i) in {k-2,k-1}.

Thus matching `a_i=k-1-d_A(i)-d_b(i)` makes the full odd tournament regular. This is finite structural evidence, not a necessary-condition theorem for all maximizers. That regular completion need not have trivial matching stabilizer and should not be confused with the completion attaining maximum blue mass.

## Elementary all-order block bound

For any positive block sizes `s_1+...+s_r=k`, discarding the transitions between consecutive index-word blocks gives

    M_k/k! <= (M_(s_r)/s_r!) product_(i<r) (2 M_(s_i)/s_i!).

Indeed an ordered partition into unordered blocks contributes `k!/product s_i!` choices. A nonterminal block has at most `2 M_s` signed directed words: those ending at `x` are bounded by `M_s`, and swapping every `x/y` gives the same bound for those ending at `y`. The final block must end at `x` and has at most `M_s` choices.

Inserting the exact value `M_7=8295` yields

    M_k/k! = O((79/24)^(k/7)).

This improves the previous elementary exponential bound but does not give a constant asymptotic upper bound. The sharper asymptotics remain open in this packet.
