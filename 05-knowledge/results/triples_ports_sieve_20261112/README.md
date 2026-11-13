# Triples, cycle defects, and exact return interfaces

**PROVED IN NOTES / FINITE-EXACT, with independent internal reviews. Chan–Pak and RH remain OPEN.** This packet recovers older Sierpinski, five-core, square-sum, arithmetic and lattice work and follows the connections through explicit maps, tests, and proofs. It classifies several infinite *families*, not every aspect of the tournament metagraph. The new results and the inherited mechanisms are distinguished below.

## The strongest constructive result

[Square-sum endpoint relocation](squares/PROOF.md#7-a-uniform-cycle-opening-endpoint-relocation) gives an exact bijection: deleting an old selected edge and adjoining a two-edge ear can turn a path plus one cycle into a Hamiltonian path. The observed Q25→Q26 birth does exactly this. An explicit two-parameter family places the same operation in infinitely many actual square-sum graphs. Existence of a suitable spanning source support is an additional condition, not inferred from the arithmetic host.

The complete endpoint-conditioned reduction of the first 25 square-sum graphs leaves only two theta kernels, with arm lengths `(1,4,4)` and `(2,3,6)`. Their degree-completion spaces are each K2. The all-size multi-arm theorem classifies those spaces as products of Johnson graphs. The same literal kernel and residual demands give one Hamiltonian completion at Q25 and two at Q26: the outside forced-forest partition changes. This supplies the missing state that explains the birth, rather than merely exhibiting different counts.

The [directed cycle-opening operator](squares/PROOF.md#9-an-actual-directed-cycle-opening-operator) is an actual all-order tournament transfer. Remove a selected cycle arc `a→v` and insert `a→z→b`, joining the opened cycle to the old path from b to t. Its exact inverse gives a bijection with Hamiltonian paths starting at v, ending at t, and using `a→z→b`. In the smallest host the conditioned count stays1 across all8 new-vertex orientations, while the full H ranges over3,5,9,17,23. This identifies precisely what a target-preserving return still needs to control.

## Exact structural bridges

1. **Theta → Cartan → triples.** [The rank-one network identity](triples/PROOF.md) gives
   `det C=pq+pr+qr-pqr` and the complete inertia from `1/p+1/q+1/r-1`.
   The two square kernels have A7 and affine-E8 Cartan shadows. The same family contains E8 at `(2,3,5)` and E10 at `(2,3,7)`. Acyclic star orientations give the cited finite-type cluster criterion. Outside Hamiltonian connectivity is not preserved by this map.
2. **Sierpinski → exact boundary automata.** [Connected dihedral transport](dihedral/PROOF.md) has exactly `1+m^3+4m^2` three-boundary states for odd m>=3 and `1+m^3/2+4m^2` for even m>=4. This recovers 64 at m3 and extends to 97 at m4, 226 at m5 and 540 at m7. Exact triangular gluing closes at every depth. Even moduli expose a circle-lift failure; CRT requires one shared reflection sign.
3. **Gasket → scale triples.** The q-corner electrical refinement has exact conductance factor `q/(q+2)`. Its geometric/cell-count/mass-resistance scales are `(2,q,q+2)`, giving `(2,3,5)` and `(2,5,7)`. These have a specified network meaning; they are not interchangeable with triangle-group or prime-sieve triples.
4. **Tournament → cluster-return component.** [The native classification](ports/PROOF.md) and [four-step return theorem](ports/FOUR_RETURNS.md) retain weighted matrix mutation. Every return of length at most four preserves H, at every order. Thus a count-changing return would need at least five mutations. The inherited SCC-hull component classification yields the new ordinary-class quotient `c_n=A000568(n)-A000568(n-1)+1`. It classifies a quotient of known counts, rather than extending the sequence by itself.
5. **Twelve paired charts → measured code.** Eleven observations modulo13 distinguish all4096 binary port choices with Hamming distance5. Exactly225 five-time subsets identify every choice; no four-time subset does. This is a nonlinear `(11,4096,5)_13` observation code, not a single tournament addition operation. The natural paired binary code has distance2, so it is not the Golay code or the Leech lattice.

## Arithmetic: a smaller exact clock and a joint branch test

[The first-25-prime sieve](arithmetic/WHEEL.md) for `A_k=99*2^k+23` has **minimal period 61,209,720**, seven times smaller than the inherited valid period. Nine irredundant exclusions suffice. Exactly32,353,278 phases survive, with density5392213/10201620. All2300 three-prime subsets collapse to213 exact exponent predicates;364 triples impose no exclusion. In particular `{2,3,5}` permits three quarters of exponents while `{2,3,7}` permits every exponent.

The stronger [joint formal-branch sieve](arithmetic/JOINT_BRANCH_SIEVE.md) for `(A_k,3A_k+2,3A_k+4)` forces `k=6+12t`, then six independent conditional CRT clocks. Its minimal period is1,642,460,820 and its exact density is5346/136735. At k6 the scalar triple is the prime triple `(6359,19079,19081)`. On the explicit infinite progression `k=6+239978082932247660*u`, u>=1, all three entries pass the25 tests and are all composite. Positive finite-wheel density therefore does not establish prime occurrence.

[Nicolas's RH criterion](arithmetic/NICOLAS.md) is CITED and verified against the primary source: the user's strict inequality, threshold `n>=p_120569#`, and coefficient3.6444150964... are correct. Its unbounded Euler-product estimate is a different predicate from the exponent wheel. The entire A_k family is outside A002496: for k>=2 it is3 mod4, and k1 is checked separately. No new A002496 term range is claimed.

For A135244, the inherited exact equivalence `M(n)=(n-1)^2 iff n-1 is prime` transports the prime seed and the all-composite progression into exact equality and strict-inequality families. The seed yields `M(6360)=40436881`, `M(19080)=364008241`, `M(19082)=364084561`; these lie inside the prior certified prefix, so they are applications rather than new sequence terms.

## What this changes for Chan–Pak

The earlier [absolute-return packet](../berggren_target_return_20261111/SYNTHESIS.md) supplies actual target-preserving returns at221 and665 and finite defect-growing branches from1997. Those constructions remain the anchor. This packet adds a precise reason to carry a cycle-defect sector through a birth, a complete infinite family of local completion spaces, a classifier for several compatibility interfaces, and arithmetic diagnostics for formal ternary branches. None supplies the missing all-target, source-valid return selector.

There is a positive exact root-return family inside the cluster theorem: `R→C3` maps to `converse(C3)→R`, preserving `H=3H(R)` and all vertex names while creating three indegree-one roots. It toggles back repeatably, adds no vertices, and changes old arcs; the intermediate quivers need not be tournaments. This requires an extreme ordinal C3 factor. Both desired next targets3H+2 and3H+4 lose divisibility by3, so this carrier cannot by itself return after either branch. That is its exact closure failure.

At checkpoint, incoming Session292 supplied a stronger hidden-factor theorem. [Its joint extension](arithmetic/JOINT_REFINEMENT.md) proves that every accepted joint phase has subprogressions with arbitrarily many distinct prime factors in each of A,B,C; the three factor lists are automatically disjoint. The entire triple also misses A002496, with C excluded by8 mod16 at k>=4. These are complete statements about this family, with no new global prime or sequence record.

The synthesis is operational: retain the variable shared by the next move. For the wheel it is the exponent phase; for CRT transport it is the sign; for square kernels it is the outside component partition; for ports it is the continuation profile; for quiver mutation it is the full weighted matrix. These are related ways of repairing a lossy summary, not a proved identification of all five objects.

- [Inheritance, corrections, and the live concept board](INHERITANCE.md)
- [Independent review of arithmetic](AUDIT-arithmetic.md)
- [Independent review of the root proofs](AUDIT-root.md)
- [Independent square review](ports/AUDIT_SQUARES.md) and [mutation review](squares/AUDIT-cluster.md)
- [Independent four-step review](ports/AUDIT_FOUR.md)
- Reproduce with `python -B 05-knowledge/results/triples_ports_sieve_20261112/reproduce.py` (and `--optimized`).

The packet's receipts declare finite universes, positive and hostile controls, source hashes, and exact stopping boundaries. No new theorem namespace is reserved and no unreviewed canon promotion is made.

For cross-platform reproduction, SHA-256 values in the packet manifest, source pins, and replay summaries normalize CRLF text to LF; compressed binary inputs are hashed without alteration. Mathematical receipts retain their declared source-hash conventions.
