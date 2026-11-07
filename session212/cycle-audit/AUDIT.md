# Session212 independent cycle-spectrum audit

**Verdict: PASS.** The all-order constructions, their stated structural
refinements, the fresh finite spectra through order nine, and the exact
inverse-size tails have passed independent mathematical review and executable
checks. Order ten remains explicitly inherited finite-exact evidence, not a
new full census. No literature-priority conclusion is drawn.

## 1. Scope and conventions

A Hamiltonian cycle is a directed spanning cycle modulo cyclic rotation,
without an automorphism quotient or a reversal quotient. Fixing a named first
vertex counts each cycle once. Tournaments are nonempty; the cycle count is
zero at orders one and two. Strong-cycle equivalence is used only at order
at least three, where it is the classical Camion theorem.

Audited packets are the sibling `cycle-bridge`, `cycle-spectrum`, and
`cycle-levels` directories. Their exact approved manifests are recorded in
`APPROVED_PACKETS.json`. The mathematical review also covers the session
index. The existing dirty checkout was read, not modified or staged.

## 2. All-order proof audit

### Binary comb and sole-exit construction

For a transitive chain and two apices, every cycle normalized at the
sole-exit apex first takes its designated root. The vertices before the
coding apex form an increasing subset. The sentinel vertex zero must lie
after that apex, so the complement is nonempty and has a forced first
vertex. Conditioning on the final pre-apex vertex gives exactly one unit
term and selected powers of two. This proves the formulas, not merely their
finite initial cases:

    HC(B-) = q + tau,       HC(B+) = q + tau + 1  (m >= 1).

The additional forward cycle is the full increasing chain after the two
apices. Every integer is therefore realized. At exact order N>=5 the
unrestricted family covers 0 through 2^(N-4)+1. For k>=2 it supplies strong
witnesses at every N>=max(5,4+ceil(log2(k-1))). This is an upper bound, not an
exact minimum-order formula.

The reversed family has no degree-defined gate: checking all possible
outdegree-one tails rules out a target of indegree one. Every positive k
has such a strong sole-exit witness, with the slightly weaker threshold
max(5,4+ceil(log2(k))). Removing the two apices leaves a transitive chain.
A tournament with feedback vertex number at most one has at most one
Hamiltonian cycle, since deleting its feedback vertex leaves a unique
Hamiltonian path. Thus every constructed k>=2 witness has feedback vertex
number exactly two. These properties hold simultaneously.

The marked insertion operation preserves all old arcs and changes
q to 2q+epsilon. Hence HC changes to 2HC-tau+epsilon. At tau=0 this is the
literal binary grammar HC -> 2HC or 2HC+1. No unproved interpolation or
closure of arbitrary tournaments is used.

### Endpoint supply and forced-entry construction

The almost-transitive tournament U_m, m>=3, has endpoint vector

    (1,1,2,4,...,2^(m-3),1).

A Hamiltonian path uses either no reversed arc or its unique reversed arc,
which partitions the remaining vertices into two increasing runs. This
proves the vector directly. Its subset sums fill the full interval from
zero through 2^(m-2)+1. Adding a forced-entry vertex and its predecessor
converts a selected endpoint sum into the cycle count, with no rotational
factor. This independently gives the same fixed-order interval and bound.
For positive outputs the return arc is a strong bridge. The zero output is
nonstrong, so that bridge terminology is not applied there.

### Strong endpoint cores and two prescribed endpoints

For a minimally wide nonzero binary code the top bit is one, so the coding
core contains the cycle apex,0,1,...,last,apex. The h=1, zero-bit core is a
triangle. Consequently every positive integer is already a specified-start
path count of a strong tournament. In its sole-exit closure, cutting the
forced arc also gives exactly that many paths with both endpoints specified.
None of this identifies a marked endpoint count with the total path count H.

## 3. Boundary and overclaim checks

The final packets correctly retain these restrictions:

1. The forward comb increment fails at m=0: its proposed direct cycle would
   require the nonexistent arc x1 -> a. The reversed zero-bit case is valid.
2. Zero is realized by nonstrong tournaments. The positive spectrum on strong
   tournaments excludes the exceptional singleton convention.
3. A sole exit is not a degree gate. The forward top-value member does have
   a gate when tau=1 and all digits are one; the entire enlarged interval is
   not labeled gate-free.
4. Neither logarithmic order bound is asserted to be minimal.
5. Transitive padding of an arbitrary vertex need not preserve cycles.
   The proved padding uses a vertex of indegree one or, dually, outdegree one.
6. A first-order cycle witness need not itself be paddable. For example
   14,17,24 belong to C_7 but not E_6. The inverse-size proof uses separate
   endpoint witnesses and explicitly covers preceding sizes.
7. The total Hamiltonian-path omissions 7 and 21 remain valid. Complete
   cycle support and complete marked endpoint support do not settle the
   unrestricted path-spectrum conjecture.
8. Fixed-order gaps are not global gaps. No all-order nesting, next-order
   filling, or minimum-order upper-ray theorem is inferred from a census.

## 4. Fresh independent computations

The audit implementation does not import either researcher's builder or
counter. `audit_binary_cycles.py` constructs adjacency matrices by ordered
pairs and counts with a top-down remaining-set recursion, separate from the
research packets' forward subset DPs. It checks:

- 511 binary constructions, all digit words of widths zero through eight;
- 63 complete-permutation cycle controls through order nine;
- all anchors on 12 selected targets, 101 anchor checks;
- 14 zero-digit padded cases and seven transitive zero controls;
- five rejected invalid parameter choices.

`audit_packet.py` independently checks:

- all 2,040 comb parameter cases for widths one through eight and both apex
  orientations, including strongness and exact degree-gate boundaries;
- 120 full-permutation comb controls, plus all four zero-bit boundary cases;
- 504 selected-endpoint constructions, every subset at core orders 3..8;
- all 1,166 saved cycle and endpoint witnesses from the level census;
- 847 endpoint closures and their physical three-vertex-chain replacements;
- C_n subset E_n at every order 3..9, the contiguous E_9 prefix 0..593,
  and the missing endpoint value 594;
- agreement of both separate census implementations through order nine.

`audit_refinements.py` checks 504 marked growth operations with unchanged
old comparisons, 256 strong minimal-width cores, and 256 counts with both
endpoints prescribed.

`brute_census.cpp` uses literal permutations, with no subset DP, on every
class of orders 3..8. It reproduces the entire class histogram, not only
its support, on all 7,410 classes. The comparison is recorded in
`brute_histograms.json`.

All independent Python audit programs have explicit checks that survive
optimization. Normal and `python3 -O` outputs agree byte for byte.

## 5. Coverage and exact inverse-size audit

The class generator's color refinement is isomorphism-invariant. It examines
every permutation within the resulting ordered cells, so its canonical key
is exact. Starting from the singleton and taking every vertex-extension cut
therefore proves coverage inductively; agreement with known class counts
is a control, not the sole reason for completeness.

The fixed-cycle enumeration through order eight also proves complete
positive support: every strong tournament can be labeled along a Hamiltonian
cycle. At order nine, all 6,880*256=1,761,280 extensions of the complete
order-eight bank are examined. Duplicate presentations do not affect
support coverage, but their frequencies are not class multiplicities.
The integer types suffice: at these orders total Hamiltonian paths are at
most 9!, well within signed int. The final level counter uses safe per-row
array initialization.

The finite inclusions C_n subset E_n, 3<=n<=9, and [0,593] subset E_9,
combined with the proved endpoint closure and degree-one padding, establish
exact all-order upper rays for every cycle count 1..593. The finite smaller
spectra establish minimality. The 25 additional listed endpoint values above
593 also have exact minimum order ten and full tails. This conclusion is
stronger than eventual existence but remains explicitly restricted to those
certified values.

## 6. Independent replay and inherited evidence

The audit separately reran:

- generation of every order-eight class and the full order-nine extension
  census, matching every generated row, histogram and witness file;
- the level verifier normally and with optimization, the 4,052 endpoint-cut
  controls, and its entire labeled-tournament census through order seven;
- both unified-bridge verifiers normally and with optimization, checking
  the regenerated result-file hashes rather than confusing summary stdout
  with the full result artifact;
- the spectrum verifier normally and with optimization, and re-extraction
  of all inherited order-eight, nine and ten data from pinned git objects.

The inherited extraction matches committed source SHA-256 values, total
class counts, strong class counts, and its original receipts. The complete
order-ten spectrum [0,961] union {971} is still inherited evidence. Its
(H,HC)=(15745,971) witness is independently recounted by literal permutations.
No independent full order-ten enumeration was performed or claimed.

All pinned source copies and citation hashes in the two main construction
packets were verified against their exact git revisions. The final session
index's six relative destinations and hashes were independently checked.
A pre-existing navigation limit failure was repaired by deleting one empty
line only; all nonwhitespace content and the initial warning were preserved.
That formatting repair has no mathematical consequence.

## 7. Remaining limits

External priority has not been established. General minimum order, global
nesting of exact-order spectra, the first-order upper-ray conjecture beyond
the certified values, stronger-connectivity spectra, and unrestricted
Hamiltonian-path spectrum completeness remain outside these results.
