# Fixed twists, the reverse potential, Antihydra, and least common ancestors: second half of the session

klein, 2026-10-25 (worktree `codex/session-tournament-twist-unify`). Reflection: provenance of ideas,
stopping reasons, and the shape of what was learned. Truth lives in THM-4772/4773 and the packet
`05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/`.

## 1. How the directive mapped onto the repository

The owner's second prompt was a compressed itinerary. Reading it against the newest packets:
"fixed-twist rational lifts and conductor ideals over non-UFD bases" are verbatim the two "next
precise targets" at the end of `normseven_boundary_20261025/README.md`; "the one-step compatibility
test" is CONTINUATION D2; "the reverse valuation potential" is `port_residual_integral_20261025/
counting/REVERSE_POTENTIAL.md`, proved there with implementation declared open; "the six order-12/14
repair failures" are the LOW_ORDERS cells; "Antihydra's odd orbit matches the twisted Syracuse map"
is the lane opus-S56 reserved THM-4770/4771 for the same day. The Petersen/K5/K3,3 and LCA requests
pointed at THM-4529, THM-261, THM-071, procgen_kuratowski and Session56. The method lesson repeats
the first half's: search the statement before deriving, and read the phrases as the other sessions'
vocabulary.

## 2. What was learned, lane by lane

**Conductors.** The UFD proof in BOUNDARY B1 uses unique factorization only to say that divisibility
by a product of distinct primes is detected prime by prime; what changes over a general normal base
is not that proof but the object: a twist with a square ideal factor, or two twists sharing a prime
with odd exponents, makes the root algebra non-normal, and the normalization is the direct sum of
*fractional ideals* `c_z^{-1}` indexed by characters. The class of `c_z` is the invariant the
owner's phrase "conductor ideal classes" was asking for; over `Z[sqrt(-21)][1/2]` it is nonzero for a
single twist and the total module can still be free. The compatibility test survives, with one extra
test per reduced root; the mixed ones are a genuinely new bank. The prime 2 is the one place where
cross-character integral elements appear, found by the computer before the proof said so.

**Fixed twists.** The pointwise criterion from THM-4765's proof turns the residual torsor of BOUNDARY
§5 into explicit squareclass equalities, and on the central line those equalities are three port
conditions plus one conic. The unit circle appearing for the twist `(1,1,-1,1,-1,1,-1)` ties the
hourglass to the square-sum theme of Session57 through Pythagorean triples. The dichotomy (conic on
the central line, genus five outside) is the kind of exact boundary the repository prefers to a
"rational points exist" heuristic.

**Reverse potential.** Implementing the construction on the production graph took one page because
the packet's formulas are exact; the surprise was how small the fixed-precision live subgraph is
(parity of `a(N)/2^{floor(N/2)}`: 550 states at `N = 100`). This is a smaller unrestricted recurrence
for normalized digits, in exactly the sense the counting packets use the word; it is not a smaller
enumeration of tournaments.

**Antihydra.** Everything structural follows from canon plus one induction: the map is the negation
mirror of THM-2228's Mahler carry map, the odd first-return map is the `3n-1` Syracuse map with the
discarded twos replaced by threes, valuation sequences form a full shift whose cylinders are 2-adic
balls, and halting is a ballot barrier of rational slope `3/2`. The golden ratio in the model
probability `phi^{-6}` comes from the step law `2v - 3`, not from any Fibonacci structure. The
contrast with Collatz is clean: rational versus irrational slope, prefix barrier versus asymptotic
condition, growth versus contraction, same valuation law.

**LCA.** The honest merges are small: valuation partial sums are trie depths; the Gale tournament of a
linked `K6` is the clone hostile with depth triple `(1,1,3)`; ten pairs of a 5-set need five
hierarchies because the crossing graph is `L(K5)`. The large claims ("Petersen explains X") were
declined for lack of a map.

## 3. Unifying reading, refined

The first half proposed "expand at the critical point". This half adds the dual move: **keep the
character index.** Every exact result here is a statement indexed by characters (ports): ideals
`c_z`, conductor components `d_z`, Walsh coefficients, pair conditions on twists. Scalar summaries
(a rank, a count, a mean over twists) lose exactly the per-character data, and the obstructions
found (non-free ports, mixed reduced roots, the central conic) live in that data. The counting lane
is the same principle in a different coordinate: the reverse potential is a per-state (not global)
valuation, and the naive global subtraction fails on the real graph.

## 4. Stopping reasons and next obligations

- The six order-12/14 cells were not attacked: the detour framework (marked paths, K moves,
  anti-reversal equations) would have to be re-implemented to search replacement templates; the
  cells remain template failures per LOW_ORDERS.
- Antihydra: opus-S56 holds the reserved IDs; this note was sent to opus and takes none.
- Open: HYP-9374 collisions (untouched); the prime-2 version of THM-4772 (wild case); interior
  rational points of twisted hourglasses; an asymptotic bound on the fixed-precision live graph;
  the signed clone-excess version of the reverse potential.
