# Independent audit of Session 47 wiring theory and suspension

Status: **PASS at the stated conditional and finite scopes.** No producer files
were edited, and no Git operation was performed in this audit lane.

The eight original inputs are frozen in `INPUT_PINS.json`. The later boundary
series note and its checker are separately frozen in
`GENERATING_INPUT_PINS.json`. Both pin files include source provenance,
byte lengths, and SHA-256 values; local snapshots make this audit portable.
Copied Markdown inputs use .txt archive names so their original relative links
are preserved as source text rather than misrepresented as live relocated links.

## 1. Findings suitable for the session synthesis

1. The minimax/cut/MST theorem is correct for the declared finite undirected
   weighted state graph. It is classical weighted connectivity, not an
   independent proof that a native state graph or helper supply exists.
2. The neutral-interval quotient is exact at each budget under the complete
   pairwise-record/additive-guard and simultaneous-arrival packet contract.
   Retaining every heavy neutral coordinate is a sufficient uniform rule,
   not a necessary rule for every special instance.
3. The interval boundary gap is an exact minimum crossing cost but is not a
   metric. Internal routing is what permits successive quotient edges to use
   different endpoint witnesses. Treating it as free below its budget is wrong.
4. Helper packing proves survival of arrivals after forced-U exclusions. It
   does not by itself prove an adjacent route, its boundary cost, original
   current guards, or a native supply.
5. The explicit order-five pair has isomorphic old targets and nonisomorphic,
   non-converse order-seven targets, despite intrinsically unique outer ports.
   Thus forgetting the source/path incidence prevents a forward map on bare
   target classes. This does not invalidate the reverse inner readout.
6. Both triangle identities and the two-port determinant identity are correct
   at all orders covered by their hypotheses. The fourth-order expression
   counts directed simple four-cycles, including their tournament chords;
   the higher trace coefficients count rooted closed walks, not simple cycles.

No mathematical correction to the reviewed producer statements is required.
The scope qualifications above are essential, not optional disclaimers.

## 2. Bottleneck proof audit

### Minimax, cuts, and minimum spanning trees

For any threshold b, a path whose edge maximum is at most b is exactly a path
in the threshold graph. Concatenation proves the strong triangle inequality.
The finite state set ensures an attained simple path whenever the minimum is
finite. Nonnegative weights allow zero-distance distinct states, so
“pseudoultrametric” is the correct general term.

Every path crosses each separating cut, giving the lower bound by its least
crossing edge. If the pair bottleneck is t>0, the component formed using
strictly smaller weights is a separating cut with no edge below t. The t=0
case follows from nonnegativity, and the unreachable case uses the full
finite-edge component. This handles zero and infinite values without an
unjustified perturbation argument.

Removing a heaviest edge on an MST path produces a cut. Any lighter crossing
edge would decrease the total tree weight. Thus the tree path maximum equals
the minimax value for every pair, including tied weights. A singleton has
budget zero; a disconnected non-singleton graph has infinite universal budget.
A minimum spanning forest recovers the finite within-component bottlenecks.

The maximal dominated pseudoultrametric characterization follows by applying
the strong triangle inequality along each legal path. Supply monotonicity
requires an actual cost-preserving state/event embedding; adding new endpoint
obligations need not lower the worst case. The fixed-partition two-sided bound
in equation (8) follows by projection and by routing internally with bound r.
The audit checks that bound directly on 2,324 graph partitions.

### Exact neutral fibres and exact boundary costs

The pinned inherited contract explicitly permits replacement of any set of
records between two feasible arrivals, with the actual additive changed-record
cost. This is necessary for realizing a quotient crossing as one event. The
proof does not silently assert a microscopic route through intermediate
owner changes, and does not apply unchanged to an arbitrary graph that bans
some feasible-to-feasible packet replacements.

For a fixed visible choice, upward closure of required choices and downward
closure of forbidden choices give exact bounds L,R. The charge equation remains
complete because every omitted component is individually neutral. The induced
neutral upset description is exact: implications into a selected visible
successor are already forced by L, while implications into a forbidden visible
successor exclude their predecessors through R.

Removing a minimal element of W\W' preserves upward closure. Any predecessor
remaining in W' would force that element into W', a contradiction. The dual
argument justifies adding a maximal missing element. The resulting route
changes each differing component once, stays within the bounds, preserves
charge, and uses only events of weight at most b.

For two intervals, every member pair must change all coordinates in
l_alpha\r_beta and l_beta\r_alpha. The proposed witnesses

    l_alpha union (l_beta intersect r_alpha),
    l_beta  union (l_alpha intersect r_beta)

are upsets in the respective bounds. Their set differences are exactly the
two forced sets. Consequently their weighted distance is the lower bound,
including when weights vanish. Adding the fixed visible and helper costs
therefore gives the exact minimum crossing cost, not merely a relaxation.

An internal route reaches the witness needed for each crossing. This is the
step that establishes sufficiency and the stated event-count upper bound.
Heavy neutral visibility is sufficient to make this construction work at all
budgets, but a forced singleton fibre may safely omit a heavy coordinate. The
producer now states this distinction explicitly.

The off/free/on example correctly violates the ordinary triangle inequality
for the interval gap. Disjoint intervals may have a legal simultaneous bridge,
and zero weights permit zero gap even for disjoint intervals. Neither fact
contradicts the pseudoultrametric obtained after legal minimax closure.

### Supply-to-route corollary

A forced-U slot meets at most one member of a pairwise disjoint packing. Thus
p candidate sets survive any t forbidden slots if p>t, under the stated
assumption that avoiding those slots is the only remaining selection
restriction. Every nonselected helper remains at its actual U record.

The route inequality K+H+D<=b follows directly from the exact crossing formula.
Two complete r-selected helper configurations differ in at most 2r slots,
giving H<=2r in slot units and H<=2rW in clock units. These bounds cannot erase
mandatory neutral changes. The two singleton intervals at neutral cost m and
helper cost one give the claimed counterexample at budget m.

The theorem does not prove any unprovided packing, infer ordinary per-port
phase preservation from nodewise-current preservation, or permit consumptive
no-return policies to be treated as reversible graphs. The integer label 11
is not a structural prime: the displayed linear parity recoding is valid.
Likewise a slot bottleneck 11 implies clock bottleneck 22 only under suitable
actual cost information, not by factorization alone.

## 3. Suspension proof audit

### Explicit class-descent obstruction

The independently constructed sources have labelled outdegrees

    T2: (2,1,2,2,3),     T3: (3,1,1,2,3).

The permutation (1,2,3,0,4) carries every old target arc from K2 to K3.
Full permutation enumeration gives canonical old code 8 for both, using the
audit's big-endian upper-triangle encoding. The enlarged canonical codes are
34881 and 34889; their converse codes are 34946 and 34818. Thus the merged
canonical codes are also different. Each seven-owner canonicalization uses
all 5,040 permutations, not a degree partition shortcut.

The enlarged score multisets and square sums are exactly those stated in the
producer note: (1,2,3,3,3,4,5) with 73, and (1,2,2,3,4,4,5) with 75. For a
seven-owner tournament, replacing each degree d by 6-d preserves the sum of
squares, since the total degree is 21. This gives an independent converse
obstruction without requiring canonicalization.

In each example the new ports are the unique degree-one and degree-five
owners. Their unique outgoing and incoming neighbours recover f and l.
The inner restriction is literally K. The failure is therefore a failure of
single-valued forward class descent, not a failure of port recognition or
of the reverse readout. No metagraph harmonicity conclusion follows.

### All-size triangle identities

Reversing every arc except the path changes the inner degree vector to

    d_K(v) = n-1-d_T(v) + 1[v=f] - 1[v=l].

After suspension the inner degrees are n-d_T(v), while the new port degrees
are 1 and n. Substitution into

    c(A) = choose(|A|,3) - sum_v choose(d_A(v),2)

proves both identities in the note. The formula for c(A) follows because every
transitive triple has exactly one vertex of induced outdegree two.
Subtracting the inner count from the enlarged count counts exactly the
triangles that touch a new port.

As an additional direct control, the two separate boundary contributions are
n-d_T(f) at the first outer port and 1+d_T(l) at the last. No cyclic triangle
contains both new ports. Their sum is the claimed boundary correction.
The original fixed-forward-path universe is exactly
1+2+8+64+1024=1,099. It is not a full tournament or metagraph census.

## 4. Added two-port generating-function audit

The literal adjacency construction agrees with the claimed block matrix.
For the old block A=I-zK, the Schur complement of A in I-zH is

    [ 1-z^2 a,       -z^2 b       ]
    [ -z-z^2 c,      1-z^2 d     ].

Its determinant is exactly

    Q=1-z^2(a+d)-z^3 b+z^4(ad-bc).

The constant term of A is the identity, so the calculation is valid in formal
power series without a spectral or convergence assumption. Multiplying Q by
det(I-zK) gives the actual polynomial det(I-zH).

Formal logarithmic differentiation gives

    -z Q'/Q = sum_(r>=1) (tr(H^r)-tr(K^r)) z^r.

Since K is induced, the difference counts precisely rooted closed walks that
visit the new ports. At lengths three and four every closed walk in a
tournament is simple: a repeated vertex would force a loop or a directed
2-cycle. The trace coefficients can therefore be divided by three and four,
respectively. This reasoning does not extend to arbitrary lengths.

The constants a(0), b(0), d(0) vanish because f and l are distinct. The fourth
coefficient of -log Q is consequently

    e_f^t K^2(1-e_f) + (1-e_l)^t K^2 e_l + K_(f,l).

The first two terms count four-cycles through exactly one new port. The last
counts the only possible pattern through both: a->f->l->b->a. The independent
checker enumerates actual cycles and verifies this split, not only traces.
It also constructs determinant coefficients by the Leibniz permutation
formula, independently of the producer's Newton-trace recurrence, then
compares formal series through degree n+8.

For the two explicit order-five witnesses, the old determinant polynomial is

    1 - 2z^3 - z^4

in both cases. The new polynomials are respectively

    1 - 9z^3 - 13z^4 - 17z^5 - 10z^6 - 2z^7,
    1 - 8z^3 - 11z^4 - 15z^5 - 9z^6 - 2z^7.

Their boundary four-cycle counts are 12 and 10. These are additional explicit
records of the missing endpoint incidence, not a spectral-transfer theorem.

## 5. Reproduction and coverage

Both independent checkers use explicit exceptions, no optimization-sensitive
assertions, and import no producer checker.

- `independent_check.py`: **2,340,879 checks in each mode**.
  - All 4,165 weighted graphs of orders 1..4 over {0,2,7,absent}, plus 120
    seeded graphs of orders five and six. Literal simple paths, all cuts,
    exhaustive trees, and every minimum-total tree are compared.
  - Every labelled poset through four elements: counts 1,3,19,219, enumerated
    by direct transitivity testing rather than closure of forward DAGs.
    There are 425,880 weighted interval pairs and 68,174 shortest-route pairs.
  - 968 three-helper quotient instances, tested at 7,744 budgets, with
    Z x C3 charges, zero and heavy weights, impossible helper fibres, and
    symmetric nonmetric helper costs. These are abstract controls, not claims
    that every helper cost table is a native Hamming metric.
  - The two explicit canonicalized witnesses; all 1,099 original fixed-path
    sources; all 32,768 additional fixed-path order-seven sources; and two
    full 40,320-permutation eight-owner canonicalizations under renaming.
  - 1,286 forbidden-U cases, 5,340 equal-margin helper pairs, and the stated
    hostile scope controls.
- `independent_generating.py`: **6,740 checks in each mode**.
  - All 1,099 sources at orders 2..6 for literal four-cycle boundary counts.
  - All 75 sources at orders 2..5 for independent permutation determinants
    and exact walk/response series through degree n+8.
  - Three additional order-six series cases, including eight-owner
    determinant expansions.
- Producer reruns reproduce **2,386,018**, **12,118**, and **3,019** checks per
  mode for the bottleneck, original suspension, and generating-function
  checkers. Their own packet manifests verify.

The normal and optimized results agree after excluding only elapsed time and
Python's optimization flag. Together, the two independent checkers execute
2,347,619 checks per mode. Finite checks corroborate, but do not replace, the
all-size proofs and the native authentication hypotheses.

Run from this audit directory:

    python independent_check.py --output receipt-normal.json
    python -O independent_check.py --output receipt-optimized.json
    python independent_generating.py --output generating-normal.json
    python -O independent_generating.py --output generating-optimized.json
    python verify_manifest.py
    python -O verify_manifest.py

The last two commands authenticate the frozen delivery. Rerunning a checker
with `--output` overwrites a timed receipt, so its hash will no longer match the
original manifest. Run without `--output` for nonmutating reproduction.

## 6. Final publication portability refresh

The producer renamed its three inherited Markdown snapshots to .txt, preserving
every byte and SHA-256 value. The mathematical THEOREM.md, SUPPLY_TO_ROUTE.md,
check.py, and verify_manifest.py remained unchanged. Its SOURCEPINS.json now
resolves the archive names and records original source locations; its README
and timed receipts were refreshed. This was an archival/provenance repair, not
a mathematical revision.

This audit adopted .txt names for all copied Markdown inputs, updated its
source provenance and the copied SOURCEPINS.json, and regenerated both modes
of the main and boundary-series independent checks and bottleneck producer
replays. Check counts and all mathematical outcomes are unchanged.
PORTABILITY_UPDATE.json records before/after hashes and identifies the one
intentional metadata-content change. The separate global-observer audit was
not modified.
