# Session 135 independent machine-semantics audit

**Verdict: PASS in the exact scopes below.** All-input conclusions come from
proofs/invariants and the explicitly inherited S126/S129 mechanisms, not from
extrapolating finite tests. There is no literature-priority claim, general
FRACTRAN halting decider, LRC(14) proof, or universal odd-H-spectrum theorem.

This directory was produced by a separate auditor. Producer executors are not
used to obtain the reference trajectories. Compiler APIs are called only to
obtain machines under test; the independent interpreters apply actual finite
transition rows or ordered rational lists. All checks raise explicit exceptions
and remain active under `python -O`. Paths are sibling-relative. No external
package, network access, Git checkout, or binary table file is required.

## Evidence inventory

| Area | Independent executed universe | Boundary |
|---|---|---|
| PRIMEGAME basins | 59,049 full exponent vectors in `{0,1,2}^10`; 4,375 one-hot states with data counters 0..4; 166,536 independently traversed states; 2,000 literal integer transition controls | Reference detects arbitrary repeats or canonical candidate arrival, not the producer's six-step/Boolean shortcut |
| Basin density | All 1,024 cap-one valuation cylinders with independent basin labels and exact weights; three analytic lower bounds; precision-tail controls m=0..64 | Density limits and polynomial precision complexity use proved tail bounds/totality, not frequency extrapolation |
| Square FRACTRAN | Literal integers `7*2^s` for every s=0..1024, with complete first-eligible scans | Exact terminal shapes and every stated step count checked; all-s result follows odd-subtraction invariant |
| Square TMs | 43,690 all-word machine cases, lengths 0..7; 11,826 permutation cases, orders 1..7; 400 isolated square-subroutine inputs; path15 both directions and cycle32 every rotation plus reversal | Executes full 50,625/72,900-row delta tables; checks exact runtime/scratch cleanup rather than only acceptance |
| Machine translations | 1,764 serialized TM instruction cases; 1,764 serialized reverse-FRACTRAN cases; 210 original-fraction boundary cases; 1,000 general-list cases; eight self-loop controls | All four published tables; source-fraction steps, IR steps, TM transitions and reverse-fraction steps remain separate |
| LRC comparator/cells | 2,401 five-fraction literal cases; 1,681 literal 12-transition TM cases; all four thirteen-speed V26 examples including input hashes; independent Farey construction k=1..120 | Rational arithmetic preprocessing is an explicit sidecar; no full LRC search TM is claimed |
| Bounded periodic points | 39,204 literal-integer horizon/return checks; PRIMEGAME complete 59,049-state cap2 density; 81,648 positive-polynomial cap checks; affine-cost example | Exactly first L rule observations, not a post-L halt observation; point sets are not eventual basins |
| Word acceleration | 341,775 guard-membership checks; 7,130 certified words on 81 two-row programs and 25 inputs each; 3,066 exact finite maximality checks; 4,064 forty-round infinite-prefix controls | Infinity is justified by the all-t guard proof, not the finite prefix controls; this is not a termination decider |

All seven normal/optimized receipt pairs are byte-identical. A relocated copy
with no repository and **no `.gz` files** reproduces the same receipts. The
square-TM loader independently expands finite disjoint clauses when needed,
checks the complete Cartesian state/read-key set, rejects duplicates, and
verifies the frozen full-stream SHA-256/count before execution. The translation
loader authenticates every literal-row shard and aggregate count before use.

## Structural proof findings

### Fixed PRIMEGAME: a positive-phase budget closes the basin classification

For the exact fourteen-row list, the coordinates
`T=e2+e5+e23`, `Q=e3+e7+e17+e19+e23+e29` and the six-control count k have the
claimed deltas. If k>0, a lone-control denominator enables some A..K row, so
L/M/N cannot fire. Thus T is constant and Q,k are nonincreasing. Conditioning
on a complete control vector fixes e23; consequently a+c=T-e23 has at most
T+1 possibilities. The valid bound is

```
B=(T+1)*binom(Q+2,2)*binom(k+6,6).
```

This is an actual full-state bound, not a finite quotient of the instruction
trace. Before k<=1, every new state lies in this box. A repeat is a genuine
periodic orbit; the inherited S129 classification permits only periods 2/6.
The seven-row one-hot escape bit is preserved by each priority-guarded move.
Its false branch reaches a JK cycle; its true branch excludes that only
one-hot cycle and therefore reaches I and then zero control. At that exit
b=c=0 and a=T. The exact reset `L^T M^d N (EF)^T K` takes 3T+d+2 steps to
`C(T+1,T)`. S126 proves that seed2 visits every candidate state of that form,
not every pure power `2^n`. The latter distinction is essential.

Hence basin membership has the stated polynomial decision bound in binary
input length. The actual entry bound B+3T+Q+2 is also valid for initially
positive control; it is not a next-prime-output waiting-time bound. The
original program itself never halts, since its last denominator is one.

Exact finite-prime valuation cylinders have weights
`product_p (1-1/p)*p^(-e_p)`. Finite-box tails have mass
`1-product_p(1-p^(-h_p-1))`, so every valuation-defined basin has a natural
density. The proved decider makes it computable; H=m+4 gives polynomial-time
precision in m requested bits, with the explicitly large degree19 bound.
The three stated rational lower bounds independently recompute as

```
prime:  50136286961664/86285158710179
cycle2: 753557212/2002524095
cycle6: 37396469504/151504041213525
```

These are lower bounds on eventual basins, not exact basin densities.
The last is the inherited periodic six-cycle density, recomputed from its
six reciprocal base states and scale density.

### Finite machine translations retain every guard and cost model

The general FRACTRAN compiler reduces fractions before denominator tests;
failed partial tests restore exactly the consumed prime exponents before
trying the next row. This proves priority preservation for arbitrary lists,
including cancellation and repeated-prime guards. The optimized one-hot
compiler is a narrower theorem, with explicit data/control separation.

INC/DECJZ unary TM primitives cost exactly 2u+3 at old register value u;
GOTO costs one. Packed-track writes preserve other bits, and inactive
multitape heads explicitly remain at `#`; no wildcard arithmetic is hidden
in the translation tables. Reverse compilation uses fresh trampoline primes
and takes exactly two fractions per nonhalting IR instruction, including
self-loops. This is not two fractions per original source row.

The pure-blue drain-loop cost and positive-matrix norm identity pass algebra:
for y=Ax, the matrix round takes
`||x||²+3||y||²+5|x|+9|y|+6d+2n+4` TM transitions. The scratch increments
telescope, explaining dependence on endpoint norms, not matrix sparsity.
The migration face preserves `8C-9t`, yields `F'=80F+3I`, and retains the
integer-lattice condition F=1 mod26 on the stated ray. Count streams do not
reconstruct graph arcs or discharge general spectrum membership.

### Named square witnesses and closed LRC cells need their sidecars

The square program retains distinct source/target controls and ordered zero
fallbacks. Its odd-subtraction invariant proves all-input termination and the
exact terminal forms. The TMs' named U/X cells reject repeated labels; checking
local square edges or root moments does not suffice, as `(1,8)^4` demonstrates.
Closing the last edge is a separate cycle obligation. The exact runtime
formulas and parity/mod12 corollaries pass. All-input O(L+1) refers to actual
expanded unary word length, not binary input or Hamiltonian search complexity.

For LRC, the ordered pair-drain fractions preserve differences and terminate
before any verdict row. Rounded `K=ceil(cb/d)` gives exactly the same residue
inequalities with smaller unary exponents. The V26 rounded examples take
53/28 steps with 2,440/2,431 input bits. Preprocessing residues/prime powers is
not performed by the comparator itself. The independent Farey-gap proof
retains denominator-sum strata k+1 (positive intervals) and k+2 (singletons),
which yields the exact harmonic completion. A measure-zero safe set can be
nonempty. Any global finite-row bound keeps its stated lower-dimensional LRC
assumption; the finite smoke census is not a global theorem.

### Priority-aware acceleration: maximal full repetitions, not guessed loops

For base b, net direction D and denominator vector a, eligibility is the
integer interval defined coordinatewise by b+tD>=a. Positive D gives a ceiling
lower bound, negative D a floor upper bound, and a failing constant coordinate
makes the interval empty. Since the first occurrence is valid, the first
future own-guard failure or earlier-row eligibility is exactly the first
invalid repetition. Taking the minimum over every phase gives the maximal
number R of **full** repetitions. A prefix of the next repetition may still
execute. Every earlier complete repetition is valid, so the endpoint x+RD
and primitive cost R|w| are exact. If no finite obstruction exists, all phases
and their priorities remain valid for every t>=0, certifying the infinite word.
The hostile `[1/10,5/2]` sends 8→20→2→5 and halts in three steps; blind drain
acceleration would incorrectly return 125. Budget exhaustion stays inconclusive.

Fixed-horizon saturation has a different, exact scope. If M_p is the largest
denominator valuation at p in a reduced list, capping the initial p exponent
at T*M_p preserves the first T rule/halting observations. At phase j<T an
initially saturated coordinate stays at least (T-j)M_p, so every denominator
comparison has the same truth value; unsaturated coordinates remain equal by
induction. This does not preserve final valuations or an unbounded trace.
Likewise nonnegative integer matrix updates commute with coordinatewise capK
for K>=1, so threshold observables have an exact finite quotient. Subtraction
does not: exponents12 and18 agree under cap5 plus residue mod6, but eight
decrements give capped values4 and5. Fixed residues do not repair countdown.

For any fixed program and L>=1, this finite-horizon cap also preserves whether
there is a return within L steps: a completed prefix returns iff its net
valuation increment is zero, independently of the starting exponent vector.
Thus least-period-at-most-L POINT membership factors through the finite cap
box. Each capped coordinate y<C has cylinder weight (1-1/p)p^(-y), and the
top cell y=C has weight p^(-C), including C=0. Summing a finite union gives an
exact rational natural density. No universal period bound or eventual-basin
conclusion follows. First L observations mean times0 through L-1; for [1/2],
inputs4 and8 share their first two decrements, but have different halt status
immediately afterward.

Independent whole-integer cap2 enumeration recovers PRIMEGAME's inherited
period-two point density 2654208/36409529. The toy [3/2,2/3] has periodic-point
density5/12 and eventual two-cycle basin density2/3; reversing the two rows
changes the point density to5/9. These values were also derived directly from
the valuation guards, not only from the finite cap census. The added positive-
polynomial saturation proof requires nonnegative INTEGER coefficients; the
affine cost progression is for finite repetition counts and requires its
explicitly stated phase-domain hypothesis.
REPLAY.initial.json retains the twelve pre-extension script/mode replays.
The current complete fourteen script/mode replays are in REPLAY.json and the
byte-identical REPLAY.normal.json, REPLAY.optimized.json and REPLAY.portable.json.
FOLLOWUP_REPLAY.json records the earlier targeted bounded-period replay.

Finite prime permutations conjugate literal first-eligible execution, including
halt status, least periods and designated output valuations when their prime
labels are renamed. The independent literal-integer controls include the swap
2<->3 and the non-involutive cycle2->3->5->2, with unchanged spectator factors.
The partial relabeling must be injective; its extension is completed on the
union of source and target primes, fixing only primes outside that union.
The conjugacy is not a natural-density invariance: the same toy two-cycle
system has point density5/12 before the2/3 swap and5/9 afterward. Arithmetic
input weights and integer bit sizes are not supplied by the abstract counter
transition graph. This strengthens the inherited prime-label gauge boundary,
not a prime/group numerical identification or a measure-preservation theorem.

## Unpublished wrapper correction and final end-to-end verification

The root's isolated final replay exposed a wrapper bug introduced when
`--only` was added: the argparse namespace `args` was overwritten by the
inner `for mode,args in ...` flag list. A full invocation failed at the next
outer audit with `AttributeError: list has no attribute only`. Running only
the last audit had masked that path. This was an unfrozen verification-harness
error, not a mathematical program change; no first-checkpoint mathematical
claim or individual component receipt was invalidated.

The repair separates `options` from `mode_flags`, selects the requested audit
labels before iterating, and checks the completed audit/mode count. The
regression `--only density periodic_points` now exercises the second selected
audit; its four completed runs are recorded in HARNESS_REGRESSION.json.
Most importantly, the FULL seven-audit wrapper has now reached exit0 normally,
under parent Python `-O`, and in a fresh UTF8-only relocated packet with no
repository or binary `.gz` tables. Each authenticated all41 producer inputs and
completed all14 normal/optimized component runs. The three complete wrapper
receipts are byte-identical. These are terminal results, not pending runs.

Only the wrapper, report and receipt/manifests changed in this repair. The
producer files, mathematical audit scripts, and their expected component
receipts remained unchanged. The older pre-extension full replay is retained
as REPLAY.initial.json rather than represented as a test of the new wrapper.

## Inheritance and remaining scope

The closest recovered mechanisms are S126's divisor loop, S129's complete
cycle classification and no-finite-trace-observer theorem, S132/S133 guarded
compilation and valuation-output cost separation, and S134's authenticated
migration. The new basin Boolean quotient preserves eventual classification,
not the full block trace, so it does not contradict S129's infinite-memory
obstruction. Fixed-list decidability likewise does not contradict language
universality. External novelty and fresh literature verification are not
conclusions of this audit.

No mathematical blocker remains within these scoped claims. Hash pins attest
which bytes were inspected; they do not replace proofs. The parent owns final
repository publication and remote-checkpoint verification.
