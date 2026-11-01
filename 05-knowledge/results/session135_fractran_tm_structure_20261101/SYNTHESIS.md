# Exact bridges, and what each simplification forgets

## 1. A fixed program can be much simpler than its language

FRACTRAN can simulate general computation. That fact does not make every
particular list universal or its particular orbit questions undecidable.
For Conway's exact PRIMEGAME list, the positive-control budgets T,Q,k bound
an entire phase in a finite box. A onehot escape predicate and a canonical
candidate return then classify every positive input: eventual period2,
eventual period6, or exact coalescence with a spectator-scaled seed2 orbit.
The original list never halts; the new basin classifier always does.

This closes the basin gap left by the inherited complete classification of
periodic points. The rational periodic-point densities and the newly
computable eventual-basin densities are different quantities. A bounded
control graph, by itself, would not give this result; the counter budget and
actual reset/coalescence proof are essential.

## 2. Register primes are labels unless an arithmetic observer fixes them

Renaming the finitely many register primes conjugates a FRACTRAN execution
exactly, provided output designations are renamed too. It preserves the
ordered rule trace and period, but need not preserve natural density of
integer inputs. Swapping2 and3 conjugates [3/2,2/3] with [2/3,3/2], while
their periodic-point densities are5/12 and5/9. Thus a numerical prime-group
analogy needs an explicitly retained arithmetic observer.

This is a multiplicative conjugacy, not a ring map: the2/3 swap sends4 to9,
whereas the sum of the images of2 and2 is6. Additive or square-sum arithmetic
is not preserved merely by renaming register primes.

The spectator factor is another such observer. The positive standard-orbit
basin density refers to coalescence AFTER removing that fixed factor. If
S>1, the literal trajectory need not emit pure powers of2. Indeed any input
whose PRIMEGAME orbit ever becomes a pure power of2 must be smooth over the
ten program primes. There are at most O((log X)^10) such inputs below X, so
that literal-observation set has natural density zero. This does not conflict
with the >0.58105 lower density for the spectator-normalized standard basin.

The finite prime alphabet also does not bound the values encoded in its
exponents. PRIMEGAME's unbounded primes appear as exponents, not new prime
factors of the state integer.

## 3. Three different finite abstractions

- The positive-polynomial cap is exact at certified macro checkpoints for
  bounded-output observables. Decrement microsteps do not respect it.
- The horizon-dependent cap at L times the largest denominator exponent
  preserves the first L instruction/halting observations. With net-update
  testing it yields exact rational densities for points of bounded period.
- PRIMEGAME's three basin labels retain eventual behavior, not the complete
  future trace or elapsed time.

The word accelerator instead retains an actual trace word, exact repetition
count, higher-priority guards and primitive cost. Its infinite-word certificate
is sufficient when found; it is not guaranteed for arbitrary nonhalting
computations. These are complementary interfaces, not interchangeable
claims that an unbounded counter system is finite-state.

## 4. Actual Turing realizations with explicit costs

The pure-blue tournament stream becomes a fully explicit76-state single-tape
TM with257 symbols and10,926 ordinary transitions. The compact multiplicity
stream has a15-tape machine with293 fully specified ordinary rows: inactive
heads explicitly read/write the origin marker. No wildcard expansion or
hidden arithmetic is used there.

On the actual radix80 invariant face, removing two identically zero state
coordinates reduces25 fractions to13 and nine data tapes to five. It saves
exactly four fraction steps and twelve TM transitions per round. This is a
model-fixed comparison, not a claim that state count measures a predicate's
intrinsic difficulty. The quadratic unary costs still dominate large runs.
The complete general FRACTRAN/register compiler retains failed-test rollback
and avoids cancelled control guards; its inverse has declared initialization
and observation boundaries.

## 5. Square-sum graphs need their named vertices

The14-fraction square predicate and explicit path/cycle TMs certify supplied
witnesses without building the whole graph. Named U/X marks enforce a
permutation. The word(1,8)^4 at n8 passes local square edges, length, label sum
and the corresponding squared-root moment, but is rejected for repeated
vertices. This refutes those particular lossy summaries, not every imaginable
arithmetic encoding of a visited set.

For a Hamiltonian cycle, sum_edges(u+v)=n(n+1). Independently running the
square FRACTRAN test on each edge therefore costs exactly5n(n+1)+2n steps,
independent of cycle order. That is the edge-testing cost only, not a whole
permutation verifier. The full TM's cost, after removing the start-vertex
term, instead records6 sum_edges sqrt(u+v) plus a fixed polynomial in n.
This is an invariant of the integer-labelled cycle under rotation/reversal,
not a claim about arbitrary unlabelled graph isomorphisms or equal-work speed.
The retained32-cycle has root sum180, arithmetic-edge cost5344 and full-TM
cost6662 for its start at1. `interface_checks.py` checks these exactly.

## 6. The odd eighths acquire different certificate roles

For the harmonic core {1,...,6} at threshold1/8, the four familiar odd eighths
are all safe but have two different types:

- 1/8 and7/8 are endpoints of positive safe intervals, with one tight owner.
- 3/8 and5/8 are isolated safe points, pinned by speeds3 and5 on opposite sides.

This is an exact pointwise connection, not an identification of unrelated
mod8 structures. The general harmonic completion retains all such isolated
unit residues. At the LRC14 deep well it supplies3/14,5/14,9/14,11/14 in
addition to twelve positive intervals. Measure alone loses them.

Replacing cross-multiplied witness counters by the rounded threshold
K=ceil(cb/d) keeps the exact predicate. For the V26 witness2/27 the same
53-fraction program drops from703 to53 steps and33888 to2440 input bits.
Modular preprocessing remains explicit. LRC14 is still OPEN; the full-row
safe-set algorithm is finite but magnitude-dependent, and the global-height
reduction retains its stated lower-dimensional assumption.

## 7. Circulant words share a certificate principle, not a proof by analogy

For a circulant tournament of odd order q>=3 on Z/q with connection set S, a path
starting at0 corresponds exactly to a word of q-1 letters in S whose partial
sums are all distinct modulo q. Translation gives H=q times the number of
such words. A constant step a is Hamiltonian iff gcd(a,q)=1. Hence every
such tournament has at least q phi(q)/2 constant-step Hamiltonian paths.
For q9, the locally permitted word3^8 visits only three residues and fails.
The role of distinct partial sums is the same kind of global named-state
constraint as square-sum U/X marks; local edge permission alone does not
supply it. This correspondence is recovered elementary/concurrent work,
not a newly claimed existence theorem. Root controls cover every connection
set at q3,5,7 and the composite hostile.

## 8. Conjecture search is another interface

Concurrent Goldbach/Fermat/Legendre programs reuse a terminating primality
kernel but change the outer controller and exit. Our theorem for the original
PRIMEGAME list does NOT classify those modified whole programs. Their
halts-iff claims, original control failures and subsequent source audit are pinned separately
in CONCURRENT_BOUNDARY.md; they are not used as proved dependencies here.
Encoding a conjecture as nonhalting preserves a question; it does not answer
it. A test that consumes a counter must also declare whether its observed
value is read before or after that consumption.

No general halting, Collatz, Goldbach, sixth-Fermat-prime, Legendre, LRC14 or
unrestricted tournament-spectrum resolution is asserted by these bridges.
