# Positive linear streaming of fixed-H tournament multiplicities

**PROVED:** the generic matrix-to-FRACTRAN compiler and its exact microstep law.
**PROVED, with a cited strong-order bound and finite-exact atom inventory:**
the emitted 397-fraction stream simultaneously counts every ordinary tournament
isomorphism class with each odd Hamiltonian-path count from 1 through 21.
No novelty is claimed for the inherited SCC generating functions or the
all-order exclusions of 7 and 21.

## 1. Mathematical input and lost information

Let N_n(h) count n-vertex tournament isomorphism classes with exactly h directed
Hamiltonian paths, with N_0(1)=1 for the empty word. Every tournament has a
unique ordered sequence of strongly connected components (SCCs), and its
Hamiltonian-path count is the product of their counts. A singleton has count 1.
A nontrivial strong component has at least 3 paths.

The complete strong inventory for h<=21 is the list (order,h,multiplicity)

    (3,3,1), (4,5,1), (5,9,2), (5,11,1), (5,13,1),
    (5,15,2), (6,15,2), (6,17,2), (6,19,2).

The sibling spectrum packet regenerates this list by exhaustive canonical
augmentation through order six, with direct permutation and odd-cycle checks.
Busch's cited inequality H(T)^3 >= 5^(n-1) for strong tournaments excludes every
order n>=7 when H<=21, because 5^6=15625 >21^3=9261. Thus this is a complete
inventory in path count, not an extrapolation from small order.

Taking the first SCC gives, for n>=1,

    N_n(h) = N_(n-1)(h)
             + sum_(k,q,c in inventory; q divides h) c N_(n-k)(h/q),

with negative-order terms zero. The recurrence counts ordered SCC words, not
labelled tournaments, paths modulo automorphisms, or converse-merged classes.
Different orders of strong components are distinct words. It forgets endpoint
incidence and metagraph wiring. No adjacency recovery is asserted.

These identities and their rational, eventual-polynomial consequences are
inherited from the 2026-10-21 packet
`actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`.
The new deliverable here is their simultaneous, positive, literal FRACTRAN
realization, with a general compilation lemma and exact execution cost.

## 2. Linear state

Order h as 1,3,5,...,21. For each h keep six registers

    x_(h,k)(n)=N_(n-k)(h),  k=0,...,5.

The recurrence above and shifts x'_(h,k)=x_(h,k-1) give a 66-by-66
nonnegative integer matrix A. The initial vector is e_0. Thus every future
state is A^n e_0. No divisions or signed subtraction are required.

The zero registers for h=7 and h=21 are mathematical consequences of the
certified atom inventory and divisor recurrence. Merely omitting transitions
would not prove these exclusions without the inventory-completeness theorem.
The checker deliberately adds a fake H=7 atom and sees a nonzero channel.

## 3. Generic compiler and proof

For any d>=1 and any nonnegative integer d-by-d matrix A, assign disjoint
primes x_i (state), y_i (scratch), C_i,D_i (drain control), E_i,F_i (return
control), Q (checkpoint), and w (round clock). Exactly one control prime has
exponent one; all other control exponents are zero.

At each drain state C_i, the first fraction consumes one x_i and replaces
control by D_i while multiplying by product_j y_j^A[j,i]. A trampoline returns
D_i to C_i. The fallback fraction advances to C_(i+1) only when x_i is empty.
After all inputs are consumed, scratch equals A times the old state, and
all state registers are empty. A second pass transfers each y_j unit to x_j,
again using a trampoline and zero fallback. Its final fallback returns to Q
and increments w. The initial Q step enters C_0.

Each fraction has a control prime in its denominator and a *different*
control prime in its numerator. Hence reduction cannot cancel the guard.
Fractions for inactive controls are ineligible. At an active drain, the
consuming row precedes its fallback, giving exactly the intended zero test.
Every pass terminates because it consumes a finite register while filling
only the other bank. Induction proves that checkpoint n contains x=A^n v,
y=0, clock=n, and the single control Q, for any initial nonnegative vector v.
The program runs forever through checkpoints; it is not a halting counter.

The construction has exactly 6d+1 fractions and 6d+2 distinct register primes.
These are counts for this compiler, not minimality claims.
A round starting at v and ending at Av takes exactly

    1 + 2d + 2 sum_i v_i + 2 sum_i (Av)_i

literal FRACTRAN steps: one entry step, two steps per consumed unit, and one
fallback per register in each pass. Thus through n rounds the cumulative cost
is

    n(1+2d) + 2|v_0|_1 + 4 sum_(j=1)^(n-1)|v_j|_1 + 2|v_n|_1.

All claims assume the stated one-hot-control initialization and zero scratch.
Arbitrary integers with multiple control factors are outside the contract.

## 4. Concrete interface

`bounded21.fractran` is the actual ordered list of 397 reduced fractions.
`program.json` lists every register prime. Start at the ordinary integer 5438
(2719 times 2). Checkpoints are precisely the orbit states divisible by 2719.
At checkpoint n, the exponent of 2729 is n. In `program.json`, the prime
`x[6*i]` has exponent N_n(2*i+1), for i=0,...,10. Other x slots are lagged
values; all y primes have exponent zero at checkpoints.

This program counts infinitely many orders using a finite certified strong
inventory. It does not enumerate individual tournaments at runtime. It does
not prove that every other odd h is attained, or discover strong inventories
for arbitrary h by itself. General bounded-h compilation first needs a
certified inventory, and hence may be expensive.

## 5. Verification

Run `python3 stream.py` to regenerate the list and metadata, then
`python3 check.py` and `python3 -O check.py`. The checker includes full literal
integer execution scanning the entire ordered fraction list through n=10,
valuation execution through n=32, closed-form controls, exact step counts,
small arbitrary matrix programs with full integer scans, and negative inputs.
Independent audit is in the sibling independent packet. Execution uses an
inactive-control index only in the valuation implementation; the literal
integer implementation does not skip any earlier fraction.

## 6. Refinement: 67 fractions after the finite prefix

The inherited closed forms allow removing all lag registers once n>=6.
Keep exactly eleven states v_h=N_n(h). Their simultaneous updates are

    v_1'=v_1; v_7'=v_7; v_21'=v_21;
    v_h'=v_h+v_1 for h=3,5,11,13;
    v_h'=v_h+2v_1 for h=17,19;
    v_9'=v_9+v_3; v_15'=v_15+2v_5.

Initialize at n=6 with (1,4,3,0,5,2,2,6,2,2,0). Subtracting consecutive
values in the inherited exact polynomials proves each update; alternatively
it follows from the SCC recurrence and the specified initial state. The
67-fraction `tail21.fractran` is the generic compiler applied to this
11-state matrix. Its exact integer seed, checkpoint, clock and output primes
are in `tail_program.json`. The clock starts at 6, so output n is literal.

This is an 83.1 percent reduction in fraction count relative to the
397-fraction stream, at the cost of starting after the finite prefix.
It is not a claim of a minimum program or of an equivalent arbitrary-input
map. Both count the same eleven sequences on their stated common domain.
The zero channels are initialized at zero and stay there, with their
mathematical interpretation certified by the graph-theoretic proof.

## 7. Alias refinement: 43 fractions with the same valuation readout

Two pairs of entire sequences coincide: N_n(11)=N_n(13) and
N_n(17)=N_n(19). The inventory gives the same single-atom formulas in each
pair. Both forbidden channels are identically zero. Consequently only seven
distinct nonzero sequences need state registers: h=1,3,5,9,11,15,17.
Applying the same matrix compiler to those seven states gives the actual
**43-fraction** `compact_tail21.fractran`, again starting at order six.

`compact_tail_program.json` lists a prime for each of the eleven requested
outputs. The equal channels deliberately share a prime; the two zero channels
use one extra prime absent from the seed and every fraction. Thus readout is
still just a prime valuation, not an external polynomial calculation. No
claim that eleven independent registers remain is made. Distinct-sequence
aliases are justified mathematically, rather than guessed from a finite
prefix. This is 35.8 percent fewer fractions than the 67-fraction tail.
The exact generic compiler proof and runtime law apply with d=7.
