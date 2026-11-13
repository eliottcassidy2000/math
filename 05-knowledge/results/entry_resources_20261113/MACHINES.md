# Exact machine translations and their construction boundaries

**PROVED scoped / inherited mechanisms identified below.** A representation
of a finite construction is not a new tournament realizing its output.

## M1. The exact first25-prime sieve as a local clock machine

For F(k)=99·2^k+23, k≥1, the first25 primes are exactly those below100.
The inherited active divisibility classes (prime: residue modulo period) are

    5:3 mod4; 13:1 mod12; 17:1 mod8; 19:9 mod18;
    29:27 mod28; 37:23 mod36; 53:25 mod52; 59:35 mod58;
    61:0 mod60; 67:56 mod66; 83:15 mod82.

The other14 primes never divide F(k) for k≥1. In particular11 divides99
but not23, while7 sees the repeating nonzero residues3,4,6. These are
actual congruence mechanisms, not consequences of being small primes.
At k=0 the separate even-prime exception exists; the combined acceptance
still rejects because61 also divides122.

The29 and37 clauses are already excluded by5. Removing them gives least
period L=61,209,720 and32,353,278 surviving residues, with density
5,392,213/10,201,620. This minimization was independently obtained here
and is already proved in incoming
[Session292, arithmetic §1](../session292_triple_bridges_20261113/arithmetic/PROOF.md).
Credit belongs to that parallel refinement as well as the inherited
[endpoint compiler E7](../endpoint_coverage_compiler_20261111/CONNECTIONS.md).
The old428,468,040 wheel is valid and contains seven copies of this one.

Here is a smaller *local-state representation*, different from Session292's
360-cell conditional CRT decoder. Use disjoint directed cycles of lengths

    8,9,5,11,13,29,41  (total116 binary cells).

Place one1 at position0 in each cycle; every cell copies its predecessor
at each step. At time k the markers record k modulo each prime power
in L. For each of the nine essential clauses, a fixed Boolean conjunction
tests the appropriate projected residues. Reject if any conjunction holds.
The receipt specifies every bitmask of this readout.

This is a binary cellular automaton on a finite directed graph, with a
fixed global Boolean readout. The116 count excludes the readout circuitry.
It is not claimed to be a one-dimensional elementary cellular automaton,
a universal computer, or a minimal cell implementation.

For each prime q|L an explicit pair k,k+L/q has opposite output. These
seven certificates prove the least period: every proper divisor period
would divide one of the rejected candidate periods L/q. They also prove
that a complete unary DFA for the same acceptance sequence needs exactly
L states. Distinct phases have distinct continuation languages; otherwise
their difference would give a smaller period. Factoring its state into
seven clocks preserves all L phases while making updates and readout small.

The machine is checked against direct modular exponentiation for50,001
successive exponents. Generalized-CRT inclusion/exclusion independently
recovers the exact density. This finite check supports the symbolic CRT
proof; it does not replace the infinite periodicity argument.

The hostile remains explicit:

    k=2+1,162,984,680t, t≥1

passes every small-prime test but has proper divisor419, since
ord_419(2)=418 and the step is lcm(L,418). F(2)=419; later terms are
larger multiples. Session292 proves the stronger result that every
surviving cylinder contains subcylinders forcing arbitrarily many extra
prime factors. These refinements concern divisibility of already realized
tournament counts F(k), not a new arbitrary-target compiler.

## M2. Counts as finite automata with consumable vertices

A directed Hamilton path is a word accepted by adjacency transitions,
using each vertex exactly once. Keep the used-vertex subset as state, or
work over the commuting square-zero algebra

    R=Z[x_1,...,x_n]/(x_1²,...,x_n²).

With adjacency matrix A and X=diag(x_i), the full-support coefficient in
1ᵀX(I−AX)^(−1)1 is H. The inverse is a finite sum because each nonzero
monomial consumes a new vertex. Products of endpoint polynomials count
disjoint covers; reused vertices vanish. A two-phase automaton introduces
one selected reset and counts E_T(O), independently of endpoint-table
convolution. This is the engine used to verify all410676 rejected gadgets.

The inherited [Session291 exact carrier](../session291_cluster_ports_20261113/SYNTHESIS.md)
already explains why elimination must retain original occupancy marks and
addressable ports. The algebra has2^n monomial coordinates: calling the
transition matrix n-state has not removed the exponential resource memory.
The min-plus endpoint recurrence in S1 is another valid projection, but
it retains optimum support and loses path multiplicities.

This yields several faithful encodings of a finite construction: a
straight-line graph program, a subset-state weighted automaton, a finite
arithmetic circuit, or a constraint system with occupancy and endpoints.
A Turing machine can execute each exact finite algorithm, and the standard
TM-to-CA simulation preserves its execution. The repository's audited
[THM-4858, exact machine translations](../../../01-canon/theorems/THM-4858-exact-translations-of-the-goldbach-machines-a-53-symbol-radius-one-cellular-automaton-and-a-568-fraction-fractran-program-from-lengs-25-state-turing-machine.md)
supplies the generic simulation mechanism. Its particular53-symbol count
belongs to Leng's Goldbach machine, not to a Chan–Pak machine here.

The reverse direction is the hard one. A general arithmetic or counting
circuit need not have a tournament realization with the same output.
Even a proposed universal compiler with output2N+1 fails at N=3 and10
because it would realize7 and21. A different offset, scaling, or restricted
source language needs its own physical theorem and initialization guard.

## M3. Chan–Pak as a precise nonhalting statement

The inherited
[Session289 finite witness bound](../session289_joint_births_20261112/absence/PROOF.md)
implies that a realizable h>1 has a witness with at most floor(3log_3 h)
vertices after singleton SCCs are removed. For a purely integral algorithm
use the weaker bound3 ceil(log_3 h), computed by powers of3. It follows
that membership of each fixed h in the global Hamilton-count spectrum is
decidable by a finite tournament search. This may be extremely expensive.

Therefore a fully specified conceptual program can enumerate the positive
odd targets other than7,21, run that finite decision procedure for each,
and halt at the first rejected target. Together with the inherited
impossibility of7,21, Chan–Pak is equivalent to this program never halting.
It is a universal decidable-predicate (Pi^0_1) formulation. No new lower
bound on computability or independence follows from that classification.
No concrete minimal-state TM or CA for this program is claimed.

The program isolates the missing invariant: prove that every target's
finite rejection branch is impossible. The frozen-band version is
equivalent: realize each odd h in [f(n),f(n+1)) from the monoid of strong
cores through n, except7,21. This is broader than closing the473 listed
repair families and may admit a structurally different induction.

## M4. Why all binary residues can still be sparse

Ordinal sums of a directed triangle (H=3) and the strong four-vertex
tournament (H=5) give actual counts3^a5^b. For every r these counts cover
all odd residues modulo2^r: for r≥3,5 generates the subgroup1 mod4,
and multiplying by3 reaches the other odd coset. The cases r=1,2 are
immediate. The script checks r≤12; the group argument proves all r.

Yet the number of distinct such counts below X is O(log² X), and they
miss7 and21 and all primes other than3,5. Complete residue coverage,
finite-state periodicity and target-decreasing arithmetic are therefore
insufficient without a physical entry domain. This is a direct hostile
to substituting a local clock argument for an exact integer realization.

## Connection map and next precise hypotheses

|Source|Target/map|Preserved|Lost or still owed|
|---|---|---|---|
|Actual marked tournament|Quadratic energy polynomial Z_T|delta, selected slope, multiplicity|A polynomial need not be realizable by a tournament|
|Selected-cut gadget|Two-phase occupancy automaton|Exact H−E on every input|A successful arithmetic transition still needs a physical constructor|
|Stable clone pools J1|Two Boolean cubes meeting at empty set|All-length feasible locus|Fixed-order exchange connectivity fails|
|Frozen graceful forest|Consecutive star buffer|Labels, colors, frozen edges|Abstract tree shape may change|
|Square graph Q24|Degree/endpoint constraints|Hamilton-path obstruction|Component cuts alone lose the obstruction|
|Exponent k|Seven local phase clocks|Exact first25-prime sieve|Prime factors outside the sieve|
|Finite target decision|TM then CA execution|Halting iff a target is missing|No invariant proving nonhalting|

**OPEN hypothesis B-FLAG:** B(T) is flag iff T is transitive. The complete
through8 census supports it; the Paley controls are consistent. This is a
question about one exact endpoint budget, not a characterization of all
metagraph structure.

**OPEN hypothesis B-ARITY:** minimal nonface sizes of B(T) are unbounded.
The first4-way obstruction occurs at9; Paley11 realizes a uniform4-way
obstruction. T6 proves unbounded face dimension, which does NOT prove
unbounded minimal-nonface size. These are distinct invariants.

**OPEN constructive route:** replace a frozen endpoint observer by a small
renewable guard, allow changes in its omitted responses, and pay the exact
change with a physical compensator. Incoming
[Session292's sparse comparison-poset grammar](../session292_triple_bridges_20261113/SYNTHESIS.md)
gives actual core choices on which to test this. A successful rule must
provide its changed adjacency, transported selected marking, full target
equation, preserved readiness predicate and decreasing inverse measure.
The graceful buffer suggests proving existence after controlled enlargement;
T1,T2 and J1 specify which simple enlargement/exchange models are inadequate.

The later incoming
[Session293 uniform-module ideal](../session293_cycle_returns_20261113/control/MODULE_IDEAL.md)
gives a complementary all-context gate: a fixed uniform module with
nonprimitive ideal cannot support a count difference2 in any such context.
Its scaled path-cover profile detects arithmetic obstruction even for
asymmetric modules. Our T1 universe permits different incidences to selected
and unselected old vertices and is not identical to that uniform-module
setting. A candidate compensator must state which uniform modules remain
protected, which selective incidences change, and how its selected budget
is transported. Primitive arithmetic content still does not give a positive
physical construction or arbitrary-target entry.
