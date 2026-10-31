# Session 132 independent audit

**Result: PASS within the stated mathematical and finite-test scopes.** This
packet was written separately from the producing implementations. Its checkers
read the emitted artifacts and use explicit exceptions, so `python -O` does not
remove checks. No claim of fraction minimality, efficient execution, universal
arbitrary-seed classification, or proof of an unprovided prime-group relation is
made.

## 1. The actual A000568 FRACTRAN artifact

Final guarded-source artifact:

- 5,543 ordered reduced fractions; 3,979 control states; 36 data registers.
- Input `157 * 5^n`, for every integer `n >= 0`.
- Successful halt control 37991; output `v_3`, with error exponent `v_2=0`.
- Fraction-file SHA256:
  `bcbc9421a05646214d7c7428715df9882bc959204cd6ce7c7e9e9e1b2ac95cfa`.
- IR SHA256:
  `1ad6b359feeaa43963edc4a90df4fc762d7abb1627117ba73e5dac4a901d06a5`.

The earlier 5,471-fraction draft was superseded during the audit by a guarded
same-cycle triangular term. The final artifacts, not that draft, were rechecked.

`audit_fractions.py` independently parses the actual ordered numerator/denominator
text. It verifies primality and disjointness of every data/control prime,
reducedness, one source and one target control per fraction, exclusion of control
self-transitions, and equality with the exact ordered translation of every IR
instruction. In particular, every nonzero DECJZ row precedes its zero row.
The runtime exponent interpreter factors the *actual fractions* and does not
execute or consult IR instructions. Its source-control index removes only rows
whose required control prime is absent, while preserving all retained row order.

Final independent exponent replay:

| n | Output | Actual fraction steps |
|---:|---:|---:|
| 0 | 1 | 784 |
| 1 | 1 | 4355 |
| 2 | 1 | 7802 |
| 3 | 2 | 42606 |
| 4 | 4 | 77994 |
| 5 | 12 | 699757 |

A full literal-integer scan of all 5,543 rows at every step independently verifies
n=0: output 1 after 784 steps, maximum integer size 79 bits. This is an additional
check, not a substitute for the exact exponent semantics or the proof below.
The producer's broader n=6 IR and n=0..3 literal replays are separately recorded
in its own receipt; this independent packet does not mislabel them as its own.

### Why this proves actual first-eligible semantics

For a configuration `q_s * product p_i^r_i`, exactly one control exponent is one.
Every fraction denominator requires its own control prime. Thus inactive rows
are ineligible. An INC or JUMP has exactly its intended effect. At DECJZ with
positive r_i, both local rows are eligible but the nonzero row is earlier; at
zero r_i only the zero row is eligible. Each transition preserves the one-hot
invariant. HALT has no outgoing row and therefore literally halts.

The guard-cancellation trap is explicitly excluded: a self-transition would
cancel the source control between numerator and denominator and could make a
row eligible in unrelated states. Every actual generated row is reduced and
every source/target control is distinct. This is a proof for all reachable
configurations, not an inference from matching sequence values.

### Arithmetic compiler and source proof

`audit_compiler.py` independently dispatches the four IR instruction forms and
checks 427 arithmetic/branch/nested-loop cases, including truncated subtraction,
zero divisor traps, zero powers, comparisons, assignment reuse, and nested
scratch reuse. It recompiles both source files and checks that their complete
IR exactly reproduces the stored artifacts.

The macro contracts were also reviewed directly. CLEAR and MOVE drain a
nonnegative source; COPY restores the source from a distinct scratch; MULTIPLY
adds a preserved operand while decrementing a copy of the other; SUBTRACT
consumes paired copies to produce truncated difference; comparison consumes
copies in lockstep; DIVMOD repeatedly subtracts a positive divisor; POWER
performs a bounded number of multiplications. All named registers are allocated
before scratch, all live expression temporaries remain allocated, and each
reused scratch is initialized before its value is read. These contracts compose
by syntax induction.

For the A000568 source, the mixed-radix code enumerates every multiplicity vector
`0 <= m_d <= floor(n/d)` for odd d exactly once. The weight test retains exactly
odd partitions of n. The exponent accumulated is

    sum_d m_d(d-1)/2 + sum_d d*binom(m_d,2)
    + sum_(e<d) gcd(d,e)*m_d*m_e,

and the denominator is `product d^m_d * m_d!`. Thus each summand is the integral
conjugacy-class size `n!/z` times the exact number of fixed tournaments. The final
average is integral by Burnside. The explicit positive-multiplicity guard keeps
all intermediate source subtraction nonnegative, including before multiplication.
For n=0 the sole empty vector contributes 1.

Every source loop is bounded by n, a multiplicity, or the finite mixed-radix
product, except Euclid's loop whose positive second argument strictly decreases.
Every actual divisor is positive. Each arithmetic macro terminates by its own
decreasing counter. Therefore the source, IR, and actual fixed FRACTRAN list
halt correctly for *every* permitted n. No literal integer replay at infeasible
orders is needed for this mathematical statement. The producer's deliberately
loose `3979*(H+1)^8` bound was reviewed: its H covers all source and scratch
registers, and the stated structured loop depth is sufficient.

The added valuation-readout lower bound was also checked. Every step increases
a designated prime valuation by at most the fixed maximum numerator-minus-
denominator exponent c. Summing gives `T >= (Y-v_p(N_0))/c`. Here c=1, so the
large program requires at least A_n steps and literal16 at least ell_n-1. The
latter is asymptotically optimal in this fixed-program, bounded-initial-output,
valuation-readout model only; it says nothing about binary-output complexity.

## 2. The small programs

`audit_small_programs.py` reads the fixed 16-fraction tuple via AST literal
extraction and the 22-fraction file directly, importing neither producer.

### Sixteen fractions: literal metagraph-line count

Independent literal ordered scans for n=3..7 validate the exact terminal integer
`59 * 2^ell_n * 3^(n-1)` and exact step formula

    ell_n = 2^(binom(n-1,2)-1),
    steps = 6*ell_n - 5 + (7t^2+29t)/2,  t=n-3.

The all-input phase proof also checks: one outer step decreases the remaining
input counter, restores growth G, performs G doublings, and increases G by one.
Every transfer drains a source not replenished in that phase. Starting from
Y=1, G=2 gives the stated closed form. These are literal line instances, not
unlabelled tournaments or converse-merged nodes.

### Twenty-two fractions: pure-blue recurrence stream

Independent literal execution validates all 25 complete checkpoint integers,
ending at n=24, output 60, after 10,686 fraction steps. All 64 zero/one initial
six-register windows were also independently checked. The seven phase transfers
send `(a,b,c,d,e,f)` to `(b,c,d,e,f,a+e)` and clear scratch. The exact time between
checkpoints is `8+2*(2a+b+c+d+2e+f)`. Each transfer terminates for arbitrary
nonnegative register values; the stream itself intentionally never halts.

The sequence interpretation follows from the stated inherited pure-blue normal
form, rather than being proved merely by the finite tournament census. The
producer's current README explicitly preserves that dependency and supersedes
the invalid historical affine-label shortcut. No independent re-proof of the
entire inherited atom-classification theorem is claimed in this audit.

`audit_geometry.py` independently factors the 22 fractions, reconstructs the
23-by-22 valuation matrix, and evaluates its deleted-prime-23-row determinant
using fraction-free Bareiss elimination. The determinant is +1. Consequently
rank=22, right kernel is zero, and the left kernel is generated by total control
exponent. Thus **no nonempty periodic orbit is possible from any positive
integer**, even an arbitrary multi-control seed. This does not imply or classify
arbitrary-seed termination. The final microstep clock was independently checked:
all 22 weighted valuation increments equal one and all 25 cumulative checkpoint
times match. The exact inequalities `p^12 >= 2^w_p` for positive weights also
prove `N_t^12 >= 2^(G(N_0)+t)` on every trajectory. This is a lower magnitude
envelope, not monotonicity of N. Both semantic control and transfer graphs are
strongly connected; neither is being identified with a tournament.

## 3. Tournament counting and folding

`audit_counts.py` independently traverses **all** integer partitions through
n=25, builds an explicit representative permutation, and computes the signed
permutation on unordered vertex pairs. It does not use the odd-partition/gcd
fixed-point formula. Across 9,296 conjugacy classes, its direct signed-orbit
Burnside computation verifies every reported A_n, S_n and M_n through n=25,
including `M_n=(A_n+S_n)/2`.

The allowed support was independently checked and proved by edge-orbit parity:
ordinary fixed tournaments require all vertex cycles odd; converse-twisted
fixed tournaments require every nonsingleton cycle to have length 2 modulo 4
and at most one singleton. Grouping doubled odd parts then gives
`S_(2m)=F_m(4,1)` and `S_(2m+1)=F_m(4,2)`. Independent pair-layer Burnside sums
verify those identities through the same range.

The interaction-state insertion recurrence is combinatorially sound: a root
cycle of odd length d contributes `(k-1)!/(k-d)!` choices, its internal exponent,
and the stored interactions with prior cycles. Updating every residual h_j by
`gcd(d,j)` accounts exactly once for each cross-cycle orbit contribution.
The strong-component identities `A=1/(1-C)` and `S=A(z^2)*(1+D)` follow from the
unique ordered strong-component decomposition and reversal of that order under
converse. Weighted unpointing by the number of child vertex orbits is necessary;
plain division by n+1 is not generally correct.

### The fold is structural, not just equal cardinality

Choose an anti-involution tau of a self-converse tournament. Its 2-cycles can be
written uniquely as `(x_i,y_i)` with `x_i -> y_i`; at most one vertex is fixed.
The induced tournament A on the winners x_i, together with the symmetric cross
relation `x_i -> y_j iff x_j -> y_i`, gives a tournament and an undirected graph.
Encode that graph equivariantly as edge disagreement between A and a second
ordered tournament. A fixed vertex contributes one binary vertex color.
The construction has an explicit inverse: copy A oppositely onto the y_i,
insert the prescribed symmetric cross relations and internal winner arcs, and
insert fixed-vertex arcs from the colors.

Choice independence holds on **isomorphism classes**, not as a canonical choice
of labelled anti-involution. Tournament automorphism groups have odd order:
an involution would swap two vertices and reverse their connecting arc, and
Cauchy's theorem excludes even group order. The automorphism/anti-automorphism
extension has order twice an odd number, so its anti-involutions are Sylow
2-subgroups and are conjugate. A conjugator can be chosen in Aut(T): if it is
anti, multiply it by the source anti-involution without changing the conjugation.
Hence every chosen fold gives the same pair-layer isomorphism class. This proves
the bijection itself, not merely the numerical F_m identity.

The final `converse_folding/THEORY.md` and checker were subsequently reviewed
directly. Their marked groupoid equivalence and unmarked class-bijection
distinction are correct, including the C3 automorphism-loss control. The final
checker is rerun under `-O` in `folding_checker_audit.json`; this is a rerun of the
producer's finite checker, not represented as separately written audit code.

## 4. Square-sum graphs and supplied prime groups

`audit_prime_graphs.py` independently constructs all unordered pairs in G_n,
for every n=1..97. It checks the complete edge lists, degree distributions,
components, articulation points and bridges (by deletion, not the producer's
Tarjan implementation), girth (edge-deletion shortest paths), odd girth (the
parity-lift shortest odd closed walk), odd-root components, and arrival roots.
All positive finite Hamiltonian witnesses and all 14 prime-order Hamiltonian
cycle certificates check exactly. The G_31 spanning open-ear decomposition has
old distinct endpoints, fresh interior vertices, valid edges, and all 31 vertices.

The all-order threshold proofs in `prime_graphs/PROOFS.md` were checked:
connectedness iff n=1 or n>=14; girth absent through 14, then 5 through 29, then 3;
2-vertex-connectivity iff n>=31. The key tail step, at least two earlier
neighbors for n>=19, follows from its finite 19..24 base and the stated square
inequality for n>=25. None implies Hamiltonicity. The producer's separate
unpruned subset-DP checker exhausts path/cycle counts through32 and is rerun
under optimization in this audit; its receipt is retained separately.

The raw groups have 25 entries, only 24 distinct primes, duplicate 71 and missing 67.
No correction is silently made. Both prime-order `p_i -> G_(p_i)` and rank-order
`p_i -> G_i` readings are preserved. Only under the first does supplied group A
match exactly the disconnected regime. The tested coarse graph flags do not
classify B/C/D/E; this is not a refutation of every possible arithmetic relation.

## Reproduction

From this directory, preserving the sibling research directories:

```sh
python3 audit_fractions.py ../fractran/a000568 --through 5
python3 audit_fractions.py ../fractran/literal_lines --through 3
python3 audit_compiler.py
python3 audit_small_programs.py
python3 audit_geometry.py
python3 audit_counts.py
python3 audit_prime_graphs.py
python3 -O ../prime_graphs/check.py
```

All scripts use only Python's standard library. The non-fraction scripts accept
an optional first argument pointing to the research root. No network, repository
write, publication, or external solver is required.
