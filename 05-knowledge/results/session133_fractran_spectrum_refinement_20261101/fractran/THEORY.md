# More finite FRACTRAN tournament counters, and a guarded compiler refinement

**Status: PROVED arithmetic/program semantics and all-input termination;
FINITE-EXACT replays.** The core-floor and extremizer identifications use the
separately audited [core reduction theorem](../core_reduction/THEORY.md).
No universality slogan, sequence lookup table, input-dependent program list,
or conjectured spectrum completeness is used.

## 1. Fixed programs and exact interfaces

All eight compiled sources have one natural input `n`. Except for the last row,
this is the tournament order. The actual input integer is `entry_prime*5^n`.
At successful halt the exponent of 3 is the answer, the exponent of error prime
2 is zero, and the stated halt-control prime is present. Other work registers
may retain values; a pure-power output is **not** claimed for compiled lists.

| Artifact stem / count | Meaning | Fractions | Entry | Halt |
|---|---|---:|---:|---:|
| `a000568` / A_n | tournament isomorphism classes | 4294 | 149 | 26879 |
| `self_converse` / S_n | self-converse isomorphism classes | 4744 | 163 | 30113 |
| `converse_merged` / M_n | classes also identifying converse | 5152 | 179 | 32993 |
| `strong` / C_n | strongly connected tournament classes | 5663 | 199 | 36677 |
| `pair_tournaments` / P_n | ordered pairs of tournament layers, same vertices | 4296 | 149 | 26891 |
| `colored_pair_tournaments` / Q_n | ordered pairs plus a 2-coloring of vertices | 4437 | 157 | 27941 |
| `colored_tournaments` / B_n | tournament plus a 2-coloring of vertices | 4436 | 157 | 27919 |
| `extremizer_tail` / e(n+6) | minimum-H reduced-core extremizers | 829 | 53 | 4229 |

In the last row, the input is `t=N-6>=0`; it returns e(N), not e(t).
The eight exact dictionaries and sizes are in `VERIFIED.json`, not estimates.
Each `*.fractions.00.txt` is the entire fixed ordered fraction list for that
program. Metadata explicitly lists the ordered instruction/state shards;
all files are below 100KB, including the strong-program IR, which has two
instruction shards. `build_sources.py` assembles sources **before** compilation;
no source assembly, procedure call, array, recursion, or hidden oracle occurs
inside any of the emitted programs.

Boundary conventions: A_0=S_0=M_0=P_0=Q_0=B_0=1, C_0=0; a singleton is strongly
connected, so C_1=1. Coloring is a genuine vertex mark, not a distinguished
root. Pair layers are ordered and their two color names are fixed, not quotiented.

## 2. Shared bivariate Burnside kernel

For an odd-part partition with multiplicities m_d, put

```
z = product_d d^m_d m_d!,
r = sum_d m_d,
E = sum_d m_d(d-1)/2 + sum_d d*binom(m_d,2)
    + sum_{d<e} gcd(d,e)m_d*m_e,
F_n(q,x) = sum_{sum d*m_d=n} q^E*x^r/z.
```

The general semantic interpretation needed here is q=2^k, k>=1, and positive
integer x: F counts isomorphism classes of ordered k-tuples of tournaments on
one vertex set together with an x-coloring of its vertices. The action relabels
all layers and colors simultaneously.

A permutation with an even vertex cycle cannot fix a tournament. For odd cycles,
there are E unordered-pair orbits, with no endpoint reversal: each layer supplies
2^E independent orientation choices. Each vertex cycle has one of x colors.
Burnside gives F_n(2^k,x). This proves integrality for every specialization used,
not for arbitrary real q,x. It gives A=F(2,1), P=F(4,1), Q=F(4,2), B=F(2,2).

The unchanged mixed-radix kernel enumerates each tuple
`0<=m_d<=floor(n/d)` exactly once, using radices `floor(n/d)+1` at odd d.
It retains exactly weight n. It computes E, r and z using only fixed loops,
Euclid remainder, addition, multiplication, powers and natural subtraction.
It accumulates `(n!//z)*2^(kE+c*r)` for x=2^c, then divides by n!.

There are two different exactness obligations:

1. `n!/z` is an integer for **each** accepted partition: it is a permutation
   conjugacy class size. Thus each accepted summand is computed exactly.
2. The **final** division by n! is exact by Burnside. This does not require
   each rational partition contribution to be an integer.

All n=0 empty-product cases are included. Denominators used for radix decoding,
Euclid on positive second argument, z and n! are positive. Every multiplicity-1
expression is guarded by positive multiplicity. The code makes no even restart.

## 3. Self-converse and merged classes

Let reversal rho complement all tournament arcs. On labeled tournaments the
number of self-converse **unlabelled** classes equals
`(1/n!)*sum_pi |Fix(pi*rho)|`. For a self-converse orbit, anti-isomorphisms form
a torsor of size |Aut(T)|; multiplying by its n!/|Aut(T)| labelings contributes
exactly n!. Non-self-converse orbits contribute zero. This is the fixed-orbit
Burnside identity, not the number of individually self-reversing labeled graphs.

Admissible permutations for pi*rho have only cycles of length 2d with d odd,
and at most one singleton. Odd cycles of length at least three fail on internal
pair orbits; lengths divisible by four fail on the opposite-pair orbit; two
singletons fail on the pair joining them. Conversely the stated cycle types
satisfy every signed pair-orbit parity requirement.

Write n=2m+epsilon, epsilon in {0,1}, and let d_1,...,d_r be the odd half-lengths.
Their centralizer factor is 2^r*z. Internal pair-orbit counts sum to m;
cross-cycle counts are `2*sum_{i<j}gcd(d_i,d_j)`; the singleton contributes
one orbit for each even cycle. Hence the fixed-count term divided by the
centralizer is

```
2^(m + 2*sum gcd + epsilon*r)/(2^r*z)
  = 4^E * 2^(epsilon*r)/z.
```

Therefore `S_(2m)=F_m(4,1)` and `S_(2m+1)=F_m(4,2)` for every m>=0.
`self_converse_source.py` retains the original parity, halves n, and evaluates
this one kernel. No anti-automorphism of even order has been mistaken for an
automorphism. Dropping the parity mark already fails at n=3 (1 instead of 2).

Reversal is an involution on A_n classes with S_n fixed points, so
`M_n=(A_n+S_n)/2`; this proves both the parity and the final exact division.
`converse_merged_source.py` uses a two-pass phase loop around one static kernel,
then divides the accumulated answer by two. Phase zero uses original n and
exponent E; phase one uses floor(n/2) and exponent 2E+(n mod 2)r.

## 4. Strongly connected classes by explicit finite composition inversion

Every tournament has a unique ordered sequence of nonempty SCCs, with all arcs
between components oriented forward. Conversely every ordered sequence of
strong tournaments has exactly that SCC sequence; its isomorphism type retains
the order and each component type. Thus, in ordinary formal power series,

```
A(z)=1/(1-C(z)),   C(z)=1-1/A(z),
C_n=sum_{k=1}^n (-1)^(k+1)
       sum_{n_1+...+n_k=n; n_i>=1} A_(n_1)*...*A_(n_k),   n>=1.
```

`strong_source.py` enumerates the `2^(n-1)` cut masks between the n positions.
The successive runs are exactly the composition parts. At each cut it sets the
kernel's input to the run length, evaluates the same fixed A kernel, and
multiplies its result into `component_product`. The product is accumulated in
`positive` for odd k and `negative` for even k. The final answer is their
natural subtraction. This is exact, since the SCC identity proves
`positive-negative=C_n>=0`; truncation cannot hide a negative intermediate
answer. There is no recursive call to A and no stored table of earlier C values.
At n=0 the program returns zero before forming n-1.

The outer mask limit, position loop, finite run sizes, and finite kernel loops
prove termination for every n. Unary costs can be large: this is a fixed
transparent enumerator, not an efficiency claim.

## 5. Compiler refinement and guard-preserving jump threading

`baseline_compiler.py` is the exact Session132 compiler, copied with a source
hash. Its register-machine semantics use INC, DECJZ, JUMP and HALT and were
already proved in Session132. The new subclass changes only:

* constant assignment: clear the named target, then emit its fixed increment
  count, without a copied temporary;
* `x += c`: emit c increments (c a nonnegative source constant);
* `x += expression`: evaluate the RHS into a fresh scratch first, then consume
  that scratch into x. This preserves RHS-before-write even for `x += x`;
* thread chains of silent JUMPs and delete unreachable control states.

Jump threading needs a nontrivial guard repair. A data transition redirected
to its own source would cancel its control prime after fraction reduction.
Whenever this happens, the optimizer inserts a private one-step JUMP bridge.
For example a CLEAR loop cannot become the uncontrolled fraction `1/p`.
A pure-JUMP cycle is explicitly rejected by the optimizer, never silently
converted to a halt. No such cycle occurs in these sources.

For the jump-threading pass on the already-shortcut IR, every retained data
operation has exactly the same register mutation and branch as in its
pre-threading program. Only zero-register-effect JUMPs disappear, and bridges
restore a separate control state when needed. Thus this pass preserves full
register trajectories up to silent steps. The earlier arithmetic shortcuts
preserve named-variable values at source-statement boundaries, not their
internal scratch traces. Composing these contracts preserves source semantics
and halting/error results. Eliminated macros also need fewer scratch registers.

For data primes p_i and distinct state primes q_s, emit

```
INC(i,t):        p_i*q_t/q_s
DECJZ(i,t,u):    q_t/(p_i*q_s), then q_u/q_s
JUMP(t):        q_t/q_s
HALT:           no row.
```

Exactly one state prime has exponent one. Inactive rows are ineligible.
For a positive counter both DECJZ rows are integral, and the nonzero row comes
first; at zero only the second row is integral. All transitions have different
source and target states, and every fraction is checked reduced. This proves
one IR step equals one literal first-eligible FRACTRAN step, on every input.

**Measured refinement on the unchanged A source:** 5543 -> 4294 fractions
(1249 fewer, about 22.53%), 3979 -> 2914 states, 36 -> 34 registers; fraction
text shrinks from 67611 to 51468 bytes. Its input changes from 157*5^n to
149*5^n because of static prime reallocation. These are exact artifact counts,
not fraction-minimality or optimal-complexity claims. The report also compares
baseline/refined named registers and step counts for every new source.

## 6. Two direct small core-floor counters

Here the mathematical input is an order N of the singleton-SCC-deleted core,
and H counts directed Hamiltonian vertex words, not tournament isomorphism
classes. The [core reduction theorem](../core_reduction/THEORY.md) proves

```
g(3)=3, g(4)=5, g(5)=9,
g(6+3k+r)=b_r*3^k,   (b_0,b_1,b_2)=(9,15,25), k>=0, r in {0,1,2}.
```

The exact TAIL14 list is

```
13/55, 9728/11,
17/65, 622592/13,
77/85, 637534208/17,
23/133, 41/19,
783/46, 31/23,
23/29,
74/93, 19/31,
31/37.
```

**Input:** `11*5^(N-6)`, N>=6. **Terminal integer:** exactly `41*2^g(N)`.
To cover every feasible positive order N>=3, prepend six rows:

```
47/215, 328/43,
53/235, 1312/47,
11/265, 20992/53.
```

The resulting FULL20 list has **input** `43*5^(N-3)` and the same pure terminal
output `41*2^g(N)`. Neither interface accepts infeasible core orders 1 and 2;
empty core order 0 has g(0)=1 by convention but lies outside this input family.

### Invariant and termination

Data registers are Y=v_2, T=v_3, N_remaining=v_5 and K=v_7. In TAIL14, states
11,13,17 consume successive groups of three input units; the third increment
adds one to K and returns to 11. If the input is exhausted at these respective
states, the zero row initializes Y to 9,15,25 and enters growth state 19.
Consequently K=k, Y=b_r, and N_remaining=0.

At state 19, a positive K is decremented and enters 23; zero K enters halt41.
The 23/29 loop consumes Y units into **three** T units each. Its zero exit
enters31. The31/37 loop consumes each T into one Y and returns to19 when empty.
Thus one growth iteration triples Y and decreases K. The invariant at19 is
`Y=b_r*3^j, K=k-j, T=N_remaining=0`. Every transfer loop decreases its source.
This proves termination and the exact terminal integer, not just its valuation.
The six prefix rows either halt with Y=3,5,9 at N=3,4,5, respectively, or
consume three input units and enter the tail program at N>=6.

All branch denominators have their state prime and the tested row precedes its
zero row. The repeated factor 3^3 in 783/46 is essential: its valuation is
three, not one. Every listed fraction is reduced and its control guard remains.

### Exact runtime and relevance to {7,21}

For N=6+3k+r, tail initialization takes 3k+r+1 steps; a tripling from Y to3Y
costs 8Y+3 steps; final halt costs one. Therefore

```
TAIL14: T = 4*(g(N)-b_r) + 6k+r+2,
FULL20: T = 4*(g(N)-b_r) + 6k+r+5,   N>=6.
```

FULL20 takes N-2 steps for N=3,4,5. Runtime is Theta(g(N)), the best possible
asymptotic order for a fixed FRACTRAN list that must write g(N) as a prime
valuation from a bounded initial valuation: each row has a fixed bounded
valuation increase. This is not a binary-computation or fraction-minimality claim.

The floor theorem and monotonicity give `H=7 => core order<=4` and
`H=21 => core order<=7`, since g(5)=9 and g(8)=25. The 20-row counter faithfully
computes the quantity used in these thresholds; it does **not** establish the
remaining small-core H spectra, and a scalar floor never proves spectral
exclusion by itself. Permanent omissions 7 and21 are inherited known results,
not newly resolved by these counters. The separate atlas supplies local spectra.

## 7. Extremizer-tail counter

The equality classification in the core theorem gives, for N>=6, ordered words
of a triangle SCCs and b strong order-four SCCs, where b=N mod3 is 0,1,2 and
`a=(N-4b)/3`. These two component types are each unique up to isomorphism, and
the SCC order is intrinsic. Thus e(N)=binomial(a+b,b), equivalently

```
e(3k)=1, e(3k+1)=k, e(3k+2)=k*(k-1)/2,  k>=2.
```

`extremizer_tail_source.py` takes t=N-6, computes residue=t mod3 and k=t//3+2,
and selects these expressions. The only final division is exact because
consecutive integers have even product. Every input t>=0 terminates. This
counts unlabelled minimizing **core tournaments**, not their Hamiltonian paths,
not all tournaments of size N, and not minimum strong tournaments in isolation.

## 8. Verification, independence and losses

Run `python3 verify.py`, save `VERIFIED.json` as `VERIFIED.normal.json`, then
`python3 -O verify.py` and compare them. All checks use explicit exceptions,
never assertions; there are no external packages or network requirements.
`run_program.py` reads the delivered fraction files and literally scans the
ordered list on the actual integer. It does not evaluate the source, regenerate
fractions, or use a sequence formula. The verifier also checks every on-disk
row and IR/state shard against deterministic compilation.

The report distinguishes source-native calculations through input12, primitive
IR replays, independent fraction-exponent replays, and literal whole-integer
replays. Large native values are **not** advertised as completed FRACTRAN runs.
Independent signed-pair permutation Burnside checks A and S through order7;
a full isomorphism census checks A,S,M,C through5; pair-orbit enumeration checks
P,Q through3. The direct core programs are replayed in both whole integers and
exponents through N=24, including exact terminal shape and exact runtime.

Hostiles include negative input, explicit zero division, a nontermination bound,
canceled control, pure-jump cycles, mandatory clear bridges, reversed DECJZ
ordering, the historical n=2/n=5 odd-cycle errors, lost S parity, naive converse
halving, unsigned strong compositions, infeasible core order, and starting the
core tripling recurrence at5 (which would predict27 instead of g(8)=25).
Independent agent tests separately stress jump guards and actual integer scans.

Connection/loss contract: tournament carriers map to scalar orbit counts;
layer/color flags are retained exactly when required by Burnside. The counts
forget individual incidence, H values, automorphism groups and extension maps.
Core replacement preserves size/H **inequalities**, not the exact H fibre.
The arithmetic refinements preserve source-boundary named states; jump threading
preserves its input IR trajectories up to silent steps. Compilation renames
arbitrary data/control primes. Prime labels themselves have no
mathematical reason to single out7 or21. No count has been used as a surrogate
for a spectrum membership predicate or for the open unbounded existence claim.
