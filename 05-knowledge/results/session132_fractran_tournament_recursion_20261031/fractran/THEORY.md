# Actual finite FRACTRAN counters for tournaments and literal metagraph lines

**Status: PROVED mathematical/compiler semantics; FINITE-EXACT executable replays.**
There are two principal *fixed*, finite ordered lists, independent of input size. No oracle,
input-dependent instruction list, lookup table of sequence values, or appeal to
universality is being substituted for a program. The large list is deliberately
unoptimized; no claim of fraction minimality or efficient unary computation is made.

## 1. The 16-fraction line counter

In the inherited fixed-Hamiltonian-path model, let

- `X_n = {0,1}^{S_n}`, where `S_n={(x,y):1<=y<x<=n, x-y>=2}`;
- `a(t)=t+1` be all-tile complement;
- `mathcal L_n = X_n/<a>` be **literal line instances**, retaining multiplicity.

This is the carrier of THM-796,
`01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`,
section 1. It is not the set of tournament isomorphism classes, not the set of
converse-merged nodes, and not the set of distinct unordered endpoint pairs.
For `n>=3`, `|S_n|=binom(n-1,2)>0`, so all-tile complement is fixed-point-free.
Consequently the number of literal lines is

    ell_n = 2^(binom(n-1,2)-1),    ell_3=1,
    ell_(n+1) = 2^(n-1) ell_n.

This recurrence has an actual object map: endpoint deletion restricts tilings,
commutes with all-tile complement, and induces `mathcal L_(n+1)->mathcal L_n`.
There are `n-1` new free bits, so every line fibre has exactly `2^(n-1)` elements.
It counts the fibres but does not recover their tournament endpoints or colours.
At `n=2`, the cube has no free coordinate and complement is the identity; the
half-cardinality formula is deliberately **not** extended to this boundary.

Here is the entire ordered FRACTRAN program, also in `literal16.py`:

    19/221, 59/17,
    1265/57, 29/19,
    19/23,
    93/319, 37/29,
    29/31,
    41/185, 51/37,
    2107/82, 47/41,
    41/43,
    106/329, 37/47,
    47/53.

**Input:** `306*13^(n-3)`, for any integer `n>=3`.
**Output:** it halts at exactly `59*2^ell_n*3^(n-1)`.
Thus the sequence value is the exponent of 2, not the final integer itself.

### Register/phase dictionary

Data exponents are `Y=v_2`, `G=v_3`, `C=v_5`, `T=v_7`, `B=v_11`, `N=v_13`.
The one active control prime is in
`17,19,23,29,31,37,41,43,47,53,59`; 59 is the halting control.
At control 17, the invariant is

    Y=ell_(j+3), G=j+2, N=n-3-j, C=T=B=0.

- At 17, if `N>0`, decrement N and enter 19; otherwise enter 59 and halt.
- The 19/23 pair transfers each unit of G into one unit each of C and B.
- The 29/31 pair transfers B back to G. Now `C=G` and B is zero.
- At 37, if C is positive, decrement it and enter the doubling routine.
  If C is zero, increment G and return to 17.
- The 41/43 pair moves each Y unit into **two** T units.
- The 47/53 pair moves T back to Y and returns to 37.

Each data-transfer loop has its own decreasing nonnegative counter. Each outer
iteration decrements N. Hence this is a halting program for every allowed input,
not just for the tested orders. There is always exactly one control prime, so only
fractions in the active phase can apply. Each conditional's tested fraction occurs
before its unconditional zero branch. This proves literal first-eligible semantics.

After an iteration starting with `G=k`, Y has been doubled exactly k times, giving
`Y_new=2^k Y_old`, and G becomes k+1. The invariant follows by induction.

### Exact step count

For `t=n-3`, the number of fraction applications is

    6 ell_n - 5 + (7t^2+29t)/2.

At growth k, copying/restoring G costs `4k+2` steps; the k doublings cost
`6(Y_new-Y_old)+3k`; entry and final C-zero exit add two. Sum
`6(Y_new-Y_old)+7k+4` for `k=2,...,t+1`, then add the final halting transition.
The formula is verified independently against literal and exponent-vector replays.

## 2. A fixed 5,543-fraction program computing A000568

`a000568.fractions.00.txt` is the complete ordered list: **5,543 fractions**, in
67,611 bytes. `a000568.ir.json` is its readable register/control metadata and
3,979-instruction IR, in 92,927 bytes. There are 36 data registers. Both files are
fixed independently of n, and each is below 100KB.

**Input:** `157*5^n`, for any integer `n>=0`.
**Successful terminal output:** the exponent of 3 is `A000568(n)`; the control prime
is 37991, and the error-prime 2 has exponent zero. Other data registers may retain
work values. The output integer is not claimed to be a pure power of 3.

The full arithmetic WHILE source is `a000568_source.py`. It contains no calls,
arrays, recursion, external inputs, or precomputed sequence values. The following
finite enumeration is directly compiled, not merely described as computable.

### Integer-only Burnside algorithm

For odd cycle lengths d, write m_d for multiplicity and define

    E(m) = sum_d m_d(d-1)/2
           + sum_d d*binom(m_d,2)
           + sum_(d<e) gcd(d,e)m_d m_e,
    z(m) = product_d d^m_d m_d!.

The source evaluates

    F=n!,
    numerator = sum_(sum d*m_d=n) (F/z(m))*2^E(m),
    A000568(n) = numerator/F.

Every division displayed here is exact. `F/z(m)` is the number of permutations
with the specified cycle type. The final quotient is integral by Burnside.
The implementation's quotient opcode is floor division, but these specific calls
therefore return the required exact integer; radix/GCD calls legitimately use floor
quotients and remainders.

The program enumerates multiplicities by the mixed radices

    r_d = floor(n/d)+1     (d odd, 1<=d<=n),
    R = product_d r_d,
    code = 0,1,...,R-1.

Successive remainder/quotient operations decode exactly one vector
`0<=m_d<=floor(n/d)`, each vector once. Keeping only `sum d*m_d=n` gives precisely
all odd-part partitions. For `n=0`, the product is empty, R=1, and the empty vector
contributes 1. No even cycle accidentally enters the enumeration. This explicitly avoids the
2026-10-28 `01-canon/MISTAKES.md` correction of earlier A000568 prototypes:
an even restart had produced 10 instead of 12 at n=5, and a missing odd-cycle
filter had produced 2 instead of 1 at n=2. Both orders are mandatory controls here.

Inside an accepted code, the outer d loop accumulates internal/same-length orbit
terms and factors `d^m_d m_d!`. The inner e loop re-decodes the earlier digits and
adds the cross terms `gcd(d,e)m_dm_e`. Euclid's algorithm computes gcd using only
remainder and assignments. Finally the source forms the summand and adds it.

### Why the formula counts tournaments

A permutation with any even cycle cannot fix a tournament: its half-cycle power
interchanges two vertices, which would reverse their unique directed edge. If all
cycles have odd length, its order is odd and no pair-orbit can reverse its endpoints.
A cycle of odd length d contributes `(d-1)/2` internal unordered-edge orbits; two
cycles of lengths d,e contribute `gcd(d,e)` cross-edge orbits. Every edge orbit can
be oriented in exactly two ways, independently. Thus the number of fixed
tournaments is `2^E`, and grouping Burnside's permutation sum by cycle type gives
the source formula. This is genuine tournament enumeration, not a tournament
attached cosmetically to program labels.

### Source invariants and all-input termination

- `factorial` is the factorial accumulated through d-1.
- `limit` is the product of the already visited odd radices.
- Each `scan`/`previous` equals the quotient after consuming precisely its visited
  digits; multiplicity is the next digit.
- `weight` is the weighted sum of visited digits.
- In an accepted code, `denominator` and `exponent` are the partial product and
  orbit sum just described. The explicit `multiplicity>0` guard prevents even the
  intermediate expression `multiplicity-1` from becoming negative.
- `total` is the sum of all previously accepted code contributions.

The factorial, radix, d, e and multiplicity loops have explicit finite integer
bounds. There are exactly R code iterations. In each Euclidean loop the positive
second argument strictly decreases. All denominators are positive. Primitive
arithmetic loops decrease a copied counter, and powers use a decreasing exponent.
Therefore all valid inputs halt with the stated output.

For an intentionally coarse explicit bound, put

    M=binom(n,2),
    H=2+R+(n+1)^2+n!*2^M.

All reachable named and temporary registers are at most H: invalid weights are
at most n*ceil(n/2), accepted denominators are at most n!, accepted orbit exponents
at most M, and the entire Burnside numerator at most `n!*2^M`. Each compiler-loop
counter is a copy of one of these bounded values. The source has at most four
nested while loops; arithmetic expansion adds at most three nested transfer loops.
An elementary structured-loop bound therefore gives fewer than

    3979*(H+1)^8

primitive fraction steps (a deliberately loose extra factor is retained). This is
not a practical runtime estimate. For the small direct program, the previous exact
step formula is much sharper.

## 3. Compiler semantics, not a universality slogan

`compiler.py` lowers the finite source into just these instructions:

    INC(i,next)
    DECJZ(i,nonzero,zero)
    JUMP(next)
    HALT

There is a distinct data prime p_i for each register and a distinct control prime
q_s for each instruction, with the two prime sets disjoint. A configuration is

    q_s * product_i p_i^r_i.

Each instruction is replaced by the following ordered fractions:

    INC:    p_i*q_next/q_s
    DECJZ:  q_nonzero/(p_i*q_s), then q_zero/q_s
    JUMP:   q_next/q_s
    HALT:   no fraction.

**Guard condition:** no transition targets its own source state. Otherwise reduction
could cancel q_s and leave an unguarded arithmetic fraction. The compiler checks
this condition, and generated fractions are checked reduced. Its loop macros use
separate transfer and return states rather than unsafe self-loop fractions.

Exactly one control exponent is 1 and the others are 0. Every denominator requires
its own source control. Therefore no inactive instruction can apply. At DECJZ with
a positive counter, both local rows would be integral but the first is selected;
with a zero counter, only the second is integral. Each step preserves the one-hot
control invariant, so induction proves one IR step equals one literal FRACTRAN
step. HALT has no eligible denominator. Prime factorization/exponent simulation
is exact because divisibility and multiplication are coordinatewise conditions,
not an approximation to the giant integer.

### Arithmetic macro contracts

- CLEAR repeatedly decrements its target to zero.
- MOVE-ADD consumes a source into a distinct destination.
- COPY consumes source into destination plus a temporary, then restores source;
  destination equals the old source and the restoration temporary is zero.
- MULTIPLY preserves operands, iterates a copy of one operand, and adds copies of
  the other to a cleared destination.
- SUBTRACT preserves operands and writes `max(a-b,0)` by paired decrements.
- LESS-THAN preserves operands and consumes copies in tandem, writing exactly 0/1.
- DIVMOD guards divisor zero (explicit ERROR halt), then repeatedly subtracts the
  positive divisor while the current remainder is at least it.
- POWER starts at one and performs the required number of multiplications.

All names are allocated before scratch registers, and active expression temporaries
are not recycled until their values are no longer needed. Each scratch-reading macro
initializes that scratch before reading it. The static reused-register allocator is
not dynamic memory and does not create new primes during execution. The error
branch is unreachable on valid A000568 inputs but is tested on a zero-divisor control.

## 4. Reproduction and exact scope of evidence

From this directory:

    python3 verify.py
    cp VERIFIED.json VERIFIED.normal.json
    python3 -O verify.py
    cmp VERIFIED.normal.json VERIFIED.json
    python3 run_program.py a000568 3
    python3 run_program.py literal16 6

The last two commands read/apply actual fractions, with ordinary integer
modulus/division/multiplication and full first-eligible scans. They do not execute
the arithmetic source or call the compiler. The A000568 fraction file is parsed
from disk; the small literal list is a fixed tuple.

The verifier uses explicit exceptions, not Python `assert`, so optimization cannot
silently remove its checks. Its exact coverage is recorded in `VERIFIED.json`:

- 560 primitive arithmetic/branch cases, including independent integer checks;
- five hostile controls: negative input, zero divisor, a nontermination bound,
  canceled control, and incorrect zero-first conditional ordering;
- A000568 primitive replay for n=0..6: `1,1,1,2,4,12,56`;
- exact fraction-exponent replay for n=0..5, comparing every final register and step;
- literal whole-integer FRACTRAN replay for n=0..3;
- a separate signed-edge-orbit Burnside implementation over **all permutations**,
  with no cycle-partition/gcd formula, checks the same values through n=6;
- the 16-fraction program is checked with literal integers for n=3..7 and exact
  exponent vectors for n=3..7, including its exact terminal shape and step formula.

The tested n=6 large-program run takes 32,725,255 fraction/IR steps. This honestly
shows the cost of a simple unary compiler; higher-n values in the separate counting
packet are **not** advertised as replayed FRACTRAN outputs. The all-input result is
proved by the invariants/termination argument, not inferred from these finite tests.

## 5. Preserved information, losses, and prime relabeling

The large program preserves the number of tournament isomorphism classes, but
forgets each tournament and its path/line incidence. The small program preserves
the size and uniform extension-fibre count of the literal-line carrier, but forgets
endpoint classes, reflection colours, cycles, and Hamiltonian information. Neither
scalar output gives a canonical extension map on bare tournament isomorphism nodes.
The richer interaction-state recursion is developed in the accompanying counting
packet, rather than asserted as a scalar recurrence for A000568.

Any injective renaming of the finitely many data/control primes conjugates the
entire FRACTRAN orbit: exponent configurations, eligible row indices, halting and
outputs are unchanged under the new dictionary. Consequently these particular
prime labels carry no intrinsic moonshine or group-order classification. The
invariant structure is the guarded register dependency/control pattern, not the
numerical identities of 2,3,5,... . No arbitrary instruction graph is being called
a tournament; the actual tournament relation belongs to the counted input objects.


## 6. A genuine time lower bound — for this valuation output model only

Let a fixed FRACTRAN list have fractions `a_j/b_j` and designate output prime p.
Put `c=max(0,max_j(v_p(a_j)-v_p(b_j)))`. After T steps from initial integer N_0,

    v_p(N_T) <= v_p(N_0)+c*T.

This follows simply by summing the per-step valuation increments, whether or not
some steps decrease the output register. Thus if `c>0` and the required terminal
output is Y,

    T >= ceil((Y-v_p(N_0))/c).

If c=0 no larger output can be produced at all. This is a statement about a fixed
fraction list and a specified exponent readout, not a general complexity bound
for computing the same number in binary.

Direct inspection of every emitted fraction, also checked in `verify.py`, gives
`c=1` in both principal programs. Each has exactly one fraction that increments
its output exponent and one that decrements it:

- A000568 uses output prime 3, initially absent. Therefore **T>=A000568(n)**.
- Literal16 uses output prime 2, initially at exponent 1. Therefore **T>=ell_n-1**.
  Its proved runtime is `6ell_n+O(n^2)`, and the ratio to this elementary bound
  tends to 6. Consequently it has the optimal *asymptotic order* among fixed
  FRACTRAN programs with valuation readout and bounded initial output exponent.
  This does not say 16 fractions are minimal or the factor six is optimal.

The contrast matters: `ell_n` can be described by the small integer exponent
`binom(n-1,2)-1`, but requesting `ell_n` itself as a prime valuation makes the
terminal integer's bit length at least `ell_n*log_2(p)`. The constant-factor
statement is about that intentionally unary-like output contract, not about
fast binary calculation of the metagraph cardinality.
