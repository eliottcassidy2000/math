---
id: THM-3480
title: "Rule 30 staircase-transducer entropy and nonrectangular macroblock compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The physical zero-right-
  boundary h-step Rule 30 map has an exact ternary staircase transducer.  An
  exact integer SFT certificate bounds its reachable-state entropy by
  log_2(117/64)<7/8.  A charged diagonal-state chunk compiler consequently
  computes the nth center bit in O(n^2/log^2 n) word-RAM operations with a
  strictly smaller one-lookup cone tariff than the raw rectangular compiler.
  This claims no Rule 30 prize and no literature novelty.
source: root-rule30-next-targets-20260815
audit: >
  An independent hostile audit rederived the elementary quotient and serial
  cascade, height cocycle, factor-language envelope, exact integer SFT
  certificate (including the sharp certificate constant 41), charged table
  construction, padding/alignment, and asymptotic query tariff.  Independent
  whole-row controls passed; ordinary and optimized replays match the stored
  output byte-for-byte: ACCEPT.
depends_on:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
related:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
script: 04-computation/rule30_staircase_entropy_compiler_thm3480.py
output: 05-knowledge/results/rule30_staircase_entropy_compiler_thm3480.out
script_sha256: c6b58d4f5f0b62a70d4f4f2fe1cc0b9b54295c176efced290d19aaea8c3af59d
output_sha256: e2060bffeeba608760905946e33e2fecb6cce26f519aaa9a8fab5d0d3edbfde6
hash_basis: raw bytes
---

# THM-3480 -- Rule 30 staircase-transducer entropy and nonrectangular macroblock compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The raw-cone compiler of THM-3471 spends `s+2h` address bits to move an
`s`-cell block through `h` Rule 30 rows.  The two `h`-bit halos are necessary
if every lookup is presented as an independent rectangle.  They are not the
right boundary object when a finite row is scanned from its physical zero
exterior.  Left permutivity turns the shared halo into a diagonal staircase
state, and the states which are actually reachable from the zero exterior
have entropy strictly below one bit per simulated row.

The gain below is structural but not a new exponent.  It improves the exact
one-lookup cone tariff and exposes a new prize-facing object, the staircase
entropy.  It does not prove an asymptotic lower bound, a center-density result,
or a Rule 30 prize.

## 1. Model and theorem

Write a Rule 30 row as `x=(x_j)` and

```text
(F x)_j=x_(j-1) xor (x_j or x_(j+1)).                (1)
```

Rows in the algorithm have finite support and are explicitly padded by zeros.
Use the uniform random-access word-RAM of THM-3471: on input `n`,

```text
w=ceil(log_2(n+2)),                                  (2)
```

Boolean word operations, shifts, addition, and one addressed word load/store
have unit cost; the program is constant; there is no advice, oracle, or
uncharged preprocessing.

There is an exact three-state right-to-left transducer `E` for one Rule 30
step.  Its height-`h` serial cascade has a physical reachable state set `S_h`.
Let

```text
N_h=|S_h|,       eta=lim_(h->infinity) log_2(N_h)/h. (3)
```

Then:

1. the limit in (3) exists;
2. an exact finite graph certificate proves

```text
N_h <= 41(117/64)^h       for h>=12,
eta <= log_2(117/64) < 7/8;                           (4)
```

3. for `h=8m`, a reachable staircase state therefore needs at most
   `7m+6` binary address bits; and
4. a charged table on one such state and `b` fresh row bits advances the row
   by `h` times while emitting `b` new bottom-row bits.

For a choice specified in Section 6, the number of main chunk queries is

```text
Q(n)=(7/2+o(1)) n^2/w^2.                             (5)
```

Including construction of every table and the final short remainder gives

```text
T(n)=O(n^2/log^2 n) word operations,
S(n)=O(n/log n) words                                (6)
```

for the single center bit `c_n`.  Small inputs use direct packed simulation.

## 2. The elementary ternary staircase

Scan a top row from right to left.  Immediately before consuming `a=x_j`,
the raw state for one Rule 30 step is the pair

```text
(x_(j+1),x_(j+2)).                                   (7)
```

The two raw states `10` and `11` have the same output and the same future
response.  The exact quotient has three states

```text
A=00,       B=01,       C=1*.                        (8)
```

Put `e(A)=0` and `e(B)=e(C)=1`.  The emitted bit belongs at site `j+1` and is

```text
z=a xor e(q).                                        (9)
```

The next state is

```text
delta(A,0)=A,  delta(B,0)=A,  delta(C,0)=B,
delta(A,1)=delta(B,1)=delta(C,1)=C.                 (10)
```

These equations are just (1): the emitted bit is
`a xor (x_(j+1) or x_(j+2))`, and the new raw pair is
`(a,x_(j+1))`.  Thus no dynamical information has been guessed or averaged.

For `h` steps, cascade `h` copies of `E`.  If

```text
q=(q_1,...,q_h),       z_0=a,
```

one external transition is

```text
z_i       = z_(i-1) xor e(q_i),
q_i'      = delta(q_i,z_(i-1))       (1<=i<=h).      (11)
```

The final carrier `z_h` is the emitted `F^h` bit.  Induction over the stages
proves this because each copy of `E` consumes the output stream of the
preceding copy.  The offset grows by one at each stage: after consuming a
length-`2h+1` cone from right to left, the last carrier is the direct depth-`h`
center.

There is also a useful two-mask form.  Let

```text
E_i=[q_i != A],       C_i=[q_i=C],
p_i=a xor E_1 xor ... xor E_(i-1).                   (12)
```

Then (10)--(11) are exactly

```text
C_i'=p_i,       E_i'=p_i or C_i,
z_h=a xor E_1 xor ... xor E_h.                       (13)
```

So the diagonal update is a prefix-parity transport, not a `4^h` raw-context
black box.  The companion checks (13) independently against the ternary
transition on every state used by the finite certificate.

## 3. Physical reachability and the height cocycle

The zero exterior initializes every stage at `A`.  Define

```text
S_h=Orb_{Phi_0,Phi_1}(A^h),                          (14)
```

where `Phi_a` is (11).  This is precisely the set encountered while scanning
an arbitrary finite binary row from its zero right tail.  It is not the full
universal state cube `{A,B,C}^h`.

Height addition is serial concatenation.  If a height-`h+k` state is split as
`(p,q)` with lengths `h,k`, one transition first applies `Phi_a` to `p` and
then uses its final carrier as the input to the `k`-cascade on `q`.  Therefore

```text
S_(h+k) -> S_h x S_k,
r |-> (first h coordinates, last k coordinates)     (15)
```

is an injection.  The second projection is in `S_k` because it starts at
`A^k` and is driven by a binary output word from the first cascade.  Hence

```text
N_(h+k) <= N_h N_k.                                  (16)
```

Fekete's lemma applied to `log_2 N_h` proves the existence of (3), with the
usual infimum formula.  Equation (15) is the exact nonrectangular height
cocycle; it preserves the ordered diagonal boundary rather than widening a
top rectangle as in THM-3471.

More locally, let `R_k=S_k` as a language of ternary words.  Every contiguous
length-`k` factor of a word in `S_h` lies in `R_k`: the corresponding subcascade
starts at `A^k` and is driven by the carrier word entering its first stage.
This factorial-language observation makes a fixed finite orbit into an
all-height entropy certificate.

## 4. An exact finite SFT certificate

Build the directed graph `G_13` as follows:

```text
vertices: R_12,
edge u->v: a word r in R_13 has prefix u and suffix v. (17)
```

Exact breadth-first construction gives

```text
|R_12|=5873,       |R_13|=10483.                    (18)
```

Let `A_13` be the `0/1` adjacency matrix of this graph.  The companion
constructs a positive integer vector `v`, in lexicographic `R_12` order, with

```text
min(v)=15,       max(v)=457,       sum(v)=854695,
64 A_13 v <= 117 v.                                 (19)
```

This is an exact Collatz certificate, not a floating-point eigenvalue.  Its
canonical `state:weight` payload has SHA-256

```text
52d2b23844f52fb6af9666de6034440887aafeb67f07ec2ffc9a95202a576a94. (20)
```

The vector is generated by a committed monotone integer recurrence and every
coordinate of (19) is then checked.  The `R_12` and `R_13` languages are also
committed by separate hashes in the script.

For `h>=12`, every word of `S_h` determines a path of length `h-12` in
`G_13`.  Since the all-one vector is at most `v/15`, (19) gives

```text
N_h <= 1^T A_13^(h-12) 1
    <= (854695/15)(117/64)^(h-12)
    <= 41(117/64)^h.                                 (21)
```

The last inequality is checked after clearing denominators.  Consequently
`eta<=log_2(117/64)`.  The entirely integral comparison

```text
117^8=35114532758015841 < 2^55=36028797018963968     (22)
```

proves the strict `7/8` bound in (4).  If `h=8m>=16`, then

```text
N_h < 41*2^(7m) < 2^(7m+6).                         (23)
```

Thus `7m+6` bits suffice for a dense ID of every physically reachable
staircase state.

The finite certificate is an upper envelope: a path in `G_13` need not be a
globally reachable staircase.  No equality or lower entropy estimate is being
inferred from (19).

## 5. The diagonal-state chunk table

For `h,b>=1`, define the scan-order table

```text
D_(h,b): S_h x {0,1}^b -> S_h x {0,1}^b.            (24)
```

It consumes `b` adjacent top-row bits from right to left through (11), returns
the new staircase state, and emits the corresponding `b` adjacent `F^h` bits.
This is a nonrectangular tile: `S_h` is the shared diagonal boundary and the
binary word is the fresh top boundary.  Neighboring calls pass the staircase
state instead of rereading a `2h`-cell raw halo.

For each fixed initial state, the input-word to output-word part of (24) is a
permutation.  Indeed, at one scan position (11) gives

```text
z_h=a xor e(q_1) xor ... xor e(q_h),                 (25)
```

so `a` is recovered from `z_h` and the current state.  Then the next state is
known, and induction recovers an entire chunk.  The diagonal quotient therefore
does not evade left-permutive rank; it compresses the physical boundary states,
not the `b` fresh data bits.

All preprocessing is charged.  One deterministic construction is:

1. allocate a direct marker array indexed by the base-three code of a ternary
   height-`h` state;
2. breadth-first search from `A^h`, computing both binary successors, and
   assign dense IDs to `S_h`; and
3. build (24) by extending chunk words one bit at a time from the two stored
   state transitions.

This costs

```text
P(h,b)=O(3^h+hN_h+N_h 2^b) word operations,
M(h,b)=O(3^h+N_h 2^b) words.                        (26)
```

The first term initializes the deterministic membership array.  The last term
is linear in the final table size: tables of all shorter word lengths form a
geometric sum.  Hashing, advice, and an uncharged minimal automaton are not
being assumed.

If a state ID and a `b`-bit chunk together occupy at most one word, one call to
(24) uses `O(1)` word operations and loads.  Its returned state ID and output
chunk obey the same one-word budget.

## 6. Charged all-n optimization

For sufficiently large `w`, put

```text
d=ceil(log_2 w),       L=w-8d-8,
m=floor((L-6)/14),     h=8m,
b=floor((L-6)/2).                                    (27)
```

Then

```text
7m+b+6 <= L.                                         (28)
```

By (23), the number of entries and the bits returned by one entry satisfy

```text
N_h 2^b < 2^(7m+b+6) <= 2^L.                        (29)
```

The deterministic ternary marker array is also sublinear.  The exact integer
inequality `3^8=6561<8192=2^13` gives

```text
3^h=3^(8m)<2^(13m)=2^((13/14+o(1))w)=o(n/w).        (30)
```

Moreover `2^L=O(n/w^8)`.  Equations (23), (26), and (29)--(30) show that all
table construction and marker-array costs are `o(n/w)` words and lower order
than the main simulation.

At each macro time, scan the explicitly padded active row from right to left,
starting at `A^h`, in `b`-bit chunks.  Induction using (11) proves that the
emitted row is exactly `h` Rule 30 times later.  Over `floor(n/h)` macrosteps,
the active widths are `2rh+O(h)`.  Hence the number of table calls is

```text
Q(n) = n^2/(hb)+O(n/b+n/h)
     = (7/2+o(1)) n^2/w^2,                           (31)
```

because

```text
h=(4/7+o(1))w,       b=(1/2+o(1))w,
hb=(2/7+o(1))w^2.                                  (32)
```

The fewer-than-`h` final rows use ordinary packed Rule 30 updates in `O(n)`
word operations.  Extracting the center cell then proves (6).

For comparison, a one-word raw rectangular address has `s+2h<=w+o(w)` and
therefore `sh<=(1/8+o(1))w^2`.  Its cone-query tariff is at least
`(8-o(1))n^2/w^2` within that declared model.  Equation (31) improves the
lookup throughput by the factor

```text
8/(7/2)=16/7.                                        (33)
```

This is not a lower bound for general Rule 30 algorithms, and (33) is a table-
query comparison rather than a machine-independent leading constant.

## 7. Preservation and loss ledger

| map | preserves | destroys / admits | needed boundary |
|---|---|---|---|
| raw right context -> `A,B,C` | exact response to every future left bit | distinction `10` versus `11` | scan orientation |
| `h` stages -> staircase state | exact `F^h` stream and height cocycle | internal spacetime away from the cut | ordered stage coordinate |
| universal staircase -> `S_h` | every finite zero-right-tail row | arbitrary bi-infinite right phases | physical zero exterior |
| `S_h` -> `G_13` path | a necessary factor condition and entropy upper bound | spurious locally compatible paths | exact vector (19) |
| `(state,chunk)` -> table output | bottom chunk and successor state | intermediate rows | `h,b` and chunk alignment |

Three hostiles delimit the result.

First, fixed-state chunk maps are full permutations, so the universal-input
Walsh/rank obstruction of THM-3471 survives.  Second, exact Myhill--Nerode
refinement finds every state of `S_h` distinguishable for `1<=h<=13`, with

```text
N_h=(3,7,16,35,71,141,272,517,971,1792,3263,5873,10483). (34)
```

This blocks naive further state merging at those finite depths but is not an
asymptotic lower bound.  Third, the zero-tail restriction is essential: the
algorithm says nothing about the state distribution of a universal bi-infinite
row or about the marked center orbit's density.

The next exact target is the opposite side of (4): prove a positive lower
bound for `eta`, or find a still smaller sufficient boundary carrier.  A
positive lower entropy would turn the persistent `log^2 n` ceiling into a
restricted diagonal-lookup obstruction; zero entropy would point toward a
genuinely stronger time compression.  Neither direction is claimed here.

## 8. Verification

The exact replay commands are

```bash
python3 04-computation/rule30_staircase_entropy_compiler_thm3480.py
python3 -O 04-computation/rule30_staircase_entropy_compiler_thm3480.py
```

The companion has no assertion-dependent gates.  It checks the elementary
transducer against every direct cone through `h=6`, the ternary/mask recurrence,
the reachable languages through `h=13`, factor closure, exact Myhill--Nerode
refinement, the committed `R_12/R_13` graph and integer certificate, fixed-state
chunk permutivity, and the integer address inequalities in (27)--(30).  The
graph proof uses only exact integer arithmetic; no numerical spectral radius is
a proof dependency.
